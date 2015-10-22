#!/usr/bin/env python2.7
"""
smartSam2Fastq.py: turn sorted-by-names SAM input into a properly deduplicated
FASTQ. Also work around the bug in bwa mem where some even-length alignments to
reverse strands of alts will have incorrect bases for one of the middle two
bases.

accounts for both secondary and supplementary alignments of the same read

"""

import argparse, sys, os, os.path, random, itertools, string, re
import doctest

import pysam

# We need to do reverse complements
RC_TABLE=string.maketrans("ACGT", "TGCA")

def reverse_complement(dna):
    # Simplest way to do it according to
    # <http://stackoverflow.com/a/26615937/402891>
    return dna.translate(RC_TABLE)[::-1]

# Global SAM constants
BAM_FREVERSE=16
BAM_FREAD1=64
BAM_FREAD2=128

def parse_args(args):
    """
    Takes in the command-line arguments list (args), and returns a nice argparse
    result with fields for all the options.
    
    Borrows heavily from the argparse documentation examples:
    <http://docs.python.org/library/argparse.html>
    """
    
    # Construct the parser (which is stored in parser)
    # Module docstring lives in __doc__
    # See http://python-forum.com/pythonforum/viewtopic.php?f=3&t=36847
    # And a formatter class so our examples in the docstring look good. Isn't it
    # convenient how we already wrapped it to 80 characters?
    # See http://docs.python.org/library/argparse.html#formatter-class
    parser = argparse.ArgumentParser(description=__doc__, 
        formatter_class=argparse.RawDescriptionHelpFormatter)
    
    # General options
    parser.add_argument("--input_sam", type=argparse.FileType("r"),
        default=sys.stdin,
        help="input SAM in name-sorted order. Cannot be standard input.")
    parser.add_argument("--fq1", type=argparse.FileType("w"),
        default=sys.stdout,
        help="FASTQ file to save all the READ1 reads in")
    parser.add_argument("--fq2", type=argparse.FileType("w"),
        default=sys.stdout,
        help="FASTQ file to save all the READ2 reads in")
    parser.add_argument("--interleaved", action="store_true",
        help="write interleaved FASTQ to fq1")
    
    # The command line arguments start with the program name, which we don't
    # want to treat as an argument for argparse. So we remove it.
    args = args[1:]
        
    return parser.parse_args(args)

class Read(object):
    """
    Represent a Read as reconstructed from an alignment.
    
    """
    def __init__(self, sam_line):
        """
        Parse the given SAM line and construct a read.
        
        """
        
        # Save the line
        self.line = sam_line
        
        # Parse out the fields
        parts = sam_line.split("\t")
        # Get the template name
        self.template = parts[0]
        
        # Grab the flags
        self.flags = int(parts[1])
        
        # What end are we (1, 2, or 0 for unpaired)
        if self.flags & BAM_FREAD1:
            if self.flags & BAM_FREAD2:
                # This shouldn't happen
                raise RuntimeError("Alignment flagged as both READ1 and READ2")
            
            # Otherwise we're READ1
            self.end = 1
        elif self.flags & BAM_FREAD2:
            # We're READ2
            self.end = 2
        else:
            # We're unpaired
            self.end = 0
            
        # Grab sequence and qualities
        self.sequence = parts[9]
        self.qualities = parts[10]
        
        # Count edits
        self.edits = float("inf")
        for tag in parts[11:]:
            if tag.startswith("NM:i:"):
                self.edits = int(tag[5:])
        
        
        
        if self.flags & BAM_FREVERSE:
            # Flip to the other strand by RCing sequence and reversing qualities
            self.sequence = reverse_complement(self.sequence)
            self.qualities = self.qualities[::-1]
            self.is_reverse = True
        else:
            self.is_reverse = False
            
        # Grab the contig we mapped to
        self.contig = parts[2]
        
        # Say we are suspect if we're on an alt.
        self.is_suspect = self.contig.endswith("_alt")
            
    def get_name(self):
        """
        Produce a name for the read based on the template and the end.
        
        """
        
        if self.end == 0:
            # Unpaired reads are named for the template
            return self.template
        else:
            # Paired reads get /1 and /2
            return "{}/{}".format(self.template, self.end)
            
    def __eq__(self, other):
        """
        Two Reads are equal when they are the same read, but other meta-info
        about the original alignment may differ.
        
        """
        
        if not isinstance(other, self.__class__):
            return False
        
        
        if self.template != other.template:
            return False
        if self.sequence != other.sequence:
            return False
        if self.qualities != other.qualities:
            return Fasle
        if self.end != other.end:
            return False
            
        return True
        
    def __ne__(self, other):
        return not (self.__eq__(other))
        
    def __gt__(self, other):
        """
        A Read is greater than another read if it should replace the other read
        (they are the same end of the same template, and the other read is
        suspect while this one is not, or this one has a longer sequence and
        isn't suspect).
        
        """
        
        # Also just take ones with less edits if they aren't bad-looking
        return (self.template == other.template and self.end == other.end and 
            not (self.is_suspect and not other.is_suspect) and
            (len(self.sequence) >= len(other.sequence) or 
            self.edits <= other.edits))
            
    def __str__(self):
        """
        Turn this Read into a string for displaying.
        
        """
        
        mark = "?" if self.is_suspect else ""
        return "{} end {} on {}: {}{}".format(self.template, self.end,
            self.contig, self.sequence, mark)
            
def parse_MD_tag(md_string, offset, length):
    """
    Given the string value of an MD tag, the offset into the read at which it
    starts (due to e.g. hard clipping), and the total length of the read, return
    a dict from mismatch positions to replaced bases.
    
    Delete is not present since those obviously are not in the read.
    """
    
    # Fill this in with mismatches.
    to_return = {}
    
    # Split the MD tag between indels/substitutions and digits
    parts = re.split("([ACGTN^]+)", md_string)
    
    # Keep a cursor for what we need to update next
    cursor = offset
    
    for part in parts:
        if re.match("[0-9]+", part):
            # This is a number of matches or inserts
            part = int(part)
            
            # All these bases are matches or inserts
            
            # Move along to the base after that
            cursor += part
        else:
            if part.startswith("^"):
                # Nobody here cares. This was just a deletion.
                pass
            else:
                for i in xrange(len(part)):
                    # All the bases we have letters for here were replaced. Put
                    # in the bases that were replaced.
                    to_return[cursor + i] = part[i]
                
                cursor += len(part)
            
    
    print md_string
    print parts
    print to_return
    
    return to_return
    

def pysam_parse_reads(sam_file):
    """
    Parse and deduplicate with pysam. Takes a filename or File object.
    
    Yields Reads with their is_suspect flag modified so we only suspect alt
    alignments which have central mismatches.
    
    """
    
    sam = pysam.AlignmentFile(sam_file)
    
    for read in sam:
    
        print("Handling read on template {}".format(read.query_name))
    
        is_suspect = False
        
        # Calculate the length of the original input read, including bases that
        # were hard clipped.
        input_length = 0
        
        if read.cigartuples is None:
            # Just use the read itself.
            # TODO: why do we have no CIGAR?
            input_length = len(read.query_sequence)
        else:
        
            for (op, count) in read.cigartuples:
                if op == 2 or op == 3 or op == 6:
                    # Skip deletions, reference skips, and padding
                    continue
                else:
                    input_length += count
        
        reference_name = sam.getrname(read.reference_id)
        
        if input_length is None:
            print read
        
        if input_length % 2 == 0 and reference_name.endswith("_alt"):
            # It could be corrupted.

            # Work out the clipping offset from the CIGAR (4 = soft clip, 5 =
            # hard clip)
            first_cigar = read.cigartuples[0]
            if (first_cigar[0] == 4 or
                first_cigar[0] == 5):
                offset = first_cigar[1]
            else:
                offset = 0

            # Work out what's up with every query base from the MD tag
            status_per_base = parse_MD_tag(read.get_tag("MD"), offset,
                input_length)

            
                
            print("Errors: {}: {} {}: {}".format(input_length / 2, 
                status_per_base.has_key(input_length / 2),
                input_length / 2 + 1, 
                status_per_base.has_key(input_length / 2 + 1)))
                
            # Use truncating division to check around the center
            if (status_per_base.has_key(input_length / 2) or
                status_per_base.has_key(input_length / 2 + 1)):
                
                # There's a mismatch near the center
                is_suspect = True
        
        # Re-parse with our code after fixing up the reference name.
        # TODO ugly hack in need of redesign
        parts = str(read).split("\t")
        parts[2] = reference_name
        our_read = Read("\t".join(parts))
        
        # Fix up the suspect flag
        our_read.is_suspect = is_suspect
        
        yield our_read
    
        
def parse_and_deduplicate_sam(sam_input):
    """
    Given a source of input SAM lines, parses lines into Read objects, and
    deduplicates them, discarding suspect ones when non-suspect ones are
    available.
    
    Yields dicts form end number to Read object for each template.
    
    """
    
    # What was the template for the last read
    last_template = None
    
    # For this template, we keep the best read for each end we find.
    reads_by_end = {}
    
    for line in sam_input:
        if line.startswith("@"):
            continue
            
        read = Read(line)
        
        # Work on the reads
        
        if read.template != last_template:
            
            # It's OK if we spit out suspect stuff as long as we got the best
            # suspect alignment
            yield reads_by_end
        
            # Start a new state
            last_template = read.template
            reads_by_end = {}
            
        
        if not reads_by_end.has_key(read.end):
            # This is the only read for this end so far
            reads_by_end[read.end] = read
        else:
            if reads_by_end[read.end] < read:
                # Replace the existing read
                reads_by_end[read.end] = read
            elif (not read.is_suspect and 
                len(read.sequence) >= len(reads_by_end[read.end].sequence) and 
                read != reads_by_end[read.end]):
                # We aren't suspect, we differ, and we can't replace the other
                # read.
                raise RuntimeError("Non-suspect alignments don't agree on end "
                    "{} of template {}:\n{}\n{}".format(read.end, read.template,
                    read.line, reads_by_end[read.end].line))
    
    # It's OK if we spit out suspect stuff as long as we got the best
    # suspect alignment
    yield reads_by_end
            
def write_fastq(stream, read):
    """
    Write the given record as FASTQ to the given stream
    """
    
    # Unpack and format the record
    stream.write("@{}\n{}\n+\n{}\n".format(read.get_name(), read.sequence,
        read.qualities))
    
def main(args):
    """
    Parses command line arguments and do the work of the program.
    "args" specifies the program arguments, with args[0] being the executable
    name. The return value should be used as the program's exit code.
    """
    
    if len(args) == 2 and args[1] == "--test":
        # Run the tests
        return doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
    
    options = parse_args(args) # This holds the nicely-parsed options object
    
    for reads_by_end in parse_and_deduplicate_sam(options.input_sam):
        if not (reads_by_end.has_key(1) and reads_by_end.has_key(2)):
            # Skip unpaired reads
            continue
            
        # Split up the reads to their files
        write_fastq(options.fq1, reads_by_end[1])
        
        if options.interleaved:
            # Both go to the same file
            write_fastq(options.fq1, reads_by_end[2])
        else:
            write_fastq(options.fq2, reads_by_end[2])
        
    
    
    
if __name__ == "__main__" :
    sys.exit(main(sys.argv))
        
        
        
        
        
        
        
        
        
        

