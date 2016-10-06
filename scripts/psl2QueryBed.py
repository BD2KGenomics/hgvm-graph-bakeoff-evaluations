#!/usr/bin/env python2.7
"""
psl2QueryBed.py: take a PSL and convert it into a BED in the space of the query

Allows overriding contig names and providing offsets.

"""

import argparse, sys, os, os.path, random, itertools, string, re
import doctest

import Bio.SearchIO
import tsv

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
    parser.add_argument("--input_psl", type=argparse.FileType("r"),
        default=sys.stdin,
        help="PSL file to read as input")
    parser.add_argument("--output_bed", type=argparse.FileType("w"),
        default=sys.stdout,
        help="BED file to write as output")
    parser.add_argument("--contig", default=None,
        help="replace all contig names with this name")
    parser.add_argument("--offset", type=int, default=0,
        help="add this offset to all coordinates")
    
    # The command line arguments start with the program name, which we don't
    # want to treat as an argument for argparse. So we remove it.
    args = args[1:]
        
    return parser.parse_args(args)


    
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
    
    # Actually do the work. We structure it like this so we can use it as a
    # script or a module.
    run(options)
    
def run(options):
    """
    Do the actual work of the program.
    """
    
    # Set up our BED output
    writer = tsv.TsvWriter(options.output_bed)
    
    # Read the PSL
    for result in Bio.SearchIO.parse(options.input_psl, "blat-psl"):
        for hit in result:
            for hsp in hit:
                # A hit is the equivalent of a PSL line; we're going to make it a BED record
                
                # Pull out the block sizes and starts
                block_sizes = list(hsp.hit_span_all)
                block_starts = list(hsp.query_start_all)
                
                # We need to find out the strand
                strand = "+"
                for fragment in hsp:
                    # Loop over all the fragments. We know they'll all be on the
                    # same strand, because that's how PSL can articulate them.
                    if fragment.query_strand == -1:
                        # We ought to be on the - strand. This means our blocks
                        # are going to be in backwards order.
                        strand = "-"
                        block_sizes.reverse()
                        block_starts.reverse()
                        break
                        
                # Convert starts to be relative to aligned region
                block_starts = [x - hsp.query_start for x in block_starts]
                
                # The first should always be 0
                assert(block_starts[0] == 0)
                
                # BED is: chrom, chromStart, chromEnd, name, score, strand,
                # thickStart, thickEnd, itemRgb, blockCount, blockSizes,
                # blockStarts
                bed_record = [
                    # chrom
                    options.contig or result.id,
                    # chromStart
                    hsp.query_start + options.offset,
                    # chromEnd
                    hsp.query_end + options.offset,
                    # name
                    hit.id,
                    # score
                    0,
                    # strand
                    strand,
                    # thickStart
                    hsp.query_start + options.offset,
                    # thickEnd
                    hsp.query_start + options.offset,
                    # itemRgb
                    "0,0,0", 
                    # blockCount
                    len(hsp),
                    # blockSizes
                    ",".join((str(x) for x in block_sizes)),
                    # blockStarts
                    ",".join((str(x) for x in block_starts))
                ]
                
                writer.list_line(bed_record)
        
    
    
    
if __name__ == "__main__" :
    sys.exit(main(sys.argv))
        
        
        
        
        
        
        
        
        
        

