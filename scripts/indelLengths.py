#!/usr/bin/env python2.7

"""
Script to extract the lengths of indels from a VCF file.

Given a VCF file as input on standard in, produces, for each variant, an indel
length on standard out. Indel lengths are absolute, and represent the largest
number of bases that can be gained or lost by changing from the shortest allele
of the variant to the longest. Lengths are emitted one per line.

Symbolic alleles and breakends are not supported. Symbolic alleles are not
detected, and will be interpreted as sequence rather than allele names.

"""


import argparse
import sys
import os
import doctest

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
    
    
    parser.add_argument("--in_file", type=argparse.FileType("r"), default=sys.stdin,
        help="VCF file to read")
    parser.add_argument("--out_file", type=argparse.FileType("w"),
        default=sys.stdout,
        help="file of lengths to write")
    parser.add_argument("--indels_only", action="store_true",
        help="ignore variants with no length change")
    
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
    
    # Read the VCF in as a TSV
    reader = tsv.TsvReader(options.in_file)
    
    # Open a writer to spit out length differences (as a 1-column TSV)
    writer = tsv.TsvWriter(options.out_file)
    
    for parts in reader:
        # For every VCF line
        
        if len(parts) < 4:
            # Skip things that aren't valid records.
            continue
    
        # Ref allele is field 3
        ref_allele = parts[3]
        # Alt alleles are field 4, separated by commas
        alt_alleles = parts[4].split(",")
        
        # Squish them all together
        all_alleles = alt_alleles + [ref_allele]
        
        # Turn them into a sorted list of deduplicated lengths
        allele_lengths = sorted({len(allele) for allele in all_alleles})
        
        # We know there must be at least one allele, so we can subtract the
        # shortest length from the longest length.
        length_difference = allele_lengths[-1] - allele_lengths[0]
        
        if not options.indels_only or length_difference > 0:
            # Emit the length difference
            writer.line(length_difference)
    
    
    
if __name__ == "__main__" :
    sys.exit(main(sys.argv))
    
    
    
