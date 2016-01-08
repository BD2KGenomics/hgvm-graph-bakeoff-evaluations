#!/usr/bin/env python2.7

"""
Tool to take a region of a VCF, cut it out, and re-assign all its coordinates to
a new contig name and starting offset.

Only rewrites schom and start columns; doesn't deal with any references to
positions that may occur in INFO or FORMAT values.

Variants only partially overlapping the specified range will be dropped.

"""

import argparse
import sys
import os
import doctest

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
    
    parser.add_argument("--vcf_in", default=sys.stdin,
        type=argparse.FileType("r"),
        help="VCF file to read")
    parser.add_argument("--vcf_out", default=sys.stdout,
        type=argparse.FileType("w"),
        help="VCF file to write")
        
    parser.add_argument("--source_contig", required=True,
        help="contig name to extract variants from")
    parser.add_argument("--source_start", default=1, type=int,
        help="first base (inclusive) at which to start collecting variants")
    parser.add_argument("--source_end", default=float("+inf"), type=int,
        help="last base (exclusive) at which to stop collecting variants")
    
    parser.add_argument("--dest_contig", default="ref",
        help="contig name to place variants on")
    parser.add_argument("--dest_start", default=1, type=int,
        help="base on the destination contig corresponding to --source_start")
    
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
    
    for line in options.vcf_in:
        # Loop through the lines in the input VCF. We need to copy the headers,
        # and possibly copy and rewrite the records.
        
        # Drop the newline
        line = line.rstrip("\n")
        
        if len(line) == 0:
            # Skip blank lines
            continue
            
        if line[0] == "#":
            # It's a header. Keep it
            options.vcf_out.write("{}\n".format(line))
            continue
            
        # Otherwise it's a record
        
        # Split on tabs
        parts = line.split("\t")
        
        if parts[0] != options.source_contig:
            # It's not on the right contig
            continue
            
        # Parse where the variant starts
        variant_start = int(parts[1])
            
        if (variant_start < options.source_start or
            variant_start >= options.source_end):
            # It's out of range
            return
          
        # Rewrite position and contig  
        variant_start = (variant_start - options.source_start +
            options.dest_start)
        parts[1] = str(variant_start)
        parts[0] = options.dest_contig
        
        # Spit out the fixed variant
        options.vcf_out.write("{}\n".format("\t".join(parts)))
    
if __name__ == "__main__" :
    sys.exit(main(sys.argv))
    
    
    
