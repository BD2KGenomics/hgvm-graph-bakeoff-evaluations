#!/usr/bin/env python2.7
"""
fasta2reads.py: Turn a FASTA into chunks of at most the given number of bases,
one per line.

These will be compatible with `vg map`'s "read file" input option.

"""

import argparse, sys, os, os.path, random, subprocess, shutil, itertools
import collections, urllib2, shutil, subprocess, glob, doctest
import math

from Bio import AlignIO, SeqIO, Align, Entrez
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

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
    parser.add_argument("--fasta_file", type=argparse.FileType('r'),
        default=sys.stdin,
        help="FASTA to parse")
    parser.add_argument("--reads_file", type=argparse.FileType('w'),
        default=sys.stdout,
        help="Reads file to write")
    parser.add_argument("--max_chunk_size", type=int, default=1000000,
        help="maximum chunk size to generate")
    parser.add_argument("--uppercase", action="store_true",
        help="make sequence upper-case")

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
    
    for record in SeqIO.parse(options.fasta_file, "fasta"):
        # For every FASTA record...
        
        # Get its total length
        length = len(record.seq)
        
        # Determine number of chunks needed
        chunks = int(math.ceil(float(length) / options.max_chunk_size))
        
        # Split into that many equal chunks (so we have no tiny chunks if we can
        # help it).
        chunk_size = length / chunks
        
        for i in xrange(chunks):
            # Grab out each chunk
            chunk = record.seq[chunk_size * i : chunk_size * (i + 1)]
            
            if options.uppercase:
                # Fix up the case
                chunk = chunk.upper()
            
            # Spit it out
            options.reads_file.write(str(chunk))
            options.reads_file.write("\n")
        
                
        
            
if __name__ == "__main__" :
    sys.exit(main(sys.argv))
