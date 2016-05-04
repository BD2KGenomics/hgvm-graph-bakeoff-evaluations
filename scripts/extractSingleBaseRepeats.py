#!/usr/bin/env python

#Copyright (C) 2014 by Glenn Hickey
#
#Released under the MIT license, see LICENSE.txt
import sys
import os
import argparse
import copy
import array

"""
Extract runs of a single nucleotide out of a FASTA file and into a BED file
"""
# Copied from bioio.py from sonLib (https://github.com/benedictpaten/sonLib):
# Copyright (C) 2006-2012 by Benedict Paten (benedictpaten@gmail.com)
# Released under the MIT license, see LICENSE.txt
def fastaRead(fileHandle):
    """iteratively a sequence for each '>' it encounters, ignores '#' lines
    """
    line = fileHandle.readline()
    while line != '':
        if line[0] == '>':
            name = line[1:-1]
            line = fileHandle.readline()
            seq = array.array('c')
            while line != '' and line[0] != '>':
                if line[0] != '#':
                    seq.extend([ i for i in line[:-1] if i != '\t' and i != ' ' ])
                line = fileHandle.readline()
            for i in seq:
                #For safety and sanity I only allows roman alphabet characters in fasta sequences. 
                if not ((i >= 'A' and i <= 'Z') or (i >= 'a' and i <= 'z') or i == '-'):
                    raise RuntimeError("Invalid FASTA character, ASCII code = \'%d\', found in input sequence %s" % (ord(i), name))
            yield name, seq.tostring()
        else:
            line = fileHandle.readline()

            
def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Extract runs of a single nucleotide from a "
        "FASTA file into a BED file.")
    parser.add_argument("inputFa", help="Input FASTA file")
    parser.add_argument("minLength", help="Minimum interval length", type=int)
    parser.add_argument("--caseSensitive", help="Case sensitive comparison",
                        action="store_true", default=False)
    parser.add_argument("--nuc", help="Commas-separated list of nucleotides "
                        " to consider",
                        default="N")
    
    args = parser.parse_args()
    assert os.path.isfile(args.inputFa)
    faFile = open(args.inputFa, "r")
    nucs = set()
    for nuc in args.nuc.split(","):
        assert len(nuc) == 1
        if args.caseSensitive is False:
            nucs.add(nuc.upper())
        else:
            nucs.add(nuc)
        
    for seqName, seqString in fastaRead(faFile):
        curInterval = [None, -2, -2, None]
        for i in xrange(len(seqString)):
            
            if seqString[i] in nucs or (args.caseSensitive is False and
                                        seqString[i].upper() in nucs):
                # print then re-init curInterval
                if i != curInterval[2]:
                    if curInterval[0] is not None:
                        printInterval(sys.stdout, curInterval, args)
                    curInterval = [seqName, i, i + 1, seqString[i]]
                # extend curInterval
                else:
                    assert seqName == curInterval[0]
                    assert i > curInterval[1]
                    curInterval[2] = i + 1

            # print last interval if exists
            if i == len(seqString) - 1 and curInterval[0] is not None:
                printInterval(sys.stdout, curInterval, args)
                
    faFile.close()


def printInterval(ofile, curInterval, args):
    l = curInterval[2] - curInterval[1]
    if l >= args.minLength:
        ofile.write("%s\t%d\t%d\t%s\t%d\n" % (
            curInterval[0],
            curInterval[1],
            curInterval[2],
            curInterval[3],
            l))
    
if __name__ == "__main__":
    sys.exit(main())
