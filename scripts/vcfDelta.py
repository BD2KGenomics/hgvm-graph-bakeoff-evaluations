#!/usr/bin/env python2.7

"""
5-minute vcf compare script to help debug calls. return vc1 - vcf2
"""


import argparse, sys, os, os.path, random, subprocess, shutil, itertools

def parse_args(args):
    parser = argparse.ArgumentParser(description=__doc__, 
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("vcf1", type=str,
                        help="Input vcf file 1"),
    parser.add_argument("vcf2", type=str,
                        help="Input vcf file 2")
    parser.add_argument("-i", action="store_true", default=False,
                        help="Do vcf2 - vcf1 instead of vcf1 - vcf2")
    parser.add_argument("-c", action="store_true", default=False,
                        help="Ignore sequence name (only look at start)")
    parser.add_argument("-g", action="store_true", default=False,
                        help="Take into account GT (10th col)")
    parser.add_argument("-a", action="store_true", default=False,
                        help="Do intersection of vcf1 AND vcf2")
    
    
    args = args[1:]
    options = parser.parse_args(args)
    return options
    
def main(args):
    options = parse_args(args)

    if options.i is True:
        options.vcf1, options.vcf2 = options.vcf2, options.vcf1

    # load all the coordinates in vcf2 (dont mind alleles or gt)
    coord_set = set()
        
    vcf2 = open(options.vcf2, "r")
    
    while True:
        line = vcf2.readline()
        if line == "":
            break
        elif line[0] != "#":
            toks = line.split()
            assert len(toks) > 3
            if options.c is True:
                toks[0] = None
            gt = None
            if options.g is True and len(toks) >= 10:
                gt = toks[9].replace("|","/")
                if gt == "1/0":
                    gt = "0/1"
            coord_set.add((toks[0], int(toks[1]), gt))
    vcf2.close()

    # read vcf 1, output stuff not found above
    vcf1 = open(options.vcf1, "r")
    
    while True:
        line = vcf1.readline()
        if line == "":
            break
        elif line[0] == "#":
            sys.stdout.write(line)
        else:
            toks = line.split()
            assert len(toks) > 3
            if options.c is True:
                toks[0] = None
            gt = None
            if options.g is True and len(toks) >= 10:
                gt = toks[9].replace("|","/")                
                if gt == "1/0":
                    gt = "0/1"

            found = (toks[0], int(toks[1]), gt) in coord_set
            if (options.a and found) or (not options.a and not found):
                sys.stdout.write(line)
    vcf1.close()

if __name__ == "__main__" :
    sys.exit(main(sys.argv))
