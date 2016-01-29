#!/usr/bin/env python2.7

"""
 Filter out records that don't meet a particular quality threshold.
Quality assumed to be 6th column
"""


import argparse, sys, os, os.path, random, subprocess, shutil, itertools

def parse_args(args):
    parser = argparse.ArgumentParser(description=__doc__, 
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("in_vcf", type=str,
                        help="Input vcf file (- for stdin)")
    parser.add_argument("min_qual", type=float,
                        help="Mininum quality value to keep")
    parser.add_argument("--pct", action="store_true",
                        help="Interpret min_qual as a precentile.  So min_qual == 0.7 will return the top 30% scores")  
    
    args = args[1:]
    options = parser.parse_args(args)
    return options

def get_qual_from_line(line, options):
    return float(line.split()[5])

def compute_cutoff(options):
    if options.pct is False:
        return options.min_qual
    else:
        assert options.min_qual >= 0. and options.min_qual <= 1.
        if options.in_vcf == "-":
            vcf_file = sys.stdin
        else:
            vcf_file = open(options.in_vcf)
        quals = []
        for line in vcf_file:
            if line[0] != "#":
                quals.append(get_qual_from_line(line, options))
        if options.in_vcf != "-":
            vcf_file.close()
        return sorted(quals)[int(options.min_qual * len(quals))]
        

def main(args):
    options = parse_args(args)

    cutoff = compute_cutoff(options)

    if options.in_vcf == "-":
        vcf_file = sys.stdin
    else:
        vcf_file = open(options.in_vcf)

    for line in vcf_file:
        if line[0] == "#" or get_qual_from_line(line, options) >= cutoff:
            sys.stdout.write(line)

    if options.in_vcf != "-":
        vcf_file.close()
	 
if __name__ == "__main__" :
    sys.exit(main(sys.argv))
