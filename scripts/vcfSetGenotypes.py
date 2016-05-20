#!/usr/bin/env python2.7

"""
 Set the genotype of every record.  Assume it's just the last column (ie only one sample)
"""


import argparse, sys, os, os.path, random, subprocess, shutil, itertools

def parse_args(args):
    parser = argparse.ArgumentParser(description=__doc__, 
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("in_vcf", type=str,
                        help="Input vcf file (- for stdin)")
    parser.add_argument("--gt", type=str, default="0/1",
                        help="Genotype to apply to each record")
                        
    
    args = args[1:]
    options = parser.parse_args(args)
    return options
        

def main(args):
    options = parse_args(args)

    if options.in_vcf == "-":
        vcf_file = sys.stdin
    else:
        vcf_file = open(options.in_vcf)

    for line in vcf_file:
        if line[0] == "#":
            sys.stdout.write(line)
        else:
            toks = line.split("\t")
            toks[-1] = options.gt
            sys.stdout.write("\t".join(toks) + "\n")

    if options.in_vcf != "-":
        vcf_file.close()
	 
if __name__ == "__main__" :
    sys.exit(main(sys.argv))
