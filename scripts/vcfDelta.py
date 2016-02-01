#!/usr/bin/env python2.7

"""
5-minute vcf compare script to help debug calls. return vc1 - vcf2
"""


import argparse, sys, os, os.path, random, subprocess, shutil, itertools

from vcfCompare import make_vcf_dict, find_alt

def parse_args(args):
    parser = argparse.ArgumentParser(description=__doc__, 
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("vcf1", type=str,
                        help="Input vcf file 1"),
    parser.add_argument("vcf2", type=str,
                        help="Input vcf file 2")
    parser.add_argument("-r", action="store_true", default=False,
                        help="Do vcf2 - vcf1 instead of vcf1 - vcf2")
    parser.add_argument("-c", action="store_true", default=False,
                        help="Ignore sequence name (only look at start)")
    parser.add_argument("-g", action="store_true", default=False,
                        help="Take into account GT (10th col)")
    parser.add_argument("-a", action="store_true", default=False,
                        help="Do intersection of vcf1 AND vcf2")
    parser.add_argument("-s", action="store_true", default=False,
                        help="Only snps")
    parser.add_argument("-i", action="append", default=[],
                        help="Ignore lines contaning keyword")    
    
    
    args = args[1:]
    options = parser.parse_args(args)
    return options
    
def main(args):
    options = parse_args(args)

    if options.r is True:
        options.vcf1, options.vcf2 = options.vcf2, options.vcf1

    # load vcf2 into dictionary lookup
    vcf_dict2 = make_vcf_dict(options.vcf2, options)

    # scan vcf1, checking dictionary to find elements not in vcf2
    with open(options.vcf1) as vcf1:
        for line in vcf1:
            if line[0] != "#":
                toks = line.split()
                chrom, pos, ref, alts = toks[0], int(toks[1]), toks[3], toks[4].split(",")
                for alt in alts:
                    if not options.s or (len(ref) == 1 and len(alt) == 1):
                        if not find_alt(chrom, pos, ref, alt, vcf_dict2):
                            sys.stdout.write("\t".join(toks[0:4]) + "\t{}\t".format(alt) + "\t".join(toks[5:]) + "\n")
    
if __name__ == "__main__" :
    sys.exit(main(sys.argv))
