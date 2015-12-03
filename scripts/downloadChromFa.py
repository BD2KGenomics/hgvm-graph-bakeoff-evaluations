#!/usr/bin/env python2.7
"""
Baseline comparison of calliing needs 1000 genomes sample graphs
which are cosntructed from 1000 genomes vcf data.  To build these,
we also need the entire chromosome fasta sequence for each region.
This script downloads these from ucsc into one file that can be 
used by vg construct.

mostly cribbed from:
https://github.com/glennhickey/vg2sg/blob/master/data/fetch1kgpRegion.py

"""

import argparse, sys, os, os.path, random, subprocess, shutil, itertools, glob

def parse_args(args):
    parser = argparse.ArgumentParser(description=__doc__, 
        formatter_class=argparse.RawDescriptionHelpFormatter)
        
    # General options
    parser.add_argument("--out_fa", type=str, default="data/g1kvcf/chrom.fa",
                        help="output directory")
    parser.add_argument("--chroms", type=str, default="5,6,13,17,19",
                        help="comma-separated list of chroms to download")
    args = args[1:]
        
    return parser.parse_args(args)

def downloadChrom(chrom, options, assembly = "hg38"):
    """ concatenate whole chromosome onto output file
    """
    # This is the URL to get it from
    fa_url = ("http://hgdownload.cse.ucsc.edu/goldenPath/{}/"
              "chromosomes/chr{}.fa.gz").format(assembly, chrom)

    # download and uncompress 
    cmd = "curl {} | zcat".format(fa_url)
    # take chr out of sequence names
    cmd += " | sed \"s/chr{}/{}/\" ".format(chrom, chrom)
    # convert to upper case
    cmd += " | sed \"s/a/A/g;s/c/C/g;s/g/G/g;s/t/T/g;s/n/N/g\" "
    # tack onto output
    cmd += " >> {}".format(options.out_fa)
    os.system(cmd)
    

def main(args):
    
    options = parse_args(args)

    # make sure output directory exists and file doesn't
    if not os.path.isdir(os.path.dirname(options.out_fa)):
        os.makedirs(os.path.dirname(options.out_fa))
    os.system("rm -f {}".format(options.out_fa))

    # append each chrom to file
    for chrom in options.chroms.split(","):
        downloadChrom(chrom, options)

    #make index
    os.system("samtools faidx {}".format(options.out_fa))
    
        
if __name__ == "__main__" :
    sys.exit(main(sys.argv))
