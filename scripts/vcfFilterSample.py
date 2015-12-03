#!/usr/bin/env python2.7

"""
 Filter a sample out of a vcf, keeping only snps that are present in 
 genotype information for that sample.  

 https://github.com/ekg/vcflib must be installed (with vcflib/bin) in PATH. 
"""


import argparse, sys, os, os.path, random, subprocess, shutil, itertools

def parse_args(args):
    parser = argparse.ArgumentParser(description=__doc__, 
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("vcf", type=str,
                        help="Input vcf file"),
    parser.add_argument("sample", type=str,
                        help="Sample ID to extract")
    parser.add_argument("--dot", action="store_true", default=False,
                        help="Accept GT ., GT ./. and GT 0|. genotypes")    
    args = args[1:]
    options = parser.parse_args(args)
    return options

def main(args):
    options = parse_args(args)

    filterCmd = "vcfkeepsamples {} {}".format(options.vcf, options.sample)

    filterProc = subprocess.Popen(filterCmd, shell=True, stdout=subprocess.PIPE)

    first = True
    while True:
        line = filterProc.stdout.readline()
        if line != "":
	    if line.find("GT\t0|0\n") < 0 and (
                # not too sure what to make of these genotypes,
                # so toggle with "dot" option for now. 
	        options.dot or (line.find("GT\t.\n") < 0 and
				line.find("GT\t.|.\n") < 0 and
                                line.find("GT\t.|0\n") < 0 and
                                line.find("GT\t0|.\n") < 0)):
                sys.stdout.write(line)
            elif first is True and len(line.lstrip()) > 0 and line.lstrip()[0] != "#":
                # vg construct doesn't like empty vcf.  make sure that there's always
                # a trivial reference position if nothing else
                toks = line.split("\t")
                toks[4] = toks[3]
                sys.stdout.write("\t".join(toks))
                first = False
        else:
            break
	 
if __name__ == "__main__" :
    sys.exit(main(sys.argv))
