#!/usr/bin/env python2.7
"""
Computes Mendelian concordance between trio of sample graphs, as output to vcf. 

"""
# todo: this was written before vcf conversion, then modified to run on vcf.
#       should look into using off-the-shelf tool on vcfs...

import argparse, sys, os, os.path, random, subprocess, shutil, itertools, glob
import doctest, re, json, collections, time, timeit, string
from collections import defaultdict
from vcfCompare import parse_alts, parse_ref

def parse_args(args):
    parser = argparse.ArgumentParser(description=__doc__, 
        formatter_class=argparse.RawDescriptionHelpFormatter)
        
    # General options
    parser.add_argument("child", type=str,
                        help="child calls (vcf)")    
    parser.add_argument("parent1", type=str,
                        help="parent1 calls (vcf)")
    parser.add_argument("parent2", type=str,
                        help="parent2 calls (vcf)")
    parser.add_argument("-c", action="store_true", default=False,
                        help="Ignore sequence name (only look at start)")
    parser.add_argument("-i", action="append", default=[],
                        help="Ignore lines contaning keyword")
    parser.add_argument("-g", action="store_true", default=False,
                        help="Ignore genotype (which isnt reliable for graph calls) assume all alts 0/1")
    
    
    
    
    args = args[1:]
        
    return parser.parse_args(args)


def parse_alleles(toks, options):
    """ return the alleles (first sample /last column) todo: more general?"""
    if toks[-2].split(":")[0] == "GT": 
        gtcol = len(toks) - 1
        if options.g:
            gttok = "0|1"
        else:
            gttok = toks[gtcol].split(":")[0]
        gts = "|".join(gttok.replace(".", "0").split("/")).split("|")
        vals = [parse_ref(toks)] + parse_alts(toks)
        alleles = [vals[int(x)] for x in gts]
        return alleles
    return []

def make_allele_dict(vcf_path, options):
    """ load up all variants by their coordinates
    map (chrom, pos) -> [allele1, allele2]
    """
    vcf_dict = dict()
    ref_dict = dict()
    with open(vcf_path) as f:
        for line in f:
            skip = line[0] == "#"
            for ignore_keyword in options.i:
                if ignore_keyword in line:
                    skip = True
            if not skip:
                toks = line.split()
                chrom = toks[0] if options.c is False else None
                pos = int(toks[1])
                alleles = parse_alleles(toks, options)
                vcf_dict[(chrom, pos)] = alleles
                ref_dict[(chrom, pos)] = [parse_ref(toks)]
    return vcf_dict, ref_dict

def score_call(child_alleles, parent1_alleles, parent2_alleles, options):
    """ basic consistency.  score 1 out of 1 if child alleles can come from combination of
parents, 0 / 1 otherwise 
    """
    if len(child_alleles) == 2:
        # only handle diploid calls
        if child_alleles[0] in parent1_alleles and child_alleles[1] in parent2_alleles or\
           child_alleles[0] in parent2_alleles and child_alleles[1] in parent1_alleles:
            return 1, 1
        else:
            return 0, 1
    else:

        return 0, 0    
    
def main(args):

    options = parse_args(args)

    parent1_alleles, parent1_refs = make_allele_dict(options.parent1, options)
    parent2_alleles, parent2_refs = make_allele_dict(options.parent2, options)
    child_alleles, child_refs = make_allele_dict(options.child, options)

    # keep track of all reference positions in one place (rather than
    # bother with fasta)
    ref_alleles = parent1_refs
    ref_alleles.update(parent2_refs)
    ref_alleles.update(child_refs)

    score = 0, 0
    for pos, alleles in child_alleles.items():
        p1 = parent1_alleles[pos] if pos in parent1_alleles else ref_alleles[pos]
        p2 = parent2_alleles[pos] if pos in parent2_alleles else ref_alleles[pos]
        a = score_call(alleles, p1, p2, options)
        #sys.stderr.write("{} {} -> {}\n".format(pos, alleles, a))
        score = score[0] + a[0], score[1] + a[1]
        
    print "{}\t{}\t{}".format(score[0], score[1] - score[0],
                              float(score[0]) / max(1, (score[1])))
    return 0
    
if __name__ == "__main__" :
    sys.exit(main(sys.argv))
