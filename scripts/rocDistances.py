#!/usr/bin/env python2.7
"""
Stitch some outputs of computeVariatnsDistances.py into ROC plots.  These can 
be plotted with plotVariantsDistances.py
"""

import argparse, sys, os, os.path, random, subprocess, shutil, itertools, glob
import doctest, re, json, collections, time, timeit, string, math, copy
from collections import defaultdict
from Bio.Phylo.TreeConstruction import _DistanceMatrix, DistanceTreeConstructor
from Bio import Phylo
import matplotlib
matplotlib.use('Agg')
import pylab
import networkx as nx
from collections import defaultdict
from toil.job import Job
from toillib import RealTimeLogger, robust_makedirs
from callVariants import alignment_sample_tag, alignment_region_tag, alignment_graph_tag, run
from callVariants import graph_path, sample_vg_path, g1k_vg_path, graph_path
from callStats import vg_length
from evaluateVariantCalls import defaultdict_set
from computeVariantsDistances import vcf_dist_header



def parse_args(args):
    parser = argparse.ArgumentParser(description=__doc__, 
        formatter_class=argparse.RawDescriptionHelpFormatter)
    
    # General options
    parser.add_argument("comp_dirs", nargs="+",
                        help="directories of comparison output written by computeVariantsDistances.py")
    parser.add_argument("out_dir",
                        help="output directory where /comp_tables will be written to.")
    parser.add_argument("--skip_first", type=int, default=0,
                        help="skip (alphabetically) first N input dirs")
    parser.add_argument("--skip_last", type=int, default=0,
                        help="skip (alphabetically) last N input dirs")
                            
    args = args[1:]

    return parser.parse_args(args)

    
def main(args):
    
    options = parse_args(args)

    robust_makedirs(os.path.join(options.out_dir, "comp_tables"))

    first = True
    # sort the directories, assuming their names give info on their order
    # in the roc
    for comp_dir in sorted(options.comp_dirs)[options.skip_first:len(options.comp_dirs) - options.skip_last]:
        # this can happen easily using wildcards in input
        if comp_dir == options.out_dir:
            continue
        # look through tsvs in comp_tables.
        for tsv in glob.glob(os.path.join(comp_dir, "comp_tables", "*.tsv")):
            # overwrite if first
            # strip header and append if second
            c = "cp {} ".format(tsv) if first is True else "tail -n +2 {} >> ".format(tsv)
            # just cat into the output directory
            os.system("{} {}".format(c, os.path.join(options.out_dir, "comp_tables",
                                                     os.path.basename(tsv))))
        first = False
    
if __name__ == "__main__" :
    sys.exit(main(sys.argv))
        
        

        
