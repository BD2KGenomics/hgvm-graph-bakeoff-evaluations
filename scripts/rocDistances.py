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


def smooth_table(linetoks, threshold = 0.001):
    """ precisions dont seem to always be roc-like.  remove outliers to smooth
    into curves until I figure out what's going on
    """

    # compute biggest "dip" in entire table.  remove it. then repeat.
    # (n^2 is slower than necessary, but good enough for our tables)
    spike_idx = -1
    spike_val = 0
    for i in range(1, len(linetoks)):
        if i == 0:
            continue
        # we expect precision to increase and recall to decrease as i increases
        # spike measures deviation of this from i-1 to i and i to i+1
        spike = 0.
        # same graph as previous
        if linetoks[i-1][0] == linetoks[i][0]:
            # compute how much precision we *lose*
            spike += max(0., float(linetoks[i-1][1]) - float(linetoks[i][1]))
            # compute how much recall we *gain*
            spike += max(0., float(linetoks[i][2]) - float(linetoks[i-1][2]))

        # same graph as next
        if i < len(linetoks) - 1 and linetoks[i][0] == linetoks[i+1][0]:
            # compute how much precision we *lose*
            spike += max(0., float(linetoks[i][1]) - float(linetoks[i+1][1]))
            # compute how much recall we *gain*
            spike += max(0., float(linetoks[i+1][2]) - float(linetoks[i][2]))

        if spike > spike_val:
            spike_val, spike_idx = spike, i

    if spike_val > threshold:
        return smooth_table(linetoks[0:spike_idx] + linetoks[spike_idx+1:])
    else:
        return ["\t".join(x) + "\n" for x in linetoks]
    
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

    # let's sort the output to make it easier to remove dead points
    for tsv in glob.glob(os.path.join(options.out_dir, "comp_tables", "*.tsv")):
        with open(tsv) as f:
            lines = [line for line in f]
            lines = [lines[0]] + sorted(lines[1:], key = lambda x : (x.split()[0], float(x.split()[1]), 1 - float(x.split()[2])))
            # precisions can be bumpy (need to change to sensitivy?)
            # use simple smoother in the meantime
            lines = smooth_table([x.split() for x in lines])
        with open(tsv, "w") as f:
            for line in lines:
                f.write(line)

if __name__ == "__main__" :
    sys.exit(main(sys.argv))
        
        

        
