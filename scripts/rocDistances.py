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
    parser.add_argument("--pcol", type=int, default=1,
                        help="column to look for for precision")
    parser.add_argument("--rcol", type=int, default=2,
                        help="column to look for for recall")
    parser.add_argument("--smooth", action="store_true",
                        help="apply greedy smoothing to remove outliers")
                            
    args = args[1:]

    return parser.parse_args(args)


def smooth_table(linetoks, options, threshold = 0.001):
    """ precisions dont seem to always be roc-like.  remove outliers to smooth
    into curves until I figure out what's going on
    """

    # compute biggest "dip" in entire table.  remove it. then repeat.
    # (n^2 is slower than necessary, but good enough for our tables)
    keep_going = True
    linetoks = copy.deepcopy(linetoks)
    while keep_going:
        spike_idx = -1
        spike_val = 0
        for i in range(1, len(linetoks)):
            # we expect precision to increase and recall to decrease as i increases
            # spike measures deviation of this from i-1 to i and i to i+1
            spike = 0.
            # same graph as previous
            if linetoks[i-1][0] == linetoks[i][0]:
                # compute how much precision we *lose*
                spike += max(0., float(linetoks[i-1][options.pcol]) - float(linetoks[i][options.pcol]))
                # compute how much recall we *gain*
                spike += max(0., float(linetoks[i][options.rcol]) - float(linetoks[i-1][options.rcol]))

            # same graph as next
            if i < len(linetoks) - 1 and linetoks[i][0] == linetoks[i+1][0]:
                # compute how much precision we *lose*
                spike += max(0., float(linetoks[i][options.pcol]) - float(linetoks[i+1][options.pcol]))
                # compute how much recall we *gain*
                spike += max(0., float(linetoks[i+1][options.rcol]) - float(linetoks[i][options.rcol]))

            if spike > spike_val:
                spike_val, spike_idx = spike, i

        if spike_val > threshold:
            del linetoks[spike_idx]
        else:
            keep_going = False
            
    return ["\t".join(x) + "\n" for x in linetoks]

def avg_acc(tsv, options):
    """ expects a vcf comp table.  takes average of 2nd two columns skipping gatk and platypus and g1kvcf """
    precs = []
    recs = []
    with open(tsv) as f:
        for line in f:
            toks = line.split()
            if toks[0] == "Graph":
                continue
            if toks[0] not in ["gatk3", "platypus", "g1kvcf"]:
                precs.append(float(toks[options.pcol]))
                recs.append(float(toks[options.rcol]))
    avg_prec = float(sum(precs)) / len(precs)
    avg_rec = float(sum(recs)) / len(recs)
    return ((avg_prec + avg_rec) / 2, avg_prec, avg_rec, tsv)

def update_best_table(best_table, tsv, options):
    """ keep track of best comp dir for each graph for each region.  will use to make a set of best calls """
    tb = os.path.splitext(os.path.basename(tsv))[0].split("-")
    region = tb[-1]
    with open(tsv) as f:
        table = [line for line in f]
    f1_lists = defaultdict(list)
    
    for line in table[1:]:
        toks = line.split("\t")
        # hardcode prec and recall columns
        graph, precision, recall = toks[0], float(toks[options.pcol]), float(toks[options.rcol])
        if precision + recall != 0:
            f1 = 2. * (precision * recall) / (precision + recall)
        else:
            f1 = 0.
        f1_lists[graph].append(f1)

    for graph, f1s in f1_lists.items():
        avg_f1 = sum(f1s) / len(f1s)
        cur_val = best_table[region][graph]
        if avg_f1 > cur_val[1]:
            best_table[region][graph] = (tsv, avg_f1)

def make_best_calls(best_table, options):
    """ using softlinks, make a call set with best f1s from the roc.  this is dependent on the 
call directories being obtainable from the comparison directory by dropping extension """
    best_dir = options.out_dir.strip("/") + ".best"
    for region in best_table.keys():
        for graph in best_table[region].keys():

            comp_tsv_path = best_table[region][graph][0]
            comp_tsv_path = comp_tsv_path[:comp_tsv_path.find("/comp_tables")]
            call_base_path = os.path.splitext(comp_tsv_path)[0]
            call_path = os.path.join(call_base_path, region, graph)
            # gatk3 and platypus: we just link in their vcf since they don't have call directory
            if graph in ["gatk3", "platypus"]:
                robust_makedirs(os.path.join(best_dir, region, graph))
                
            else:
                robust_makedirs(os.path.join(best_dir, region))
                os.system("ln -fs {} {}".format(os.path.abspath(call_path),
                                                os.path.abspath(os.path.join(best_dir, region))))
            # link in the preprocessed vcf from the comp dir to the same directory
            comp_path = os.path.join(call_base_path +".comp")
            for pvcf in glob.glob(os.path.join(comp_path, "preprocessed_vcfs", region, "*_{}.vcf".format(graph))):
                os.system("ln -fs {} {}".format(os.path.abspath(pvcf),
                                                os.path.abspath(os.path.join(best_dir, region, graph, os.path.basename(pvcf).replace(graph, "sample_preprocessed")))))
            # link in the truth while we're at it
            for pvcf in glob.glob(os.path.join(comp_path, "preprocessed_vcfs", region, "*_platvcf*.vcf")):
                os.system("ln -fs {} {}".format(os.path.abspath(pvcf),
                                                os.path.abspath(os.path.join(best_dir, region, graph))))
            
            
def main(args):
    
    options = parse_args(args)

    robust_makedirs(os.path.join(options.out_dir, "comp_tables"))

    # compute average score for each roc dir in this table
    avg_table = []
    # [region][method] --> (path, f1)
    best_table = defaultdict(lambda : defaultdict(lambda : (None, -1)))
    first = True
    # sort the directories, assuming their names give info on their order
    # in the roc
    for comp_dir in sorted(options.comp_dirs)[options.skip_first:len(options.comp_dirs) - options.skip_last]:
        print comp_dir
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
            print "{} {}".format(c, os.path.join(options.out_dir, "comp_tables",
                                                     os.path.basename(tsv)))
            tb = os.path.basename(tsv).split("-")
            if len(tb) > 2 and tb[0] == "platvcf" and tb[1] == "sompy":
                avg_table.append(avg_acc(tsv, options))
                update_best_table(best_table, tsv, options)
                
        first = False
        
    # make a call directory of links to the best in roc points for each graph
    make_best_calls(best_table, options)

    # write out our sompy vcf snp accutacy
    with open(os.path.join(options.out_dir, "platvcf-sompy-SNP-avg.tsv"), "w") as f:
        lines = sorted(avg_table)
        for line in lines:
            for tok in line:
                f.write(str(tok) + "\t")
            f.write("\n")

    # let's sort the output to make it easier to remove dead points
    if options.smooth is True:
        for tsv in glob.glob(os.path.join(options.out_dir, "comp_tables", "*.tsv")):
            print "smoothing {}".format(tsv)
            with open(tsv) as f:
                lines = [line for line in f]
                lines = [lines[0]] + sorted(lines[1:], key = lambda x : (x.split()[0], float(x.split()[options.pcol]), 1 - float(x.split()[options.rcol])))
                # precisions can be bumpy (need to change to sensitivy?)
                # use simple smoother in the meantime
                lines = smooth_table([x.split() for x in lines], options)
            with open(tsv, "w") as f:
                for line in lines:
                    f.write(line)

if __name__ == "__main__" :
    sys.exit(main(sys.argv))
        
        

        
