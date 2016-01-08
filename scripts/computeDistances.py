#!/usr/bin/env python2.7
"""
Do an all-to-all comparison of all input graphs.  Two distance measures are used:
1)kmer set (jaccard and recall)
2)corg overlap
"""

import argparse, sys, os, os.path, random, subprocess, shutil, itertools, glob
import doctest, re, json, collections, time, timeit, string, math, copy
from Bio.Phylo.TreeConstruction import _DistanceMatrix, DistanceTreeConstructor
from Bio import Phylo
import matplotlib
matplotlib.use('Agg')
import pylab
import networkx as nx
from collections import defaultdict
from toil.job import Job
from toillib import RealTimeLogger, robust_makedirs
from callVariants import alignment_sample_tag, run
from callStats import vg_length

def parse_args(args):
    parser = argparse.ArgumentParser(description=__doc__, 
        formatter_class=argparse.RawDescriptionHelpFormatter)
    
    # Add the Toil options so the job store is the first argument
    Job.Runner.addToilOptions(parser)
    
    # General options
    parser.add_argument("graphs", nargs="+",
                        help="graphs to compare to each other")
    parser.add_argument("out_dir", type=str,
                        help="directory to store intermeidate comparison output")
    parser.add_argument("out_file", type=str,
                        help="path of tsv distance matrix to output")
    parser.add_argument("comp_type", type=str,
                        help="comparison type from {jaccard,recall,precision,corg}")
    parser.add_argument("--kmer", type=int, default=27,
                        help="kmer size for indexing")
    parser.add_argument("--edge_max", type=int, default=5,
                        help="edge-max parameter for vg kmer index")    
    parser.add_argument("--overwrite", action="store_true", default=False,
                        help="overwrite existing files (indexes and comparison output)")
    parser.add_argument("--vg_cores", type=int, default=1,
                        help="number of cores to give to vg commands")
    parser.add_argument("--avg_samples", action="store_true", default=False,
                        help="Average samples into mean value")
    parser.add_argument("--dir_tag", action="store_true", default=False,
                         help="Use directory of graph as name tag")
    parser.add_argument("--orig_tag", type=str, default="graphs",
                        help="When dir_tag used, change this tag to original")
    parser.add_argument("--skip", type=str, default="",
                        help="comma-separated list of keywords"
                        " such that input file will ignored if it contains one")
    parser.add_argument("--ignore_ns", action="store_true", default=False,
                        help="Dont include ns in kmer comparison")
    parser.add_argument("--timeout", type=int, default=sys.maxint,
                        help="timeout in seconds for long jobs (vg index and corg in this case)")
                            
    args = args[1:]

    return parser.parse_args(args)

def index_path(graph, options):
    """ get the path of the index given the graph
    """
    return graph + ".index"

def compute_kmer_index(job, graph, options):
    """ run vg index (if necessary) and vg compare on the input
    vg indexes are just created in place, ie same dir as graph,
    so need to have write permission there
    """
    out_index_path = index_path(graph, options)
    do_index = options.overwrite or not os.path.exists(out_index_path)

    index_opts = "-s -k {} -t {}".format(options.kmer, options.vg_cores)
    if options.edge_max > 0:
        index_opts += " -e {}".format(options.edge_max)

    if do_index:
        os.system("rm -rf {}".format(out_index_path))
        run("vg index {} {}".format(index_opts, graph), timeout_sec=options.timeout,
            timeout_dep=out_index_path)

def dir_tag(graph, options):
    """ optionally use directory for unique prefix
    """
    if not options.dir_tag:
        return ""
    tag = graph.split("/")[-2] + "_"
    if tag == options.orig_tag + "_":
        tag = "original_"
    return tag

def comp_path(graph1, graph2, options):
    """ get the path for json output of vg compare
    """
    name1 = dir_tag(graph1, options) + os.path.splitext(os.path.basename(graph1))[0]
    name2 = dir_tag(graph2, options) + os.path.splitext(os.path.basename(graph2))[0]
    
    return os.path.join(options.out_dir, name1 + "_vs_" + name2 + ".json")

def corg_path(graph1, graph2, options):
    """ get the path for the distance computed via corg lengths
    """
    name1 = dir_tag(graph1, options) + os.path.splitext(os.path.basename(graph1))[0]
    name2 = dir_tag(graph2, options) + os.path.splitext(os.path.basename(graph2))[0]
    
    return os.path.join(options.out_dir, name1 + "_vs_" + name2 + ".txt")

def corg_graph_path(graph1, graph2, options):
    """ get the path for vg output of corg
    """
    name1 = dir_tag(graph1, options) + os.path.splitext(os.path.basename(graph1))[0]
    name2 = dir_tag(graph2, options) + os.path.splitext(os.path.basename(graph2))[0]
    
    return os.path.join(options.out_dir, name1 + "_vs_" + name2 + ".vg")

def jaccard_dist_fn(graph1, graph2, options):
    """ scrape jaccard dist from vg compare output
    """
    jpath = comp_path(min(graph1, graph2), max(graph1, graph2), options)    
    with open(jpath) as f:
        j = json.loads(f.read())
        if float(j["union"]) == 0:
            jaccard = 2.
        else:
            jaccard = float(j["intersection"]) / float(j["union"])
        return 1. - jaccard

def recall_dist_fn(graph1, graph2, options):
    """ assymmetric version of above to compute recall of graph1 on graph2
    return 1 - recall to be consistent with other functions where similar is smaller. 
    """
    jpath = comp_path(min(graph1, graph2), max(graph1, graph2), options)
    with open(jpath) as f:
        j = json.loads(f.read())
        if graph1 <= graph2:
            db2_total = float(j["db2_total"])
        else:
            # we go the other direction if graphs flipped
            db2_total = float(j["db1_total"])
        intersection = float(j["intersection"])
        recall = intersection / db2_total
        return 1. - recall

def precision_dist_fn(graph1, graph2, options):
    """ get 1 - precision of graph1 on graph2
    """
    return recall_dist_fn(graph2, graph1, options)

def corg_dist_fn(graph1, graph2, options):
    """ scrape corg dist from corg output 
    """
    cpath = corg_path(min(graph1, graph2), max(graph1, graph2), options)
        
    with open(cpath) as f:
        c = f.readline().strip()
        dist = float(c)
        return dist

def sample_name(graph):
    """ Try to extract a sample name from the graph path.  Will succeed
    only in cases like .../<SAMPLE_ID>_[linear|augmented|sample].vg
    """
    toks = os.path.splitext(os.path.basename(graph))[0].split("_")
    if len(toks) == 2 and toks[1] in ["sample", "linear", "augmented"]:
        return toks[0]
    return ""

def different_sample(graph1, graph2):
    """ Are two graphs from different samples?  If either doesn't 
    have a sample name, then return false
    """
    sample1 = sample_name(graph1)
    sample2 = sample_name(graph2)
    return len(sample1) > 0 and len(sample2) > 0 and sample1 != sample2
    
def compute_matrix(options, dist_fn):
    """ make a distance matrix (dictionary), also write it to file
     dist_fn takes (graph1, graph2, options) and returns a float
    """
    def label_fn(graph):
        if options.avg_samples:
            # ex: NA3453456_agumented.vg -> augmented
            label = "".join(os.path.splitext(os.path.basename(graph))[0].split("_")[1:])
            if label == "":
                label = "".join(os.path.splitext(os.path.basename(graph))[0].split("_")[0])
                toks = label.split("-")
                label = toks[0]
                if label == "debruijn":
                    label += "-{}".format(toks[-1])
            label = dir_tag(graph, options) + label

            # hack (original_cactus -> cactus_original)
            if label.split("_")[0] == "original" and len(label.split("_")) > 1:
                label = "_".join(label.split("_")[1:]) + "_" + label.split("_")[0]

            assert len(label) > 0
            return label
        else:
            return dir_tag(graph, options) + os.path.splitext(os.path.basename(graph))[0]

    # make empty distance matrix and counts table (for mean)
    mat = dict()
    counts = dict()
    for graph in options.graphs:
        mat[label_fn(graph)] = defaultdict(float)
        counts[label_fn(graph)] = defaultdict(float)
            
    # fill the matrix, summing if two graphs map to same label 
    for graph1 in options.graphs:
        for graph2 in options.graphs:
            if not options.avg_samples or not different_sample(graph1, graph2):
                try:
                    val = dist_fn(graph1, graph2, options)
                except Exception as e:
                    RealTimeLogger.get().warning("comp {} vs {}: {}".format(graph1, graph2, str(e)))
                    val = None
                if None in [mat[label_fn(graph1)][label_fn(graph2)], val]:
                    mat[label_fn(graph1)][label_fn(graph2)] = None
                else:
                    mat[label_fn(graph1)][label_fn(graph2)] += val
                counts[label_fn(graph1)][label_fn(graph2)] += 1.

    # divide by counts to get mean
    for graph1 in set(map(label_fn, options.graphs)):
        for graph2 in set(map(label_fn, options.graphs)):
            if mat[graph1][graph2] is not None and counts[graph1][graph2] > 0:
                mat[graph1][graph2] /= counts[graph1][graph2]
            else:
                mat[graph1][graph2] = None

    return mat, list(set(map(label_fn, options.graphs)))

def write_tsv(out_path, mat, col_names, row_names, row_label):
    """ write tsv distance matrx
    """
    with open(out_path, "w") as f:
        # header
        f.write("{}\t".format(row_label) + "\t".join(col_names) + "\n")
        for graph1 in col_names:
            f.write(graph1)
            for graph2 in row_names:
                f.write("\t{}".format(mat[graph1][graph2]))
            f.write("\n")

def read_tsv(in_path):
    """ opposite of above
    """
    with open(in_path) as f:
        # header
        line = f.readline()
        toks = line.split()
        row_label = toks[0]
        col_names = toks[1:]
        row_names = []
        # body
        mat = []
        for line in f:
            toks = line.split()
            row_names.append(toks[0])
            toks = map(lambda x : None if x == "None" else float(x), toks[1:])
            mat.append(toks)
    return mat, col_names, row_names, row_label

def remove_nones(mat, col_names, row_names):
    """ Naive greedy remove of rows and columns with None elements.  
    idea find row or column with most Nones.  remove it.  repeat.
    haven't given this too much thought.
    """
    keep_going = True
    while keep_going is True:
        row_counts = [0 for x in range(len(row_names))]
        col_counts = [0 for x in range(len(col_names))]
        # could be moved outside loop but that'd be too clever
        for i in range(len(row_names)):
            for j in range(len(col_names)):
                if mat[row_names[i]][col_names[j]] == None:
                    row_counts[i] += 1
                    col_counts[j] += 1

        row_max = max(row_counts)
        col_max = max(col_counts)
        if row_max > 0 and row_max >= col_max:
            idx = row_counts.index(row_max)
            del mat[row_names[idx]]
            del row_names[idx]
        elif col_max > row_max:
            idx = col_counts.index(col_max)
            for i in range(len(row_names)):
                del mat[row_names[i]][col_names[idx]]
            del col_names[idx]
        else:
            keep_going = False

    return mat, col_names, row_names

def compute_kmer_comparison(job, graph1, graph2, options):
    """ run vg compare between two graphs
    """
    out_path = comp_path(graph1, graph2, options)
    graph1_index_path = index_path(graph1, options)
    assert os.path.exists(graph1_index_path)
    graph2_index_path = index_path(graph2, options)
    assert os.path.exists(graph2_index_path)

    do_comp = options.overwrite or not os.path.exists(out_path)
    ignore = ""
    if options.ignore_ns is True:
        ignore = "-i"

    if do_comp:
        if os.path.isfile(out_path):
            os.remove(out_path)
        robust_makedirs(os.path.dirname(out_path))        
        run("vg compare {} {} {} -t {} > {}".format(graph1, graph2, ignore,
                                                    min(options.vg_cores, 2), out_path))

def compute_corg_comparison(job, graph1, graph2, options):
    """ run corg on the graphs.  store the output in a text file
    """
    out_path = corg_path(graph1, graph2, options)
    corg_vg = corg_graph_path(graph1, graph2, options)
    do_comp = options.overwrite or not os.path.exists(out_path)
    if do_comp:
        if os.path.isfile(out_path):
            os.remove(out_path)
        robust_makedirs(os.path.dirname(out_path))
        run("corg {} {} -e {} -k {} -t {} > {} 2> {}".format(graph1, graph2, options.edge_max,
                                                             options.kmer, options.vg_cores, corg_vg,
                                                             out_path.replace(".txt", ".log")),
            timeout_sec=options.timeout,
            timeout_dep=out_path)
        len1 = vg_length(graph1, options)
        len2 = vg_length(graph2, options)
        lenC = vg_length(corg_vg, options)
        # corg screwing up will leave an empty vg which gives length 0
        if lenC == 0:
            corg_val = "error: corg graph not computed. see .log"
        else:
            corg_val = abs((2. * lenC) / float(len1 + len2) -1.) 

        with open(out_path, "w") as f:
            f.write("{}\n".format(corg_val))

def compute_kmer_comparisons(job, options):
    """ run vg compare in parallel on all the graphs,
    outputting a json file for each
    """
    RealTimeLogger.get().info("Running vg compare on pairs of input graphs")
    for graph1 in options.graphs:
        for graph2 in options.graphs:
            if graph1 <= graph2:
                if not options.avg_samples or not different_sample(graph1, graph2):
                    out_path = comp_path(graph1, graph2, options)
                    if options.overwrite or not os.path.exists(out_path):
                        job.addChildJobFn(compute_kmer_comparison, graph1, graph2, options,
                                          cores=min(options.vg_cores, 2))
        
def compute_corg_comparisons(job, options):
    """ run corg compare on all corg-ablegraphs. 
    """
    RealTimeLogger.get().info("Running corg on pairs of input graphs")
    for graph1 in options.graphs:
        for graph2 in options.graphs:
            if graph1 <= graph2:
                if not options.avg_samples or not different_sample(graph1, graph2):
                    out_path = corg_path(graph1, graph2, options)
                    if options.overwrite or not os.path.exists(out_path):
                        job.addChildJobFn(compute_corg_comparison, graph1, graph2,
                                          options, cores=options.vg_cores)

    
def compute_kmer_indexes(job, options):
    """ run everything (root toil job)
    first all indexes are computed,
    then all comparisons (follow on)
    then summary (follow on of that)
    """
    # do all the indexes
    for graph in options.graphs:
        job.addChildJobFn(compute_kmer_index, graph, options, cores=options.vg_cores)

    # do the comparisons
    if options.comp_type != "corg":
        job.addFollowOnJobFn(compute_kmer_comparisons, options, cores=1)
    else:
        job.addFollowOnJobFn(compute_corg_comparisons, options, cores=options.vg_cores)
    
def main(args):
    
    options = parse_args(args)

    assert options.comp_type in ["corg", "jaccard", "recall", "precision"]
    
    RealTimeLogger.start_master()

    skipList = options.skip.split(",")
    goodGraphs = []
    for graph in options.graphs:
        if os.path.splitext(graph)[1] != ".vg":
            raise RuntimeError("Input graphs expected to have .vg extension")
        toAdd = True
        for skip in skipList:
            if len(skip) > 0 and skip in graph:
                toAdd = False
                break
        if toAdd is True:
            goodGraphs.append(graph)

    options.graphs = sorted(goodGraphs)

    # Make a root job
    root_job = Job.wrapJobFn(compute_kmer_indexes, options,
        cores=1, memory="2G", disk=0)
    
    # Run it and see how many jobs fail
    toil_exception = None
    try:
        Job.Runner.startToil(root_job,  options)
    except Exception as e:
        toil_exception = e
                                   
    RealTimeLogger.stop_master()

    # munge through results to make a matrix (saved to tsv file)
    if options.comp_type == "jaccard":
        dist_fn = jaccard_dist_fn
    elif options.comp_type == "recall":
        dist_fn = recall_dist_fn
    elif options.comp_type == "precision":
        dist_fn = precision_dist_fn
    elif options.comp_type == "corg":
        dist_fn = corg_dist_fn

    mat, names = compute_matrix(options, dist_fn)

    # write a "_raw" copy with nones included
    opath, oext = os.path.splitext(options.out_file)
    write_tsv(opath + "_raw" + oext, mat, names, names, "graph")
    
    # strip out all nones and write output tsv
    mat, col_names, row_names = remove_nones(mat, names, copy.deepcopy(names))
    write_tsv(options.out_file, mat, names, names, "graph")
        
    if toil_exception is not None:
        raise Exception("{}".format(str(toil_exception)))

if __name__ == "__main__" :
    sys.exit(main(sys.argv))
        
        
        
        
        
        
        
        
        
        

