#!/usr/bin/env python2.7
"""
Do an all-to-all comparison of all input graphs.  Two distance measures are used:
1)kmer set jaccard 
2)corg overlap
Heatmaps and trees are created for each metric in the output directory.  
Several heatmap versions are written with various combinations of scaling and ymax values
"""

import argparse, sys, os, os.path, random, subprocess, shutil, itertools, glob
import doctest, re, json, collections, time, timeit, string, math
from Bio.Phylo.TreeConstruction import _DistanceMatrix, DistanceTreeConstructor
from Bio import Phylo
import matplotlib
matplotlib.use('Agg')
import pylab
import networkx as nx
from collections import defaultdict
from toil.job import Job
from toillib import RealTimeLogger, robust_makedirs
from heatmap import plotHeatMap
from callVariants import alignment_sample_tag
from callStats import vg_length

def parse_args(args):
    parser = argparse.ArgumentParser(description=__doc__, 
        formatter_class=argparse.RawDescriptionHelpFormatter)
    
    # Add the Toil options so the job store is the first argument
    Job.Runner.addToilOptions(parser)
    
    # General options
    parser.add_argument("graphs", nargs="+",
                        help="other graph(s) to compare to baseline")
    parser.add_argument("out_dir", type=str,
                        help="directory to which results will be written.")
    parser.add_argument("--kmer", type=int, default=27,
                        help="kmer size for comparison")
    parser.add_argument("--edge_max", type=int, default=5,
                        help="edge-max parameter for vg kmer index")    
    parser.add_argument("--overwrite", action="store_true", default=False,
                        help="overwrite existing files")
    parser.add_argument("--vg_cores", type=int, default=1,
                        help="number of cores to give to vg commands")
    parser.add_argument("--avg_samples", action="store_true", default=False,
                        help="Average samples into mean value")
    parser.add_argument("--dir_tag", action="store_true", default=False,
                         help="Use directory of graph as name tag")
    parser.add_argument("--orig_tag", type=str, default="graphs",
                        help="When dir_tag used, change this tag to original")
    parser.add_argument("--only_summary", action="store_true", default=False,
                        help="Only generate summary output.  Do not do any"
                        " compute")
    parser.add_argument("--no_corg", action="store_true", default=False,
                        help="Dont try to do corg comparison or heatmap")
    parser.add_argument("--no_kmer", action="store_true", default=False,
                        help="dont try to do kmer comparison or heatmap")
    parser.add_argument("--skip", type=str, default="",
                        help="comma-separated list of keywords"
                        " such that input file will ignored if it contains one")
                            
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
        os.system("vg index {} {}".format(index_opts, graph))

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
    
    return os.path.join(options.out_dir, "compare", name1 + "_vs_" + name2 + ".json")

def corg_path(graph1, graph2, options):
    """ get the path for the distance computed via corg lengths
    """
    name1 = dir_tag(graph1, options) + os.path.splitext(os.path.basename(graph1))[0]
    name2 = dir_tag(graph2, options) + os.path.splitext(os.path.basename(graph2))[0]
    
    return os.path.join(options.out_dir, "corg", name1 + "_vs_" + name2 + ".txt")

def corg_graph_path(graph1, graph2, options):
    """ get the path for vg output of corg
    """
    name1 = dir_tag(graph1, options) + os.path.splitext(os.path.basename(graph1))[0]
    name2 = dir_tag(graph2, options) + os.path.splitext(os.path.basename(graph2))[0]
    
    return os.path.join(options.out_dir, "corg", name1 + "_vs_" + name2 + ".vg")

def mat_path(options):
    """ get the path of the distance matrix
    """
    return os.path.join(options.out_dir, "distmat.tsv")

def tree_path(options, tag=""):
    """ path for newick tree
    """
    return os.path.join(options.out_dir, "tree{}.newick".format(tag))

def heatmap_path(options, tag=""):
    """ path for heatmap
    """
    return os.path.join(options.out_dir, "heatmap{}.pdf".format(tag))

def draw_len(weight):
    """ actual weights are between 0 and 1 but vary by many orders of
    magnitude.  try to map them into something for graphviz edge length hint
    """
    if weight < 0.0001:
        return .5
    elif weight < 0.001:
        return .75
    elif weight < 0.01:
        return 1.
    elif weight < 0.1:
        return 1.25
    elif weight < 0.2:
        return 1.6
    else:
        return 1.6 + weight

def jaccard_dist_fn(graph1, graph2, options):
    """ scrape jaccard dist from vg compare output
    """
    jpath = comp_path(min(graph1, graph2), max(graph1, graph2), options)    
    jaccard = -1.
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
    jaccard = -1.
    with open(jpath) as f:
        j = json.loads(f.read())
        if graph1 <= graph2:
            db2_total = float(j["db2_total"])
        else:
            # we go the other direction if graphs flipped
            db2_total = float(j["db1_total"])
        intersection = float(j["intersection"])
        recall = intersection / db2_total
        print graph1, graph2, db2_total, intersection, recall, (1-recall)
        return 1. - recall
    
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
                val = dist_fn(graph1, graph2, options)
                mat[label_fn(graph1)][label_fn(graph2)] += val
                counts[label_fn(graph1)][label_fn(graph2)] += 1.

    # divide by counts to get mean
    for graph1 in set(map(label_fn, options.graphs)):
        for graph2 in set(map(label_fn, options.graphs)):
            mat[graph1][graph2] /= max(1, counts[graph1][graph2])

    return mat, list(set(map(label_fn, options.graphs)))

def compute_tree(options, mat, names, tag):
    """ make upgma hierarchical clustering and write it as png and
    graphviz dot
    """
    # oops, convert to biopython matrix
    matrix = []
    for i in xrange(len(names)):
        row = []
        for j in xrange(i + 1):
            # tree constructor writes 0-distances as 1s for some reason
            # so we hack around here
            val = float(mat[names[i]][names[j]])
            if val == 0.:
                val = 1e-10
            elif val == 1.:
                val = 1.1
            row.append(val)
        matrix.append(row)
    dm = _DistanceMatrix(names, matrix)

    # upgma tree
    constructor = DistanceTreeConstructor()
    tree = constructor.upgma(dm)
    robust_makedirs(os.path.dirname(tree_path(options)))
    Phylo.write(tree, tree_path(options, tag), "newick")

    # png tree -- note : doesn't work in toil
    def f(x):
        if "Inner" in str(x):
            return ""
        else:
            return x
    Phylo.draw_graphviz(tree, label_func = f, node_size=1000, node_shape="s", font_size=10)
    pylab.savefig(tree_path(options, tag).replace("newick", "png"))

    # graphviz
    # get networkx graph
    nxgraph = Phylo.to_networkx(tree)
    # make undirected
    nxgraph = nx.Graph(nxgraph)
    # push names to name labels
    nxgraph = nx.convert_node_labels_to_integers(nxgraph, label_attribute="label")
    for node_id in nxgraph.nodes():
        node = nxgraph.node[node_id]
        if "Inner" in str(node["label"]):
            node["label"] = "\"\""
            node["width"] = 0.001
            node["height"] = 0.001
        else:
            node["fontsize"] = 18
    for edge_id in nxgraph.edges():
        edge = nxgraph.edge[edge_id[0]][edge_id[1]]
        # in graphviz, weight means something else, so make it a label
        weight = float(edge["weight"])
        # undo hack from above
        if weight > 1:
            weight = 1.
        if weight <= 1e-10 or weight == 1.:
            weight = 0.
        edge["weight"] = None
        edge["label"] = "{0:.3g}".format(float(weight) * 100.)
        edge["fontsize"] = 14
        edge["len"] = draw_len(weight)
    nx.write_dot(nxgraph, tree_path(options, tag).replace("newick", "dot"))
        
def compute_heatmap(options, mat, names, tag):
    """ make a pdf heatmap out of the matrix
    """
    array_mat = []
    for graph1 in names:
        array_mat.append([])
        for graph2 in names:
            array_mat[-1].append(mat[graph1][graph2])

    plotHeatMap(array_mat, names, names,
                heatmap_path(options, tag),
                leftTree=True,
                topTree=True,
                logNorm=False)

    plotHeatMap(array_mat, names, names,
                heatmap_path(options, "_log" + tag),
                leftTree=True,
                topTree=True,
                logNorm=True)

    plotHeatMap(array_mat, names, names,
                heatmap_path(options, "_vm1" + tag),
                leftTree=True,
                topTree=True,
                logNorm=False,
                vmax=1.0)

    plotHeatMap(array_mat, names, names,
                heatmap_path(options, "_log_vm1" + tag),
                leftTree=True,
                topTree=True,
                logNorm=True,
                vmax=1.0)


def cluster_comparisons(options, dist_fn, tag):
    """ write a (tsv) distance matrix
              a graphviz dot upgma tree (and png)
              a heatmap (pdf)
    dist_fn is used to make the matrix, and the tag string
    is tacked on to the results files. 
    """
    mat, names = compute_matrix(options, dist_fn)

    print mat, names, dist_fn

    compute_tree(options, mat, names, tag)

    compute_heatmap(options, mat, names, tag)

    
def compute_kmer_comparison(job, graph1, graph2, options):
    """ run vg compare between two graphs
    """
    out_path = comp_path(graph1, graph2, options)
    graph1_index_path = index_path(graph1, options)
    assert os.path.exists(graph1_index_path)
    graph2_index_path = index_path(graph2, options)
    assert os.path.exists(graph2_index_path)

    do_comp = options.overwrite or not os.path.exists(out_path)
    
    if do_comp:
        robust_makedirs(os.path.dirname(out_path))        
        os.system("vg compare {} {} -t {} > {}".format(graph1, graph2,
                                                       min(options.vg_cores, 2), out_path))

def compute_corg_comparison(job, graph1, graph2, options):
    """ run corg on the graphs.  store the output in a text file
    """
    out_path = corg_path(graph1, graph2, options)
    corg_vg = corg_graph_path(graph1, graph2, options)
    do_comp = options.overwrite or not os.path.exists(out_path)
    corg_val = -1
    if do_comp:
        robust_makedirs(os.path.dirname(out_path))
        try:
            os.system("corg {} {} -e {} -k {} -t {} > {} 2> {}".format(graph1, graph2, options.edge_max,
                                                                       options.kmer, options.vg_cores, corg_vg,
                                                                 out_path.replace(".txt", ".log")))
            len1 = vg_length(graph1, options)
            len2 = vg_length(graph2, options)
            lenC = vg_length(corg_vg, options)
            corg_val = abs((2. * lenC) / float(len1 + len2) -1.) 
        except:
            pass
    with open(out_path, "w") as f:
        f.write("{}\n".format(corg_val))
    return corg_val

def compute_kmer_comparisons(job, options):
    """ run vg compare in parallel on all the graphs,
    outputting a json file for each
    """
    if not options.no_kmer:
        for graph1 in options.graphs:
            for graph2 in options.graphs:
                if graph1 <= graph2:
                    if not options.avg_samples or not different_sample(graph1, graph2):
                        out_path = comp_path(graph1, graph2, options)
                        if options.overwrite or not os.path.exists(out_path):
                            job.addChildJobFn(compute_kmer_comparison, graph1, graph2, options,
                                              cores=min(options.vg_cores, 2))

    if not options.no_corg:
        job.addFollowOnJobFn(compute_corg_self_comparisons, options, cores=1)

def compute_corg_self_comparisons(job, options):
    """ run corg compare on all graphs with themselves.  the results are used 
    to decide if a graph is corg-able. 
    """    
    for graph1 in options.graphs:
        graph2 = graph1
        out_path = corg_path(graph1, graph2, options)
        if options.overwrite or not os.path.exists(out_path):
            job.addChildJobFn(compute_corg_comparison, graph1, graph2, options,
                              cores=options.vg_cores)

    job.addFollowOnJobFn(compute_corg_comparisons, options, cores=1)

def can_corg(graph, options, threshold):
    """ return whether a corg self compare returned distance 0 for this graph
    """
    try:
        out_path = corg_path(graph, graph, options)
        with open(out_path) as f:
            val = float(f.readline().strip())
        return val <= threshold
    except:
        return False
        
def compute_corg_comparisons(job, options):
    """ run corg compare on all corg-ablegraphs. 
    """
    for graph1 in options.graphs:
        for graph2 in options.graphs:
            if graph1 <= graph2:
                if not options.avg_samples or not different_sample(graph1, graph2):
                    out_path = corg_path(graph1, graph2, options)
                    if options.overwrite or not os.path.exists(out_path):
                        job.addChildJobFn(compute_corg_comparison, graph1, graph2,
                                          options, cores=1)

    
def compute_kmer_indexes(job, options):
    """ run everything (root toil job)
    first all indexes are computed,
    then all comparisons (follow on)
    then summary (follow on of that)
    """
    # do all the indexes
    if not options.no_kmer:
        for graph in options.graphs:
            job.addChildJobFn(compute_kmer_index, graph, options, cores=options.vg_cores)

    # do the comparisons
    job.addFollowOnJobFn(compute_kmer_comparisons, options, cores=1)
    
def main(args):
    
    options = parse_args(args)
    
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
    if not options.only_summary:
        failed_jobs = Job.Runner.startToil(root_job,  options)
    else:
        failed_jobs = 0
    
    if failed_jobs > 0:
        raise Exception("{} jobs failed!".format(failed_jobs))
                               
    RealTimeLogger.stop_master()

    # Do the drawing outside toil to get around weird import problems
    cluster_comparisons(options, jaccard_dist_fn, "_kmer")

    cluster_comparisons(options, recall_dist_fn, "_recall")
    
    if not options.no_corg:
        # hack out non-corg-able graphs
        graphs = []
        input_graphs = options.graphs
        for graph in input_graphs:
            if can_corg(graph, options, 0.0001):
                graphs.append(graph)
        options.graphs = graphs
        cluster_comparisons(options, corg_dist_fn, "_corg")

        graphs = []
        for graph in input_graphs:
            if can_corg(graph, options, 0.5):
                graphs.append(graph)
        options.graphs = graphs
        cluster_comparisons(options, corg_dist_fn, "_corg1")

if __name__ == "__main__" :
    sys.exit(main(sys.argv))
        
        
        
        
        
        
        
        
        
        

