#!/usr/bin/env python2.7
"""
Compare all sample graphs to baseline graphs (platvcf and g1kvcf). 
depends on callVariants.py output directory structure. Can do:
1)kmer set (jaccard and recall)
2)corg overlap
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
from callVariants import graph_path, sample_vg_path, g1k_vg_path, graph_path, sample_txt_path
from callStats import vg_length
from evaluateVariantCalls import defaultdict_set

def parse_args(args):
    parser = argparse.ArgumentParser(description=__doc__, 
        formatter_class=argparse.RawDescriptionHelpFormatter)
    
    # Add the Toil options so the job store is the first argument
    Job.Runner.addToilOptions(parser)
    
    # General options
    parser.add_argument("in_gams", nargs="+",
                        help="input alignment files")
    parser.add_argument("var_dir", type=str,
                        help="output dir for callVariants.py")
    parser.add_argument("graph_dir", type=str,
                        help="name of input graphs directory")
    parser.add_argument("comp_type", type=str,
                        help="comparison type from {kmer,corg,vcf}")    
    parser.add_argument("comp_dir", type=str,
                        help="directory to write comparison output")    
    parser.add_argument("--kmer", type=int, default=27,
                        help="kmer size for indexing")
    parser.add_argument("--edge_max", type=int, default=5,
                        help="edge-max parameter for vg kmer index")    
    parser.add_argument("--overwrite", action="store_true", default=False,
                        help="overwrite existing files (indexes and comparison output)")
    parser.add_argument("--g1kvcf_path", type=str, default="data/g1kvcf",
                        help="path to search for 1000 genomes vcf and sequences. expects "
                        "these to be in <g1kvcf_path>/BRCA1/BRCA1.vcf. etc. ")
    parser.add_argument("--platinum_path", type=str, default="data/platinum",
                        help="path to search for platinum genomes vcf. expects "
                        "these to be in <g1kvcf_path>/BRCA1/BRCA1.vcf. etc. ")        
    parser.add_argument("--vg_cores", type=int, default=1,
                        help="number of cores to give to vg commands")
    parser.add_argument("--timeout", type=int, default=sys.maxint,
                        help="timeout in seconds for long jobs (vg index and corg in this case)")
    parser.add_argument("--orig", action="store_true",
                        help="do all vs all comparison of input graphs")
    parser.add_argument("--sample", action="store_true",
                        help="do all vs all comparison of sample graphs")
    parser.add_argument("--orig_and_sample", action="store_true",
                        help="do all vs all comparison of sample + input graphs")
                            
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
    
def comp_path(graph1, graph2, options):
    """ get the path for json output of vg compare
    """
    # canonical path for pair
    if graph1 > graph2:
        graph1, graph2 = graph2, graph1    
    region1, sample1, method1 = options.tags[graph1]
    region2, sample2, method2 = options.tags[graph2]
    assert region1 == region2
    if sample1 is not None and sample2 is not None:
        assert sample1 == sample2

    s1tag = "_" + sample1 if sample1 is not None else ""
    s2tag = "_" + sample2 if sample2 is not None else ""
    return os.path.join(options.comp_dir, "compare_data", region1,
                        method1 + s1tag + "_vs_" + method2 + s2tag + ".json")

def corg_path(graph1, graph2, options):
    """ get the path for the distance computed via corg lengths
    """
    # canonical path for pair (todo corg isnt symmetric!)
    if graph1 > graph2:
        graph1, graph2 = graph2, graph1
    region1, sample1, method1 = options.tags[graph1]
    region2, sample2, method2 = options.tags[graph2]
    assert region1 == region2
    if sample1 is not None and sample2 is not None:
        assert sample1 == sample2

    s1tag = "_" + sample1 if sample1 is not None else ""
    s2tag = "_" + sample2 if sample2 is not None else ""
    return os.path.join(options.comp_dir, "corg_data", region1,
                        method1 + s1tag + "_vs_" + method2 + s2tag + ".txt")

def comp_path_vcf(graph1, graph2, options):
    """ get the path for json output of vcf compare
    """
    region1, sample1, method1 = options.tags[graph1]
    region2, sample2, method2 = options.tags[graph2]
    assert region1 == region2
    if sample1 is not None and sample2 is not None:
        assert sample1 == sample2

    s1tag = "_" + sample1 if sample1 is not None else ""
    s2tag = "_" + sample2 if sample2 is not None else ""
    return os.path.join(options.comp_dir, "vcf_compare_data", region1,
                        method1 + s1tag + "_vs_" + method2 + s2tag + ".json")

    
def corg_graph_path(graph1, graph2, options):
    """ get the path for vg output of corg
    """
    b, e = os.path.splitext(corg_path(graph1, graph2, options))
    return b + ".vg"

def out_tsv_path(options, region, category, distance):
    """ get the output tsv path
    """
    return os.path.join(options.comp_dir, "comp_tables",
                        category + "-" + distance + "-" + region + ".tsv")

def raw_tsv_path(options, region, category, distance):
    """ get the output tsv path for "raw" tables (ie with nones for missing data)
    """
    return os.path.join(options.comp_dir, "comp_tables_raw",
                        category + "-" + distance + "-" + region + ".tsv")
                        
def jaccard_dist_fn(graph1, graph2, options):
    """ scrape jaccard dist from vg compare output
    """
    jpath = comp_path(graph1, graph2, options)
    with open(jpath) as f:
        j = json.loads(f.read())
        if float(j["union"]) == 0:
            jaccard = 2.
        else:
            jaccard = float(j["intersection"]) / float(j["union"])
        return [1. - jaccard]

def recall_dist_fn(graph1, graph2, options):
    """ assymmetric version of above to compute recall of graph1 on graph2
    return recall to be consistent with other functions where similar is smaller. 
    """
    jpath = comp_path(graph1, graph2, options)
    with open(jpath) as f:
        j = json.loads(f.read())
        if index_path(graph2, options) == j["db2_path"]:
            denom = float(j["db2_total"])
        else:
            assert index_path(graph2, options) == j["db1_path"]
            denom = float(j["db1_total"])
        intersection = float(j["intersection"])
        recall = intersection / denom
        return [recall]

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
        return [dist]

def vcf_dist_fn(graph1, graph2, options):
    """ scrape vcfCompare data"""
    jpath = vcf_comp_path(graph1, graph2, options)
    with open(jpath) as f:
        j = json.loads(f.read())
        if index_path(graph2, options) == j["db2_path"]:
            denom = float(j["db2_total"])
        else:
            assert index_path(graph2, options) == j["db1_path"]
            denom = float(j["db1_total"])
        intersection = float(j["intersection"])
        recall = intersection / denom
        return [recall]

def vcf_dist_header(options):
    """ header"""
    return ["Alt-SNP-Precision",
            "Alt-SNP-Recall",
            "Alt-MB-Precision",
            "Alt-MB-Recall",
            "Alt-INS-Precision",
            "Alt-INS-Recall",
            "Alt-DEL-Precision",
            "Alt-DEL-Recall",
            "Alt-TOT-Precision",
            "Alt-TOT-Recall",
            "Allele-REF-Precision",
            "Allele-REF-Recall",
            "Allele-SNP-Precision",
            "Allele-SNP-Recall",
            "Allele-MB-Precision",
            "Allele-MB-Recall",
            "Allele-INS-Precision",
            "Allele-INS-Recall",
            "Allele-DEL-Precision",
            "Allele-DEL-Recall",
            "Allele-TOT-Precision",
            "Allele-TOT-Recall"]

def make_mat(options, row_graphs, column_graphs, dist_fns):
    """ make a distance matix """
    mat = []
    for row in row_graphs:
        mat.append([])
        for col in column_graphs:
            for dist_fn in dist_fns:
                mat[-1].append(dist_fn(row, col))
    return mat

def make_tsvs(options):
    """ make some tsv files in the output dir
    """
    # break apart by region
    for region in options.sample_graphs.keys():

        # do the baseline tsvs.  this is one row per graph,
        # with one column per comparison type with truth
        for baseline in ["g1kvcf", "platvcf", "g1kvcf_filter", "platvcf_filter"]:
            RealTimeLogger.get().info("Making {} baseline tsv for {}".format(baseline, region))
            mat = []
            row_labels = []
            if options.comp_type == "kmer":
                header = ["Jaccard-Dist", "Precision", "Recall"]
                dist_fns = [jaccard_dist_fn, precision_dist_fn, recall_dist_fn]
            elif options.comp_type == "corg":
                header = ["Corg-Dist"]
                dist_fns = [corg_dist_fn]
            else:
                header = vcf_dist_header()
                dist_fns = [vcf_dist_fn]
            for sample in options.sample_graphs[region].keys():
                for truth in options.baseline_graphs[region][sample]:
                    if options.tags[truth][2] == baseline:
                        for graph in options.sample_graphs[region][sample]:
                            row_labels.append(options.tags[graph][2])
                            row = []
                            for d in dist_fns:
                                row += d(graph, truth, options)
                            mat.append(row)
                        break # shoud not be necessary
                
            # write the baseline matrix (with None for missing data) to file 
            tsv_path = raw_tsv_path(options, region, baseline, options.comp_type)
            write_tsv(tsv_path, mat, header, row_labels, "Graph")

            # remove Nones and write tsv again
            clean_mat, clean_header, clean_row_labels = remove_nones(mat, header, row_labels)
            tsv_path = out_tsv_path(options, region, baseline, options.comp_type)
            write_tsv(tsv_path, clean_mat, clean_header, clean_row_labels, "Graph")
                
        # do the all vs all tsvs.  with averaging over samples, for heatmaps
        # todo
        
def write_tsv(out_path, mat, col_names, row_names, row_label):
    """ write tsv distance matrx
    """
    if len(mat) == 0 or len(mat[0]) == 0:
        RealTimeLogger.get().warning("Unable to write {} because input matrix empty".format(out_path))
        return
        
    robust_makedirs(os.path.dirname(out_path))
    with open(out_path, "w") as f:
        # header
        f.write("{}\t".format(row_label) + "\t".join(col_names) + "\n")
        for i, row_name in enumerate(row_names):
            f.write(row_name)
            for j, col_name in enumerate(col_names):
                f.write("\t{}".format(mat[i][j]))
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
                if mat[i][j] == None:
                    row_counts[i] += 1
                    col_counts[j] += 1

        row_max = max(row_counts)
        col_max = max(col_counts)
        # normalize by length
        row_frac_max = float(row_max) / float(len(col_counts))
        col_frac_max = float(col_max) / float(len(row_counts))
        if row_max > 0 and row_frac_max >= col_frac_max:
            idx = row_counts.index(row_max)
            del mat[idx]
            del row_names[idx]
        elif col_frac_max > row_frac_max:
            idx = col_counts.index(col_max)
            for i in range(len(row_names)):
                del mat[i][idx]
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

    if do_comp:
        if os.path.isfile(out_path):
            os.remove(out_path)
        robust_makedirs(os.path.dirname(out_path))        
        run("vg compare {} {} -i -t {} > {}".format(graph1, graph2,
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

def compute_vcf_comparison(job, graph1, graph2, options):
    """ run vcf compare between two graphs
    """
    out_path = vcf_comp_path(graph1, graph2, options)
    # we expect graph1 to be a sample graph
    region1, sample1, method1  = options.tags[graph1]
    assert method1 != "None" and "g1kvcf" not in method1 and "platvcf" not in method1

    # we expect graph2 to be a baseline graph
    region2, sample2, method2 = options.tags[graph2]
    assert method2 in ["g1kvcf", "platvcf"]
    assert region1 == region2
    assert sample1  == sample2
    
    # get the vcf of the sample graph
    query_vcf_path = sample_txt_path(graph1, options).replace(".txt", ".vcf")

    # and the baseline
    if method2 == "g1kvcf":
        truth_vcf_path = os.path.join(options.g1kvcf_path, region2.upper() + ".vcf")
    else:
        truth_vcf_path = os.path.join(options.platinum_path, sample2, region2.upper() + ".vcf")

    do_comp = options.overwrite or not os.path.exists(out_path)

    if do_comp:
        if os.path.isfile(out_path):
            os.remove(out_path)
        robust_makedirs(os.path.dirname(out_path))        
        run("vg compare {} {} -c -t > {}".format(query_vcf_path, truth_vcf_path))
        
def compute_kmer_comparisons(job, options):
    """ run vg compare in parallel on all the graphs,
    outputting a json file for each
    """
    RealTimeLogger.get().info("Running vg compare on {} pairs of input graphs".format(
        len(options.pair_comps)))
    for pair_comp in options.pair_comps:
        graph1, graph2 = pair_comp[0], pair_comp[1]
        out_path = comp_path(graph1, graph2, options)
        if options.overwrite or not os.path.exists(out_path):
            job.addChildJobFn(compute_kmer_comparison, graph1, graph2, options,
                                          cores=min(options.vg_cores, 2))
        
def compute_corg_comparisons(job, options):
    """ run corg compare on all corg-ablegraphs. 
    """
    RealTimeLogger.get().info("Running corg on pairs of input graphs")
    
    RealTimeLogger.get().info("Running vg compare on {} pairs of input graphs".format(
        len(options.pair_comps)))
    for pair_comp in options.pair_comps:
        graph1, graph2 = pair_comp[0], pair_comp[1]
        out_path = comp_path(graph1, graph2, options)
        if options.overwrite or not os.path.exists(out_path):
            job.addChildJobFn(compute_corg_comparison, graph1, graph2, options,
                                          cores=min(options.vg_cores, 2))

def compute_vcf_comparisons(job, options):
    """ run vg compare in parallel on all the graphs,
    outputting a json file for each
    """
    RealTimeLogger.get().info("Running vcf comparison {} pairs of input graphs".format(
        len(options.pair_comps)))
    for pair_comp in options.pair_comps:
        graph1, graph2 = pair_comp[0], pair_comp[1]
        out_path = vcf_comp_path(graph1, graph2, options)
        if options.overwrite or not os.path.exists(out_path):
            job.addChildJobFn(compute_vcf_comparison, graph1, graph2, options,
                                          cores=min(options.vg_cores, 2))

def compute_kmer_indexes(job, options):
    """ run everything (root toil job)
    first all indexes are computed,
    then all comparisons (follow on)
    then summary (follow on of that)
    """
    # do all the indexes
    input_set = set()
    for pair_comp in options.pair_comps:
        input_set.add(pair_comp[0])
        input_set.add(pair_comp[1])

    RealTimeLogger.get().info("Computing indexes for {} input graphs".format(len(input_set)))

    if options.comp_type in ["kmer", "corg"]:
        for graph in input_set:
            if options.overwrite or not os.path.exists(index_path(graph, options)):
                job.addChildJobFn(compute_kmer_index, graph, options, cores=options.vg_cores)

    # do the comparisons
    if options.comp_type == "kmer":
        job.addFollowOnJobFn(compute_kmer_comparisons, options, cores=1)
    elif options.comp_type == "corg":
        job.addFollowOnJobFn(compute_corg_comparisons, options, cores=1)
    elif options.comp_type == "vcf":
        job.addFollowOnJobFn(compute_vcf_comparisons, options, cores=1)

def breakdown_gams(in_gams, orig, orig_and_sample, options):
    """ use callVariants methods to find all the relevant graphs given
    a list of input gams"""

    # sort through the input, and make some dictionaries splitting
    # up the different vg files by region, sample, algorithm. 
    orig_graphs = defaultdict(set)
    sample_graphs = defaultdict(defaultdict_set)
    baseline_graphs = defaultdict(defaultdict_set)
    #other direction
    tags = dict()

    for input_gam in in_gams:
        region = alignment_region_tag(input_gam, options)
        sample = alignment_sample_tag(input_gam, options)
        method = alignment_graph_tag(input_gam, options)

        orig_path = graph_path(input_gam, options)
        sample_path = sample_vg_path(input_gam, options)
        g1kvcf_path = g1k_vg_path(input_gam, False, False, options)
        platvcf_path = g1k_vg_path(input_gam, True, False, options)
        g1kvcf_filter_path = g1k_vg_path(input_gam, False, True, options)
        platvcf_filter_path = g1k_vg_path(input_gam, True, True, options)
        
        
        if orig or orig_and_sample:
            orig_graphs[region].add(orig_path)
            tags[orig_path] = (region, None, method)

        sample_graphs[region][sample].add(sample_path)
        # we dont expect to have baselines for every sample
        if os.path.isfile(g1kvcf_path):
            baseline_graphs[region][sample].add(g1kvcf_path)
        if os.path.isfile(platvcf_path):
            baseline_graphs[region][sample].add(platvcf_path)
        if options.comp_type != "vcf" and os.path.isfile(g1kvcf_filter_path):
            baseline_graphs[region][sample].add(g1kvcf_filter_path)
        if optoins.comp_type != "vcf" and os.path.isfile(platvcf_filter_path):
            baseline_graphs[region][sample].add(platvcf_filter_path)

        tags[sample_path] = (region, sample, method)
        tags[g1kvcf_path] = (region, sample, "g1kvcf")
        tags[platvcf_path] = (region, sample, "platvcf")
        tags[g1kvcf_filter_path] = (region, sample, "g1kvcf_filter")
        tags[platvcf_filter_path] = (region, sample, "platvcf_filter")

    return orig_graphs, sample_graphs, baseline_graphs, tags
    
def main(args):
    
    options = parse_args(args)

    assert options.comp_type in ["corg", "kmer", "vcf"]
    if options.comp_type == "vcf":
        assert not options.orig and not options.orig_and_sample
    
    RealTimeLogger.start_master()

    # since we re-use callVariants methods
    options.out_dir = options.var_dir

    # find all the graphs
    breakdown = breakdown_gams(options.in_gams, options.orig,
                               options.orig_and_sample, options)
    options.orig_graphs = breakdown[0]
    options.sample_graphs = breakdown[1]
    options.baseline_graphs = breakdown[2]
    options.tags = breakdown[3]

    options.pair_comps = []
    
    # determine all pairwise comparisons between original gaphs
    for region in options.orig_graphs.keys():
        for graph1 in options.orig_graphs[region]:
            for graph2 in options.orig_graphs[region]:
                if graph1 <= graph2:
                    options.pair_comps.append((graph1, graph2))
                    
            # optional original vs sample
            if options.sample_and_original:
                for sample in options.sample_graphs[region].keys():
                    for graph2 in options.sample_graphs[region][sample]:
                        options.pair_comps.append((graph1, graph2))
                

    # now all sample vs baseline comparisons
    for region in options.sample_graphs.keys():
        for sample in options.sample_graphs[region].keys():
            for graph1 in options.sample_graphs[region][sample]:
                for graph2 in options.baseline_graphs[region][sample]:
                    # baseline graphs are expected not to exist sometimes
                    if os.path.isfile(graph2):
                        options.pair_comps.append((graph1, graph2))

                # optional smaple vs sample
                if options.sample:
                    for graph2 in options.sample_graphs[region][sample]:
                        if graph1 <= graph2:
                            options.pair_comps.append((graph1, graph2))
                        
    # Make a root job
    root_job = Job.wrapJobFn(compute_kmer_indexes, options,
        cores=1, memory="2G", disk=0)
    
    # Run it and see how many jobs fail
    Job.Runner.startToil(root_job,  options)
                                   
    # munge through results to make a matrix (saved to tsv file)
    make_tsvs(options)

    RealTimeLogger.stop_master()
        

if __name__ == "__main__" :
    sys.exit(main(sys.argv))
        
        
        
        
        
        
        
        
        
        

