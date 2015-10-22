#!/usr/bin/env python2.7
"""
Run vg compare on a bunch of graphs and make a table on how there kmer sets match up. 
"""

import argparse, sys, os, os.path, random, subprocess, shutil, itertools, glob
import doctest, re, json, collections, time, timeit, string
from collections import defaultdict
from operator import sub
from toil.job import Job
from toillib import RealTimeLogger, robust_makedirs
from callVariants import sample_vg_path, augmented_vg_path, alignment_region_tag, alignment_graph_tag
from callVariants import graph_path, index_path, augmented_vg_path, linear_vg_path, linear_vcf_path, sample_vg_path
from callVariants import alignment_region_tag
from clusterGraphs import comp_path

def parse_args(args):
    parser = argparse.ArgumentParser(description=__doc__, 
        formatter_class=argparse.RawDescriptionHelpFormatter)
    
    # Add the Toil options so the job store is the first argument
    Job.Runner.addToilOptions(parser)
    
    # General options
    parser.add_argument("baseline", type=str,
                        help="baseline graph to compare all other graphs to. "
                        "(just the name part, ex snp1000g. graph_dir/<baseline>-<region>.vg"
                        " must exist for each region)")
    parser.add_argument("in_gams", nargs="+",
                        help="input alignment files. (must have been run through callVariants.py!)")
    parser.add_argument("--out_dir", type=str, default="variants",
                        help="output directory (used for callVariants.py)")
    parser.add_argument("--graph_dir", type=str, default="graphs",
                        help="name of input graphs directory")
    parser.add_argument("--index_ext", type=str, default=".index",
                        help="extension to find input grpah index")
    parser.add_argument("--overwrite", action="store_true", default=False,
                        help="overwrite files if dont exist")
    parser.add_argument("--kmer", type=int, default=27,
                        help="kmer size for comparison")
    parser.add_argument("--edge_max", type=int, default=7,
                        help="edge-max parameter for vg kmer index")
    parser.add_argument("--vg_cores", type=int, default=1,
                        help="number of cores to give to vg commands")
    parser.add_argument("--avg_samples", action="store_true", default=False,
                        help="Average samples into mean value")
    parser.add_argument("--dir_tag", action="store_true", default=False,
                         help="Use directory of graph as name tag")
    parser.add_argument("--orig_tag", type=str, default="graphs",
                        help="When dir_tag used, change this tag to original")
    parser.add_argument("--out_sub", type=str, default="",
                        help="make a subfolder with this name for output")
    parser.add_argument("--only_summary", action="store_true", default=False,
                        help="Only generate summary output.  Do not do any"
                        " compute")    
                            
    args = args[1:]
        
    return parser.parse_args(args)


def index_path(graph, options):
    """ get the path of the index given the graph
    """
    return graph + ".index"

def compare_out_path(options):
    """ get root output dir for comparison output
    """
    tag = "compare"
    if len(options.out_sub) > 0:
        tag += "/" + options.out_sub
    return os.path.join(options.out_dir, tag)

def json_out_path(options):
    """ where to put the json compareison output
    """
    return os.path.join(compare_out_path(options), "json")

def dist_tsv_path(options):
    """ path where the tsv for tdistance goes
    """
    tag = ""
    if len(options.out_sub) > 0:
        tag = options.out_sub + "_"
    return os.path.join(compare_out_path(options),
                        "{}call_dist_{}.tsv".format(tag, options.baseline))

def acc_tsv_path(options):
    """ path where the tsv for prec/rec/f1 goes
    """
    tag = ""
    if len(options.out_sub) > 0:
        tag = options.out_sub + "_"    
    return os.path.join(compare_out_path(options),
                        "{}call_acc_{}.tsv".format(tag, options.baseline))

def count_tsv_path(options):
    tag = ""
    if len(options.out_sub) > 0:
        tag = options.out_sub + "_"    
    return os.path.join(compare_out_path(options),
                        "{}call_count_{}.tsv".format(tag, options.baseline))

def size_tsv_path(options):
    tag = ""
    if len(options.out_sub) > 0:
        tag = options.out_sub + "_"    
    return os.path.join(compare_out_path(options),
                        "{}call_size_{}.tsv".format(tag, options.baseline))

def baseline_path(gam, options):
    """ put together path for baseline graph 
    """
    region = alignment_region_tag(gam, options)
    if options.baseline.find("debruijn") == 0:
        dbtag = options.baseline[8:]
    else:
        dbtag = ""
    return os.path.join(options.graph_dir,
                        "{}-{}{}.vg".format(options.baseline, region, dbtag))
    
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

def compute_comparison(job, baseline, graph, options):
    """ run vg compare between two graphs
    """
    graph_index_path = index_path(graph, options)
    assert os.path.exists(graph_index_path)
    baseline_index_path = index_path(baseline, options)
    assert os.path.exists(baseline_index_path)

    out_path = comp_path(baseline, graph, options)
    do_comp = options.overwrite or not os.path.exists(out_path)
    
    if do_comp:        
        os.system("vg compare {} {} -t {} > {}".format(baseline, graph,
                                                       min(2, options.vg_cores),
                                                       out_path))
           
def count_vg_paths(vg, options):
    """ assuming output of vg call here, where one path written per snp 
    """
    if not os.path.exists(vg):
        return -1
    cmd = "vg view -j {} | jq .path | jq length".format(vg)
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                         stderr=sys.stderr, bufsize=-1)
    output, _ = p.communicate()
    assert p.wait() == 0
    num_paths = int(output.strip())
    return num_paths

def vg_length(vg, options):
    """ get sequence length out of vg stats
    """
    if not os.path.exists(vg):
        return -1
    cmd = "vg stats -l {}".format(vg)
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                         stderr=sys.stderr, bufsize=-1)
    output, _ = p.communicate()
    assert p.wait() == 0
    length = int(output.split()[1])
    return length    

def count_vcf_snps(vcf, options):
    """ get number of snps from bcftools
    """
    if not os.path.exists(vcf):
        return -1
    cmd = "./vcfCountSnps.sh {}".format(vcf)
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                         stderr=sys.stderr, bufsize=-1)
    output, _ = p.communicate()
    assert p.wait() == 0
    num_paths = int(output.strip())
    return num_paths

def jaccard_dist(comp_json_path):
    """ get a distance from the kmer set comparison json output
    """
    if not os.path.isfile(comp_json_path):
        return -1.
    with open(comp_json_path) as f:
        comp_json = json.loads(f.read())

    intersection_size = float(comp_json["intersection"])
    union_size = float(comp_json["union"])
    if union_size == 0.:
        return -1.
    else:
        return 1. - intersection_size / union_size

def accuracy(comp_json_path):
    """ compute precision, recall, f1 from the kmer set comparison json output
    (assuming db1 is the "truth")
    """
    if not os.path.isfile(comp_json_path):
        return -1., -1., -1.
    with open(comp_json_path) as f:
        comp_json = json.loads(f.read())
        
    intersection_size = float(comp_json["intersection"])
    db1_size = float(comp_json["db1_total"])
    db2_size = float(comp_json["db2_total"])
    if db2_size > 0:
        precision = intersection_size / db2_size
    else:
        precision = -1.
    if db1_size > 0:
        recall = intersection_size / db1_size
    else:
        recall = -1.
    if recall >= 0 and precision >= 0 and precision + recall > 0:
        f1 = 2. * ((precision * recall) / (precision + recall))
    else:
        f1 = -1.
    return precision, recall, f1

def dist_table(options):
    """ make the jaccard distance table by scraping together all the comparison
    json files
    """
    # tsv header
    dist_table =  "#\t{}\t\t\t\t\\t\tn".format(options.baseline)
    dist_table += "#graph\tgraph_dist\tlinear_dist\taugmented_dist\tsample_dist\tdelta_linear\tdelta_augmented\tdelta_sample\n"
    
    for gam in options.in_gams:
        baseline = baseline_path(gam, options)
        graph_comp_path = comp_path(baseline, graph_path(gam, options), options)
        graph_dist = jaccard_dist(graph_comp_path)
        aug_comp_path = comp_path(baseline, augmented_vg_path(gam, options), options)
        aug_graph_dist = jaccard_dist(aug_comp_path)
        lin_comp_path = comp_path(baseline, linear_vg_path(gam, options), options)
        lin_graph_dist = jaccard_dist(lin_comp_path)
        sam_comp_path = comp_path(baseline, sample_vg_path(gam, options), options)
        sam_graph_dist = jaccard_dist(sam_comp_path)        

        delta_lin = lin_graph_dist - graph_dist
        delta_aug = aug_graph_dist - graph_dist
        delta_sam = sam_graph_dist - graph_dist

        dist_table += "{}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\n".format(
            os.path.splitext(os.path.basename(graph_path(gam, options)))[0],
            graph_dist,
            lin_graph_dist,
            aug_graph_dist,
            sam_graph_dist,
            delta_lin,
            delta_aug,
            delta_sam)

    with open(dist_tsv_path(options), "w") as ofile:
        ofile.write(dist_table)


def acc_table(options):
    """ make the accuracy table by scraping together all the comparison
    json files
    """
    # tsv header
    acc_table =  "#\t{}\t\t\t\t\t\t\t\t\t\n".format(options.baseline)
    acc_table += "#graph\tgraph_prec\tgraph_rec\tgraph_f1"
    acc_table += "\tlinear_prec\tlinear_rec\tlinear_f1"
    acc_table += "\taugmented_prec\taugmented_rec\taugmented_f1"
    acc_table += "\tsample_prec\tsample_rec\tsample_f1"

    sums = defaultdict(lambda : (0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.))
    counts = defaultdict(lambda : 0.)
    
    for gam in options.in_gams:
        baseline = baseline_path(gam, options)
        graph_comp_path = comp_path(baseline, graph_path(gam, options), options)
        graph_acc = accuracy(graph_comp_path)
        aug_comp_path = comp_path(baseline, augmented_vg_path(gam, options), options)
        aug_graph_acc = accuracy(aug_comp_path)
        lin_comp_path = comp_path(baseline, linear_vg_path(gam, options), options)
        lin_graph_acc = accuracy(lin_comp_path)
        sam_comp_path = comp_path(baseline, sample_vg_path(gam, options), options)
        sam_graph_acc = accuracy(sam_comp_path)

        name = graph_path(gam, options)

        sums[name] = (sums[name][0] +  graph_acc[0],
                      sums[name][1] +  graph_acc[1],
                      sums[name][2] +  graph_acc[2],
                      sums[name][3] +  lin_graph_acc[0],
                      sums[name][4] +  lin_graph_acc[1],
                      sums[name][5] +  lin_graph_acc[2],
                      sums[name][6] +  aug_graph_acc[0],
                      sums[name][7] +  aug_graph_acc[1],
                      sums[name][8] +  aug_graph_acc[2],
                      sums[name][9] +  sam_graph_acc[0],
                      sums[name][10] + sam_graph_acc[1],
                      sums[name][11] + sam_graph_acc[2])

        counts[name] = counts[name] + 1
        
    for name in list(set(map(lambda x : graph_path(x, options), options.in_gams))):
        acc_table += "{}\t{:.4}\t{:.4}\t{:.4}\t".format(
            os.path.splitext(os.path.basename(graph_path(gam, options)))[0],
            float(sums[name][0]) / float(counts[name]),
            float(sums[name][1]) / float(counts[name]),
            float(sums[name][2]) / float(counts[name]))
        
        acc_table += "{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t".format(
            float(sums[name][3]) / float(counts[name]),
            float(sums[name][4]) / float(counts[name]),
            float(sums[name][5]) / float(counts[name]),
            float(sums[name][6]) / float(counts[name]),
            float(sums[name][7]) / float(counts[name]),
            float(sums[name][8]) / float(counts[name]))
        acc_table +="{:.4}\t{:.4}\t{:.4}\n".format(
            float(sums[name][9]) / float(counts[name]),
            float(sums[name][10]) / float(counts[name]),
            float(sums[name][11]) / float(counts[name]))

    with open(acc_tsv_path(options), "w") as ofile:
        ofile.write(acc_table)

def snp_count_table(options):
    """ make a table of snp counts.  there are serious problems with this now:
    1) don't have snp count for baseline (as it's not gam or vcf)
    2) snps counted differenty for gam/vcf (multiple alternates at same site
    counted in former but not latter)
    """
    # tsv header
    count_table =  "#\t{}\t\n".format(options.baseline)
    count_table += "#graph\tlinear_snp_count\tsample_snp_count\taugmented_snp_count\n"

    sums = defaultdict(lambda : (0,0,0))
    counts = defaultdict(lambda : 0)

    for gam in options.in_gams:
        linear_vcf = linear_vcf_path(gam, options) + ".gz"
        vg_sample = sample_vg_path(gam, options)
        vg_augmented = augmented_vg_path(gam, options)
        vcf_snps = count_vcf_snps(linear_vcf, options)
        sample_snps = count_vg_paths(vg_sample, options)
        augmented_snps = count_vg_paths(vg_augmented, options)

        name = graph_path(gam, options)

        sums[name] = (sums[name][0] + vcf_snps,
                      sums[name][1] + sample_snps,
                      sums[name][2] + augmented_snps)
        counts[name] = counts[name] + 1

    for name in list(set(map(lambda x : graph_path(x, options), options.in_gams))):
        avg_vcf = float(sums[name][0]) / float(counts[name])
        avg_sam = float(sums[name][1]) / float(counts[name])
        avg_aug = float(sums[name][2]) / float(counts[name])
        count_table +="{}\t{}\t{}\t{}\n".format(
            os.path.splitext(os.path.basename(name))[0],
            avg_vcf,
            avg_sam,
            avg_aug)

    with open(count_tsv_path(options), "w") as ofile:
        ofile.write(count_table)

def graph_size_table(options):
    """ make a table of sequence lengths for the vg call outputs
    """
    # tsv header
    length_table =  "#\t{}\t\n".format(options.baseline)
    length_table += "#graph\tsample_snp_length\taugmented_snp_length\toriginal_length\n"

    sums = defaultdict(lambda : (0,0,0))
    counts = defaultdict(lambda : 0)

    for gam in options.in_gams:
        linear_vcf = linear_vcf_path(gam, options) + ".gz"
        vg_sample = sample_vg_path(gam, options)
        vg_augmented = augmented_vg_path(gam, options)
        sample_snps = vg_length(vg_sample, options)
        augmented_snps = vg_length(vg_augmented, options)
        vg_original = graph_path(gam, options)
        original_snps = vg_length(vg_original, options)

        name = graph_path(gam, options)

        sums[name] = (sums[name][0] + sample_snps,
                      sums[name][1] + augmented_snps,
                      sums[name][2] + original_snps)
        counts[name] = counts[name] + 1

    for name in list(set(map(lambda x : graph_path(x, options), options.in_gams))):
        avg_sam = float(sums[name][0]) / float(counts[name])
        avg_aug = float(sums[name][1]) / float(counts[name])
        avg_ori = float(sums[name][2]) / float(counts[name])
        length_table +="{}\t{}\t{}\t{}\n".format(
            os.path.splitext(os.path.basename(name))[0],
            avg_sam,
            avg_aug,
            avg_ori)

    with open(size_tsv_path(options), "w") as ofile:
        ofile.write(length_table)

    
def compute_all_comparisons(job, options):
    """ run vg compare in parallel on all the graphs,
    outputting a json file for each
    """
    ncores = min(2, options.vg_cores)
    for gam in options.in_gams:
        baseline = baseline_path(gam, options)
        job.addChildJobFn(compute_comparison, baseline,
                          graph_path(gam, options), options, cores=ncores)
        job.addChildJobFn(compute_comparison, baseline,
                          augmented_vg_path(gam, options), options, cores=ncores)
        job.addChildJobFn(compute_comparison, baseline,
                          linear_vg_path(gam, options), options, cores=ncores)

def compute_all_indexes(job, options):
    """ run everything (root toil job)
    first all indexes are computed,
    then all comparisons (follow on)
    then summary (follow on of that)
    """

    # do all the indexes
    baseline_set = set()
        
    for gam in options.in_gams:
        baseline = baseline_path(gam, options)
        if not os.path.isfile(baseline):
            raise RuntimeError("baseline {} for gam {} not found".format(baseline, gam))
        if baseline not in baseline_set:
            job.addChildJobFn(compute_kmer_index, baseline, options, cores=options.vg_cores)
            baseline_set.add(baseline)
        if graph_path(gam, options) != baseline:
            job.addChildJobFn(compute_kmer_index, graph_path(gam, options), options, cores=options.vg_cores)
        if augmented_vg_path(gam, options) != baseline:
            job.addChildJobFn(compute_kmer_index, augmented_vg_path(gam, options), options, cores=options.vg_cores)
        if linear_vg_path(gam, options) != baseline:
            job.addChildJobFn(compute_kmer_index, linear_vg_path(gam, options), options, cores=options.vg_cores)

    # do the comparisons
    job.addFollowOnJobFn(compute_all_comparisons, options, cores=1)

    
def main(args):
    
    options = parse_args(args) 
    
    RealTimeLogger.start_master()

    for gam in options.in_gams:
        if len(gam.split("/")) < 3 or os.path.splitext(gam)[1] != ".gam":
            raise RuntimeError("Input gam paths must be of the form "
                               ".../<alg>/<reads>/<filename>.gam")
    robust_makedirs(json_out_path(options))
    robust_makedirs(compare_out_path(options))
                    
    # Make a root job
    root_job = Job.wrapJobFn(compute_all_indexes, options,
        cores=1, memory="2G", disk=0)
    
    # Run it and see how many jobs fail
    if not options.only_summary:
        failed_jobs = Job.Runner.startToil(root_job,  options)
    else:
        failed_jobs = 0
    
    if failed_jobs > 0:
        raise Exception("{} jobs failed!".format(failed_jobs))
                               
    RealTimeLogger.stop_master()

    # make some tables from the json comparison output
    #dist_table(options)
    #acc_table(options)
    snp_count_table(options)
    graph_size_table(options)
    
if __name__ == "__main__" :
    sys.exit(main(sys.argv))
        
        
        
        
        
        
        
        
        
        

