#!/usr/bin/env python2.7
"""
Compute concordance statistics for all gams that are children in trios.  This is 
determined by checking to see if there are two other gams with same path but 
with number at end of names incremented by 1 and 2.  
"""

import argparse, sys, os, os.path, random, subprocess, shutil, itertools, glob
import doctest, re, json, collections, time, timeit, string
from collections import defaultdict
from operator import sub
from toil.job import Job
from toillib import RealTimeLogger, robust_makedirs
from callVariants import sample_vg_path, sample_txt_path, augmented_vg_path, alignment_region_tag, alignment_graph_tag
from callVariants import graph_path, index_path, augmented_vg_path, linear_vg_path, linear_vcf_path, sample_vg_path
from callVariants import alignment_region_tag, alignment_sample_tag, alignment_graph_tag
from callStats import boxPlot

def parse_args(args):
    parser = argparse.ArgumentParser(description=__doc__, 
        formatter_class=argparse.RawDescriptionHelpFormatter)

    # Add the Toil options so the job store is the first argument
    Job.Runner.addToilOptions(parser)

    # General options
    parser.add_argument("in_gams", nargs="+",
                        help="input alignment files. (must have been run through callVariants.py!)")
    parser.add_argument("output_path", type=str,
                        help="directory to write output to.  not to be confused with --out_dir"
                        "which is the output directory used for callVariants.py")    
    parser.add_argument("--out_dir", type=str, default="variants",
                        help="output directory given to callVariants.py")
    parser.add_argument("--overwrite", action="store_true", default=False,
                        help="Only compute trios if output not found on disk")
    parser.add_argument("--tag", type=str, default=None,
                        help="Tag to add to ouput files")
    parser.add_argument("--vg_cores", type=int, default=1,
                        help="number of cores to give to vg commands")

                            
    args = args[1:]
        
    return parser.parse_args(args)

def concordance_path(gam, options):
    """ output path for call to trioConcordance.py """
    calls_txt = sample_txt_path(gam, options)
    return calls_txt.replace("_sample.txt", "_concord.txt")

def trio_tsv_path(options):
    """ output table path """
    name = "trio_concordance"
    if options.tag is not None:
        name = name + "_" + options.tag
    return os.path.join(options.output_path, name + ".tsv")

def get_relative(gam, offset, options):
    """ search gam list for relative by offset"""
    try:
        base, ext = os.path.splitext(gam)
        idx = int(base[-3:])
        rel_name = "{}{}{}".format(base[:-3], idx + offset, ext)
        if rel_name in options.in_gams:
            return rel_name
    except:
        return None
    return None

def compute_concordance(job, child_txt_path, mom_txt_path, dad_txt_path,
                        concord_path, options):
    """ run trioConcordance on a trio of text call files
    """
    cmd = "scripts/trioConcordance.py {} {} {} > {}".format(child_txt_path,
                                                            mom_txt_path,
                                                            dad_txt_path,
                                                            concord_path)
    RealTimeLogger.get().info("RUN: {}".format(cmd))
    os.system(cmd)

def compute_all_trios(job, options):
    """ run trioConcordance.py in parallel on all trios found in the input
    """

    for gam in options.in_gams:
        mom_gam = get_relative(gam, -1, options)
        dad_gam = get_relative(gam, -2, options)
        if mom_gam is not None and dad_gam is not None:
            concord_path = concordance_path(gam, options)
            if options.overwrite or not os.path.isfile(concord_path):
                child_txt_path = sample_txt_path(gam, options)
                mom_txt_path = sample_txt_path(mom_gam, options)
                dad_txt_path = sample_txt_path(dad_gam, options)
                job.addChildJobFn(compute_concordance, child_txt_path, mom_txt_path,
                                  dad_txt_path, concord_path, options)

def concordance_tsv(options):
    """ scrape up all the concordance data and make a table
    """

    table = "#region\tgraph\tsample\tconcordant\tdiscordant\tconcordance\n"

    # map child gam path to pct concordance
    for gam in options.in_gams:
        mom_gam = get_relative(gam, -1, options)
        dad_gam = get_relative(gam, -2, options)
        if mom_gam is not None and dad_gam is not None:
            concord_path = concordance_path(gam, options)
            results = (-1, -1, -1)
            with open(concord_path) as f:
                line = f.readline()
                toks = line.split()
                results = (int(toks[0]), int(toks[1]), float(toks[2]))
                
                sample = alignment_sample_tag(gam, options)
                region = alignment_region_tag(gam, options)
                graph = alignment_graph_tag(gam, options)

                table += "{}\t{}\t{}\t{}\t{}\t{}\n".format(region,
                                                           graph,
                                                           sample,
                                                           results[0],
                                                           results[1],
                                                           results[2])

    with open(trio_tsv_path(options), "w") as out_file:
        out_file.write(table)
            
def main(args):
    
    options = parse_args(args)

    robust_makedirs(options.output_path)
    
    RealTimeLogger.start_master()

    # Make a root job
    root_job = Job.wrapJobFn(compute_all_trios, options,
                             cores=1, memory="2G", disk=0)
    
    # Run it and see how many jobs fail
    failed_jobs = Job.Runner.startToil(root_job,  options)
    
    if failed_jobs > 0:
        raise Exception("{} jobs failed!".format(failed_jobs))
                               
    RealTimeLogger.stop_master()

    # write the table
    concordance_tsv(options)

    # write boxplot
    boxPlot(trio_tsv_path(options), trio_tsv_path(options).replace(".tsv", ".pdf"),
            1, 5, "Trio\\ Concordance")

    
if __name__ == "__main__" :
    sys.exit(main(sys.argv))
