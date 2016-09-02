#!/usr/bin/env python2.7
"""
Run as second pass after all call(s) of callVariants.  Makes a total calls directory by merging VCFs
Other outputs (vgs and gams) are *not* merged -- just vcfs.
Needs rtg tools
"""

import argparse, sys, os, os.path, random, subprocess, shutil, itertools, glob
import doctest, re, json, collections, time, timeit, string
import signal
from collections import defaultdict
from callVariants import run
from toillib import RealTimeLogger, robust_makedirs

def parse_args(args):
    parser = argparse.ArgumentParser(description=__doc__, 
        formatter_class=argparse.RawDescriptionHelpFormatter)
        
    # General options
    parser.add_argument("call_dir", type=str,
                        help="input alignment files")
    parser.add_argument("--name", type=str, default="total",
                        help="name of output directory")

    args = args[1:]
        
    return parser.parse_args(args)

def main(args):
    
    options = parse_args(args)

    RealTimeLogger.start_master()
    # expect call_dir/<REGION>/<GRAPH>/<SAMPLE>_sample.vcf

    # count up regions
    regions = set()
    for regiondir in glob.glob(os.path.join(options.call_dir, "*")):
        if os.path.isdir(regiondir):
            region = os.path.basename(regiondir)
            # avoid crufty directories (including outputs of previous runs of this script)
            if region in ["brca1", "brca2", "mhc", "lrc_kir", "sma"]:
                regions.add(region)

    print regions

    # count up graphs (that are present in every region)
    graphs = set()
    gcount = defaultdict(int)
    for region in regions:
        for graphdir in glob.glob(os.path.join(options.call_dir, region, "*")):
            if os.path.isdir(graphdir):
                graph = os.path.basename(graphdir)
                gcount[graph] = gcount[graph] + 1
    
    for graph, count in gcount.items():
        if count == len(regions):
            graphs.add(graph)

    print graphs

    # count up samples (that are present in every graph in every region)
    samples = set()
    scount = defaultdict(int)
    for region in regions:
        for graph in graphs:
            for vcf in glob.glob(os.path.join(options.call_dir, region, graph, "*_sample.vcf")):
                sample = os.path.basename(vcf).split("_")[0]
                scount[sample] = scount[sample] + 1

    for sample, count in scount.items():
        if count == len(regions) * len(graphs):
            samples.add(sample)

    print samples

    # make our output directory
    out_dir = os.path.join(options.call_dir, options.name)
    robust_makedirs(out_dir)

    for graph in graphs:
        g_out_dir = os.path.join(out_dir, graph)

        for sample in samples:
            vcf_files = []

            for region in regions:
                vcf = os.path.join(options.call_dir, region, graph, "{}_sample.vcf".format(sample))
                assert os.path.isfile(vcf)
                vcf_files.append((region, vcf))

            # output vcf
            merge_vcf_path = os.path.join(out_dir, graph, "{}_sample.vcf".format(sample))

            # working directory for intermediates / debugging
            work_path = os.path.join(out_dir, graph, "input", sample)
            robust_makedirs(work_path)

            # preprocess all the vcfs and leave in input dir
            input_files = []
            for region, vcf in vcf_files:
                outbase = os.path.join(work_path, region)
                run("vcfsort {} > {}.vcf".format(vcf, outbase), fail_hard = True)
                run("bgzip -f {}.vcf".format(outbase))
                run("tabix -f -p vcf {}.vcf.gz".format(outbase))
                input_files.append("{}.vcf.gz".format(outbase))
            
            # run the merge
            run("rtg vcfmerge {} -o {}.gz".format(" ".join(input_files), merge_vcf_path), fail_hard = True)

            # decompress
            run("gzip -dc {}.gz > {}".format(merge_vcf_path, merge_vcf_path))
        
    return 0
    
if __name__ == "__main__" :
    sys.exit(main(sys.argv))
        
        
