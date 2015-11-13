#!/usr/bin/env python2.7
"""
Get some basic statistics from the vg call output
"""

import argparse, sys, os, os.path, random, subprocess, shutil, itertools, glob
import doctest, re, json, collections, time, timeit, string
from collections import defaultdict
from operator import sub
from callVariants import sample_vg_path, augmented_vg_path, alignment_region_tag, alignment_graph_tag
from callVariants import graph_path, index_path, augmented_vg_path, linear_vg_path, linear_vcf_path, sample_vg_path
from callVariants import alignment_region_tag

def parse_args(args):
    parser = argparse.ArgumentParser(description=__doc__, 
        formatter_class=argparse.RawDescriptionHelpFormatter)
        
    # General options
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


def count_tsv_path(options):
    tag = ""
    if len(options.out_sub) > 0:
        tag = options.out_sub + "_"    
    return os.path.join(compare_out_path(options),
                        "{}call_count.tsv".format(tag))

def size_tsv_path(options):
    tag = ""
    if len(options.out_sub) > 0:
        tag = options.out_sub + "_"    
    return os.path.join(compare_out_path(options),
                        "{}call_size.tsv".format(tag))
               
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

def snp_count_table(options):
    """ make a table of snp counts.  there are serious problems with this now:
    1) don't have snp count for baseline (as it's not gam or vcf)
    2) snps counted differenty for gam/vcf (multiple alternates at same site
    counted in former but not latter)
    """
    # tsv header
    count_table = "#graph\tlinear_snp_count\tsample_snp_count\taugmented_snp_count\n"

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
    length_table = "#graph\tsample_snp_length\taugmented_snp_length\toriginal_length\n"

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
    
def main(args):
    
    options = parse_args(args) 
                        
    # make some tables from the json comparison output
    snp_count_table(options)
    graph_size_table(options)
    
if __name__ == "__main__" :
    sys.exit(main(sys.argv))
        
        
        
        
        
        
        
        
        
        

