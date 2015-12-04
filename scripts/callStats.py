#!/usr/bin/env python2.7
"""
Get some basic statistics from the vg call output
"""

import argparse, sys, os, os.path, random, subprocess, shutil, itertools, glob
import doctest, re, json, collections, time, timeit, string
from collections import defaultdict
from operator import sub
from callVariants import sample_vg_path, sample_txt_path, augmented_vg_path, alignment_region_tag, alignment_graph_tag
from callVariants import graph_path, index_path, augmented_vg_path, linear_vg_path, linear_vcf_path, sample_vg_path
from callVariants import alignment_region_tag, alignment_sample_tag, alignment_graph_tag

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
    parser.add_argument("--avg_samples", action="store_true", default=False,
                        help="Average samples into mean value")
    parser.add_argument("--dir_tag", action="store_true", default=False,
                         help="Use directory of graph as name tag")
    parser.add_argument("--out_sub", type=str, default="",
                        help="make a subfolder with this name for output")
                            
    args = args[1:]
        
    return parser.parse_args(args)

def compare_out_path(options):
    """ get root output dir for comparison output
    """
    tag = ""
    if len(options.out_sub) > 0:
        tag += os.path.join(tag, options.out_sub)
    return os.path.join(options.out_dir, tag)
    
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
    cmd = "scripts/vcfCountSnps.sh {}".format(vcf)
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
        # ugly hack to get g1kvcf sample from vcf
        if vg_sample.split("/")[-2] == "g1kvcf":
            sample_snps = count_vcf_snps(vg_sample.replace(".vg", ".vcf.gz"), options)
        else:
            sample_snps = count_vg_paths(vg_sample, options)
        augmented_snps = count_vg_paths(vg_augmented, options)

        if options.avg_samples:
            # reduce to reference graph
            name_fn = lambda x : graph_path(x, options)
        else:
            # keep split apart by sample / graph corresponding to gam
            name_fn = lambda x : alignment_graph_tag(x, options)  \
                   + "_" + alignment_sample_tag(x, options)
            
        name = name_fn(gam)
        
        sums[name] = (sums[name][0] + vcf_snps,
                      sums[name][1] + sample_snps,
                      sums[name][2] + augmented_snps)
        counts[name] = counts[name] + 1

    for name in sorted(list(set(map(name_fn, options.in_gams)))):
        avg_vcf = float(sums[name][0]) / float(counts[name])
        avg_sam = float(sums[name][1]) / float(counts[name])
        avg_aug = float(sums[name][2]) / float(counts[name])
        count_table +="{}\t{}\t{}\t{}\n".format(
            name,
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

        if options.avg_samples:
            # reduce to reference graph
            name_fn = lambda x : graph_path(x, options)
        else:
            # keep split apart by sample / graph corresponding to gam
            name_fn = lambda x : alignment_graph_tag(x, options)  \
                   + "_" + alignment_sample_tag(x, options)
            
        name = name_fn(gam)

        sums[name] = (sums[name][0] + sample_snps,
                      sums[name][1] + augmented_snps,
                      sums[name][2] + original_snps)
        counts[name] = counts[name] + 1

    for name in list(set(map(name_fn, options.in_gams))):
        avg_sam = float(sums[name][0]) / float(counts[name])
        avg_aug = float(sums[name][1]) / float(counts[name])
        avg_ori = float(sums[name][2]) / float(counts[name])
        length_table +="{}\t{}\t{}\t{}\n".format(
            name,
            avg_sam,
            avg_aug,
            avg_ori)

    with open(size_tsv_path(options), "w") as ofile:
        ofile.write(length_table)

def boxPlot(inFile, outFile, x_column = 0, y_column = 1, title = None, x_label = None, y_label = None):
    """ make a box plot out of one of the tsv's generated above
    """

    # input format not quite right.  We strip off sample names from row
    # headers and put the column we want 2nd, and also skip -1 data points
    tempTsv = inFile.replace(".tsv", "_plot_{}_v_{}.tsv".format(x_column, y_column))
    with open(inFile) as f, open(tempTsv, "w") as o:
        for line in f:
            if len(line) > 0 and line[0] == "#":
                continue
            toks = line.split()
            # strip last _ and beyond
            name = "_".join(toks[x_column].split("_")[:-1])
            if name == "":
                name = toks[x_column]
            val = toks[y_column]
            missing = False
            try:
                missing = float(val) == -1.
            except:
                pass
            if not missing:
                o.write("{}\t{}\n".format(name, val))
    cmd = "scripts/boxplot.py {} --save {} --x_sideways".format(tempTsv, outFile)
    if title is not None:
        cmd += " --title {}".format(title)
    if x_label is not None:
        cmd += " --x_label {}".format(x_label)
    if y_label is not None:
        cmd += " --y_label {}".format(y_label)
    print cmd
    os.system(cmd)
                
def main(args):
    
    options = parse_args(args)

    # hack in some fake alignments so that the g1kvcf sample graphs
    # get added to the tables
    added = []
    for gam in options.in_gams:
        if os.path.dirname(gam).split("/")[-1] == "trivial":
            added.append(gam.replace("trivial/" + os.path.basename(gam),
                                     "g1kvcf/" + os.path.basename(gam)))
    options.in_gams = options.in_gams + added
        
                        
    # make some tables from the json comparison output
    snp_count_table(options)
    graph_size_table(options)

    # make some boxplots with adams super script
    boxPlot(count_tsv_path(options), count_tsv_path(options).replace(".tsv", ".pdf"),
            0, 2, "SNP\\ Count")
    boxPlot(size_tsv_path(options), size_tsv_path(options).replace(".tsv", ".pdf"),
            0, 1, "Sample\\ Graph\\ Size")
    
if __name__ == "__main__" :
    sys.exit(main(sys.argv))
        
        
        
        
        
        
        
        
        
        

