#!/usr/bin/env python2.7
"""
Get some basic statistics from the vg call output (via rocDistances.py)
"""

import argparse, sys, os, os.path, random, subprocess, shutil, itertools, glob
import doctest, re, json, collections, time, timeit, string, tempfile
from collections import defaultdict
from operator import sub
from toillib import robust_makedirs
from callVariants import sample_vg_path, sample_txt_path, augmented_vg_path, alignment_region_tag, alignment_graph_tag
from callVariants import graph_path, index_path, augmented_vg_path, linear_vg_path, linear_vcf_path, sample_vg_path
from callVariants import alignment_region_tag, alignment_sample_tag, alignment_graph_tag, run
from plotVariantsDistances import name_map

def parse_args(args):
    parser = argparse.ArgumentParser(description=__doc__, 
        formatter_class=argparse.RawDescriptionHelpFormatter)
        
    # General options
    parser.add_argument("roc_dir", type=str,
                        help="links directory _comp.best output by rocDistances.py")
    parser.add_argument("--sample", type=str, default="NA12878",
                        help="sample name")
                            
    args = args[1:]
        
    return parser.parse_args(args)


def sompy_stats(sample_vcf, truth_vcf, filter_xref, options):
    """ run sompy (copied from computeVariantsDistances, mostly) """

    out_base = tempfile.mkdtemp(prefix = "callStats_", dir = ".")

    if filter_xref is True:
        filter_vcf = os.path.join(out_base, "filter.vcf")
        os.system("grep -v XREF {} > {}".format(sample_vcf, filter_vcf))
        sample_vcf = filter_vcf

    run("som.py {} {} -P -o {} > /dev/null".format(truth_vcf, sample_vcf, os.path.join(out_base, "sp_out")), fail_hard=True)


    indels, snps = None, None
    with open(os.path.join(out_base, "sp_out.stats.csv")) as sp_result:
        for line in sp_result:
            toks = line.split(",")
            if len(toks) < 2:
                continue
            if toks[1] == "type":
                header = toks
                tp_idx = toks.index("tp")
                fp_idx = toks.index("fp")
            elif toks[1] == "indels":
                indels = toks
            elif toks[1] == "SNVs":
                snps = toks
            elif toks[1] == "records":
                total = toks

    os.system("rm -rf {}".format(out_base))

    # indels optional
    if indels is None:
        indels = [0] * 100
    if snps is None:
        snps = [0] * 100

    ret = dict()
    ret["SNP-TP"] = int(snps[tp_idx])
    ret["SNP-FP"] = int(snps[fp_idx])
    ret["INDEL-TP"] = int(indels[tp_idx])
    ret["INDEL-FP"] = int(indels[fp_idx])
    ret["TOTAL-TP"] = int(total[tp_idx])
    ret["TOTAL-FP"] = int(total[fp_idx])
    
    return ret

def vcf_stats(sample_vcf, truth_vcf, options):
    """ compute breakdown of snps indels reference non reference false positives true positives
    """    
    # compare everything
    tot_stats = sompy_stats(sample_vcf, truth_vcf, False, options)

    # compare without xref in sample
    aug_stats = sompy_stats(sample_vcf, truth_vcf, True, options)

    return tot_stats, aug_stats

def stats_row_header(options):
    """ header """
    return "\t".join(["Region", "Graph"] +
                     ["SNPS", "SNPS-TP", "SNPS-FP", "SNPS-FPR", "REF-SNPS", "REF-SNPS-TP", "REF-SNPS-FP", "REF-SNPS-FPR", "AUG-SNPS", "AUG-SNPS-TP", "AUG-SNPS-FP", "AUG-SNPS-FPR", "AUG-SNPS-FPR-DELTA"] +
                     ["INDELS", "INDELS-TP", "INDELS-FP", "INDELS-FPR", "REF-INDELS", "REF-INDELS-TP", "REF-INDELS-FP", "REF-INDELS-FPR", "AUG-INDELS", "AUG-INDELS-TP", "AUG-INDELS-FP", "AUG-INDELS-FPR", "AUG-INDELS-FPR-DELTA"] +
                     ["TOTAL", "TOTAL-TP", "TOTAL-FP", "TOTAL-FPR", "REF-TOTAL", "REF-TOTAL-TP", "REF-TOTAL-FP", "REF-TOTAL-FPR", "AUG-TOTAL", "AUG-TOTAL-TP", "AUG-TOTAL-FP", "AUG-TOTAL-FPR", "AUG-TOTAL-FPR-DELTA"])

def stats_row(region, graph, query_vcf, truth_vcf, options):
    """ make a stats line """
    tot, aug = vcf_stats(query_vcf, truth_vcf, options)

    def make_row(s):
        count = tot["{}-FP".format(s)] + tot["{}-TP".format(s)]
        tp = tot["{}-TP".format(s)]
        fp = tot["{}-FP".format(s)]
        aug_count = aug["{}-FP".format(s)] + aug["{}-TP".format(s)]
        aug_fp = aug["{}-FP".format(s)]
        aug_tp = aug["{}-TP".format(s)]
        fpr = float(fp) / count if count > 0 else 0
        aug_fpr = float(aug_fp) / aug_count if aug_count > 0 else 0
        ref_count = count - aug_count
        ref_fp = fp - aug_fp
        ref_tp = tp - aug_tp
        ref_fpr = float(ref_fp) / ref_count if ref_count > 0 else 0
        fpr_delta = aug_fpr - ref_fpr

        return [count, tp, fp, fpr, ref_count, ref_tp, ref_fp, ref_fpr, aug_count, aug_tp, aug_fp, aug_fpr, fpr_delta]

    row = [region, graph]
    row += make_row("SNP")
    row += make_row("INDEL")
    row += make_row("TOTAL")
    
    return "\t".join([str(x) for x in row])
                     
def best_stats(options):
    """ get the vcf stats from the "best" call dir as previously created by rocDistances.py"""
    # get mapping from file names to names we want to show
    names = name_map()
    sys.stdout.write(stats_row_header(options) + "\n")
    for region_path in glob.glob(os.path.join(options.roc_dir, "*")):
        region = os.path.basename(region_path)
        for graph_path in glob.glob(os.path.join(region_path, "*")):
            graph = os.path.basename(graph_path)
            if graph in names:
                graph = names[graph]
            else:
                sys.stderr.write("Unable to find name for {}".format(graph))
            truth_vcf = os.path.join(graph_path, "{}_platvcf_basline.vcf".format(options.sample))
            query_vcf = os.path.join(graph_path, "{}_sample_preprocessed.vcf".format(options.sample))
            if os.path.exists(truth_vcf) and os.path.exists(query_vcf):
                row = stats_row(region, graph, query_vcf, truth_vcf, options)
                sys.stdout.write(row + "\n")
            else:
                sys.stderr.write("Cannot find VCF input for {} {}\n".format(region, graph))
            
    
        
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

    best_stats(options)

    
if __name__ == "__main__" :
    sys.exit(main(sys.argv))
        
        
        
        
        
        
        
        
        
        

