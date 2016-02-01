#!/usr/bin/env python2.7
"""
Make some figures for the .tsv output of computeVariantsDistances.py
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
from callStats import vg_length
from evaluateVariantCalls import defaultdict_set
from computeVariantsDistances import vcf_dist_header

# Set up the plot parameters
# Include both versions of the 1kg SNPs graph name
# copied from plotVariantComparison.sh
PLOT_PARAMS = [
    "--categories",
    "snp1kg",
    "snp1000g",
    "haplo1kg",
    "sbg",
    "cactus",
    "camel",
    "curoverse",
    "debruijn-k31",
    "debruijn-k63",
    "level1",
    "level2",
    "level3",
    "prg",
    "refonly",
    "simons",
    "trivial",
    "vglr",
    "haplo1kg30",
    "haplo1kg50",
    "shifted1kg",
    "gatk3",
    "platypus",
    "--category_labels ",
    "1KG",
    "1KG",
    "\"1KG Haplo\"",
    "7BG",
    "Cactus",
    "Camel",
    "Curoverse",
    "\"De Bruijn 31\"",
    "\"De Bruijn 63\"",
    "Level1",
    "Level2",
    "Level3",
    "PRG",
    "Primary",
    "SGDP",
    "Unmerged",
    "VGLR",
    "\"1KG Haplo 30\"",
    "\"1KG Haplo 50\"",
    "Control",
    "GATK3",
    "Platypus",
    "--colors",
    "\"#fb9a99\"",
    "\"#fb9a99\"",
    "\"#fdbf6f\"",
    "\"#b15928\"",
    "\"#1f78b4\"",
    "\"#33a02c\"",
    "\"#a6cee3\"",
    "\"#e31a1c\"",
    "\"#ff7f00\"",
    "\"#FF0000\"",
    "\"#00FF00\"",
    "\"#0000FF\"",
    "\"#6a3d9a\"",
    "\"#000000\"",
    "\"#b2df8a\"",
    "\"#b1b300\"",
    "\"#cab2d6\"",
    "\"#00FF00\"",
    "\"#0000FF\"",
    "\"#FF0000\"",
    "\"#25BBD4\"",
    "\"#9E7C72\"",
    "--font_size 20 --dpi 90"]



def parse_args(args):
    parser = argparse.ArgumentParser(description=__doc__, 
        formatter_class=argparse.RawDescriptionHelpFormatter)
    
    # General options
    parser.add_argument("comp_dir", type=str,
                        help="directory of comparison output written by computeVariantsDistances.py")    
                            
    args = args[1:]

    return parser.parse_args(args)

def plot_kmer_comp(tsv_path, options):
    """ take a kmer compare table and make a 
    jaccard boxplot for the first column and a 
    recall / precision ploot for the 2nd and third column
    """
    out_dir = os.path.join(options.comp_dir, "comp_plots")
    robust_makedirs(out_dir)
    out_name = os.path.basename(os.path.splitext(tsv_path)[0])
    out_base_path = os.path.join(out_dir, out_name)
    region = out_name.split("-")[-1].upper()

    params = " ".join(PLOT_PARAMS)
    # jaccard boxplot
    jac_tsv = out_base_path + "_jac.tsv"
    awkstr = '''awk '{if (NR!=1) print $1 "\t" $2}' '''
    run("{} {} > {}".format(awkstr, tsv_path, jac_tsv))
    jac_png = out_base_path + "_jac.png"
    run("scripts/boxplot.py {} --save {} --title \"{} KMER Set Jaccard\" --x_label \"Graph\" --y_label \"Jaccard Index\" --x_sideways {}".format(jac_tsv, jac_png, region, params))

    # precision recall scatter plot
    acc_tsv = out_base_path + "_acc.tsv"
    awkstr = '''awk '{if (NR!=1) print $1 "\t" $4 "\t" $3}' '''
    run("{} {} > {}".format(awkstr, tsv_path, acc_tsv))
    acc_png = out_base_path + "_acc.png"
    run("scripts/scatter.py {} --save {} --title \"{} KMER Set Accuracy\" --x_label \"Recall\" --y_label \"Precision\" --width 12 --height 9 --lines {}".format(acc_tsv, acc_png, region, params))
    
def plot_vcf_comp(tsv_path, options):
    """ take the big vcf compare table and make precision_recall plots for all the categories"""
    out_dir = os.path.join(options.comp_dir, "comp_plots")
    robust_makedirs(out_dir)
    out_name = os.path.basename(os.path.splitext(tsv_path)[0])
    out_base_path = os.path.join(out_dir, out_name)
    region = out_name.split("-")[-1].upper()

    params = " ".join(PLOT_PARAMS)

    # precision recall scatter plot
    header = vcf_dist_header(options)
    for i in range(len(header) / 2):
        prec_idx = 2 * i
        rec_idx = prec_idx + 1
        print prec_idx, header[prec_idx], rec_idx, header[rec_idx]
        ptoks = header[prec_idx].split("-")
        rtoks = header[rec_idx].split("-")
        assert ptoks[1] == "Precision"
        assert rtoks[1] == "Recall"
        assert ptoks[:1] == rtoks[:1]
        comp_cat  = ptoks[0]
        if comp_cat not in ["TOT", "SNP", "INDEL"]:
            continue
        label = header[prec_idx].replace("Precision", "acc")
        acc_tsv = out_base_path + "_" + label + ".tsv"
        print "Make {} tsv with cols {} {}".format(label, rec_idx, prec_idx)
        # +1 to convert to awk 1-base coordinates. +1 again since header doesnt include row_label col
        awkcmd = '''if (NR!=1) print $1 "\t" 1-${} "\t" ${}'''.format(prec_idx + 2, rec_idx + 2)
        awkstr = "awk \'{" + awkcmd + "}\'"
        run("{} {} > {}".format(awkstr, tsv_path, acc_tsv))
        acc_png = out_base_path + "_" + label + ".png"
        title = "VCF"
        if comp_cat == "TOT":
            title += " Total Accuracy"
        else:
            title += " {} Accuracy".format(comp_cat)
        title += " for {}".format(region)
        cmd = "scripts/scatter.py {} --save {} --title \"{}\" --x_label \"1-Precision\" --y_label \"Recall\" --width 18 --height 9 {} --lines --no_n --line_width 1.5 --marker_size 5 --min_x -0.01 --max_x 1 --min_y 0 --max_y 1.01".format(acc_tsv, acc_png, title, params)
        print cmd
        os.system(cmd)
        # top 20
        cmd = "scripts/scatter.py {} --save {} --title \"{}\" --x_label \"1-Precision\" --y_label \"Recall\" --width 18 --height 9 {} --lines --no_n --line_width 1.5--marker_size 5 --min_x 0 --max_x 0.2 --min_y 0.8 --max_y 1".format(acc_tsv, acc_png.replace(".png", "_top20.png"), title, params)
        print cmd
        os.system(cmd)
        # top 40
        cmd = "scripts/scatter.py {} --save {} --title \"{}\" --x_label \"1-Precision\" --y_label \"Recall\" --width 18 --height 9 {} --lines --no_n --line_width 1.5 --marker_size 5 --min_x 0 --max_x 0.4 --min_y 0.6 --max_y 1".format(acc_tsv, acc_png.replace(".png", "_top40.png"), title, params)
        print cmd
        os.system(cmd)
        
    
def main(args):
    
    options = parse_args(args)

    # look through tsvs in comp_tables
    for tsv in glob.glob(os.path.join(options.comp_dir, "comp_tables", "*.tsv")):
        if "kmer" in tsv.split("-"):
            plot_kmer_comp(tsv, options)
        elif "vcf" in tsv.split("-") or "sompy" in tsv.split("-") or "happy" in tsv.split("-"):
            plot_vcf_comp(tsv, options)
                                

if __name__ == "__main__" :
    sys.exit(main(sys.argv))
        
        
