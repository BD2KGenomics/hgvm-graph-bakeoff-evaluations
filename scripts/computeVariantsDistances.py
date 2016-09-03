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
from evaluateVariantCalls import defaultdict_set
from vcfQualStats import vcf_qual_stats, balance_tables

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
                        help="comparison type from {kmer,corg,vcf,sompy,happy,vcfeval}")    
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
                        "these to be in <g1kvcf_path>BRCA1.vcf. etc. ")
    parser.add_argument("--platinum_path", type=str, default="data/platinum",
                        help="path to search for platinum genomes vcf. expects "
                        "these to be in <platinum_path>/<sample>/BRCA1.vcf. etc. ")
    parser.add_argument("--chrom_fa_path", type=str, default="data/g1kvcf/chrom.fa",
                        help="fasta file with entire chromosome info for all regions")
    parser.add_argument("--happy_fa_path", type=str, default="data/g1kvcf/chrom2.fa",
                        help="fasta file with chrXX names for chromosomes. todo- use for above")
    parser.add_argument("--gatk3_path", type=str, default="data/gatk3",
                        help="path to search for gatk3 vcf. expects "
                        " these to bein <gatk3_path>/<sample>/BRCA1.vcf etc.")
    parser.add_argument("--platypus_path", type=str, default="data/platypus",
                        help="path to search for platypus vcf. expects "
                        " these to bein <platypus_path>/<sample>/BRCA1.vcf etc.")
    parser.add_argument("--freebayes_path", type=str, default="data/freebayes",
                        help="path to search for freebayes vcf. expects "
                        " these to bein <freebayes_path>/<sample>/BRCA1.vcf etc.")
    parser.add_argument("--samtools_path", type=str, default="data/samtools",
                        help="path to search for samtools vcf. expects "
                        " these to bein <samtools_path>/<sample>/BRCA1.vcf etc.")
    parser.add_argument("--vg_cores", type=int, default=1,
                        help="number of cores to give to vg commands (and hap.py)")
    parser.add_argument("--timeout", type=int, default=sys.maxint,
                        help="timeout in seconds for long jobs (vg index and corg in this case)")
    parser.add_argument("--orig", action="store_true",
                        help="do all vs all comparison of input graphs")
    parser.add_argument("--sample", action="store_true",
                        help="do all vs all comparison of sample graphs")
    parser.add_argument("--orig_and_sample", action="store_true",
                        help="do all vs all comparison of sample + input graphs")
    parser.add_argument("--ignore", action="append", default=[],
                        help="keyword to ignore in vcf comparison")
    parser.add_argument("--normalize", action="store_true", default =False,
                        help="run vt normalization on all input vcfs")
    parser.add_argument("--clip", type=str, default=None,
                        help="clip vcf using specified bed file before call comparisons")
    parser.add_argument("--clip_fp", type=str, default=None,
                        help="false positives outside region will be called unknown (hap.py or som.py)")
    parser.add_argument("--roc", action="store_true", default=False,
                        help="generate happy rocs for gatk3 and platypus")
    parser.add_argument("--qpct", type=float, default=None,
                        help="apply quality percentile filter for gatk and platypus and freebayes")
    parser.add_argument("--qgraph", action="store_true", default=False,
                        help="apply quality percentile filter to graphs too")
    parser.add_argument("--baseline", type=str, default="platvcf",
                        help="baseline to use (platvcf or g1kvcf) for vcf comparisons")
    parser.add_argument("--gt", action="store_true",
                        help="take into account genotype information (sompy or vcfeval)")
    parser.add_argument("--new", action="store_true",
                        help="use new caller (vg genotype)")
    parser.add_argument("--min_ll", type=float, default=None,
                        help="apply likelihood filter to vg call vcfs")
    parser.add_argument("--filter_type", type=str, default="xaad",
                        help="used for vcfFilter for curves (when not --new): {xaad, ad, ll, xl}")
    parser.add_argument("--dedupe", action="store_true", default=False,
                        help="use --dedupe option in vcfFilterQuality.py")
    parser.add_argument("--vroc", action="store_true", default=False,
                        help="use vcfevals roc logic (only gives total, not indel snp breakdown) and wont work with clipping")
    parser.add_argument("--cwd", default=os.getcwd(),
                        help="set Toil job working directory")
                            
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
    
    # Move to the appropriate working directory from wherever Toil dropped us
    os.chdir(options.cwd)
    
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
    return os.path.join(options.comp_dir, "kmer_compare_data", region1,
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
    return os.path.join(options.comp_dir, "corg_compare_data", region1,
                        method1 + s1tag + "_vs_" + method2 + s2tag + ".txt")

def input_vcf_path(graph, options, region = None, sample = None, method = None):
    """ translate a gam to a vcf, with hack ot get at the baseline graphs"""
    if region is None or sample is None or method is None:
        region, sample, method = options.tags[graph]
    # path names are a bit hacky here, but work for platinum samples
    if method == "g1kvcf":
        return os.path.join(options.platinum_path, sample + ".g1k", region.upper() + ".vcf")
    elif method == "platvcf":
        return os.path.join(options.platinum_path, sample, region.upper() + ".vcf")
    elif method == "gatk3":
        return os.path.join(options.gatk3_path, sample, region.upper() + ".vcf")
    elif method == "platypus":
        return os.path.join(options.platypus_path, sample, region.upper() + ".vcf")
    elif method == "freebayes":
        return os.path.join(options.freebayes_path, sample, region.upper() + ".vcf")
    elif method == "samtools":
        return os.path.join(options.samtools_path, sample, region.upper() + ".vcf")    
    else:
        return graph.replace(".vg", ".vcf")    
 
def preprocessed_vcf_path(graph, options):
    """ get the path of the sorted vcf (normalized and/or clipped)
    """
    region, sample, method = options.tags[graph]
    # little hack because preprocessing different depending if graph is
    # baseline or not (note it should never be both!)
    if method == options.baseline:
        method = "{}_basline".format(method)
    return os.path.join(options.comp_dir, "preprocessed_vcfs", region,
                        sample + "_" + method + ".vcf")

def clip_bed_path(graph, options):
    """ get a bed to clip vcf by """
    return options.clip
    
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
    
def comp_path_sompy(graph1, graph2, options):
    """ get the path for json output of vcf compare
    """
    region1, sample1, method1 = options.tags[graph1]
    region2, sample2, method2 = options.tags[graph2]
    assert region1 == region2
    if sample1 is not None and sample2 is not None:
        assert sample1 == sample2

    s1tag = "_" + sample1 if sample1 is not None else ""
    s2tag = "_" + sample2 if sample2 is not None else ""
    return os.path.join(options.comp_dir, "sompy_compare_data", region1,
                            method1 + s1tag + "_vs_" + method2 + s2tag + ".stats.csv")

def comp_path_happy(graph1, graph2, options):
    """ get the path for json output of vcf compare
    """
    region1, sample1, method1 = options.tags[graph1]
    region2, sample2, method2 = options.tags[graph2]
    assert region1 == region2
    if sample1 is not None and sample2 is not None:
        assert sample1 == sample2

    s1tag = "_" + sample1 if sample1 is not None else ""
    s2tag = "_" + sample2 if sample2 is not None else ""
    return os.path.join(options.comp_dir, "happy_compare_data", region1,
                            method1 + s1tag + "_vs_" + method2 + s2tag + ".summary.csv")

def comp_path_vcfeval(graph1, graph2, options):
    """ get the path for json output of vcf compare
    """
    region1, sample1, method1 = options.tags[graph1]
    region2, sample2, method2 = options.tags[graph2]
    assert region1 == region2
    if sample1 is not None and sample2 is not None:
        assert sample1 == sample2

    s1tag = "_" + sample1 if sample1 is not None else ""
    s2tag = "_" + sample2 if sample2 is not None else ""
    return os.path.join(options.comp_dir, "vcfeval_compare_data", region1,
                            method1 + s1tag + "_vs_" + method2 + s2tag + ".vcfeval")
    
def corg_graph_path(graph1, graph2, options):
    """ get the path for vg output of corg
    """
    b, e = os.path.splitext(corg_path(graph1, graph2, options))
    return b + ".vg"

def out_tsv_path(options, region, category, distance, sample = None):
    """ get the output tsv path
    """
    rtag = region if sample is None else region + "-" + sample
    return os.path.join(options.comp_dir, "comp_tables",
                        category + "-" + distance + "-" + rtag + ".tsv")

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

def raw_tsv_path(options, region, category, distance, sample = None):
    """ get the output tsv path for "raw" tables (ie with nones for missing data)
    """
    rtag = region if sample is None else region + "-" + sample
    return os.path.join(options.comp_dir, "comp_tables_raw",
                        category + "-" + distance + "-" + rtag + ".tsv")
                        
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
        return [[1. - jaccard]]

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
        return [[recall]]

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
        return [[dist]]

def vcf_dist_fn(graph1, graph2, options):
    """ scrape vcfCompare data"""
    jpath = comp_path_vcf(graph1, graph2, options)
    with open(jpath) as f:
        j = json.loads(f.read())
        path1 = j["Path1"]
        path2 = j["Path2"]

        query_vcf_path = preprocessed_vcf_path(graph1, options)

        # we expect graph2 to be a baseline graph
        region2, sample2, method2 = options.tags[graph2]
        assert method2 in ["g1kvcf", "platvcf"]
        truth_vcf_path = preprocessed_vcf_path(graph2, options)

        # do we need to flip ever?
        assert path1 == query_vcf_path
        assert path2 == truth_vcf_path
        return [[j["Alts"]["SNP"]["Precision"],
                j["Alts"]["SNP"]["Recall"],
                j["Alts"]["MULTIBASE_SNP"]["Precision"],
                j["Alts"]["MULTIBASE_SNP"]["Recall"],
                j["Alts"]["INDEL"]["Precision"],
                j["Alts"]["INDEL"]["Recall"],
                j["Alts"]["TOTAL"]["Precision"],
                 j["Alts"]["TOTAL"]["Recall"],
                 0]]

def vcf_dist_header(options):
    """ header"""
    return ["SNP-Precision",
            "SNP-Recall",
            "MB-Precision",
            "MB-Recall",
            "INDEL-Precision",
            "INDEL-Recall",
            "TOT-Precision",
            "TOT-Recall",
            "QUAL"]

def sompy_dist_fn(graph1, graph2, options):
    jpath = comp_path_sompy(graph1, graph2, options)
    header = None
    snps = None
    indels = None
    total = None

    # read sompy csv output.  be a little flexible in terms of row column order (but not names!)
    RealTimeLogger.get().warning(jpath)
    with open(jpath) as f:
        for line in f:
            toks = line.split(",")
            if len(toks) < 2:
                continue
            if toks[1] == "type":
                header = toks
                rec_idx = toks.index("recall")
                prec_idx = toks.index("precision")
            elif toks[1] == "indels":
                indels = toks
            elif toks[1] == "SNVs":
                snps = toks
            elif toks[1] == "records":
                total = toks

    # indels optional
    if indels is None:
        indels = [0] * 100

    # shoehorn into vcfCompre style output (todo, revise this)
    return [[snps[prec_idx],
            snps[rec_idx],
            0,
            0,
            indels[prec_idx],
            indels[rec_idx],
            total[prec_idx],
            total[rec_idx],
             0]]

def happy_dist_fn(graph1, graph2, options):
    jpath = comp_path_happy(graph1, graph2, options)
    header = None
    snps = None
    indels = None
    total = None

    if options.roc is True and options.tags[graph1][2] in ["gatk3", "platypus", "g1kvcf", "freebayes", "samtools"]:
        # read happy roc output.
        # todo : indels and total:  problem= rocs have differen numbers of lines wwhich doesnt fit interface as is
        snp_roc_path = jpath.replace("summary.csv", "roc.snp.all.tsv")
        rows = []
        with open(snp_roc_path) as f:
            for line in f:
                toks = line.split()
                if "precision" in toks:
                    prec_idx = toks.index("precision")
                    rec_idx = toks.index("recall")
                else:
                    rows.append([toks[prec_idx], toks[rec_idx],
                                 0, 0,
                                 0, 0,
                                 0, 0])
        return rows
    else:
        # read happy csv output.  be a little flexible in terms of row column order (but not names!)
        with open(jpath) as f:
            for line in f:
                toks = line.strip().split(",")
                if len(toks) < 2:
                    continue
                if "METRIC.Recall" in toks:
                    header = toks
                    rec_idx = toks.index("METRIC.Recall")
                    prec_idx = toks.index("METRIC.Precision")
                elif toks[0] == "INDEL":
                    indels = toks
                elif toks[0] == "SNP":
                    snps = toks
                elif toks[0] == "Locations":
                    total = toks

        # indels optional
        if indels is None:
            indels = [0] * 100
        # toal optioal
        if total is None:
            total = [0] * 100
        
        # shoehorn into vcfCompre style output (todo, revise this)
        return [[snps[prec_idx],
                 snps[rec_idx],
                 0,
                 0,
                 indels[prec_idx],
                 indels[rec_idx],
                 total[prec_idx],
                 total[rec_idx],
                 0]]

def save_vcfeval_stats(out_path, fn_table, fp_table, tp_table):
    """ write some summary counts from the vceval vcf output """
    def write_table(t, name):
        with open(os.path.join(out_path, "comp_counts_{}.tsv".format(name)), "w") as f:
            for line in t:
                f.write("{}\t{}\t{}\t{}\n".format(line[0], line[1], line[2], line[3]))

    write_table(fn_table, "fn")
    write_table(fp_table, "fp")
    write_table(tp_table, "tp")
            
def load_vcfeval_stats(out_path):
    """ read the counts back from file as list """
    def read_table(name):
        t = []
        with open(os.path.join(out_path, "comp_counts_{}.tsv".format(name))) as f:
            for line in f:
                if len(line) > 0:
                    row = line.split("\t")
                    t += [[float(row[0]), int(row[1]), int(row[2]), int(row[3])]]
        return t

    fn_table = read_table("fn")
    fp_table = read_table("fp")
    tp_table = read_table("tp")
    balance_tables(fn_table, fp_table, tp_table)

    return fn_table, fp_table, tp_table
                
def vcf_num_records(vcf_path, bed_path = None):
    """ use bcftools stats to get the number of snps indels other in vcf """
    if not os.path.exists(vcf_path):
        return -1
    cmd = "bcftools stats {}".format(vcf_path)
    if bed_path is not None:
        cmd += " -R {}".format(bed_path)
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                         stderr=sys.stderr, bufsize=-1)
    output, _ = p.communicate()
    assert p.wait() == 0
    hits = 0
    for line in output.split("\n"):
        toks = line.split("\t")
        if len(toks) == 4 and toks[0] == "SN":
            if toks[2] == "number of SNPs:":
                num_snps = int(toks[3])
                hits += 1
            if toks[2] == "number of MNPs:":
                num_mnps = int(toks[3])
                hits += 1
            elif toks[2] == "number of indels:":
                num_indels = int(toks[3])
                hits += 1
            elif toks[2] == "number of records:":
                hits += 1
                num_records = int(toks[3])
            elif toks[2] == "number of others:":
                hits += 1
                num_other = int(toks[3])
        if hits == 5:
            break
                
    # generally, mnps are snps for our purposes (TODO: can separate?)
    num_snps += num_mnps
    return num_snps, num_indels, num_other

def vcfeval_dist_fn(graph1, graph2, options):
    out_dir = comp_path_vcfeval(graph1, graph2, options)
    summary_path = os.path.join(out_dir, "summary.txt")

    # read in the text summary, looks like this
    #2016-04-25 08:59:09 Threshold  True-pos  False-pos  False-neg  Precision  Sensitivity  F-measure
    #----------------------------------------------------------------------------
    #     None       985       1398       1075     0.4133       0.4782     0.4434

    with open(summary_path) as f:
        lines = [line for line in f]
        header = lines[0].split()
        prec_idx = header.index("Precision")
        rec_idx = header.index("Sensitivity")        
        data = lines[-1].split()
        total_precision = float(data[prec_idx])
        total_recall = float(data[rec_idx])

    ## use vcfevals roc
    if options.vroc is True:
        roc_path = os.path.join(out_dir, "weighted_roc.tsv")
        run("gzip -dc {} > {}".format(roc_path + ".gz", roc_path))
        ret = []
        with open(roc_path) as f:
            for line in f:
                if line[0] != "#" and len(line) > 0:
                    toks = line.split("\t")
                    precision, sensitivity = float(toks[4]), float(toks[5])
                    ret += [[0, 0, 0, 0, 0, 0, precision, sensitivity, toks[0]]]
        return ret

    # it's actually better to compute precision and recall from the vcf output
    # because we can do clipping here and it won't break normalization and
    # we can do a snp indel breakdown

    fn_table, fp_table, tp_table  = load_vcfeval_stats(out_dir)
    
    assert len(fn_table) == len(fp_table) and len(fp_table) == len(tp_table)
    rows = (0, len(fn_table)) if options.roc else (len(fn_table)-1, len(fn_table))
    ret = []
    for rownum in range(rows[0], rows[1]):
        fns = fn_table[rownum][1]
        fni = fn_table[rownum][2]
        fno = fn_table[rownum][3]
        fps = fp_table[rownum][1]
        fpi = fp_table[rownum][2]
        fpo = fp_table[rownum][3]
        tps = tp_table[rownum][1]
        tpi = tp_table[rownum][2]
        tpo = tp_table[rownum][3]

        fnt = fns + fni + fno
        fpt = fps + fpi + fpo
        tpt = tps + tpi + tpo

        precs = 0. if tps + fps == 0 else float(tps) / float(tps + fps)
        recs = 0. if tps + fns  == 0 else float(tps) / float(tps + fns)
        preci = 0. if tpi + fpi == 0 else float(tpi) / float(tpi + fpi)
        reci = 0. if tpi + fni  == 0 else float(tpi) / float(tpi + fni)
        prect = 0. if tpt + fpt == 0 else float(tpt) / float(tpt + fpt)
        rect = 0. if tpt + fnt  == 0 else float(tpt) / float(tpt + fnt)

        qual = fn_table[rownum][0]
        assert qual == fp_table[rownum][0]
        assert qual == tp_table[rownum][0]

        if options.clip is None and options.clip_fp is None and rownum == rows[1] - 1:
            # keep old way of computing as sanity check
            #assert abs(rect - total_recall) < 0.05 and abs(prect - total_precision) < 0.05
            pass
        
        # shoehorn into vcfCompre style output (todo, revise this) + qual
        ret += [[precs, recs, 0, 0, preci, reci, prect, rect, qual]]

    return ret
    
def make_mat(options, row_graphs, column_graphs, dist_fns):
    """ make a distance matix """
    mat = []
    for row in row_graphs:
        mat.append([])
        for col in column_graphs:
            for dist_fn in dist_fns:
                mat[-1].append(dist_fn(row, col))
    return mat

def n_avg(row):
    """ get average skipping Nones """
    if not type(row) is list:
        return row
    a = [x for x in row if x is not None]
    if len(a) == 0:
        return None
    else:
        return float(sum(a)) / len(a)

def make_tsvs(options):
    """ make some tsv files in the output dir
    """

    if options.comp_type == "kmer":
        dist_names = ["Jaccard-Dist", "Precision", "Recall"]
        dist_fns = [jaccard_dist_fn, precision_dist_fn, recall_dist_fn]
    elif options.comp_type == "corg":
        dist_names = ["Corg-Dist"]
        dist_fns = [corg_dist_fn]
    elif options.comp_type == "vcf":
        dist_names = vcf_dist_header(options)
        dist_fns = [vcf_dist_fn]
    elif options.comp_type == "sompy":
        dist_names = vcf_dist_header(options)
        dist_fns = [sompy_dist_fn]
    elif options.comp_type == "happy":
        dist_names = vcf_dist_header(options)
        dist_fns = [happy_dist_fn]
    elif options.comp_type == "vcfeval":
        dist_names = vcf_dist_header(options)
        dist_fns = [vcfeval_dist_fn]
    else:
        assert False
        
    # break apart by region
    for region in options.sample_graphs.keys():

        # do the baseline tsvs.  this is one row per graph,
        # with one column per comparison type with truth
        #for baseline in ["g1kvcf", "platvcf"]:
        for baseline in ["platvcf"]:
            for sample in options.sample_graphs[region].keys():
                header = dist_names
                RealTimeLogger.get().info("Making {} baseline tsv for {} {}".format(baseline, region, sample))
                mat = []
                row_labels = []
                for truth in options.baseline_graphs[region][sample]:
                    if options.tags[truth][2] == baseline:
                        for graph in options.sample_graphs[region][sample]:
                            rows = []
                            for d in dist_fns:
                                try:
                                    dist_res = d(graph, truth, options)
                                except Exception as e:
                                    RealTimeLogger.get().warning("Unable to retrieve distance between {} and {} because {}".format(graph, truth, e))
                                    dist_res = [[None] * len(header)]
                                for ri, r in enumerate(dist_res):
                                    if d == dist_fns[0]:
                                        rows.append(r)
                                        row_labels.append(options.tags[graph][2])
                                    else:
                                        rows[ri] += r
                            for row in rows:
                                mat.append(row)
                        break # shoud not be necessary
                # write the baseline matrix (with None for missing data) to file 
                tsv_path = raw_tsv_path(options, region, baseline, options.comp_type, sample)
                write_tsv(tsv_path, mat, header, row_labels, "Graph")

                # remove Nones and write tsv again
                clean_mat, clean_header, clean_row_labels = remove_nones(mat, header, row_labels)
                tsv_path = out_tsv_path(options, region, baseline, options.comp_type, sample)
                write_tsv(tsv_path, clean_mat, clean_header, clean_row_labels, "Graph")

        # sample vs sample heatmap
        if options.sample:
            for di, d in enumerate(dist_fns):
                RealTimeLogger.get().info("Making {} sample tsv for {}".format(dist_names[di], region))
                sample_vs_sample = defaultdict(lambda : defaultdict(list))
                for sample in options.sample_graphs[region].keys():
                    for graph1 in options.sample_graphs[region][sample]:
                        for graph2 in options.sample_graphs[region][sample]:
                            try:
                                dist_res = d(graph1, graph2, options)
                            except Exception as e:
                                RealTimeLogger.get().warning("Unable to retrieve distance between {} and {} because {}".format(graph1, graph2, e))
                                dist_res = [[None] * len(header)]
                                # only ever deal with 1 row here
                            assert len(dist_res) == 1
                            sample_vs_sample[options.tags[graph1][2]][options.tags[graph2][2]].append(dist_res[0][0])
                #average over samples
                for x in sample_vs_sample.keys():
                    for y in sample_vs_sample.keys():
                        sample_vs_sample[x][y] = n_avg(sample_vs_sample[x][y])
                # make matrix
                mat = []
                header = sample_vs_sample.keys()
                row_labels = sample_vs_sample.keys()
                for x in row_labels:
                    row = []
                    for y in header:
                        row.append(sample_vs_sample[x][y])
                    mat.append(row)
                # write the heatmap matrix (with None for missing data) to file 
                tsv_path = raw_tsv_path(options, region, "hm-sample", dist_names[di])
                write_tsv(tsv_path, mat, header, row_labels, "Graph")

                # remove Nones and write tsv again
                clean_mat, clean_header, clean_row_labels = remove_nones(mat, header, row_labels)
                tsv_path = out_tsv_path(options, region, "hm-sample", dist_names[di])
                write_tsv(tsv_path, clean_mat, clean_header, clean_row_labels, "Graph")

        # original vs original heatmap
        if options.orig:
            RealTimeLogger.get().info("Making {} orig tsv for {}".format(dist_names[di], region))
            for di, d in enumerate(dist_fns):
                orig_vs_orig = defaultdict(lambda : defaultdict(list))
                for graph1 in options.orig_graphs[region]:
                    for graph2 in options.orig_graphs[region]:
                        try:
                            dist_res = d(graph1, graph2, options)
                        except Exception as e:
                            RealTimeLogger.get().warning("Unable to retrieve distance between {} and {} because {}".format(graph1, graph2, e))
                            dist_res = [[None] * len(header)]
                            # only ever deal with 1 row here
                        assert len(dist_res) == 1
                        orig_vs_orig[options.tags[graph1][2]][options.tags[graph2][2]] = dist_res[0][0]
                # make matrix
                mat = []
                header = orig_vs_orig.keys()
                row_labels = orig_vs_orig.keys()
                for x in row_labels:
                    row = []
                    for y in header:
                        row.append(orig_vs_orig[x][y])
                    mat.append(row)
                # write the heatmap matrix (with None for missing data) to file 
                tsv_path = raw_tsv_path(options, region, "hm-orig", dist_names[di])
                write_tsv(tsv_path, mat, header, row_labels, "Graph")

                # remove Nones and write tsv again
                clean_mat, clean_header, clean_row_labels = remove_nones(mat, header, row_labels)
                tsv_path = out_tsv_path(options, region, "hm-orig", dist_names[di])
                write_tsv(tsv_path, clean_mat, clean_header, clean_row_labels, "Graph")

        # sample vs original heatmap
        if options.orig_and_sample:
            for di, d in enumerate(dist_fns):
                RealTimeLogger.get().info("Making {} sample tsv for {}".format(dist_names[di], region))
                sample_vs_orig = defaultdict(lambda : defaultdict(list))
                for sample in options.sample_graphs[region].keys():
                    for graph1 in options.sample_graphs[region][sample].union(options.orig_graphs[region]):
                        for graph2 in options.sample_graphs[region][sample].union(options.orig_graphs[region]):
                            try:
                                dist_res = d(graph1, graph2, options)
                            except Exception as e:
                                RealTimeLogger.get().warning("Unable to retrieve distance between {} and {} because {}".format(graph1, graph2, e))
                                dist_res = [[None] * len(header)]
                                # only ever deal with 1 row here
                            assert len(dist_res) == 1
                            name_1 = options.tags[graph1][2]
                            if graph1 in options.sample_graphs[region][sample]:
                                name_1 = "sample-" + name_1
                            else:
                                name_1 = "base-" + name_1
                            name_2 = options.tags[graph2][2]
                            if graph2 in options.sample_graphs[region][sample]:
                                name_2 = "sample-" + name_2
                            else:
                                name_2 = "base-" + name_2
                            sample_vs_orig[name_1][name_2].append(dist_res[0][0])
                #average over samples                            
                for x in sample_vs_orig.keys():
                    for y in sample_vs_orig.keys():
                        sample_vs_orig[x][y] = n_avg(sample_vs_orig[x][y])
                # make matrix
                mat = []
                header = sample_vs_orig.keys()
                row_labels = sample_vs_orig.keys()
                for x in row_labels:
                    row = []
                    for y in header:
                        row.append(sample_vs_orig[x][y])
                    mat.append(row)
                # write the heatmap matrix (with None for missing data) to file 
                tsv_path = raw_tsv_path(options, region, "hm-orig_and_sample", dist_names[di])
                write_tsv(tsv_path, mat, header, row_labels, "Graph")

                # remove Nones and write tsv again
                clean_mat, clean_header, clean_row_labels = remove_nones(mat, header, row_labels)
                tsv_path = out_tsv_path(options, region, "hm-orig_and_sample", dist_names[di])
                write_tsv(tsv_path, clean_mat, clean_header, clean_row_labels, "Graph")
                
        
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
        toks = line[:-1].split("\t")
        row_label = toks[0]
        col_names = toks[1:]
        row_names = []
        # body
        mat = []
        for line in f:
            toks = line[:-1].split("\t")
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
        if len(row_names) == 0:
            break
        assert len(mat) == len(row_names) 
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
    
    # Move to the appropriate working directory from wherever Toil dropped us
    os.chdir(options.cwd)
    
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
    
    # Move to the appropriate working directory from wherever Toil dropped us
    os.chdir(options.cwd)
    
    out_path = corg_path(graph1, graph2, options)
    corg_vg = corg_graph_path(graph1, graph2, options)
    do_comp = options.overwrite or not os.path.exists(out_path)
    if do_comp:
        if os.path.isfile(out_path):
            os.remove(out_path)
        robust_makedirs(os.path.dirname(out_path))
        run("corg {} {} --kmers-only -e {} -k {} -t {} > {} 2> {}".format(graph1, graph2, options.edge_max,
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
            corg_val = abs(1. - (2. * lenC) / float(len1 + len2)) 

        with open(out_path, "w") as f:
            f.write("{}\n".format(corg_val))

def preprocess_vcf(job, graph, options):
    """ run vt normalize and bed clip"""
    
    # Move to the appropriate working directory from wherever Toil dropped us
    os.chdir(options.cwd)
    
    input_vcf = input_vcf_path(graph, options)
    output_vcf = preprocessed_vcf_path(graph, options)
    robust_makedirs(os.path.dirname(output_vcf))

    run("scripts/vcfsort {} > {}".format(input_vcf, output_vcf), fail_hard=True)

    if options.qpct is not None and (options.tags[graph][2] in ["gatk3", "platypus", "freebayes", "samtools"] or
                                     (options.tags[graph][2] == "g1kvcf" and options.baseline != "g1kvcf") or
                                     options.qgraph is True):
        # g1kvcf has no quality info.  proxy with read depth to at least get a curve
        #filter_opts = "--info DP" if options.tags[graph][2] == "g1kvcf" else ""
        if "platvcf" not in options.tags[graph][2]:
            filter_opts = ""
            vcfFQcmd = "cat {} | ".format(output_vcf)
            if options.tags[graph][2] not in ["gatk3", "platypus", "freebayes", "samtools", "g1kvcf", "platvcf", "platvcf-baseline"]:
                #filter_opts += " --info DP"
                filter_opts += " --{}".format(options.filter_type) if not options.new else " --ad"
                if options.dedupe:
                    filter_opts += " --dedupe"

            vcfFQcmd += "scripts/vcfFilterQuality.py - {} --pct {} --set_qual > {}".format(options.qpct,
                                                                                           filter_opts,
                                                                                           output_vcf + ".qpct")
            run(vcfFQcmd)
            run("cp {} {}".format(output_vcf + ".qpct", output_vcf))
    

    if options.normalize is True:# and options.tags[graph][2] not in ["gatk3", "platypus", "freebayes", "samtools", "g1kvcf", "platvcf", "platvcf-baseline"]:
        # run blocksub (only on single allelic)
        sts = run("bcftools filter -e \"N_ALT > 1\" {} | vt decompose_blocksub -a - > {}".format( 
            output_vcf,
            output_vcf + ".vt"))
        # add back multiallelic variants
        run("bcftools filter -i \"N_ALT > 1\" {} | bcftools view - -H >> {}".format(
            output_vcf,
            output_vcf + ".vt"))
        if sts != 0:
            run("rm {}".format(output_vcf))
            return
        run("cp {} {}".format(output_vcf + ".vt", output_vcf))

    if options.clip is not None and options.comp_type != "vcfeval":
        clip_bed = clip_bed_path(graph, options)        
        if not os.path.isfile(clip_bed):
            RealTimeLogger.get().warning("Clip bed file not found {}".format(clip_bed))
        else:
            run("bgzip {}".format(output_vcf), fail_hard=True)
            run("tabix -f -p vcf {}".format(output_vcf + ".gz"), fail_hard=True)
            run("bcftools view {} -R {} > {}".format(output_vcf + ".gz", clip_bed, output_vcf), fail_hard=True)
            run("rm {}".format(output_vcf + ".gz*"))
            run("scripts/vcfsort {} | vcfuniq > {}".format(output_vcf, output_vcf + ".sort"), fail_hard=True)
            run("cp {} {}".format(output_vcf + ".sort", output_vcf), fail_hard=True)

    # one final sort, and strip ignored variants
    ig_opts = ""
    for ignore_keyword in options.ignore:
        ig_opts += " | grep -v {}".format(ignore_keyword)
    # also strip genotypes
    if not options.gt and options.comp_type != "vcfeval":
        ig_opts += " | scripts/vcfSetGenotypes.py -"
    run("mv {} {} ; scripts/vcfsort {} {} | vcfuniq > {}".format(output_vcf, output_vcf + ".ig",
                                                                 output_vcf + ".ig", ig_opts,
                                                                 output_vcf), fail_hard=True)

    # need compressed index for vcfeval
    run("bgzip {} -c > {}".format(output_vcf, output_vcf + ".gz"), fail_hard=True)
    run("tabix -f -p vcf {}".format(output_vcf + ".gz"), fail_hard=True)


def compute_vcf_comparison(job, graph1, graph2, options):
    """ run vcf compare between two graphs
    """
    
    # Move to the appropriate working directory from wherever Toil dropped us
    os.chdir(options.cwd)
    
    if options.comp_type == "sompy":
        out_path = comp_path_sompy(graph1, graph2, options)
    elif options.comp_type == "happy":
        out_path = comp_path_happy(graph1, graph2, options)
    elif options.comp_type == "vcfeval":
        out_path = comp_path_vcfeval(graph1, graph2, options)
    else:
        out_path = comp_path_vcf(graph1, graph2, options)

    # we expect graph1 to be a sample graph
    region1, sample1, method1  = options.tags[graph1]
    assert method1 != "None" and "platvcf" not in method1

    # we expect graph2 to be a baseline graph
    region2, sample2, method2 = options.tags[graph2]
    assert method2 in ["g1kvcf", "platvcf"]
    assert region1 == region2
    assert sample1  == sample2
    
    # get the vcf of the sample graph
    query_vcf_path = preprocessed_vcf_path(graph1, options)

    # and the baseline
    truth_vcf_path = preprocessed_vcf_path(graph2, options)

    do_comp = options.overwrite or not os.path.exists(out_path)

    if do_comp:
        if os.path.isfile(out_path):
            os.remove(out_path)
        robust_makedirs(os.path.dirname(out_path))

        if options.comp_type == "sompy":
            sp_opts = "-P"
            if options.clip_fp:
                sp_opts += " -f {}".format(options.clip_fp)
            # doesn't seem to work
            #if options.clip:
            #    sp_opts += " -R {}".format(options.clip)
                
            run("export HGREF={} ; som.py {} {} --output {} {} 2> {}".format(options.chrom_fa_path, truth_vcf_path,
                                                                             query_vcf_path, out_path.replace(".stats.csv", ""),
                                                                             sp_opts, out_path + ".stderr"))
        elif options.comp_type == "happy":
            hp_opts = ""
            # cant figure out how to get vcfeval working in happy.  leave off for now
            #hp_opts += "--engine vcfeval --engine-vcfeval-path rtg --engine-vcfeval-template {}".format(options.chrom_sdf_path)
            if options.clip_fp:
                hp_opts += " -f {}".format(options.clip_fp)
            # doesn't seem to work
            #if options.clip:
            #    hp_opts += " -R {}".format(options.clip)
            
            # make roc curves for gatk and platypus (hardcoding name check as hack for now)
            if method1 in ["gatk3", "platypus", "g1kvcf", "freebayes", "samtools"] and options.roc is True:
                hp_opts += " -V --roc Q_GQ --roc-filter LowGQX"
            # since we use just numbers for chrom names in the vcf, use options.happy_fa_path as a hack to make happy happy.
            run("export HGREF={} ; hap.py {} {} -o {} --threads {} {} 2> {}".format(options.happy_fa_path, truth_vcf_path,
                                                                                    query_vcf_path, out_path.replace(".summary.csv", ""),
                                                                                    options.vg_cores, hp_opts, out_path + ".stderr"))
        elif options.comp_type == "vcf":
            vc_opts = ""
            for ignore_keyword in options.ignore:
                vc_opts += " -i {}".format(ignore_keyword)
            run("scripts/vcfCompare.py {} {} {} > {}".format(query_vcf_path, truth_vcf_path, vc_opts, out_path))
            
        elif options.comp_type == "vcfeval":
            ve_opts = "" if options.gt else "--squash-ploidy"
            # doesn't seem to work
            #if options.clip:
            #    ve_opts += " --bed-regions={}".format(options.clip)
            # indexing and compression was done by preprocessing phase
            run("rm -rf {}".format(out_path))
            run("rtg vcfeval -b {}.gz -c {}.gz --all-records --ref-overlap --vcf-score-field QUAL -t {} {} -o {}".format(truth_vcf_path, query_vcf_path,
                                                                                                                         options.chrom_sdf_path, ve_opts,
                                                                                                                         out_path))
            # count up output
            fn_path = os.path.join(out_path, "fn.vcf.gz")
            fp_path = os.path.join(out_path, "fp.vcf.gz")
            tp_path = os.path.join(out_path, "tp.vcf.gz")

            try:
                fn_table = vcf_qual_stats(fn_path, options.clip, ignore_keywords = ["OverlapConflict"])
                fp_table = vcf_qual_stats(fp_path, options.clip_fp if options.clip_fp else options.clip)
                tp_table = vcf_qual_stats(tp_path, options.clip)
                save_vcfeval_stats(out_path, fn_table, fp_table, tp_table)
            except:
                pass
            
        
def compute_kmer_comparisons(job, options):
    """ run vg compare in parallel on all the graphs,
    outputting a json file for each
    """
    
    # Move to the appropriate working directory from wherever Toil dropped us
    os.chdir(options.cwd)
    
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
    
    # Move to the appropriate working directory from wherever Toil dropped us
    os.chdir(options.cwd)
    
    RealTimeLogger.get().info("Running corg comparison {} pairs of input graphs".format(
        len(options.pair_comps)))
    for pair_comp in options.pair_comps:
        graph1, graph2 = pair_comp[0], pair_comp[1]
        out_path = corg_path(graph1, graph2, options)
        if options.overwrite or not os.path.exists(out_path):
            job.addChildJobFn(compute_corg_comparison, graph1, graph2, options,
                              cores=options.vg_cores)

def compute_vcf_comparisons(job, options):
    """ run vg compare in parallel on all the graphs,
    outputting a json file for each
    """
    
    # Move to the appropriate working directory from wherever Toil dropped us
    os.chdir(options.cwd)
    
    RealTimeLogger.get().info("Running vcf comparison {} pairs of input graphs".format(
        len(options.pair_comps)))
    for pair_comp in options.pair_comps:
        graph1, graph2 = pair_comp[0], pair_comp[1]
        out_path = comp_path_vcf(graph1, graph2, options)
        if options.overwrite or not os.path.exists(out_path):
            cores = options.vg_cores if options.comp_type == "happy" else 1
            job.addChildJobFn(compute_vcf_comparison, graph1, graph2, options,
                                          cores=cores)

def compute_kmer_indexes(job, options):
    """ run everything (root toil job)
    first all indexes are computed,
    then all comparisons (follow on)
    then summary (follow on of that)
    """
    
    # Move to the appropriate working directory from wherever Toil dropped us
    os.chdir(options.cwd)
    
    # do all the indexes
    input_set = set()
    for pair_comp in options.pair_comps:
        input_set.add(pair_comp[0])
        input_set.add(pair_comp[1])

    if options.comp_type in ["kmer", "corg"]:
        RealTimeLogger.get().info("Computing indexes for {} input graphs".format(len(input_set)))
        for graph in input_set:
            if options.overwrite or not os.path.exists(index_path(graph, options)):
                job.addChildJobFn(compute_kmer_index, graph, options, cores=options.vg_cores)

    if options.comp_type in ["vcf", "sompy", "happy", "vcfeval"]:
        RealTimeLogger.get().info("Preprocessing {} input vcfs".format(len(input_set)))
        for graph in input_set:
            if options.overwrite or not os.path.isfile(preprocessed_vcf_path(graph, options)):
                job.addChildJobFn(preprocess_vcf, graph, options, cores=1)

    # do the comparisons
    if options.comp_type == "kmer":
        job.addFollowOnJobFn(compute_kmer_comparisons, options, cores=1)
    elif options.comp_type == "corg":
        job.addFollowOnJobFn(compute_corg_comparisons, options, cores=1)
    elif options.comp_type in ["vcf", "sompy", "happy", "vcfeval"]:
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

    # hack in dummy files for total vcfs
    total_gams = set()
    if options.comp_type in ["vcf", "sompy", "happy", "vcfeval"]:
        for input_gam in in_gams:
            region = alignment_region_tag(input_gam, options)
            # dangerous :: fix
            dummy_gam = input_gam.replace(region, "total")
            total_gams.add(dummy_gam)

    for input_gam in in_gams + list(total_gams):
        region = alignment_region_tag(input_gam, options)
        sample = alignment_sample_tag(input_gam, options)
        method = alignment_graph_tag(input_gam, options)

        orig_path = graph_path(input_gam, options)
        sample_path = sample_vg_path(input_gam, options)
        g1kvcf_path = g1k_vg_path(input_gam, False, False, options)
        platvcf_path = g1k_vg_path(input_gam, True, False, options)
        
        if orig or orig_and_sample:
            orig_graphs[region].add(orig_path)
            tags[orig_path] = (region, None, method)

        def test_path(graph, method):
            if options.comp_type in ["vcf", "sompy", "happy", "vcfeval"]:
                return input_vcf_path(graph, options, region, sample, method)
            return graph

        if os.path.isfile(test_path(sample_path, method)):
            sample_graphs[region][sample].add(sample_path)
        else:
            sys.stderr.write("WARNING, input VCF not found: {}\n".format(
                test_path(sample_path, method)))
        # we dont expect to have baselines for every sample
        if options.baseline == "g1kvcf" and os.path.isfile(test_path(g1kvcf_path, "g1kvcf")):
            baseline_graphs[region][sample].add(g1kvcf_path)
        if options.baseline == "platvcf" and os.path.isfile(test_path(platvcf_path, "platvcf")):
            baseline_graphs[region][sample].add(platvcf_path)

        # add external vcfs as sample graphs, rely on tags to tell them apart
        if options.comp_type in ["vcf", "sompy", "happy", "vcfeval"]:
            gatk3_path = input_vcf_path(None, options, region, sample, "gatk3")
            if os.path.isfile(gatk3_path):
                sample_graphs[region][sample].add(gatk3_path)
                tags[gatk3_path] = (region, sample, "gatk3")
            platypus_path = input_vcf_path(None, options, region, sample, "platypus")
            if os.path.isfile(platypus_path):
                sample_graphs[region][sample].add(platypus_path)
                tags[platypus_path] = (region, sample, "platypus")
            freebayes_path = input_vcf_path(None, options, region, sample, "freebayes")
            if os.path.isfile(freebayes_path):
                sample_graphs[region][sample].add(freebayes_path)
                tags[freebayes_path] = (region, sample, "freebayes")
            samtools_path = input_vcf_path(None, options, region, sample, "samtools")
            if os.path.isfile(samtools_path):
                sample_graphs[region][sample].add(samtools_path)
                tags[samtools_path] = (region, sample, "samtools")
            #if options.baseline != "g1kvcf" and os.path.isfile(test_path(g1kvcf_path, "g1kvcf")):                
            #    sample_graphs[region][sample].add(g1kvcf_path)

        tags[sample_path] = (region, sample, method)
        tags[g1kvcf_path] = (region, sample, "g1kvcf")
        tags[platvcf_path] = (region, sample, "platvcf")

    return orig_graphs, sample_graphs, baseline_graphs, tags
    
def main(args):
    
    options = parse_args(args)

    assert options.comp_type in ["corg", "kmer", "vcf", "sompy", "happy", "vcfeval"]
        
    if options.comp_type in ["vcf", "sompy", "happy", "vcfeval"]:
        assert not options.orig and not options.orig_and_sample

    if options.orig_and_sample:
        options.orig = True
        options.sample = True

    if options.comp_type in ["vcfeval"]:
        options.chrom_sdf_path = os.path.join(options.comp_dir, "chrom.sdf")
        if options.overwrite or not os.path.exists(options.chrom_sdf_path):
            # need to reformat the fasta as sdf for vcfeval
            run("rm -rf {}; rtg format {} -o {}".format(
                options.chrom_sdf_path, options.chrom_fa_path, options.chrom_sdf_path))
    
    RealTimeLogger.start_master()

    # since we re-use callVariants methods
    options.out_dir = options.var_dir

    # find all the graphs
    # clean out shell artifacts
    in_gams = []
    for gam in options.in_gams:
        if os.path.isfile(gam):
            in_gams.append(gam)
        else:
            RealTimeLogger.get().warning("skipping {}".format(gam))
    options.in_gams = in_gams
    
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
            if options.orig_and_sample:
                for sample in options.sample_graphs[region].keys():
                    for graph2 in options.sample_graphs[region][sample]:
                        options.pair_comps.append((graph1, graph2))
                

    # now all sample vs baseline comparisons
    for region in options.sample_graphs.keys():
        for sample in options.sample_graphs[region].keys():
            for graph1 in options.sample_graphs[region][sample]:
                for graph2 in options.baseline_graphs[region][sample]:
                    # skip g1kvcf vs g1kvcf comparison
                    if options.tags[graph1] != options.tags[graph2]:
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
        
        
        
        
        
        
        
        
        
        

