#!/usr/bin/env python2.7
"""
Get some basic statistics from the vg call output, (via computeVariantsDistances.py vcfeval)
"""

import argparse, sys, os, os.path, random, subprocess, shutil, itertools, glob
import doctest, re, json, collections, time, timeit, string, tempfile, copy
from collections import defaultdict
from operator import sub
from toillib import RealTimeLogger, robust_makedirs
from callVariants import sample_vg_path, sample_txt_path, augmented_vg_path, alignment_region_tag, alignment_graph_tag
from callVariants import graph_path, index_path, augmented_vg_path, linear_vg_path, linear_vcf_path, sample_vg_path
from callVariants import alignment_region_tag, alignment_sample_tag, alignment_graph_tag, run
from plotVariantsDistances import name_map

def parse_args(args):
    parser = argparse.ArgumentParser(description=__doc__, 
        formatter_class=argparse.RawDescriptionHelpFormatter)
        
    # General options
    parser.add_argument("comp_dir", type=str,
                        help="output directory of computeVariantsDistances.py")
    parser.add_argument("out_dir", type=str,
                        help="directory where output of this script will go")
    parser.add_argument("--ped", type=str,
                        default="0\tNA12879\tNA12877\tNA12878\t2\t0",
                        help="PED format pedigree")
    parser.add_argument("--chrom_fa_path", type=str, default="data/g1kvcf/chrom.fa",
                        help="fasta file with entire chromosome info for all regions")
                            
    args = args[1:]
        
    return parser.parse_args(args)

def munge_vcf_results(comp_dir):
    """ make this map [REGION][SAMPLE][GRAPH] -> [PREPROCESSED VCF] """
    vcfmap = dict()
    for regiondir in glob.glob(os.path.join(comp_dir, "preprocessed_vcfs", "*")):
        if os.path.isdir(regiondir):
            region = os.path.basename(regiondir)
            for pvcf in glob.glob(os.path.join(regiondir, "*.vcf")):
                vcfname = os.path.basename(pvcf)
                toks = vcfname.split("_")
                #expect sample_graph
                sample = toks[0]
                graph = os.path.splitext("".join(toks[1:]))[0]
                if region not in vcfmap:
                    vcfmap[region] = dict()
                if sample not in vcfmap[region]:
                    vcfmap[region][sample] = dict()
                assert graph not in vcfmap[region][sample]
                vcfmap[region][sample][graph] = pvcf
    return vcfmap

def make_trio_vcfs(vcfmap, options):
    """ merge up samples into same vcf using rtg return index of merged files"""
    robust_makedirs(options.out_dir)
    ped_file = os.path.join(options.out_dir, "predigree.ped")
    with open(ped_file, "w") as f:
        f.write(options.ped + "\n")

    mergetable = dict()
    
    for region, rd in vcfmap.items():
        mergetable[region] = dict()
        region_dir = os.path.join(options.out_dir, "trio_vcfs", region)
        robust_makedirs(region_dir)
        # round up all sampels for graph
        bygraph = dict()
        for sample, sd in rd.items():
            for graph, pvcf in sd.items():
                if graph not in bygraph:
                    bygraph[graph] = dict()
                bygraph[graph][sample] = pvcf

        # make a merged vcf for each graph
        for graph, sd in bygraph.items():
            input_vcfs = { "snp" : [], "indel" : [], "all" : [] }
            for sample, pvcf in sd.items():
                work_dir = os.path.join(region_dir, "input_vcf")
                merge_dir = os.path.join(region_dir, "merged_vcf")
                robust_makedirs(work_dir)
                robust_makedirs(merge_dir)
                for kind in input_vcfs.keys():
                    filter_vcf = os.path.join(work_dir, "{}_{}_{}.vcf".format(graph, sample, kind))
                    vstr = "-v snps,mnps" if kind is "snp" else "-V snps,mnps" if kind is "indel" else ""
                    run("bcftools view {} -f PASS,. {} | bcftools norm - -f {} > {}".format(
                        pvcf, vstr, options.chrom_fa_path, filter_vcf))
                    run("bgzip -f {}".format(filter_vcf))
                    run("tabix -f -p vcf {}.gz".format(filter_vcf))
                    input_vcfs[kind].append("{}.gz".format(filter_vcf))

            if len(sd.items()) >= 3 and \
               len(input_vcfs["all"]) == len(sd.items()) and\
               len(input_vcfs["snp"]) == len(sd.items()) and\
               len(input_vcfs["indel"]) == len(sd.items()):

                mergetable[region][graph] = dict()
                for kind in input_vcfs.keys():

                    out_vcf = os.path.join(merge_dir, "{}_{}_merged.vcf.gz".format(graph, kind))
                    run("rm -f {}".format(out_vcf))
                    run("rtg vcfmerge {} -o {}".format(" ".join(input_vcfs[kind]), out_vcf), fail_hard = True)
                
                    mergetable[region][graph][kind] = out_vcf

    return mergetable

def scrape_mendel(mendel_out):
    """ get the concordance from a file """
    with open(mendel_out) as f:
        for line in f:
            toks = line.split()
            if toks[0] == "Concordance" and len(toks) == 8:
                concordance = toks[-1]
                return 0.01 * float(concordance.strip("()%"))
    return None

def do_mendel(mergetable, options):
    """ run rtg mendelian on all our merged vcfs """

    header = ["graph", "all", "snp", "indel"]
    for region, gd in mergetable.items():
        table = []
        for graph, mergefiles in gd.items():
            annot_dir = os.path.join(options.out_dir, "mendel", region, graph)
            robust_makedirs(annot_dir)
            concordance = dict()
            for kind, mergefile in mergefiles.items():
                out_vcf = os.path.join(annot_dir, "mendel_{}.vcf.gz".format(kind))
                con_vcf = os.path.join(annot_dir, "consistent_{}.vcf.gz".format(kind))
                incon_vcf = os.path.join(annot_dir, "inconsistent_{}.vcf.gz".format(kind))
                out_stdout = os.path.join(annot_dir, "mendel_{}.stdout".format(kind))

                run("rtg mendelian -l -i {} -t {} --pedigree {} --output {} --output-consistent {} --output-inconsistent {} > {}".format(
                    mergefile,
                    os.path.join(options.comp_dir, "chrom.sdf"),
                    os.path.join(options.out_dir, "predigree.ped"),
                    out_vcf,
                    con_vcf,
                    incon_vcf,
                    out_stdout))

                concordance[kind] = scrape_mendel(out_stdout)

            table.append([graph, concordance["all"], concordance["snp"], concordance["indel"]])
            
        # write the tsv for this region
        with open(os.path.join(options.out_dir, "mendel-{}.tsv".format(region)), "w") as f:
            f.write("\t".join(header) + "\n")
            for row in table:
                if None not in row:
                    line = [str(s) for s in row]
                    f.write("\t".join(line) + "\n")
        
def main(args):
                                        
    options = parse_args(args)

    RealTimeLogger.start_master()

    robust_makedirs(options.out_dir)
                                    
    vcfmap = munge_vcf_results(options.comp_dir)

    mergetable = make_trio_vcfs(vcfmap, options)

    do_mendel(mergetable, options)
    
    
if __name__ == "__main__" :
    sys.exit(main(sys.argv))
        
        
        
        
        
        
        
        
        
        

