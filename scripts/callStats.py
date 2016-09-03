#!/usr/bin/env python2.7
"""
Get some basic statistics from the vg call output, (via computeVariantsDistances.py vcfeval)
"""

import argparse, sys, os, os.path, random, subprocess, shutil, itertools, glob
import doctest, re, json, collections, time, timeit, string, tempfile, copy
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
    parser.add_argument("comp_dir", type=str,
                        help="output directory of computeVariantsDistances.py")
    parser.add_argument("out_dir", type=str,
                        help="directory where output of this script will go")
    parser.add_argument("--clip", type=str, default=None,
                        help="bed regions to subset on")
                            
    args = args[1:]
        
    return parser.parse_args(args)

def munge_vcfeval_results(comp_dir):
    """ make this map [REGION][SAMPLE][GRAPH] -> [VCFEVAL OUTPUT DIR] """
    evalmap = dict()
    for regiondir in glob.glob(os.path.join(comp_dir, "vcfeval_compare_data", "*")):
        if os.path.isdir(regiondir):
            region = os.path.basename(regiondir)
            for vcfevaldir in glob.glob(os.path.join(regiondir, "*")):
                if os.path.isdir(vcfevaldir):
                    vcfevalname = os.path.basename(vcfevaldir)
                    toks = vcfevalname.split("_")
                    # expect graph_sample_vs_graph_sample
                    vsi = toks.index("vs")
                    assert vsi > 1
                    sample = toks[vsi-1]
                    graph = "".join(toks[0:vsi-1])
                    if region not in evalmap:
                        evalmap[region] = dict()
                    if sample not in evalmap[region]:
                        evalmap[region][sample] = dict()
                    assert graph not in evalmap[region][sample]
                    evalmap[region][sample][graph] = vcfevaldir
    return evalmap

def count_variants(vcf_path, filter_string, xref, kind):
    """ use vcftools to count up lines in a vcf that meet criteria.
    kind in [indels, snps, all]
    """
    vstr = "-v snps,mnps" if kind is "snps" else "-V snps,mnps" if kind is "indels" else ""
    if options.clip is not None:
        vstr += " -R {}".format(options.clip)
    xstr = "grep XREF" if xref is True else "grep -v XREF"
    cmd = "bcftools view {} -H -f {} {} | {} | wc -l".format(vcf_path, filter_string, vstr, xstr)
    proc = subprocess.Popen(cmd, shell=True, bufsize=-1,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, errors = proc.communicate()
    sts = proc.wait()
    assert sts == 0
    return int(output)

    
def get_counts(vcfeval_dir, options):
    """ count up the true and false positives, break down by type and xrefness"""
    tpvcf = os.path.join(vcfeval_dir, "tp.vcf.gz")
    fpvcf = os.path.join(vcfeval_dir, "fp.vcf.gz")
    counts = dict()
    counts["TP-REF-SNP"] = count_variants(tpvcf, "PASS,.", True, "snps")
    counts["TP-AUG-SNP"] = count_variants(tpvcf, "PASS,.", False, "snps")
    counts["TP-REF-INDEL"] = count_variants(tpvcf, "PASS,.", True, "indels")
    counts["TP-AUG-INDEL"] = count_variants(tpvcf, "PASS,.", False, "indels")
    counts["TP-REF-TOTAL"] = counts["TP-REF-SNP"] + counts["TP-REF-INDEL"]
    counts["TP-AUG-TOTAL"] = counts["TP-AUG-SNP"] + counts["TP-AUG-INDEL"]

    counts["FP-REF-SNP"] = count_variants(fpvcf, "PASS,.", True, "snps")
    counts["FP-AUG-SNP"] = count_variants(fpvcf, "PASS,.", False, "snps")
    counts["FP-REF-INDEL"] = count_variants(fpvcf, "PASS,.", True, "indels")
    counts["FP-AUG-INDEL"] = count_variants(fpvcf, "PASS,.", False, "indels")
    counts["FP-REF-TOTAL"] = counts["FP-REF-SNP"] + counts["FP-REF-INDEL"]
    counts["FP-AUG-TOTAL"] = counts["FP-AUG-SNP"] + counts["FP-AUG-INDEL"]

    return counts
                        
def do_all_counts(evalmap, options):
    """ count up all our tp and fp stats and return in table """
    # [region][sample][graph] -> counts
    count_table = copy.deepcopy(evalmap)
    for region, rd in evalmap.items():
        for sample, sd in rd.items():
            for graph, evaldir in sd.items():
                count_table[region][sample][graph] = get_counts(evaldir, options)

    return count_table
    
def counts_tsv(graph_table, options):
    """ make a tsv for a given table (corresponding to sample and graph) """
    if len(graph_table) == 0:
        return
    table = []
    header = None
    for graph, count_table in graph_table.items():
        if len(table) == 0:
            header = [h for h in count_table.keys()]
            header = ["graph"] + sorted(header)
        row = [graph]
        for category in header[1:]:
            row.append(count_table[category])
        table.append(row)

    tsv = "#" + "\t".join(header) + "\n"
    for row in table:
        tsv_row = [str(x) for x in row]
        tsv_line = "\t".join(tsv_row)
        tsv += tsv_line + "\n"
    return tsv
            
    
def main(args):
                                        
    options = parse_args(args)

    robust_makedirs(options.out_dir)
                                    
    evalmap = munge_vcfeval_results(options.comp_dir)
    counts_table = do_all_counts(evalmap, options)
    
    for region, rd in counts_table.items():
        for sample, graph_table in rd.items():
            tsv = counts_tsv(graph_table, options)
            with open(os.path.join(options.out_dir, "call_stats_{}_{}.tsv".format(region, sample)), "w") as f:
                f.write(tsv)
    
if __name__ == "__main__" :
    sys.exit(main(sys.argv))
        
        
        
        
        
        
        
        
        
        

