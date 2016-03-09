#!/usr/bin/env python2.7
"""
Get some basic statistics from the vg call output
"""

import argparse, sys, os, os.path, random, subprocess, shutil, itertools, glob
import doctest, re, json, collections, time, timeit, string
from collections import defaultdict
from operator import sub
from toillib import robust_makedirs
from callVariants import sample_vg_path, sample_txt_path, augmented_vg_path, alignment_region_tag, alignment_graph_tag
from callVariants import graph_path, index_path, augmented_vg_path, linear_vg_path, linear_vcf_path, sample_vg_path
from callVariants import alignment_region_tag, alignment_sample_tag, alignment_graph_tag

def parse_args(args):
    parser = argparse.ArgumentParser(description=__doc__, 
        formatter_class=argparse.RawDescriptionHelpFormatter)
        
    # General options
    parser.add_argument("in_gams", nargs="+",
                        help="input alignment files. (must have been run through callVariants.py!)")
    parser.add_argument("output_path", type=str,
                        help="directory to write output to.  not to be confused with --out_dir"
                        "which is the output directory used for callVariants.py")
    parser.add_argument("--out_dir", type=str, default="variants",
                        help="output directory (used for callVariants.py)")
    parser.add_argument("--graph_dir", type=str, default="graphs",
                        help="name of input graphs directory")
    parser.add_argument("--avg_samples", action="store_true", default=False,
                        help="Average samples into mean value")
    parser.add_argument("--tag", type=str, default = "",
                        help="add tag to all output files")
    parser.add_argument("--chrom_fa_path", type=str, default="data/g1kvcf/chrom.fa",
                        help="fasta file with entire chromosome info for all regions")
    
                            
    args = args[1:]
        
    return parser.parse_args(args)

def compare_out_path(options):
    """ get root output dir for comparison output
    """
    return options.output_path

def tsv_out_path(options, name):
    """ path for output tsv"""
    tag = "_{}".format(options.tag) if len(options.tag) > 0 else ""
    return os.path.join(compare_out_path(options),
                        "{}{}.tsv".format(name, tag))
        
def count_tsv_path(options):
    return tsv_out_path(options, "call_count")

def size_tsv_path(options):
    return tsv_out_path(options, "call_size")

def detailed_call_coutns_tsv_path(options):
    return tsv_out_path(options, "call_details")
               
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
    cmd = "scripts/vcfCountSnps.sh {} {}".format(vcf, options.chrom_fa_path)
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                         stderr=sys.stderr, bufsize=-1)
    output, _ = p.communicate()
    assert p.wait() == 0
    num_paths = int(output.strip())
    return num_paths

def normalize_vcf(vcf, options):
    """ run normalizeation pipeline
    """
    cmd = "vcfsort {} | vt decompose - | vt decompose_blocksub -a - | vt normalize -r {} - | uniq".format(vcf, options.chrom_fa_path)
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE, bufsize=-1)
    output, _ = p.communicate()
    assert p.wait() == 0
    return output

def count_aug_snps(sample_vcf, sample_txt, options):
    """ how many snps in the vcf are in the augmented graph """
    if not os.path.isfile(sample_txt) or not os.path.isfile(sample_vcf):
        return -1
    augs = set()
    with open(sample_txt) as f:
        for line in f:
            if line[0] != "#":
                toks = line.split()
                node = int(toks[0])
                pos = int(toks[1]) - 1
                if toks[4] == "SNP":
                    augs.add((node, pos))
    aug_count = 0
    f = normalize_vcf(sample_vcf, options).split("\n")
    found_header = False
    for line in f:
        if len(line) == 0:
            continue
        if line[0] == "#":
            found_header = True
        elif found_header is True:
            toks = line.split()
            ref = toks[3]
            alts = toks[4].split(",")
            # assume normalized
            assert len(alts) == 1
            # check if snp
            if not (len(ref) == 1 and len(alts[0]) == 1 and alts[0] != ref):
                continue
            vid = toks[2]
            if "." in vid:
                pos = int(vid.split(".")[-1])
            else:
                pos = 0
            vid = vid.split(".")[0]
            if "_" in vid:
                first_node, last_node = vid.split("_")[0], vid.split("_")[0][-1]
            else:
                first_node = vid
                last_node = vid
            first_node, last_node = int(first_node), int(last_node)
            is_aug = False
            for node in range(first_node, last_node+1):
                if (node, pos) in augs:
                    is_aug = True
            if is_aug is True:
                aug_count += 1
    return aug_count

def vg_detailed_call_counts(gam, options):
    """ return tuple (vcf-snp, aug-snp, ref-snp, non-calls, ref/ref, ref/alt1, alt1/alt1, alt1/alt2)
    """
    original_graph = graph_path(gam, options)    
    sample_graph = sample_vg_path(gam, options)
    sample_txt = sample_txt_path(gam, options)
    sample_vcf = sample_txt.replace(".txt", ".vcf")

    non_calls = 0
    ref_calls = 0
    alt_calls = 0
    ref_ref = 0
    ref_alt = 0
    alt1_alt1 = 0
    alt1_alt2 = 0
    with open(sample_txt) as f:
        line = f.readline()
        while line:
            toks = line.split()
            vals = toks[3].split(",")
            assert len(vals) == 2
            if vals == ['-', '-']:
                pass # non-call computed above from vg files
            elif vals == ['-', '.'] or vals == ['.','-'] or vals == ['.','.']:
                ref_ref += 1 # vg call doesn't distinguish .- and .. in general atm
            elif vals[0] == '.' or vals[1] == '.':
                ref_alt += 1
            elif vals[0] == vals[1] or vals[0] == '-' or vals[1] == '-':
                alt1_alt1 += 1 # vg call doesn't distinguish A- and AA atm
            else:
                assert vals[0] != vals[1]
                assert vals[0] != toks[2] and vals[1] != toks[2]
                alt1_alt2 += 1
            line = f.readline()

    # check normalized output of glenn2vcf (number ofalts)
    sample_vcf_snps = count_vcf_snps(sample_vcf, options)

    # total alt_calls
    alt_snps = count_aug_snps(sample_vcf, sample_txt, options)
    
    # todo: use more verybose calls.txt for a lot of this
    non_calls = vg_length(original_graph, options) - vg_length(sample_graph, options) + ref_alt + alt1_alt1 + 2 * alt1_alt2
    
    
    return (sample_vcf_snps, alt_snps, sample_vcf_snps - alt_snps,
            non_calls,
            ref_ref * 2 + ref_alt,
            (alt1_alt1 + alt1_alt2) * 2 + ref_alt,
            ref_ref, ref_alt, alt1_alt1, alt1_alt2)

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
        # ugly hack to get g1kvcf / platvcf sample from vcf
        if vg_sample.split("/")[-2] == "g1kvcf" or vg_sample.split("/")[-2] == "platvcf":
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

def detailed_call_count_table(options):
    """ make a table of the detailed call count. this makes count and size tables
    obselete
    """
    # tsv header
    detailed_table = "#graph\ttot-snps\taug-snps\tref-snps\tnon-calls\tref-calls\talt-calls\tref-ref\tref-alt1\talt1-alt1\talt1-alt2\n"

    sums = defaultdict(lambda : (0,0,0,0,0,0,0,0,0,0))
    counts = defaultdict(lambda : 0)

    for gam in options.in_gams:
        result = vg_detailed_call_counts(gam, options)

        if options.avg_samples:
            # reduce to reference graph
            name_fn = lambda x : graph_path(x, options)
        else:
            # keep split apart by sample / graph corresponding to gam
            name_fn = lambda x : alignment_graph_tag(x, options)  \
                   + "_" + alignment_sample_tag(x, options)
            
        name = name_fn(gam)

        sums[name] = (sums[name][0] + result[0],
                      sums[name][1] + result[1],
                      sums[name][2] + result[2],
                      sums[name][3] + result[3],
                      sums[name][4] + result[4],
                      sums[name][5] + result[5],
                      sums[name][6] + result[6],
                      sums[name][7] + result[7],
                      sums[name][8] + result[8],
                      sums[name][9] + result[9])


        counts[name] = counts[name] + 1

    for name in list(set(map(name_fn, options.in_gams))):
        detailed_table +="{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
            name,
            float(sums[name][0]) / float(counts[name]),
            float(sums[name][1]) / float(counts[name]),
            float(sums[name][2]) / float(counts[name]),
            float(sums[name][3]) / float(counts[name]),
            float(sums[name][4]) / float(counts[name]),
            float(sums[name][5]) / float(counts[name]),
            float(sums[name][6]) / float(counts[name]),
            float(sums[name][7]) / float(counts[name]),
            float(sums[name][8]) / float(counts[name]),
            float(sums[name][9]) / float(counts[name]))


    with open(detailed_call_coutns_tsv_path(options), "w") as ofile:
        ofile.write(detailed_table)
        

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

    robust_makedirs(options.output_path)

    # hack in some fake alignments so that the g1kvcf sample graphs
    # get added to the tables (same for platvcf)
    added = []
    orig_gams = options.in_gams
    for gam in options.in_gams:
        if os.path.dirname(gam).split("/")[-1] == "trivial":
            added.append(gam.replace("trivial/" + os.path.basename(gam),
                                     "g1kvcf/" + os.path.basename(gam)))
            added.append(gam.replace("trivial/" + os.path.basename(gam),
                                     "platvcf/" + os.path.basename(gam)))            
    options.in_gams = options.in_gams + added
        
                        
    # make some tables from the json comparison output
    snp_count_table(options)
    graph_size_table(options)

    # make some better stats (to eventuall replace above)
    options.in_gams = orig_gams
    detailed_call_count_table(options)

    # make some boxplots with adams super script
    boxPlot(count_tsv_path(options), count_tsv_path(options).replace(".tsv", ".pdf"),
            0, 2, "SNP\\ Count")
    boxPlot(size_tsv_path(options), size_tsv_path(options).replace(".tsv", ".pdf"),
            0, 1, "Sample\\ Graph\\ Size")

    
if __name__ == "__main__" :
    sys.exit(main(sys.argv))
        
        
        
        
        
        
        
        
        
        

