#!/usr/bin/env python2.7

"""
quick, from scratch vcf compare script to help debug calls / sanity check gatk results
"""


import argparse, sys, os, os.path, random, subprocess, shutil, itertools, json
from collections import defaultdict
from toillib import RealTimeLogger, robust_makedirs
import tempfile

def parse_args(args):
    parser = argparse.ArgumentParser(description=__doc__, 
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("vcf1", type=str,
                        help="Input vcf file 1 (query)"),
    parser.add_argument("vcf2", type=str,
                        help="Input vcf file 2 (truth)")
    parser.add_argument("-c", action="store_true", default=False,
                        help="Ignore sequence name (only look at start)")
    parser.add_argument("-i", action="append", default=[],
                        help="Ignore lines contaning keyword")

    args = args[1:]
    options = parser.parse_args(args)
    return options

def parse_alts(toks):
    """ get the list of alts """
    return toks[4].split(",")

def parse_ref(toks):
    """ return reference """
    return toks[3]


def make_vcf_dict(vcf_path, options):
    """ load up all variants by their coordinates
    map (chrom, pos) -> [(ref, alt), (ref, alts) etc.]
    """
    vcf_dict = defaultdict(set)
    with open(vcf_path) as f:
        for line in f:
            skip = line[0] == "#"
            for ignore_keyword in options.i:
                if ignore_keyword in line:
                    skip = True
            if not skip:
                toks = line.split()
                chrom = toks[0] if options.c is False else None
                pos = int(toks[1])
                ref = parse_ref(toks)
                alts = parse_alts(toks)
                for i in range(len(alts)):
                    vcf_dict[(chrom, pos)].add((ref, alts[i]))
    return vcf_dict

def alt_cat(ref, alt):
    """ 0 ref, 1 snp, 2 multibase snp, 3 insert 4 delete """
    if ref == alt:
        return 0
    elif len(ref) == len(alt):
        return 1 if len(ref) == 1 else 2
    return 3 if len(ref) < len(alt) else 4

def cat_name(c):
    return ["REF", "SNP", "MULTIBASE_SNP", "INSERT", "DELETE", "TOTAL"][c]

def find_alt(chrom, pos, ref, alt, vcf_dict):
    """ find an alt in a dict """
    for val in vcf_dict[(chrom, pos)]:
        if val[0] == ref and val[1] == alt:
            return True
    return False

def compare_vcf_dicts(vcf_dict1, vcf_dict2):
    """ check dict1 against dict2 """
    # see alt_cat() for what 5 values mean
    total_alts = [0, 0, 0, 0, 0]
    found_alts = [0, 0, 0, 0, 0]
    total_alleles = [0, 0, 0, 0, 0]
    found_alleles = [0, 0, 0, 0, 0]
    
    for key1, val1 in vcf_dict1.items():
        chrom1, pos1 = key1
        for ref1, alt1 in val1:
            total_alts[alt_cat(ref1, alt1)] += 1
            if find_alt(chrom1, pos1, ref1, alt1, vcf_dict2):
                found_alts[alt_cat(ref1, alt1)] += 1

    total_alts += [sum(total_alts)]
    found_alts += [sum(found_alts)]
    total_alleles += [sum(total_alleles)]
    found_alleles += [sum(found_alleles)]

    return total_alts, found_alts, total_alleles, found_alleles

def json_acc(vcf1, vcf2, options):
    """ compute the accuracy """
    vcf_dict1 = make_vcf_dict(vcf1, options)
    vcf_dict2 = make_vcf_dict(vcf2, options)
    total_alts1, found_alts1, total_alleles1, found_alleles1 = compare_vcf_dicts(vcf_dict1, vcf_dict2)
    total_alts2, found_alts2, total_alleles2, found_alleles2 = compare_vcf_dicts(vcf_dict2, vcf_dict1)    

    json_data = dict()
    json_data["Path1"] = options.vcf1
    json_data["Path2"] = options.vcf2
    json_data["Alts"] = dict()
    json_data["Alleles"] = dict()
    for c in range(len(found_alts1) + 1):
        label = cat_name(c) if c < len(found_alts1) else "INDEL"
        if c > 0:
            if label == "INDEL":
                tp = found_alts1[3] + found_alts1[4]
                fp = total_alts1[3] + total_alts1[4] - tp
                fn = total_alts2[3] + total_alts2[4] - tp
            else:
                tp = found_alts1[c]
                fp = total_alts1[c] - tp
                fn = total_alts2[c] - tp
            prec = 0. if tp + fp == 0 else float(tp) / float(tp + fp)
            rec = 0. if tp + fn  == 0 else float(tp) / float(tp + fn)
            json_data["Alts"][label] = {"TP" : tp, "FP" : fp, "FN" : fn, "Precision" : prec, "Recall" : rec}

        if label == "INDEL":
                tp = found_alleles1[3] + found_alleles1[4]
                fp = total_alleles1[3] + total_alleles1[4] - tp
                fn = total_alleles2[3] + total_alleles2[4] - tp
        else:
            tp = found_alleles1[c]
            fp = total_alleles1[c] - tp
            fn = total_alleles2[c] - tp
        prec = 0. if tp + fp == 0 else float(tp) / float(tp + fp)
        rec = 0. if tp + fn  == 0 else float(tp) / float(tp + fn)
        json_data["Alleles"][label] = {"TP" : tp, "FP" : fp, "FN" : fn, "Precision" : prec, "Recall" : rec}


    return json.dumps(json_data)

def system(command):
    sts = subprocess.call(command, shell=True, bufsize=-1, stdout=sys.stdout, stderr=sys.stderr)
    if sts != 0:
        raise RuntimeError("Command: %s exited with non-zero status %i" % (command, sts))
    return sts

def main(args):
    options = parse_args(args)

    vcf1 = options.vcf1
    vcf2 = options.vcf2        
    
    print json_acc(vcf1, vcf2, options)

if __name__ == "__main__" :
    sys.exit(main(sys.argv))
