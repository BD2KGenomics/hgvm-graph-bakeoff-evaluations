#!/usr/bin/env python2.7

"""
Stick pileups onto vcf coordinates to help debugging
"""


import argparse, sys, os, os.path, random, subprocess, shutil, itertools, json
from collections import defaultdict

def parse_args(args):
    parser = argparse.ArgumentParser(description=__doc__, 
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("vcf_call", type=str,
                        help="vcf derived from vg call (- for stdin)"),
    parser.add_argument("txt_call", type=str,
                        help="corresponding -c output of vg call")
    args = args[1:]
    options = parser.parse_args(args)
    return options

def parse_id(tok):
    snp_toks = tok.split(".")
    node_id = int(snp_toks[0])
    node_offset = int(snp_toks[1]) if len(snp_toks) > 1 else 0
    return (node_id, node_offset)
    
def main(args):
    options = parse_args(args)

    # open input
    in_vcf = open(options.vcf_call) if options.vcf_call != "-" else sys.stdin
    in_txt = open(options.txt_call)

    # in_vcf not necessarily sorted
    vcf_lines = sorted([line for line in in_vcf], key = lambda x : (0,0) if x[0] == "#" else parse_id(x.split()[2]))

    # in_sample not necessarily sorted?
    txt_lines = sorted([line for line in in_txt], key = lambda x : (int((x.split()[0])), int(x.split()[1])) )
    txt_line_no = 0

    for vcf_line in vcf_lines:
        if vcf_line[0] != "#":
            toks = vcf_line.split()
            chrom, start, snp_id, ref, alts = toks[0:5]

            try:
                pileup = "PILEUP\tMissing-Coord"
                node_id, node_offset = parse_id(snp_id)

                while txt_line_no < len(txt_lines):
                    txt_line = txt_lines[txt_line_no]
                    txt_toks = txt_line.split()
                    txt_node_id = int(txt_toks[0])
                    txt_node_offset = int(txt_toks[1]) - 1
                    if txt_node_id == node_id and txt_node_offset == node_offset:
                        pileup = "\t".join(["PILEUP"] + txt_toks[0:6])
                        break
                    txt_line_no += 1
            except Exception as e:
                sys.stderr.write(str(e))
                pileup = "PILEUP\tMissing-Coord"
            sys.stdout.write("\t".join(toks[0:5] + [pileup]) + "\n")
            
    
    if options.vcf_call != "-":
        in_vcf.close()
    in_txt.close()

if __name__ == "__main__" :
    sys.exit(main(sys.argv))
