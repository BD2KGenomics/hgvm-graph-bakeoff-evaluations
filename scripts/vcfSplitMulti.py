#!/usr/bin/env python2.7

"""
 Split multiallelic genotpes of form 1|2, 2|1, 1/2, or 2/1 into two single allelic
 lines, each with 0/1 using vt decompose
"""


import argparse, sys, os, os.path, random, subprocess, shutil, itertools, tempfile

def parse_args(args):
    parser = argparse.ArgumentParser(description=__doc__, 
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("in_vcf", type=str,
                        help="Input vcf file (- for stdin)")
                        
    args = args[1:]
    options = parser.parse_args(args)
    return options

# todo, really need to centralize vcf parse code (better yet use actual api)
def get_gt(toks, options):
    """ get genotype as list of strings """
    gth = toks[-2].split(":")
    gts = toks[-1].split(":")
    gt_idx = gth.index("GT")
    gt = gts[gt_idx]
    gt = gt.split("/") if "/" in gt else gt.split("|")
    return gt

def fix_gt(toks, options):
    """ replace genotype dots with zeros """
    gth = toks[-2].split(":")
    gts = toks[-1].split(":")
    gt_idx = gth.index("GT")
    gt = gts[gt_idx].rstrip()
    gt = sorted(gt.split("/") if "/" in gt else gt.split("|"))
    # unambiguous, leave alone
    if gt in [["0", "1"], ["1", "1"]]:
        new_gt = "{}/{}".format(gt[0], gt[1])
    # completely ambiguous or trivial, don't bother
    elif gt in [[".", "0"], [".", "."], ["0", "0"]]:
        return None
    # full alt will get same haplotype count when broken apart 
    elif gt == [".", "1"]:
        new_gt = "1/1"
    else:
        assert False
    new_gts = ":".join(gts[: gt_idx] + [new_gt] + gts[gt_idx + 1:]) + "\n"
    return "\t".join(toks[:-1] + [new_gts])

def main(args):
    options = parse_args(args)

    if options.in_vcf == "-":
        vcf_file = [line for line in sys.stdin]
    else:
        with open(options.in_vcf) as f:
            vcf_file = [line for line in f]

    # split out all all 1/2 and 2/1 multiallelic variants
                
    out_base = tempfile.mkdtemp(prefix = "splitMulti_", dir = ".")
    out_mult = os.path.join(out_base, "multi.vcf")
    out_sing = os.path.join(out_base, "single.vcf")
    
    with open(out_mult, "w") as mf, open(out_sing, "w") as sf:
        for line in vcf_file:
            if line[0] != "#":
                toks = line.split()
                gtset = set(get_gt(toks, options))
                # todo: maybe we can just run decomp on whole thing?
                # idea of splitting here is to know exactly which variants
                # we want to go an fix later... 
                if any(g not in ["0", "1", "."] for g in gtset):
                    mf.write(line)
                else:
                    sf.write(line)
            else:
                mf.write(line)
                sf.write(line)

    # run decompose on our multi allelic variants

    out_dec = os.path.join(out_base, "decomp.vcf")
    out_tmp = os.path.join(out_base, "out_unsrt.vcf")
    
    os.system("vt decompose {} > {}".format(out_mult, out_dec))

    # spit out the output unsorted temp file
    header = False
    with open(out_tmp, "w") as tf:
        with open(out_dec) as df:
            for line in df:
                if line[0] != "#":
                    fixed_line = fix_gt(line.split("\t"), options)
                    if fixed_line != None:
                        tf.write(fixed_line)
                else:
                    tf.write(line)
                    header = True
        with open(out_sing) as sf:
            for line in sf:
                if line[0] != "#" or not header:
                    tf.write(line)

    # sort it
    out_srt = os.path.join(out_base, "out_sort.vcf")
    os.system("scripts/vcfsort {} > {}".format(out_tmp, out_srt))

    # output to stdout
    with open(out_srt) as sf:
        for line in sf:
            if len(line) > 1:
                sys.stdout.write(line)

    # 
    # remove temp directory
    os.system("rm -rf {}".format(out_base))
	 
if __name__ == "__main__" :
    sys.exit(main(sys.argv))
