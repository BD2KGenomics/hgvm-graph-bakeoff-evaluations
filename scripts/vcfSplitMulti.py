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
    parser.add_argument("--merge", action="store_true",
                        help="Undo split my merging pairs of lines with 1/1 genotypes (and same coord) into multiallic site")
                        
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

def merge_multi(vcf_file, options):
    """ merge consecutive lines into single multiallele site if they were split
    previously (remember: we can't do this when splitting since the vt decompose_blocksub
    needs to get run in between) """
    i = 0
    while i < len(vcf_file):
        line = vcf_file[i]
        if i < len(vcf_file) - 1 and line[0] != "#":
            next_line = vcf_file[i+1]
            toks = line.split("\t")
            next_toks = next_line.split("\t")
            gt = get_gt(toks, options)
            next_gt = get_gt(next_toks, options)
            if all(toks[x] == next_toks[x] for x in range(4)) and\
               toks[5] != next_toks[5] and gt == next_gt and gt == ['1','1']:
                qual = float(toks[5])
                next_qual = float(next_toks[5])
                merge_toks = toks
                if next_qual < qual:
                    merge_toks[5] = next_toks[5]
                    merge_toks[7] = next_toks[7]
                merge_toks[4] = "{},{}".format(toks[4], next_toks[4])
                merge_toks[8] = "GT"
                merge_toks[9] = "1/2"
                sys.stdout.write("\t".join(merge_toks) + "\n")
                i += 2
                continue
        sys.stdout.write(line)
        i += 1
    

def main(args):
    options = parse_args(args)

    if options.in_vcf == "-":
        vcf_file = [line for line in sys.stdin]
    else:
        with open(options.in_vcf) as f:
            vcf_file = [line for line in f]

    if options.merge:
        # do join instead of split
        return merge_multi(vcf_file, options)
    
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

    return 0
	 
if __name__ == "__main__" :
    sys.exit(main(sys.argv))
