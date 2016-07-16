#!/usr/bin/env python2.7

"""
 Filter out records that don't meet a particular quality threshold.
Quality assumed to be 6th column
"""


import argparse, sys, os, os.path, random, subprocess, shutil, itertools, math

def parse_args(args):
    parser = argparse.ArgumentParser(description=__doc__, 
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("in_vcf", type=str,
                        help="Input vcf file (- for stdin)")
    parser.add_argument("min_qual", type=float,
                        help="Mininum quality value to keep")
    parser.add_argument("--pct", action="store_true",
                        help="Interpret min_qual as a precentile.  So min_qual == 0.7 will return the top 30pct scores")
    parser.add_argument("--info", type=str, default=None,
                        help="Use given info field")
    parser.add_argument("--ad", action="store_true",
                        help="Use minimum AD from genotype")
    parser.add_argument("--xaad", action="store_true",
                        help="Use XAAD field")
    parser.add_argument("--al", action="store_true",
                        help="Use minimum AL from genotype")
    parser.add_argument("--gq", action="store_true",
                        help="Use GQ from genotype")    
    parser.add_argument("--xl", action="store_true",
                        help="Use AL * XAAD field")
    parser.add_argument("--dedupe", action="store_true",
                        help="Filter entries with same coordinate, choosing one with higest quality")
    parser.add_argument("--max_depth", type=int, default=None,
                        help="Filter entries with DEPTH > MAX_DEPTH - CUTOFF (also applied with --ad, --xaad)")
    parser.add_argument("--qual", action="store_true",
                        help="Use Quality field (default)")
    parser.add_argument("--keep_trivial", action="store_true",
                        help="Dont ignore 0/0 and ./. genotypes")
    parser.add_argument("--set_qual", action="store_true",
                        help="Write whatever is used as quality in the quality field")
                        
    args = args[1:]
    options = parser.parse_args(args)
    return options

def trivial_gt(line, options):
    if options.keep_trivial is True:
        return False
    # filter out ./. and 0/0
    toks = line.split("\t")
    gth = toks[-2].split(":")
    gts = toks[-1].split(":")
    gt_idx = gth.index("GT")
    gt = gts[gt_idx]
    gt = gt.split("/") if "/" in gt else gt.split("|")
    for g in gt:
        if g not in ["0", "."]:
            return False
    return True
    
def get_qual_from_line(line, options):
    # note should switch to vcf api...
    if options.info is not None:
        assert options.ad is False
        info = line.split("\t")[7]
        toks = info.split(";")
        for tok in toks:
            if tok[:len(options.info) + 1] == "{}=".format(options.info):
                return float(tok[len(options.info) + 1:])
        assert False
    elif options.ad is True:
        assert options.xaad == False
        # this block is deprecated because of xaad.  keeping it
        # around for near term in case we need sanity check. 
        toks = line.split("\t")
        gth = toks[-2].split(":")
        gts = toks[-1].split(":")
        gt_idx = gth.index("GT")
        gt = gts[gt_idx]
        gt = gt.split("/") if "/" in gt else gt.split("|")
        ad_idx = gth.index("AD")
        ads = gts[ad_idx]
        ads = ads.split(",")
        ads = [int(x) for x in ads]
        min_ad = sys.maxint
        for i, g in enumerate(gt):
            g = i if g == "." else int(g)
            min_ad = min(min_ad, ads[g])
        return 0. if all(g == '0' for g in gt) else float(min_ad)
    elif options.xaad is True:
        toks = line.split("\t")
        gth = toks[-2].split(":")
        gts = toks[-1].split(":")
        xaad_idx = gth.index("XAAD")
        xaad = int(gts[xaad_idx])
        return float(xaad)
    elif options.gq is True:
        toks = line.split("\t")
        gth = toks[-2].split(":")
        gts = toks[-1].split(":")
        gq_idx = gth.index("GQ")
        gq = int(gts[gq_idx])
        return float(gq)
    elif options.al is True:
        toks = line.split("\t")
        gth = toks[-2].split(":")
        gts = toks[-1].split(":")
        gt_idx = gth.index("GT")
        gt = gts[gt_idx]
        gt = gt.split("/") if "/" in gt else gt.split("|")
        al_idx = gth.index("AL")
        als = gts[al_idx]
        als = als.split(",")
        als = [float(x) for x in als]
        assert len(gt) <= len(als)
        min_al = float(sys.maxint)
        for i, g in enumerate(gt):
            g = i if g == "." else int(g)
            min_al = min(min_al, als[g])
        return 0. if all(g == '0' for g in gt) else float(min_al)
    elif options.xl is True:
        toks = line.split("\t")
        gth = toks[-2].split(":")
        gts = toks[-1].split(":")
        ll_idx = gth.index("AL")
        lls = gts[ll_idx].split(",")
        ll = float(lls[1])
        xaad_idx = gth.index("XAAD")
        xaad = float(gts[xaad_idx])
        return ll * xaad
    else:
        # quality 
        return float(line.split("\t")[5])
        

def compute_cutoff(vcf_file, options):
    if options.pct is False:
        return options.min_qual
    else:
        assert options.min_qual >= 0. and options.min_qual <= 1.
        quals = [] 
        for line in vcf_file:
            if line[0] != "#" and not trivial_gt(line, options):
                quals.append(get_qual_from_line(line, options))

        # do our percentile on unique values
        quals = list(set(quals))
                
        return sorted(quals)[int(options.min_qual * len(quals))]
        
def main(args):
    options = parse_args(args)

    if options.in_vcf == "-":
        vcf_file = [line for line in sys.stdin]
    else:
        with open(options.in_vcf) as f:
            vcf_file = [line for line in f]
    
    cutoff = compute_cutoff(vcf_file, options)
    max_cutoff = sys.maxint if options.max_depth is None else options.max_depth - cutoff
    sys.stderr.write("Cutoff = ({}, {})\n".format(cutoff, max_cutoff))

    buf = None, None, None, None # chrom , start, qual ,line
    for line in vcf_file:
        if line[0] == "#":
            sys.stdout.write(line)
        elif not trivial_gt(line, options):
            toks = line.split("\t")
            chrom, start = toks[0], int(toks[1])
            qual = get_qual_from_line(line, options)

            if options.set_qual is True:
                line = "\t".join(toks[:5] + [str(qual)] + toks[6:])

            # get depth
            if options.max_depth is not None:
                gth = toks[-2].split(":")
                gts = toks[-1].split(":")
                dp_idx = gth.index("DP")
                depth = float(gts[dp_idx])
            else:
                depth = -1

            # new coordinate, write and clear buffer
            if not options.dedupe or (chrom, start) != (buf[0], buf[1]):
                if buf[0] != None:
                    sys.stdout.write(buf[3])
                    buf = None, None, None, None

            # update buffer
            if qual >= cutoff and depth <= max_cutoff and (buf[2] == None or qual > buf[2]):
                if buf[3] is not None:
                    sys.stderr.write("favouring {} over\n{}\n\n".format(str([chrom, start, qual, line]), str(buf)))
                buf = chrom, start, qual, line

    # write buffer
    if buf[0] != None:
        sys.stdout.write(buf[3])
	 
if __name__ == "__main__" :
    sys.exit(main(sys.argv))
