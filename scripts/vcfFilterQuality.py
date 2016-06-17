#!/usr/bin/env python2.7

"""
 Filter out records that don't meet a particular quality threshold.
Quality assumed to be 6th column
"""


import argparse, sys, os, os.path, random, subprocess, shutil, itertools

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
                        help="Use XAAD file")
    parser.add_argument("--dedupe", action="store_true",
                        help="Filter entries with same coordinate, choosing one with higest quality")
                        
    args = args[1:]
    options = parser.parse_args(args)
    return options

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
        assert len(gt) <= len(ads)
        min_ad = sys.maxint
        for i, g in enumerate(gt):
            g = i if g == "." else int(g)
            min_ad = min(min_ad, ads[g])
        return float(min_ad)
    elif options.xaad is True:
        toks = line.split("\t")
        gth = toks[-2].split(":")
        gts = toks[-1].split(":")
        xaad_idx = gth.index("XAAD")
        xaad = int(gts[xaad_idx])
        return float(xaad)
    else:
        # quality 
        return float(line.split("\t")[5])
        

def compute_cutoff(options):
    if options.pct is False:
        return options.min_qual
    else:
        assert options.min_qual >= 0. and options.min_qual <= 1.
        if options.in_vcf == "-":
            vcf_file = sys.stdin
        else:
            vcf_file = open(options.in_vcf)
        quals = []
        for line in vcf_file:
            if line[0] != "#":
                quals.append(get_qual_from_line(line, options))
        if options.in_vcf != "-":
            vcf_file.close()
        return sorted(quals)[int(options.min_qual * len(quals))]
        

def main(args):
    options = parse_args(args)

    cutoff = compute_cutoff(options)

    if options.in_vcf == "-":
        vcf_file = sys.stdin
    else:
        vcf_file = open(options.in_vcf)

    buf = None, None, None, None # chrom , start, qual ,line
    for line in vcf_file:
        if line[0] == "#":
            sys.stdout.write(line)
        else:
            toks = line.split("\t")
            chrom, start = toks[0], int(toks[1])
            qual = get_qual_from_line(line, options)
            # new coordinate, write and clear buffer
            if not options.dedupe or (chrom, start) != (buf[0], buf[1]):
                if buf[0] != None:
                    sys.stdout.write(buf[3])
                    buf = None, None, None, None

            # update buffer
            if qual >= cutoff and (buf[2] == None or qual > buf[2]):
                buf = chrom, start, qual, line

    # write buffer
    if buf[0] != None:
        sys.stdout.write(buf[3])

    if options.in_vcf != "-":
        vcf_file.close()
	 
if __name__ == "__main__" :
    sys.exit(main(sys.argv))
