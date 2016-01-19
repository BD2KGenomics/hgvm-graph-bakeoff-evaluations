#!/usr/bin/env python2.7

"""
 Filter out indels and, optionally, overlaps and multibase snps. Goal is to make easy test
 set to compare snp caller to. 
"""


import argparse, sys, os, os.path, random, subprocess, shutil, itertools

def parse_args(args):
    parser = argparse.ArgumentParser(description=__doc__, 
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("in_vcf", type=str,
                        help="Input vcf file (- for stdin)"),
    parser.add_argument("--overlap", action="store_true",
                        help="Filter snps that overlap indels")
    parser.add_argument("--multi", action="store_true",
                        help="Filter multibase snps")
    parser.add_argument("--qual", action="store_true",
                        help="Require 7th column = PASS")
    
    args = args[1:]
    options = parser.parse_args(args)
    return options

def main(args):
    options = parse_args(args)

    if options.in_vcf == "-":
        vcf_file = sys.stdin
    else:
        vcf_file = open(options.in_vcf)

    multibase_count = 0
    insert_count = 0
    delete_count = 0
    overlap_count = 0
    pass_count = 0

    prev = 0
    
    while True:
        line = vcf_file.readline()
        if line != "":
            # copy comments
            if line[0] == "#":
                sys.stdout.write(line)
            else:
                toks = line.split()
                seq, vcf_pos, ref, alts = toks[0], int(toks[1]), toks[3], toks[4]

                # find the longest alt
                alts = alts.split(",")
                alt_lens = [x for x in enumerate(map(len, alts))]
                max_alt = alts[max(alt_lens, key=lambda x : x[1])[0]]

                skip = False

                # check pass field
                if options.qual is True and toks[6] != "PASS":
                    pass_count +=1
                    skip = True

                # check multibase
                if skip is False and (len(ref) > 1 or len(max_alt) > 1):
                    # keep track of number of each case
                    if len(ref) > len(max_alt):
                        delete_count += 1
                        skip = True
                    elif len(ref) < len(max_alt):
                        insert_count += 1
                        skip = True
                    elif options.multi is True:
                        multibase_count += 1
                        skip = True
                        
                # check overlap with previous
                if skip is False and options.overlap is True and prev >= vcf_pos:
                    skip = True
                    overlap_count += 1

                # write to output
                if skip is False:
                    sys.stdout.write(line)

                # update prev
                end_pos = vcf_pos + max(len(ref), len(max_alt)) - 1
                prev = max(prev, end_pos)

        else:
            break

    if options.in_vcf != "-":
        vcf_file.close()

    sys.stderr.write("Filter stats {}\nInserts:{}\nDeletes:{}\nMultibaseSNPS:{}\nOverlaps:{}\nTotal:{}\n".format(
        options.in_vcf,
        insert_count,
        delete_count,
        multibase_count,
        overlap_count,
        (insert_count + delete_count + multibase_count + overlap_count)))
	 
if __name__ == "__main__" :
    sys.exit(main(sys.argv))
