#!/usr/bin/env python2.7

"""
 Filter a sample out of a vcf, keeping only snps that are present in 
 genotype information for that sample.  Produce filtered vcf *and*
 reference fasta file, which will have double-alt positions
 embedded to help vg construct not write the reference.

 https://github.com/ekg/vcflib must be installed (with vcflib/bin) in PATH. 
"""


import argparse, sys, os, os.path, random, subprocess, shutil, itertools
from Bio import SeqIO
from Bio.Seq import Seq

def parse_args(args):
    parser = argparse.ArgumentParser(description=__doc__, 
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("in_vcf", type=str,
                        help="Input vcf file"),
    parser.add_argument("in_fa", type=str,
                        help="Input fa")
    parser.add_argument("sample", type=str,
                        help="Sample ID to extract")
    parser.add_argument("out_vcf", type=str,
                        help="Output vcf file")
    parser.add_argument("out_fa", type=str,
                        help="Output fasta file")
    
    args = args[1:]
    options = parser.parse_args(args)
    return options

def main(args):
    options = parse_args(args)

    out_vcf = open(options.out_vcf, "w")
    in_fa = open(options.in_fa)
    record_dict = SeqIO.to_dict(SeqIO.parse(in_fa, "fasta"))
    record = None
    sequence = None
    indel_offset = 0
    
    filterCmd = "vcfkeepsamples {} {}".format(options.in_vcf, options.sample)
    
    filterProc = subprocess.Popen(filterCmd, shell=True, stdout=subprocess.PIPE)

    first = True
    last_pos = None
    while True:
        line = filterProc.stdout.readline()
        if line != "":
            # copy comments
            if line[0] == "#":
                out_vcf.write(line)
            else:
                toks = line.split()
                seq, vcf_pos, ref, alts = toks[0], int(toks[1]), toks[3], toks[4]
                assert last_pos is None or int(vcf_pos) >= last_pos
                last_pos = vcf_pos
                gti = line.rfind("GT")
                gt = None
                if gti > 0:
                    gts = line[gti + 2:].split()
                    # only expect single sample
                    assert len(gts) == 1
                    gt = gts[0]
                if gt and ("1" in gt or "2" in gt):
                    # load up fasta for first time
                    if record is None:
                        record = record_dict[seq]
                        sequence = str(record.seq)
                    assert record.id == seq

                    # sanity check between vcf and fasta
                    assert record.seq[(vcf_pos - 1) : (vcf_pos - 1 + len(ref))] == ref

                    # we have a legitimate snp for this genotype
                    alts = alts.split(",")
                    
                    # keep only alts for this genotype
                    gen_alts = []
                    if "1" in gt:
                        gen_alts.append(alts[0])
                    if "2" in gt:
                        gen_alts.append(alts[1])
                    alts = ",".join(gen_alts)
                    
                    # apply cumulative indel offset from previous events
                    vcf_pos = vcf_pos + indel_offset
                                        
                    # change the reference to first alt if it's not in genotype
                    if "0" not in gt:
                        alt = gen_alts[0]
                        indel_len = len(alt) - len(ref)
                        if len(alt) == len(ref):
                            # snp
                            sequence = sequence[:vcf_pos - 1] + alt + sequence[vcf_pos + len(ref) - 1:]
                        elif len(alt) > len(ref):
                            # insertion
                            sequence = sequence[:vcf_pos - 1] + alt + sequence[vcf_pos + len(ref) - 1:]
                        elif len(alt) < len(ref) and len(ref) > 1:
                            # deletion
                            sequence = sequence[:vcf_pos] + sequence[(vcf_pos - indel_len):]
                        else:
                            print line
                            assert False
                        indel_offset += indel_len
                        assert len(sequence) == len(record) + indel_offset
                        
                    else:
                        # write out the snp
                        out_vcf.write("\t".join([seq] + [str(vcf_pos)] + [ref] + [alts] + toks[4:]) + "\n")

        else:
            break

    in_fa.close()
    out_vcf.close()

    out_fa = open(options.out_fa, "w")
    assert len(sequence) == len(record) + indel_offset
    record.seq = Seq(sequence)
    out_fa.write(record.format("fasta"))
    out_fa.close()
	 
if __name__ == "__main__" :
    sys.exit(main(sys.argv))
