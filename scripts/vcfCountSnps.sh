#!/bin/bash
# count snps in vcf using bcftools stats.  

set -e

VCF_PATH=$1
FA_PATH=$2

#bcftools stats $VCF_PATH  | grep "number of SNPs" | head -1 | awk '{print $6}'

vcfsort $VCF_PATH | vt decompose - | vt decompose_blocksub -a - | vt normalize -r $FA_PATH - | uniq | scripts/vcfFilterIndels.py - | grep -v "\#" | wc -l



