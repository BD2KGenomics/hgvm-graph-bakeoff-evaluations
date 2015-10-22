#!/bin/bash
# count snps in vcf using bcftools stats.  

set -e

VCF_PATH=$1

bcftools stats $VCF_PATH  | grep "number of SNPs" | head -1 | awk '{print $6}'


