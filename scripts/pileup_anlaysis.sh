#!/bin/bash

# Generate VCF files for true and false positives that are annotated with pileup columns

if [ "$#" -ne 2 ]; then
	 echo "Syntax $0 <calls_dir> <out_dir>"
	 exit 1
fi

VARIANTS=$1
OUT_DIR=$2
REGIONS=( "brca1" "brca2" "sma" "lrc_kir" "mhc" )
GRAPHS=( "refonly" "snp1kg" "cactus" "haplo1kg30" )

mkdir ${OUT_DIR}

for i in "${REGIONS[@]}"
do
	 for j in "${GRAPHS[@]}"
	 do
		  # True Positives
		  scripts/vcfDelta.py data/platinum/NA12878/${i^^}_orig.vcf ${VARIANTS}/${i}/${j}/NA12878_sample_orig.vcf -ia | scripts/vcfPileups.py - ${VARIANTS}/${i}/${j}/NA12878_sample.txt > ${OUT_DIR}/${i}_${j}_tp.vcf

		  # False Positives
		  scripts/vcfDelta.py data/platinum/NA12878/${i^^}_orig.vcf ${VARIANTS}/${i}/${j}/NA12878_sample_orig.vcf -i | scripts/vcfPileups.py - ${VARIANTS}/${i}/${j}/NA12878_sample.txt > ${OUT_DIR}/${i}_${j}_fp.vcf
	 done
done
