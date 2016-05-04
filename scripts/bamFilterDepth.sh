#!/bin/bash

# make (vcf coord style) beds of bam regions that have a maximum depth

if [ "$#" -ne 2 ]; then
	 echo "Syntax $0 <reads_dir> <out_dir>"
	 exit 1
fi

READS_DIR=$1
OUT_DIR=$2

REGIONS=( "BRCA2" "MHC" "BRCA1" "SMA" "LRC_KIR" )
SAMPLES=( "NA12878" "NA12877" "NA12879" )
DEPTHS=( 50 75 100 200 500 )

function run_filter {
	 
	 BAM=$1
	 MAX_DEPTH=$2
	 OUT_FILE=$3

	 # the +1 in the first awk is to print a real BED coordinate (where last position is exclusive)
	 # the -1 in the second awk brings this back to Inclusive coordinate which (i think) is what bcftools wants
	 samtools depth $BAM | awk -v t=$MAX_DEPTH '{ if ($3 <= t) print $1"\t"$2"\t"$2 + 1}' | mergeBed -d 10 | awk '{ print $1"\t"$2"\t"$3-1 }' | sed "s/chr//g" >> ${OUT_FILE}

}

mkdir $OUT_DIR

for SAMPLE in "${SAMPLES[@]}"
do
	 for DEPTH in "${DEPTHS[@]}"
	 do
		  BED=${OUT_DIR}/${SAMPLE}_${DEPTH}.bed
		  rm ${BED}
		  for REGION in "${REGIONS[@]}"
		  do
				run_filter ${READS_DIR}/${REGION}/${SAMPLE}/${SAMPLE}.bam ${DEPTH} ${BED}
		  done
		  sortBed -i ${BED} > ${BED}.sort ; mv ${BED}.sort ${BED}
	 done
done
 



