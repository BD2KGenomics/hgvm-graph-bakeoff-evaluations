#!/bin/bash

# Given a sorted BAM, produce a 1-based inclusive BED File of coordinates that are >= a certain depth threshold

if [ "$#" -ne 2 ]; then
	 echo "Syntax $0 <sorted_bam> <min_depth>"
	 exit 1
fi

BAM=$1
MIN_DEPTH=$2

# the +1 in the first awk is to print a real BED coordinate (where last position is exclusive)
# the -1 in the second awk brings this back to Inclusive coordinate which (i think) is what bcftools wants
samtools depth $BAM | awk -v t=$MIN_DEPTH '{ if ($3 > t) print $1"\t"$2"\t"$2 + 1}' | mergeBed -d 10 | awk '{ print $1"\t"$2"\t"$3-1 }' | sed "s/chr//g"
