#!/bin/bash

# run the variant calling pipeline on vg alignment data.  Directory structure
# of input is important, and explained in ../README.md and callVariants.py
# running is grouped by region to slightly facilitate debugging and resuming

if [ "$#" -ne 3 ]; then
	 echo "Syntax $0 <graphs_dir> <alignments_dir> <out_dir>"
	 exit 1
fi

GRAPHS=$1
ALIGNMENTS=$2
OUT_DIR=$3
TOIL_DIR=cs_toil_dir

REGIONS=( "brca1" "brca2" "sma" "lrc_kir" "mhc" )
#REGIONS=( "brca1" )
OPTS="--maxCores 36 --vg_cores 2 --vg_only"

CALL_OPTS=" -r 0.001 -d 20 -s 15"
PILEUP_OPTS=" "
FILTER_OPTS=" -e 0.0001 -f -d 1 -r 1"

mkdir $OUT_DIR

for i in "${REGIONS[@]}"
do
	 rm -rf ${TOIL_DIR} ; scripts/callVariants.py ./${TOIL_DIR} ${ALIGNMENTS}/${i}/*/*.gam --graph_dir ${GRAPHS} --out_dir ${OUT_DIR} ${OPTS} --call_opts "${CALL_OPTS}" --pileup_opts "${PILEUP_OPTS}" --filter_opts "${FILTER_OPTS}"
done
