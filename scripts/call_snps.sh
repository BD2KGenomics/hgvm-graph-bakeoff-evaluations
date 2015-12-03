#!/bin/bash

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

OPTS="--maxCores 24 --vg_cores 2 --vg_only"

CALL_OPTS=" -r 0.01 -d 40 -s 20"
PILEUP_OPTS=" -s "
#PILEUP_OPTS=" "


mkdir -f $OUT_DIR

for i in "${REGIONS[@]}"
do
	 rm -rf ${TOIL_DIR} ; scripts/callVariants.py ./${TOIL_DIR} ${ALIGNMENTS}/${i}/*/*.gam --graph_dir ${GRAPHS} --out_dir ${OUT_DIR} ${OPTS} --call_opts "${CALL_OPTS}" --pileup_opts "${PILEUP_OPTS}"
done
