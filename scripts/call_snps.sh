#!/bin/bash

if [ "$#" -ne 3 ]; then
	 echo "Syntax $0 <graphs_dir> <alignments_dir> <out_dir>"
	 exit 1
fi

GRAPHS=$1
ALIGNMENTS=$2
OUT_DIR=$3
TOIL_DIR=cs_toil_dir

OPTS="--maxCores 24 --vg_cores 2 --vg_only"

mkdir -f $OUT_DIR

for i in brca1 brca2 sma lrc_kir mhc
#for i in brca1
do
	 rm -rf ${TOIL_DIR} ; scripts/callVariants.py ./${TOIL_DIR} ${ALIGNMENTS}/${i}/*/*.gam --graph_dir ${GRAPHS} --out_dir ${OUT_DIR} ${OPTS}
done
