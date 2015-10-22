#!/bin/bash

if [ "$#" -ne  ]; then
	 echo "Syntax $0 <graphs_dir> <out_dir>"
	 exit 1
fi

GRAPHS=$1
OUT_DIR=$2

OPTS="--maxCores 20 --kmer 27 --edge_max 7"

mkdir -f $OUT_DIR

for i in brca1 brca2 sma lrc_kr mhc cenx
do
	 rm -rf zzz12 ; ./clusterGraphs.py ./zzz12 ${GRAPHS}/*${i}*.vg ${OUT_DIR} ${OPTS}
	 for j in heatmap.pdf heatmap_log.pdf tree.dot  tree.newick  tree.png
	 do
		  cp ${OUT_DIR}/${j} ${OUT_DIR}/${i}_${j}
	 done
done
