#!/bin/bash

# Generates a variety of statistics on the output of call_snps.sh
# The first three arguments should be those used for call_snps.sh
# This script is very computationally intense due to the all-to-all comparisons...

if [ "$#" -ne 4 ]; then
	 echo "Syntax $0 <graphs_dir> <alignments_dir> <calls_dir> <out_dir>"
	 exit 1
fi

GRAPHS=$1
ALIGNMENTS=$2
VARIANTS=$3
OUT_DIR=$4
TOIL_DIR=cps_toil_dir
TOIL_OPTS="--maxCores 48 --vg_cores 4"
INDEX_OPTS="--kmer 20 --edge_max 7 --ignore_ns --timeout 10000"
#INDEX_OPTS="--kmer 20 --edge_max 7 --ignore_ns --timeout 20 --overwrite"
TAG_OPTS="--orig_tag ${GRAPHS}"
COMPS=( "jaccard" "recall" "corg" )
#COMPS=( "jaccard" )
REGIONS=( "brca1" "brca2" "sma" "lrc_kir" "mhc" )
#REGIONS=( "brca1" )

# vglr lrc_kir is a bad graph.  we can censor it as input in the wildcard arguments below
# to make it disappear from the analysis
GLOBIGNORE="*/lrc_kir/vglr*:*vglr-lrc_kir*"

# These are the g1kvcf samples that don't exist.  right now heatmaps are just computed without
# them.  so squares compareing g1kvcf samples will be average of 5 datapoints whereas all
# other squares will ne average of 9 datapoints.  We can ensure only 5 datapoints ever used
# by uncommenting this, but it seems a bit drastic
#GLOBIGNORE="*HG00512*:*HG00514*:*HG00733*:*NA19240*":${GLOBIGNORE}

mkdir $OUT_DIR

for i in "${REGIONS[@]}"
do
	 mkdir ${OUT_DIR}/${i}
	 for j in "${COMPS[@]}"
	 do
		  # millions of temprorary output files go here 
		  COMPDIR="compare"
		  if [ $j == "corg" ]
		  then
				COMPDIR="corg"
		  fi
		  
		  # all vs all original graphs
		  rm ${OUT_DIR}/${i}/comp_orig_${j}_distmat_${i}.tsv
		  rm -rf ${TOIL_DIR} ; scripts/computeDistances.py ./${TOIL_DIR} ${GRAPHS}/*${i}*.vg ${OUT_DIR}/${i}/${COMPDIR} ${OUT_DIR}/${i}/comp_orig_${j}_distmat_${i}.tsv ${j} ${TOIL_OPTS} ${INDEX_OPTS} 

		  # all vs all sample graphs
		  rm ${OUT_DIR}/${i}/comp_sample_${j}_distmat_${i}.tsv
		  #rm -rf ${TOIL_DIR} ; scripts/computeDistances.py ./${TOIL_DIR} ${VARIANTS}/${i}/*/*sample*.vg ${OUT_DIR}/${i}/${COMPDIR} ${OUT_DIR}/${i}/comp_sample_${j}_distmat_${i}.tsv ${j} ${TAG_OPTS} ${TOIL_OPTS} ${INDEX_OPTS} --dir_tag --avg_sample

		  # all vs all original and sample graphs
		  rm ${OUT_DIR}/${i}/comp_sample_and_orig_${j}_distmat_${i}.tsv
		  #rm -rf ${TOIL_DIR} ; scripts/computeDistances.py ./${TOIL_DIR} ${GRAPHS}/*${i}*.vg ${VARIANTS}/${i}/*/*sample*.vg ${OUT_DIR}/${i}/${COMPDIR} ${OUT_DIR}/${i}/comp_sample_and_orig_${j}_distmat_${i}.tsv ${j} ${TAG_OPTS} ${TOIL_OPTS} ${INDEX_OPTS} --dir_tag --avg_sample 
	 done
	 
	 # tables
	 scripts/callStats.py ${ALIGNMENTS}/${i}/*/*.gam ${OUT_DIR}/${i} --out_dir ${VARIANTS}  --tag $i  --graph_dir ${GRAPHS} --avg_sample
	 rm -rf ${TOIL_DIR} ; scripts/trioStats.py ./${TOIL_DIR} ${ALIGNMENTS}/${i}/*/*.gam ${OUT_DIR}/${i} --out_dir ${VARIANTS} --tag $i ${TOIL_OPTS}
	 	 
done
