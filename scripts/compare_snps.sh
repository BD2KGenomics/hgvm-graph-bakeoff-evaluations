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
INDEX_OPTS="--kmer 20 --edge_max 7 --timeout 10000"
#COMPS=( "vcf" "kmer" "corg" )
COMPS=( "happy" )
#COMP_OPTS="--orig --orig_and_sample"
COMP_OPTS="--clip --normalize --ignore Conflict --ignore Silver"
#COMP_OPTS="--clip  --ignore Conflict --ignore Silver"
REGIONS=( "brca1" "brca2" "sma" "lrc_kir" "mhc" )
REGIONS=( "brca2" "mhc" )

# vglr lrc_kir is a bad graph.  we can censor it as input in the wildcard arguments below
# to make it disappear from the analysis
GLOBIGNORE="*/lrc_kir/vglr*:*vglr-lrc_kir*"
#GLOBIGNORE="*vglr*"

# leave simons out for now as its hg19
GLOBIGNORE="*simons*:${GLOBIGNORE}"

# These are the g1kvcf samples that don't exist.  right now heatmaps are just computed without
# them.  so squares compareing g1kvcf samples will be average of 5 datapoints whereas all
# other squares will ne average of 9 datapoints.  We can ensure only 5 datapoints ever used
# by uncommenting this, but it seems a bit drastic
#GLOBIGNORE="*HG00512*:*HG00514*:*HG00733*:*NA19240*":${GLOBIGNORE}

mkdir $OUT_DIR

for i in "${REGIONS[@]}"
do
	 for j in "${COMPS[@]}"
	 do
		  # compute distances
		  rm -rf ${TOIL_DIR} ; scripts/computeVariantsDistances.py ./${TOIL_DIR} ${ALIGNMENTS}/${i}/*/NA12878.gam ${VARIANTS} ${GRAPHS} ${j} ${OUT_DIR} ${COMP_OPTS}  ${TOIL_OPTS} ${INDEX_OPTS} 
	 done

	 # plots
	 scripts/plotVariantsDistances.py ${OUT_DIR} &
	 # tables
	 #mkdir ${OUT_DIR}/call_stats
	 #scripts/callStats.py ${ALIGNMENTS}/${i}/*/*.gam ${OUT_DIR}/call_stats --out_dir ${VARIANTS}  --tag $i  --graph_dir ${GRAPHS} --avg_sample
	 #mkdir ${OUT_DIR}/trio_stats
	 #rm -rf ${TOIL_DIR} ; scripts/trioStats.py ./${TOIL_DIR} ${ALIGNMENTS}/${i}/*/*.gam ${OUT_DIR}/trio_stats --out_dir ${VARIANTS} --tag $i ${TOIL_OPTS}
done

wait
