#!/bin/bash

# run the variant calling pipeline on vg alignment data.  Directory structure
# of input is important, and explained in ../README.md and callVariants.py
# running is grouped by region to slightly facilitate debugging and resuming

if [ "$#" -ne 3 ]; then
	 echo "Syntax $0 <graphs_dir> <alignments_dir> <out_dir>"
	 exit 1
fi


# vglr lrc_kir is a bad graph.  we can censor it as input in the wildcard arguments below
# to make it disappear from the analysis
GLOBIGNORE="*/lrc_kir/vglr*:*vglr-lrc_kir*"
#GLOBIGNORE="*vglr*"

# leave simons out for now as its hg19
GLOBIGNORE="*simons*:${GLOBIGNORE}"

GRAPHS=$1
ALIGNMENTS=$2
OUT_DIR=$3
TOIL_DIR=cs_toil_dir2
COMPS=( "sompy" )
INDEX_OPTS="--kmer 20 --edge_max 7 --timeout 10000"
COMP_OPTS="--clip --normalize --ignore Conflict --ignore Silver --vg_cores 10 --maxCores 40"
#COMP_OPTS="--clip  --ignore Conflict --ignore Silver --vg_cores 10 --maxCores 40"
COMP_TAG=comp_norm
REGIONS=( "brca2" "mhc" "brca1" "sma" "lrc_kir" )
#REGIONS=( "brca2" "mhc" )
OPTS="--maxCores 38 --vg_cores 2 --vg_only --skipBaseline --bubble"

# call variants, compute and plot baseline comparison
function run_pipeline {

	 VARIANTS_OUT_DIR=$1
	 CALL_OPTS=$2
	 PILEUP_OPTS=$3
	 FILTER_OPTS=$4
	 ROC=$5
	 mkdir  ${VARIANTS_OUT_DIR} ${VARIANTS_OUT_DIR}.${COMP_TAG}
	 
	 for i in "${REGIONS[@]}"
	 do
		  rm -rf ${TOIL_DIR}_test${TAG} ; scripts/callVariants.py ./${TOIL_DIR}_test${TAG} ${ALIGNMENTS}/${i}/*/NA12878.gam --graph_dir ${GRAPHS} --out_dir ${VARIANTS_OUT_DIR} ${OPTS} --call_opts "${CALL_OPTS}" --pileup_opts "${PILEUP_OPTS}" --filter_opts "${FILTER_OPTS}" 2> ${VARIANTS_OUT_DIR}/call_log_${i}.txt

		  for j in "${COMPS[@]}"
		  do
				# compute distances
				rm -rf ${TOIL_DIR}_testc${TAG} ; scripts/computeVariantsDistances.py ./${TOIL_DIR}_testc${TAG} ${ALIGNMENTS}/${i}/*/NA12878.gam ${VARIANTS_OUT_DIR} ${GRAPHS} ${j} ${VARIANTS_OUT_DIR}.${COMP_TAG} ${COMP_OPTS} ${INDEX_OPTS} ${ROC} 2> ${VARIANTS_OUT_DIR}.${COMP_TAG}/comp_log_${i}_${j}.txt
		  done

	 # plots
	 #scripts/plotVariantsDistances.py ${VARIANTS_OUT_DIR}.${COMP_TAG} &

	 done
}


# 1000 Genomes-like options fixed here for pileups
PILEUP_OPTS=" -w 40 -m 4 -q 10 "

mkdir ${OUT_DIR}
ident=0.90
delta=0.10
secscore=10000
depths=(01 05 10 15 20 25)
quals=(0 0.05 0.01 0.15 0.2 0.25)
for i in {0..5}
#for depth in 01
do
	 depth=${depths[${i}]}
	 qual=${quals[${i}]}
	 ROC_FLAG="--qpct ${qual}"
	 if [ "$depth" == "01" ]; then
		  ROC_FLAG="${ROC_FLAG} --roc"
	 fi
	 run_pipeline ${OUT_DIR}/primary_${depth}_i_${ident}_delt_${delta}_ss_${secscore} "-d ${depth}" "$PILEUP_OPTS" "-r ${ident} -e ${delta} -afu -s ${secscore} " "$ROC_FLAG"
done
wait
scripts/rocDistances.py ${OUT_DIR}/primary_*_i_${ident}_delt_${delta}_ss_${secscore}.${COMP_TAG} ${OUT_DIR}/primary_depth_roc_${COMP_TAG}
scripts/plotVariantsDistances.py ${OUT_DIR}/primary_depth_roc_${COMP_TAG}

secscores=(0.00 0.25 0.50 0.75 0.95 1.00)
quals=(0 0.05 0.01 0.15 0.2 0.25)
depth=10
for i in {0..5}
#for depth in 01
do
	 secscore=${secscores[${i}]}
	 qual=${quals[${i}]}
	 ROC_FLAG="--qpct ${qual}"
	 if [ "$depth" == "01" ]; then
		  ROC_FLAG="${ROC_FLAG} --roc"
	 fi
	 run_pipeline ${OUT_DIR}/primary_${depth}_i_${ident}_delt_${delta}_ss_${secscore} "-d ${depth}" "$PILEUP_OPTS" "-r ${ident} -e ${delta} -afu -s ${secscore} " "$ROC_FLAG"
done
wait
scripts/rocDistances.py ${OUT_DIR}/primary_${depth}_i_${ident}_delt_${delta}_ss_*.${COMP_TAG} ${OUT_DIR}/primary_depth_roc_${COMP_TAG}
scripts/plotVariantsDistances.py ${OUT_DIR}/primary_secscore_roc_${COMP_TAG}


# temp
exit 0

depth=10
for ident in 0.20 0.40 0.60 0.80 1.00
do
	 echo $ident
	 run_pipeline ${OUT_DIR}/primary_${depth}_i_${ident}_delt_${delta} "-d ${depth}" "$PILEUP_OPTS" "-r ${ident} -e ${delta} -afu -s 10000 "
done
scripts/rocDistances.py ${OUT_DIR}/primary_${depth}_i_*_delt_${delta}.${COMP_TAG} ${OUT_DIR}/primary_ident_roc
scripts/plotVariantsDistances.py ${OUT_DIR}/primary_ident_roc

ident=0.90
for delta in 0.00 0.05 0.15 0.25 0.50 0.75 0.99
do
	 echo $delta
	 run_pipeline ${OUT_DIR}/primary_${depth}_i_${ident}_delt_${delta} "-d ${depth}" "$PILEUP_OPTS" "-r ${ident} -e ${delta} -afu -s 10000 "
done
scripts/rocDistances.py ${OUT_DIR}/primary_${depth}_i_${ident}_delt_*.${COMP_TAG} ${OUT_DIR}/primary_delta_roc
scripts/plotVariantsDistances.py ${OUT_DIR}/primary_delta_roc
