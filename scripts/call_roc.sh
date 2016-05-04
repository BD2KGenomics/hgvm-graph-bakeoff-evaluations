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
INDEX_OPTS="--kmer 27 --edge_max 5 --timeout 5000"
COMP_OPTS="--clip --normalize --ignore Conflict --ignore Silver --vg_cores 10 --maxCores 40"
#COMP_OPTS="--clip  --ignore Conflict --ignore Silver --vg_cores 10 --maxCores 40"
COMP_TAG=comp_norm_hp
REGIONS=( "brca2" "mhc" "brca1" "sma" "lrc_kir" )
#REGIONS=( "brca1" "brca2" )
OPTS="--maxCores 38 --vg_cores 2 --vg_only --bubble"

# call variants, compute and plot baseline comparison
function run_pipeline {

	 VARIANTS_OUT_DIR=$1
	 CALL_OPTS=$2
	 PILEUP_OPTS=$3
	 FILTER_OPTS=$4
	 ROC=$5
	 	 
	 for i in "${REGIONS[@]}"
	 do
		  mkdir ${VARIANTS_OUT_DIR}
		  rm -rf ${TOIL_DIR}_test${TAG} ; scripts/callVariants.py ./${TOIL_DIR}_test${TAG} ${ALIGNMENTS}/${i}/*/*.gam --graph_dir ${GRAPHS} --out_dir ${VARIANTS_OUT_DIR} ${OPTS} --call_opts "${CALL_OPTS}" --pileup_opts "${PILEUP_OPTS}" --filter_opts "${FILTER_OPTS}" 2>> ${VARIANTS_OUT_DIR}/call_log_${i}.txt

		  for j in "${COMPS[@]}"
		  do
				for k in strict conf
				do
					 # compute distances
					 mkdir ${VARIANTS_OUT_DIR}.${COMP_TAG}.${k}
					 rm -rf ${TOIL_DIR}_testc${TAG} ; scripts/computeVariantsDistances.py ./${TOIL_DIR}_testc${TAG} ${ALIGNMENTS}/${i}/*/NA12878.gam ${VARIANTS_OUT_DIR} ${GRAPHS} ${j} ${VARIANTS_OUT_DIR}.${COMP_TAG}.${k} ${COMP_OPTS} ${INDEX_OPTS} ${ROC} --clip_suffix _${k} 2>> ${VARIANTS_OUT_DIR}.${COMP_TAG}.${k}/comp_log_${i}_${j}.txt
					 rm -rf ${TOIL_DIR}_testc${TAG} ; scripts/computeVariantsDistances.py ./${TOIL_DIR}_testc${TAG} ${ALIGNMENTS}/${i}/*/NA12878.gam ${VARIANTS_OUT_DIR} ${GRAPHS} ${j} ${VARIANTS_OUT_DIR}.${COMP_TAG}.${k} ${COMP_OPTS} ${INDEX_OPTS} ${ROC} --clip_suffix _${k} --baseline g1kvcf 2>> ${VARIANTS_OUT_DIR}.${COMP_TAG}.${k}/comp_log_${i}_${j}.txt

				done
		  done
	 # plots
	 #scripts/plotVariantsDistances.py ${VARIANTS_OUT_DIR}.${COMP_TAG} &

	 done
}

# 1000 Genomes-like options fixed here for pileups
#PILEUP_OPTS=" -w 40 -m 4 -q 10 "
PILEUP_OPTS=" -q 10 "

mkdir ${OUT_DIR}
ident=0.90
delta=0.10
secscore=10000
depths=(02 04 06 08 10 12 14 16 18 20 25 30)
quals=(0.00 0.05 0.10 0.15 0.20 0.25 0.40 0.55 0.70 0.85 0.90 0.99)
for i in {0..11}
#for depth in 01
do
	 depth=${depths[${i}]}
	 qual=${quals[${i}]}
	 ROC_FLAG="--qpct ${qual}"
	 if [ "$depth" == "01" ]; then
		  ROC_FLAG="${ROC_FLAG} --roc"
	 fi
	 run_pipeline ${OUT_DIR}/primary_${depth}_i_${ident}_delt_${delta}_ss_${secscore} "-r 0.0001 -b 0.4 -f 0.25 -d ${depth}" "$PILEUP_OPTS" "-r ${ident} -d ${delta} -e ${delta} -afu -s ${secscore} -o 15 " "$ROC_FLAG"
done
wait
for k in strict conf
do
	 scripts/rocDistances.py ${OUT_DIR}/primary_*_i_${ident}_delt_${delta}_ss_${secscore}.${COMP_TAG}.${k} ${OUT_DIR}/primary_depth_roc_${COMP_TAG}.${k}
	 scripts/plotVariantsDistances.py ${OUT_DIR}/primary_depth_roc_${COMP_TAG}.${k}
done

#temp
exit 0















secscores=(0.00 0.25 0.50 0.75 0.95 1.00)
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
scripts/rocDistances.py ${OUT_DIR}/primary_${depth}_i_${ident}_delt_${delta}_ss_*.${COMP_TAG} ${OUT_DIR}/primary_secscore_roc_${COMP_TAG}
scripts/plotVariantsDistances.py ${OUT_DIR}/primary_secscore_roc_${COMP_TAG}


idents=(0.10 0.20 0.40 0.60 0.80 1.00)
secscore=10000
for i in {0..5}
#for depth in 01
do
	 ident=${idents[${i}]}
	 qual=${quals[${i}]}
	 ROC_FLAG="--qpct ${qual}"
	 if [ "$depth" == "01" ]; then
		  ROC_FLAG="${ROC_FLAG} --roc"
	 fi
	 run_pipeline ${OUT_DIR}/primary_${depth}_i_${ident}_delt_${delta}_ss_${secscore} "-d ${depth}" "$PILEUP_OPTS" "-r ${ident} -e ${delta} -afu -s ${secscore} " "$ROC_FLAG"
done
wait
scripts/rocDistances.py ${OUT_DIR}/primary_${depth}_i_${ident}_delt_${delta}_ss_*.${COMP_TAG} ${OUT_DIR}/primary_ident_roc_${COMP_TAG}
scripts/plotVariantsDistances.py ${OUT_DIR}/primary_ident_roc_${COMP_TAG}


ident=0.90
deltas=(0.00 0.10 0.25 0.50 0.75 0.99)
for i in {0..5}
#for depth in 01
do
	 delta=${deltas[${i}]}
	 qual=${quals[${i}]}
	 ROC_FLAG="--qpct ${qual}"
	 if [ "$depth" == "01" ]; then
		  ROC_FLAG="${ROC_FLAG} --roc"
	 fi
	 run_pipeline ${OUT_DIR}/primary_${depth}_i_${ident}_delt_${delta}_ss_${secscore} "-d ${depth}" "$PILEUP_OPTS" "-r ${ident} -e ${delta} -afu -s ${secscore} " "$ROC_FLAG"
done
wait
scripts/rocDistances.py ${OUT_DIR}/primary_${depth}_i_${ident}_delt_${delta}_ss_*.${COMP_TAG} ${OUT_DIR}/primary_delta_roc_${COMP_TAG}
scripts/plotVariantsDistances.py ${OUT_DIR}/primary_delta_roc_${COMP_TAG}
