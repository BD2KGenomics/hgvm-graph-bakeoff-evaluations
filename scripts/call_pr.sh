#!/bin/bash

# create precision-recall plots for for vg variants calling using
# gold standard vcf.  Directory structure
# of input is important, and explained in ../README.md and callVariants.py
# running is grouped by region to slightly facilitate debugging and resuming

# wildcards used to work for samples, activate by wrapping in quotes

if [ "$#" -ne 7 ]; then
	 echo "Syntax $0 <toil_dir> <graphs_dir> <alignments_dir> <out_dir> <sample> <clip> <qual>"
	 exit 1
fi

if [[ "$6" -ne 0 && "$6" -ne 1 ]]; then
	 echo "clip must be 0 or 1"
	 exit 1
fi

# vglr lrc_kir is a bad graph.  we can censor it as input in the wildcard arguments below
# to make it disappear from the analysis
GLOBIGNORE="*/lrc_kir/vglr*:*vglr-lrc_kir*"
# leave out debruijn mhc for now as it makes invalid augmented graph
GLOBIGNORE="*/mhc/debruijn*:debruijn*mhc*:${GLOBIGNORE}"
GLOBIGNORE="*/brca2/debruijn*:debruijn*brca2*:${GLOBIGNORE}"
GLOBIGNORE="*/lrc_kir/debruijn*:debruijn*lrc_kir*:${GLOBIGNORE}"
#ditto camel -- vcf conversion doesn't work for crazy regions (probably due to caller making bizarre graph?)
GLOBIGNORE="*/mhc/camel*:camel*mhc*:${GLOBIGNORE}"
GLOBIGNORE="*/lrc_kir/camel*:camel*lrc_kir*:${GLOBIGNORE}"
# leave simons out for now as its hg19
GLOBIGNORE="*simons*:${GLOBIGNORE}"
GLOBIGNORE=$GLOBIGNORE
# and trivial, since its same as refonly
GLOBIGNORE="*trivial*:${GLOBIGNORE}"
GLOBIGNORE=$GLOBIGNORE

# leave gatk3 out of comparison because it's misleading
GLOBIGNORE_COMP="*gatk*:*simons*"


#temp leave out prg because happy takes forever on it
#GLOBIGNORE_COMP="*prg*:${GLOBIGNORE}"

# command line args
TOIL_DIR=$1
GRAPHS=$2
ALIGNMENTS=$3
OUT_DIR=$4
SAMPLE=$5
CLIP=$6
qual=$7

# VCF comparisons : sompy happy vcf vcfeval
#COMPS=( "sompy" "happy" "vcf" )
# graph comparison : kmer corg
#COMPS=( "kmer" "corg")
COMPS=( "vcfeval" )

# Parameters for indexing (used for corg / kmer comparison)
INDEX_OPTS="--kmer 20 --edge_max 5 --timeout 5000"

# Calling parameters
REGIONS=( "brca1" "brca2" "sma" "lrc_kir" "mhc" )
#REGIONS=( "brca2" "brca1" "sma" "lrc_kir" )
#REGIONS=( "brca2" "sma" "lrc_kir" )
#REGIONS=( "lrc_kir" )

#OPTS="--maxCores 8 --vg_cores 4 --vg_only --skipBaseline --genotype"
OPTS="--maxCores 36 --vg_cores 6 --vg_only --skipBaseline"

#hack in truth alignments as test
#TA_PATH="/cluster/home/anovak/hive/ga4gh/bake-off/hgvm-graph-bakeoff-evalutations/platinum_truth_alignments/alignments"
TA_PATH="platinum_truth_alignments"
# Comparison parameters

#CLIP_PATH="data/filters/platinum2016.bed"
CLIP_PATH="data/filters/platinum.bed"
#CLIP_PATH="data/filters/homo/homo_44_platinum.bed"
PLAT_PATH="data/platinum"

# Normalization (requires vt, recommended for som.py and vcf comparison)
#COMP_OPTS="--clip data/filters/platinum.bed --normalize --ignore Conflict --ignore Silver --vg_cores 10 --maxCores 30"

# this controls what part of vg-call vcf's get used for the roc.  will get passed as --${QF_TYPE} to vcfFilterQuality.py (which has really hacky hardcoded support for a few fields i've tried)
QF_TYPE="ad"
# see also --dedupe --genotype --vroc --clip 
COMP_OPTS="--vg_cores 6 --maxCores 30 --gt --freebayes_path platinum_classic3/freebayes --platypus_path platinum_classic3/platypus --samtools_path platinum_classic3/samtools --gatk3_path null --g1kvcf_path null --filter_type ${QF_TYPE} --normalize --platinum_path ${PLAT_PATH}"

# example how to use user-specified platypus and freebayes vcfs
#COMP_OPTS="--clip data/filters/platinum.bed --normalize --ignore Conflict --ignore Silver --vg_cores 10 --maxCores 30 --platypus_path platinum_classic/platypus --freebayes_path platinum_classic/freebayes"

COMP_TAG=comp.gt.${QF_TYPE}.norms.myroc.truth3.rtg

if [ "$CLIP" = 1 ]; then
	 COMP_OPTS="$COMP_OPTS --clip ${CLIP_PATH}"
	 COMP_TAG=${COMP_TAG}.clip
fi

# call variants, compute and plot baseline comparison
function run_pipeline {

	 local VARIANTS_OUT_DIR=$1
	 local COMP_OUT_DIR=$2
	 local CALL_OPTS=$3
	 local PILEUP_OPTS=$4
	 local FILTER_OPTS=$5
	 local ROC=$6
	 local RUN_CALLER=$7
	 	 
	 for i in "${REGIONS[@]}"
	 do
		  if [ "$RUN_CALLER" = true ]; then
				GLOBIGNORE=$GLOBIGNORE_CALL
				mkdir ${VARIANTS_OUT_DIR}				
				#rm -rf ${TOIL_DIR}_test${TAG} ; scripts/callVariants.py ./${TOIL_DIR}_test${TAG}  ${ALIGNMENTS}/${i}/*/${SAMPLE}.gam ${TA_PATH}/${i}/*/${SAMPLE}.gam --graph_dir ${GRAPHS} --out_dir ${VARIANTS_OUT_DIR} ${OPTS} --call_opts "${CALL_OPTS}" --pileup_opts "${PILEUP_OPTS}" --filter_opts "${FILTER_OPTS}" 2>> ${VARIANTS_OUT_DIR}/call_log_${i}.txt
				rm -rf ${TOIL_DIR}_test${TAG} ; scripts/callVariants.py ./${TOIL_DIR}_test${TAG}  ${ALIGNMENTS}/${i}/refonly/${SAMPLE}.gam ${ALIGNMENTS}/${i}/snp1kg/${SAMPLE}.gam ${ALIGNMENTS}/${i}/cactus/${SAMPLE}.gam  ${TA_PATH}/${i}/*/${SAMPLE}.gam --graph_dir ${GRAPHS} --out_dir ${VARIANTS_OUT_DIR} ${OPTS} --call_opts "${CALL_OPTS}" --pileup_opts "${PILEUP_OPTS}" --filter_opts "${FILTER_OPTS}" 2>> ${VARIANTS_OUT_DIR}/call_log_${i}.txt
		  fi
		  for j in "${COMPS[@]}"
		  do
				GLOBIGNORE=$GLOBIGNORE_COMP				
				# compute distances
				mkdir ${COMP_OUT_DIR}
				rm -rf ${TOIL_DIR}_testc${TAG} ; scripts/computeVariantsDistances.py ./${TOIL_DIR}_testc${TAG} ${ALIGNMENTS}/${i}/*/${SAMPLE}.gam ${TA_PATH}/${i}/*/${SAMPLE}.gam ${VARIANTS_OUT_DIR} ${GRAPHS} ${j} ${COMP_OUT_DIR} ${COMP_OPTS} ${INDEX_OPTS} ${ROC}  2>> ${COMP_OUT_DIR}/comp_log_${i}_${j}.txt
				GLOBIGNORE=$GLOBIGNORE_CALL
		  done

	 done
}

# 1000 Genomes-like options fixed here for pileups
PILEUP_OPTS=" -w 40 -m 10 -q 10 -a"
#GLENN2VCF_OPTS="--depth 10 --max_het_bias 3 --min_count 1 -C 50"
#GLENN2VCF_OPTS="--max_het_bias 3"
# filter options (no secondary, primary must have > 90% identity and > 5% improvement over secondary)
mkdir ${OUT_DIR}
ident=0.90
delta=0.00
secscore=10000

ROC_FLAG="--qpct ${qual} --roc --qgraph"
DO_CALL=true
CALL_OPTS="-b 5 -s 1 -d 1 -f 0 -D 10 -H 3 -n 1 -F 0.2 -B 250 -R 4"
#CALL_OPTS="-d 0 -e 5000 -s 3 -D 20 -n 0 -F 0.2 -B 250 -H 3 -R 4"
FILTER_OPTS="-r ${ident} -d ${delta} -e ${delta} -afu -s ${secscore} -q 15 -o 0 -E 4"
VAR_OUT_DIR=${OUT_DIR}/primary_call_i_${ident}_delt_${delta}_ss_${secscore}
COMP_OUT_DIR=${OUT_DIR}/primary_${qual}_i_${ident}_delt_${delta}_ss_${secscore}.${COMP_TAG}
run_pipeline $VAR_OUT_DIR $COMP_OUT_DIR "$CALL_OPTS" "$PILEUP_OPTS" "$FILTER_OPTS" "$ROC_FLAG" $DO_CALL

# scrape together all results into one tsv / region / comp type
scripts/rocDistances.py ${OUT_DIR}/primary_${qual}_i_${ident}_delt_${delta}_ss_${secscore}.${COMP_TAG} ${OUT_DIR}/pr_plots.${qual}.${COMP_TAG} --best_comp happy

# finally, draw out the tables created above
scripts/plotVariantsDistances.py ${OUT_DIR}/pr_plots.${qual}.${COMP_TAG}

exit 0

