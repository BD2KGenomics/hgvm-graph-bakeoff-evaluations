#!/bin/bash

# create precision-recall plots for for vg variants calling using
# gold standard vcf.  Directory structure
# of input is important, and explained in ../README.md and callVariants.py
# running is grouped by region to slightly facilitate debugging and resuming

# wildcards used to work for samples, activate by wrapping in quotes

if [ "$#" -ne 8 ]; then
	 echo "Syntax $0 <toil_dir> <graphs_dir> <alignments_dir> <out_dir> <sample> <clip> <qual> <gt>"
	 exit 1
fi

if [[ "$6" -ne 0 && "$6" -ne 1 ]]; then
	 echo "clip must be 0 or 1"
	 exit 1
fi

if [[ "$8" -ne 0 && "$8" -ne 1 ]]; then
	 echo "gt must be 0 or 1"
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
# Drop all the extra 1kg variants
GLOBIGNORE="*haplo1kg*_af*:${GLOBIGNORE}"
GLOBIGNORE="*snp1kg*_af*:${GLOBIGNORE}"
GLOBIGNORE="*snp1kg_kp*:${GLOBIGNORE}"
GLOBIGNORE="*snp1kg_norm*:${GLOBIGNORE}"
GLOBIGNORE="*snp1kg_plat*:${GLOBIGNORE}"
# Expand globs
GLOBIGNORE=$GLOBIGNORE

GLOBIGNORE_CALL=$GLOBIGNORE


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
GTCOMP=$8

# VCF comparisons : sompy happy vcf vcfeval
#COMPS=( "sompy" "happy" "vcf" )
# graph comparison : kmer corg
#COMPS=( "kmer" "corg")
COMPS=( "vcfeval" )

# Parameters for indexing (used for corg / kmer comparison)
INDEX_OPTS="--kmer 20 --edge_max 5 --timeout 5000"

# Calling parameters
REGIONS=( "brca1" "brca2" "sma" "lrc_kir" "mhc" )
#REGIONS=( "brca2" "brca1" )
#REGIONS=( "brca2" "sma" "lrc_kir" )
#REGIONS=( "lrc_kir" )

#OPTS="--maxCores 8 --vg_cores 4 --vg_only --skipBaseline --genotype"
OPTS="--maxCores 36 --vg_cores 4 --vg_only --skipBaseline"

# Comparison parameters

#CLIP_PATH="data/filters/platinum2016.bed"
CLIP_PATH="data/filters/platinum.bed"
#CLIP_PATH="data/filters/homo/homo_44_platinum.bed"
PLAT_PATH="data/platinum"

CLASSIC_PATH=/cluster/home/hickey/ga4gh/hgvm-graph-bakeoff-evalutations/platinum_classic
#CLASSIC_PATH=/cluster/home/hickey/ga4gh/hgvm-graph-bakeoff-evalutations/platinum_classic_primary

# this controls what part of vg-call vcf's get used for the roc.  will get passed as --${QF_TYPE} to vcfFilterQuality.py (which has really hacky hardcoded support for a few fields i've tried)
QF_TYPE="ad"
# see also --dedupe --genotype --vroc --clip 

COMP_OPTS="--vg_cores 6 --maxCores 30 --freebayes_path ${CLASSIC_PATH}/freebayes --platypus_path ${CLASSIC_PATH}/platypus --samtools_path ${CLASSIC_PATH}/samtools --gatk3_path null --g1kvcf_path null --filter_type ${QF_TYPE} --normalize --platinum_path ${PLAT_PATH} --combine_samples NA12877,NA12878"

STATS_OPTS=""

COMP_TAG=comp.${QF_TYPE}.norm.myroc

if [ "$GTCOMP" = 1 ]; then
	 COMP_OPTS="$COMP_OPTS --gt"
	 COMP_TAG=${COMP_TAG}.gt
fi

if [ "$CLIP" = 1 ]; then
	 COMP_OPTS="$COMP_OPTS --clip ${CLIP_PATH}"
	 COMP_TAG=${COMP_TAG}.clip
	 STATS_OPTS="$STATS_OPTS --clip ${CLIP_PATH}"
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
				rm -rf ${TOIL_DIR}_test${TAG} ; scripts/callVariants.py ./${TOIL_DIR}_test${TAG}  ${ALIGNMENTS}/${i}/*/${SAMPLE}.gam  --graph_dir ${GRAPHS} --out_dir ${VARIANTS_OUT_DIR} ${OPTS} --call_opts "${CALL_OPTS}" --pileup_opts "${PILEUP_OPTS}" --filter_opts "${FILTER_OPTS}" 2>> ${VARIANTS_OUT_DIR}/call_log_${i}.txt
		  fi
	 done

	 # make a total directory (one vcf for all regions)
	 scripts/callTotals.py ${VARIANTS_OUT_DIR} 2> ${VARIANTS_OUT_DIR}/total_log.txt

	 for j in "${COMPS[@]}"
	 do
		  GLOBIGNORE=$GLOBIGNORE_COMP				
		  # compute distances
		  mkdir ${COMP_OUT_DIR}
		  rm -rf ${TOIL_DIR}_testc${TAG} ; scripts/computeVariantsDistances.py ./${TOIL_DIR}_testc${TAG} ${ALIGNMENTS}/*/*/${SAMPLE}.gam ${TA_PATH}/*/*/${SAMPLE}.gam ${EXP_PATH}/*/*/${SAMPLE}.gam ${VARIANTS_OUT_DIR} ${GRAPHS} ${j} ${COMP_OUT_DIR} ${COMP_OPTS} ${INDEX_OPTS} ${ROC}  2>> ${COMP_OUT_DIR}/comp_log_${j}.txt
		  GLOBIGNORE=$GLOBIGNORE_CALL
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
CALL_OPTS=""  # Use the defaults that Glenn coded in as the best parameter set
#CALL_OPTS="-d 0 -e 5000 -s 3 -D 20 -n 0 -F 0.2 -B 250 -H 3 -R 4"
# Defray all the way (as reads are << 999 bases long)
FILTER_OPTS="-r ${ident} -d ${delta} -e ${delta} -afu -s ${secscore} -q 15 -o 0 --defray-ends 999"
VAR_OUT_DIR=${OUT_DIR}/primary_call_i_${ident}_delt_${delta}_ss_${secscore}
COMP_OUT_DIR=${OUT_DIR}/primary_${qual}_i_${ident}_delt_${delta}_ss_${secscore}.${COMP_TAG}
run_pipeline $VAR_OUT_DIR $COMP_OUT_DIR "$CALL_OPTS" "$PILEUP_OPTS" "$FILTER_OPTS" "$ROC_FLAG" $DO_CALL

# scrape together all results into one tsv / region / comp type
# note this script doesn't do much of anythign now that we rely on vcfeval to make rocs
# but we keep old logic of copying stuff to pr_directories
scripts/rocDistances.py ${OUT_DIR}/primary_${qual}_i_${ident}_delt_${delta}_ss_${secscore}.${COMP_TAG} ${OUT_DIR}/pr_plots.${qual}.${COMP_TAG} --best_comp happy 2>> ${OUT_DIR}/roc_cp.log


echo ""
echo "scripts/callStats.py ${OUT_DIR}/primary_${qual}_i_${ident}_delt_${delta}_ss_${secscore}.${COMP_TAG} ${OUT_DIR}/call_stats.${qual}.${COMP_TAG} ${STATS_OPTS} 2>> ${OUT_DIR}/callstats.log"
echo "scripts/trioStats.py ${OUT_DIR}/primary_${qual}_i_${ident}_delt_${delta}_ss_${secscore}.${COMP_TAG} ${OUT_DIR}/trio_stats.${qual}.${COMP_TAG} ${STATS_OPTS} 2>> ${OUT_DIR}/triostats.log"
#echo "scripts/plotVariantsDistances.py ${OUT_DIR}/pr_plots.${qual}.${COMP_TAG} --top 2>> ${OUT_DIR}/plots.log "

#scripts/callStats.py ${OUT_DIR}/primary_${qual}_i_${ident}_delt_${delta}_ss_${secscore}.${COMP_TAG} ${OUT_DIR}/call_stats.${qual}.${COMP_TAG} ${STATS_OPTS} 2>> ${OUT_DIR}/callstats.log &
# note: triostats.py only runs if there's a dummy platinum directory for NA12879 and the comparison was run on it too
# todo: fix (note the dummy comparison is never actually used, it's just the vcfs never get preprocessed without it)
#scripts/trioStats.py ${OUT_DIR}/primary_${qual}_i_${ident}_delt_${delta}_ss_${secscore}.${COMP_TAG} ${OUT_DIR}/trio_stats.${qual}.${COMP_TAG} ${STATS_OPTS} 2>> ${OUT_DIR}/triostats.log &

#scripts/plotVariantsDistances.py ${OUT_DIR}/pr_plots.${qual}.${COMP_TAG} --top >> ${OUT_DIR}/plots.log &
scripts/plotVariantsDistances.py ${OUT_DIR}/pr_plots.${qual}.${COMP_TAG} >> ${OUT_DIR}/plots.log 


exit 0

