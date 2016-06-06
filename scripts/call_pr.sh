#!/bin/bash

# create precision-recall plots for for vg variants calling using
# gold standard vcf.  Directory structure
# of input is important, and explained in ../README.md and callVariants.py
# running is grouped by region to slightly facilitate debugging and resuming

# qgraph is 0/1 flag.  0: old way running vg call with different depths. 1: new way:
# run vg call once and use vcfQalutiyFilter for plot

# wildcards used to work for samples, activate by wrapping in quotes
SAMPLE="CHM1"

if [ "$#" -ne 5 ]; then
	 echo "Syntax $0 <graphs_dir> <alignments_dir> <out_dir> <sample> <qgraph>"
	 exit 1
fi

if [[ "$4" -ne 0 && "$4" -ne 1 ]]; then
	 echo "qgraph must be 0 or 1"
	 exit 1
fi

# vglr lrc_kir is a bad graph.  we can censor it as input in the wildcard arguments below
# to make it disappear from the analysis
GLOBIGNORE="*/lrc_kir/vglr*:*vglr-lrc_kir*"
# leave out debruijn mhc for now as it makes invalid augmented graph
GLOBIGNORE="*/mhc/debruijn*:debruijn*mhc*:${GLOBIGNORE}"
GLOBIGNORE="*/lrc_kir/debruijn*:debruijn*lrc_kir*:${GLOBIGNORE}"
#ditto camel -- vcf conversion doesn't work for crazy regions (probably due to caller making bizarre graph?)
GLOBIGNORE="*/mhc/camel*:camel*mhc*:${GLOBIGNORE}"
GLOBIGNORE="*/lrc_kir/camel*:camel*lrc_kir*:${GLOBIGNORE}"
# leave simons out for now as its hg19
GLOBIGNORE_COMP="*simons*:${GLOBIGNORE}"
GLOBIGNORE_CALL=$GLOBIGNORE

#temp leave out prg because happy takes forever on it
#GLOBIGNORE_COMP="*prg*:${GLOBIGNORE}"

# command line args
GRAPHS=$1
ALIGNMENTS=$2
OUT_DIR=$3
SAMPLE=$4
QGRAPH=$5

TOIL_DIR=call_pr_toil5

# VCF comparisons : sompy happy vcf vcfeval
#COMPS=( "sompy" "happy" "vcf" )
# graph comparison : kmer corg
#COMPS=( "kmer" "corg")
COMPS=( "vcfeval" )

# Parameters for indexing (used for corg / kmer comparison)
INDEX_OPTS="--kmer 20 --edge_max 5 --timeout 5000"

# Calling parameters
REGIONS=( "brca2" "mhc" "brca1" "sma" "lrc_kir" )
#REGIONS=( "brca2" "brca1" "sma" "lrc_kir" )
REGIONS=( "brca2" )

#OPTS="--maxCores 8 --vg_cores 4 --vg_only --skipBaseline"
OPTS="--maxCores 20 --vg_cores 4 --vg_only --skipBaseline"

# Comparison parameters


# Normalization (requires vt, recommended for som.py and vcf comparison)
#COMP_OPTS="--clip data/filters/platinum.bed --normalize --ignore Conflict --ignore Silver --vg_cores 10 --maxCores 30"

# No normalization (only recommended for vcfeval and hap.py)
COMP_OPTS="--vg_cores 4 --maxCores 20 --gt --freebayes_path platinum_classic3/freebayes --platypus_path platinum_classic3/platypus --samtools_path platinum_classic3/samtools"

# example how to use user-specified platypus and freebayes vcfs
#COMP_OPTS="--clip data/filters/platinum.bed --normalize --ignore Conflict --ignore Silver --vg_cores 10 --maxCores 30 --platypus_path platinum_classic/platypus --freebayes_path platinum_classic/freebayes"

# make sure qgraph option gets into compare script
if [ "$QGRAPH" = 1 ]; then
	 COMP_OPTS="${COMP_OPTS} --qgraph"
fi

COMP_TAG=comp.gt.noclip

# call variants, compute and plot baseline comparison
function run_pipeline {

	 local VARIANTS_OUT_DIR=$1
	 local COMP_OUT_DIR=$2
	 local CALL_OPTS=$3
	 local PILEUP_OPTS=$4
	 local FILTER_OPTS=$5
	 local G2V_OPTS=$6
	 local ROC=$7
	 local RUN_CALLER=$8
	 	 
	 for i in "${REGIONS[@]}"
	 do
		  if [ "$RUN_CALLER" = true ]; then
				GLOBIGNORE=$GLOBIGNORE_CALL
				mkdir ${VARIANTS_OUT_DIR}				
				rm -rf ${TOIL_DIR}_test${TAG} ; scripts/callVariants.py ./${TOIL_DIR}_test${TAG}  ${ALIGNMENTS}/${i}/*/${SAMPLE}.gam --graph_dir ${GRAPHS} --out_dir ${VARIANTS_OUT_DIR} ${OPTS} --call_opts "${CALL_OPTS}" --pileup_opts "${PILEUP_OPTS}" --filter_opts "${FILTER_OPTS}" 2>> ${VARIANTS_OUT_DIR}/call_log_${i}.txt
		  fi
		  for j in "${COMPS[@]}"
		  do
				GLOBIGNORE=$GLOBIGNORE_COMP				
				# compute distances
				mkdir ${COMP_OUT_DIR}
				rm -rf ${TOIL_DIR}_testc${TAG} ; scripts/computeVariantsDistances.py ./${TOIL_DIR}_testc${TAG} ${ALIGNMENTS}/${i}/*/${SAMPLE}.gam ${VARIANTS_OUT_DIR} ${GRAPHS} ${j} ${COMP_OUT_DIR} ${COMP_OPTS} ${INDEX_OPTS} ${ROC}  2>> ${COMP_OUT_DIR}/comp_log_${i}_${j}.txt
				GLOBIGNORE=$GLOBIGNORE_CALL
		  done

	 done
}

# 1000 Genomes-like options fixed here for pileups
PILEUP_OPTS=" -w 20 -m 4 -q 10 "
GLENN2VCF_OPTS="--depth 10 --max_het_bias 4.2 --min_count 6 --min_fraction 0.15"
# filter options (no secondary, primary must have > 90% identity and > 5% improvement over secondary)
mkdir ${OUT_DIR}
ident=0.90
delta=0.00
secscore=10000

# points on the precision-precall chart hard-coded here. 
depths=(01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 22 24 26 28 30 35 40 45)
supports=(01 02 03 04 05 06 07 08 09 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 15 15 15 15)
quals=(0.000 0.020 0.040 0.060 0.080 0.010 0.100 0.120 0.140 0.160 0.180 0.200 0.220 0.240 0.260 0.280 0.300 0.320 0.340 0.360 0.380 0.400 0.500 0.600 0.700 0.800 0.900 1.000)
#for i in {0..27}
#for i in 1 2 3 4 5 6 7 8 1 10 20
points=( 0 1 2 5 7 10 12 15 20 21 23 25 27 )
#points=( 0 )
for i in "${points[@]}" 
do
	 depth=${depths[${i}]}
	 qual=${quals[${i}]}
	 support=${supports[${i}]}
	 ROC_FLAG="--qpct ${qual}"
	 if [ "$depth" == "01" ]; then
		  ROC_FLAG="${ROC_FLAG} --roc"
	 fi
	 if [[ "$i" = "${points[0]}" || "$QGRAPH" = 0 ]]; then
		  DO_CALL=true
		  calldepth=$depth
	 else
		  DO_CALL=false
	 fi
	 CALL_OPTS="-r 0.0001 -b 1 -f 0.05 -d 4"
	 FILTER_OPTS="-r ${ident} -d ${delta} -e ${delta} -afu -s ${secscore} -o 0 "
	 VAR_OUT_DIR=${OUT_DIR}/primary_${calldepth}_i_${ident}_delt_${delta}_ss_${secscore}
	 COMP_OUT_DIR=${OUT_DIR}/primary_${depth}_i_${ident}_delt_${delta}_ss_${secscore}.${COMP_TAG}
	 run_pipeline $VAR_OUT_DIR $COMP_OUT_DIR "$CALL_OPTS" "$PILEUP_OPTS" "$FILTER_OPTS" "$GLENN2VCF_OPTS" "$ROC_FLAG" $DO_CALL
done
wait

# scrape together all results into one tsv / region / comp type
scripts/rocDistances.py ${OUT_DIR}/primary_*_i_${ident}_delt_${delta}_ss_${secscore}.${COMP_TAG} ${OUT_DIR}/pr_plots.${COMP_TAG} --best_comp happy

# finally, draw out the tables created above
scripts/plotVariantsDistances.py ${OUT_DIR}/pr_plots.${COMP_TAG}

exit 0

