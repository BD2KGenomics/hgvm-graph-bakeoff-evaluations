#!/bin/bash

# Draw a PDF linear and log heatmap for each matrix (tsv) created with compare_snps.sh
# Should be run with exact same arguments as compare_snps.sh

if [ "$#" -ne 4 ]; then
	 echo "Syntax $0 <graphs_dir> <alignments_dir> <calls_dir> <out_dir>"
	 exit 1
fi

GRAPHS=$1
ALIGNMENTS=$2
VARIANTS=$3
OUT_DIR=$4
COMPS=( "jaccard" "recall" "corg" )
REGIONS=( "brca1" "brca2" "sma" "lrc_kir" "mhc" )
#COMPS=( "jaccard" )
#REGIONS=( "brca1" )

for i in "${REGIONS[@]}"
do
	 for j in "${COMPS[@]}"
	 do
		  # all vs all original graphs
		  rm ${OUT_DIR}/${i}/comp_orig_${j}_heatmap_${i}.pdf
        rm ${OUT_DIR}/${i}/comp_orig_${j}_heatmap_${i}.pdf		  
		  scripts/heatmap.py ${OUT_DIR}/${i}/comp_orig_${j}_distmat_${i}.tsv ${OUT_DIR}/${i}/comp_orig_${j}_heatmap_${i}.pdf
		  rm ${OUT_DIR}/${i}/comp_orig_${j}_heatmap_log_${i}.pdf
		  scripts/heatmap.py ${OUT_DIR}/${i}/comp_orig_${j}_distmat_${i}.tsv ${OUT_DIR}/${i}/comp_orig_${j}_heatmap_log_${i}.pdf --log_scale

		  # all vs all sample graphs
		  rm ${OUT_DIR}/${i}/comp_sample_${j}_heatmap_${i}.pdf
		  scripts/heatmap.py ${OUT_DIR}/${i}/comp_sample_${j}_distmat_${i}.tsv ${OUT_DIR}/${i}/comp_sample_${j}_heatmap_${i}.pdf
		  rm ${OUT_DIR}/${i}/comp_sample_${j}_heatmap_log_${i}.pdf
		  scripts/heatmap.py ${OUT_DIR}/${i}/comp_sample_${j}_distmat_${i}.tsv ${OUT_DIR}/${i}/comp_sample_${j}_heatmap_log_${i}.pdf --log_scale

		  # all vs all original and sample graphs
		  rm ${OUT_DIR}/${i}/comp_sample_and_orig_${j}_heatmap_${i}.pdf
		  scripts/heatmap.py ${OUT_DIR}/${i}/comp_sample_and_orig_${j}_distmat_${i}.tsv ${OUT_DIR}/${i}/comp_sample_and_orig_${j}_heatmap_${i}.pdf
		  rm ${OUT_DIR}/${i}/comp_sample_and_orig_${j}_heatmap_log_${i}.pdf
		  scripts/heatmap.py ${OUT_DIR}/${i}/comp_sample_and_orig_${j}_distmat_${i}.tsv ${OUT_DIR}/${i}/comp_sample_and_orig_${j}_heatmap_log_${i}.pdf --log_scale
	 done
	 	 	 
done
