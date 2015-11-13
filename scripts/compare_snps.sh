#!/bin/bash

if [ "$#" -ne 4 ]; then
	 echo "Syntax $0 <graphs_dir> <alignments_dir> <calls_dir> <out_dir>"
	 exit 1
fi

GRAPHS=$1
ALIGNMENTS=$2
VARIANTS=$3
OUT_DIR=$4
TOIL_DIR=cps_toil_dir
TOIL_OPTS="--maxCores 48 --vg_cores 8"
INDEX_OPTS="--kmer 27 --edge_max 5"
#OPTS="--dir_tag --only_summary"
OPTS="--dir_tag"

# output for clusteGraphs.py that we rename
CP_FILES=( "heatmap_kmer.pdf" "heatmap_log_kmer.pdf" "heatmap_vm1_kmer.pdf" "heatmap_log_vm1_kmer.pdf" "tree_kmer.dot" "tree_kmer.newick"  "tree_kmer.png" "heatmap_corg.pdf" "heatmap_log_corg.pdf" "heatmap_vm1_corg.pdf" "heatmap_log_vm1_corg.pdf" "tree_corg.dot" "tree_corg.newick"  "tree_corg.png" )


mkdir $OUT_DIR

#for i in brca1 brca2 sma lrc_kr mhc cenx
for i in brca1 brca2 sma lrc_kir mhc
#for i in brca1
do
	 # heatmap of everything
	 rm -rf ${TOIL_DIR} ; scripts/clusterGraphs.py ./${TOIL_DIR} ${GRAPHS}/*${i}*.vg ${VARIANTS}/${i}/*/*.vg ${OUT_DIR}/${i} ${TOIL_OPTS} ${INDEX_OPTS} ${OPTS} --avg_sample
	 for j in "${CP_FILES[@]}"
	 do
		  cp ${OUT_DIR}/${i}/${j} ${OUT_DIR}/${i}/${i}_avg_${j}
	 done

	 rm -rf ${TOIL_DIR} ; scripts/clusterGraphs.py ./${TOIL_DIR} ${GRAPHS}/*${i}*.vg ${VARIANTS}/${i}/*/*.vg ${OUT_DIR}/${i} ${TOIL_OPTS} ${INDEX_OPTS} ${OPTS}

	 for j in "${CP_FILES[@]}"
	 do
		  cp ${OUT_DIR}/${i}/${j} ${OUT_DIR}/${i}/${i}_${j}
	 done

	 #heatmap of original
	 rm -rf ${TOIL_DIR} ; scripts/clusterGraphs.py ./${TOIL_DIR} ${GRAPHS}/*${i}*.vg ${OUT_DIR}/${i} ${TOIL_OPTS} ${INDEX_OPTS} ${OPTS}

	 for j in "${CP_FILES[@]}"
	 do
		  cp ${OUT_DIR}/${i}/${j} ${OUT_DIR}/${i}/original_${i}_${j}
	 done

	 # heatmap of categories
	 for k in sample linear augmented
	 do
		  rm -rf ${TOIL_DIR} ; scripts/clusterGraphs.py ./${TOIL_DIR} ${VARIANTS}/${i}/*/*${k}*.vg ${OUT_DIR}/${i} ${TOIL_OPTS} ${INDEX_OPTS} ${OPTS} --avg_sample

		  for j in "${CP_FILES[@]}"
		  do
				cp ${OUT_DIR}/${i}/${j} ${OUT_DIR}/${i}/${k}_${i}_avg_${j}
		  done

		  rm -rf ${TOIL_DIR} ; scripts/clusterGraphs.py ./${TOIL_DIR} ${GRAPHS}/*${i}*.vg ${VARIANTS}/${i}/*/*${k}*.vg ${OUT_DIR}/${i} ${TOIL_OPTS} ${INDEX_OPTS} ${OPTS} --avg_sample

		  for j in "${CP_FILES[@]}"
		  do
				cp ${OUT_DIR}/${i}/${j} ${OUT_DIR}/${i}/${k}_and_original_${i}_avg_${j}
		  done
	 done

	 # linear sample heatmap
	 rm -rf ${TOIL_DIR} ; scripts/clusterGraphs.py ./${TOIL_DIR} ${VARIANTS}/${i}/*/*sample*.vg ${VARIANTS}/${i}/*/*linear*.vg ${OUT_DIR}/${i} ${TOIL_OPTS} ${INDEX_OPTS} ${OPTS} --avg_sample

	 for j in "${CP_FILES[@]}"
	 do
		  cp ${OUT_DIR}/${i}/${j} ${OUT_DIR}/${i}/sample_and_linear_${i}_avg_${j}
	 done

	 # tables
<<<<<<< HEAD
	 scripts/callStats.py ${ALIGNMENTS}/${i}/*/*.gam --out_dir ${VARIANTS} --avg_sample --out_sub $i ${OPTS}
=======
	 ./callStats.py ${ALIGNMENTS}/${i}/*/*.gam --out_dir ${VARIANTS} --avg_sample --out_sub $i ${OPTS}
>>>>>>> small cleanup of variants scripts. add corg to comparison logic

	 cp ${VARIANTS}/compare/${i}/*.tsv ${OUT_DIR}/${i}
	 
done
