#!/usr/bin/env bash
# Makes plots comparing indel lengths across graphs

set -ex

# What plot filetype should we produce?
PLOT_FILETYPE="png"

# Grab the input directory to look in
INPUT_DIR=${1}
# And the output directory to plot in
OUTPUT_DIR=${2}

if [[ ! -d "${INPUT_DIR}" ]]
then
    echo "Specify input directory!"
    exit 1
fi

if [[ ! -d "${OUTPUT_DIR}" ]]
then
    mkdir -p "${OUTPUT_DIR}"
fi

# Set up the plot parameters
# Include both versions of the 1kg SNPs graph name
PLOT_PARAMS=(
    --categories
    snp1kg
    snp1000g
    sbg
    cactus
    camel
    curoverse
    debruijn-k31
    debruijn-k63
    level1
    level2
    level3
    prg
    refonly
    simons
    trivial
    vglr
    haplo1kg30
    haplo1kg50
    shifted1kg
    --category_labels 
    1KG
    1KG
    7BG
    Cactus
    Camel
    Curoverse
    "De Bruijn 31"
    "De Bruijn 63"
    Level1
    Level2
    Level3
    PRG
    Primary
    SGDP
    Unmerged
    VGLR
    "1KG Haplo 30"
    "1KG Haplo 50"
    Scrambled
    --colors
    "#fb9a99"
    "#fb9a99"
    "#b15928"
    "#1f78b4"
    "#33a02c"
    "#a6cee3"
    "#e31a1c"
    "#ff7f00"
    "#FF0000"
    "#00FF00"
    "#0000FF"
    "#6a3d9a"
    "#000000"
    "#b2df8a"
    "#b1b300"
    "#cab2d6"
    "#00FF00"
    "#0000FF"
    "#FF0000"
    --font_size 20 --dpi 90
)


# `ls ${INPUT_DIR}/mapping.*.tsv | xargs -n 1 basename | sed 's/mapping.\(.*\).tsv/\1/'`



for REGION in `ls ${INPUT_DIR} | xargs -n 1 basename`
do
    # For every region we ran
    
    # We need to make a TSV to hold the plotting data
    INDEL_LENGTH_TSV="${OUTPUT_DIR}/indel_lengths_${REGION}.tsv"
    
    # And we also want mean indel length, because we really do care about the
    # long ones.
    INDEL_MEAN_TSV="${OUTPUT_DIR}/indel_means_${REGION}.tsv"
    
    # And the sum
    INDEL_SUM_TSV="${OUTPUT_DIR}/indel_sums_${REGION}.tsv"
    
    # Where do we put the plots?
    LENGTH_PLOT_FILE="${OUTPUT_DIR}/plot_indel_lengths_${REGION}.${PLOT_FILETYPE}"
    MEAN_PLOT_FILE="${OUTPUT_DIR}/plot_indel_means_${REGION}.${PLOT_FILETYPE}"
    SUM_PLOT_FILE="${OUTPUT_DIR}/plot_indel_sums_${REGION}.${PLOT_FILETYPE}"
    
    # Empty the plotting data files
    rm -f "${INDEL_LENGTH_TSV}"
    touch "${INDEL_LENGTH_TSV}"
    rm -f "${INDEL_MEAN_TSV}"
    touch "${INDEL_MEAN_TSV}"
    rm -f "${INDEL_SUM_TSV}"
    touch "${INDEL_SUM_TSV}"
    
    for GRAPH in `ls ${INPUT_DIR}/${REGION} | xargs -n 1 basename`
    do
        # For every graph we ran on this region
    
        for SAMPLE_FILE in `ls ${INPUT_DIR}/${REGION}/${GRAPH}/*_sample.vcf`
        do
    
            # For all the un-normalized sample call VCFs for the graph
            # Don't use the normalized ones because they are duplicating and weirdly left-shifting deletions.
    
            # Get all the indel lengths
            scripts/indelLengths.py --in_file ${SAMPLE_FILE} --indels_only > "${OUTPUT_DIR}/temp.tsv"
    
            # Get all the indel lengths, tack on the graph name, and save them to the file
            cat "${OUTPUT_DIR}/temp.tsv" | awk "{print \"${GRAPH}\t\", \$1}" >> "${INDEL_LENGTH_TSV}"
            
            # Also do the mean length
            cat "${OUTPUT_DIR}/temp.tsv" | scripts/mean.sh | awk "{print \"${GRAPH}\t\", \$1}" >> "${INDEL_MEAN_TSV}"
            
            # And the total length
            cat "${OUTPUT_DIR}/temp.tsv" | scripts/sum.sh | awk "{print \"${GRAPH}\t\", \$1}" >> "${INDEL_SUM_TSV}"
            
            rm "${OUTPUT_DIR}/temp.tsv"
            
        done
            
    done
    
    echo "Plotting ${REGION^^}..."
        
    # Remove underscores from region names to make them human readable
    HR_REGION=`echo ${REGION^^} | sed 's/_/ /g'`
    
    # TODO: you need to run collateStatistics.py to build the per-region-and-
    # graph stats files. We expect them to exist and only concatenate the final
    # overall files and make the plots.
    
    ./scripts/boxplot.py "${INDEL_LENGTH_TSV}" \
        --title "$(printf "Indel lengths in ${HR_REGION}")" \
        --x_label "Graph" --y_label "Indel length (bp)" --save "${LENGTH_PLOT_FILE}" \
        --x_sideways --hline_median refonly \
        --range --log_y --min 0.1 --no_n \
        "${PLOT_PARAMS[@]}"
        
    ./scripts/barchart.py "${INDEL_MEAN_TSV}" \
        --title "$(printf "Indel length means in ${HR_REGION}")" \
        --x_label "Graph" --y_label "Mean indel length (bp)" --save "${MEAN_PLOT_FILE}" \
        --x_sideways \
        --sparse_ticks --log_y \
        "${PLOT_PARAMS[@]}"
        
    ./scripts/barchart.py "${INDEL_SUM_TSV}" \
        --title "$(printf "Indel length sums in ${HR_REGION}")" \
        --x_label "Graph" --y_label "Total indel length (bp)" --save "${SUM_PLOT_FILE}" \
        --x_sideways \
        --sparse_ticks --log_y \
        "${PLOT_PARAMS[@]}"
        
    
done

