#!/usr/bin/env bash
# Makes plots comparing indel lengths across graphs

set -ex

# What plot filetype should we produce?
PLOT_FILETYPE="svg"

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
#    sbg
    cactus
#    camel
#    curoverse
#    debruijn-k31
#    debruijn-k63
#    level1
#    level2
#    level3
#    prg
    refonly
#    simons
#    trivial
#    vglr
#    haplo1kg30
#    haplo1kg50
#    shifted1kg
    --category_labels 
    1KG
    1KG
#    7BG
    Cactus
#    Camel
#    Curoverse
#    "De Bruijn 31"
#    "De Bruijn 63"
#    Level1
#    Level2
#    Level3
#    PRG
    Primary
#    SGDP
#    Unmerged
#    VGLR
#    "1KG Haplo 30"
#    "1KG Haplo 50"
#    Scrambled
    --colors
    "#fb9a99"
    "#fb9a99"
#    "#b15928"
    "#1f78b4"
#    "#33a02c"
#    "#a6cee3"
#    "#e31a1c"
#    "#ff7f00"
#    "#FF0000"
#    "#00FF00"
#    "#0000FF"
#    "#6a3d9a"
    "#000000"
#    "#b2df8a"
#    "#b1b300"
#    "#cab2d6"
#    "#00FF00"
#    "#0000FF"
#    "#FF0000"
    --dpi 90
)


# `ls ${INPUT_DIR}/mapping.*.tsv | xargs -n 1 basename | sed 's/mapping.\(.*\).tsv/\1/'`



for REGION in `ls ${INPUT_DIR} | xargs -n 1 basename`
do
    # For every region we ran
    
    # We need to make a TSV to hold the plotting data
    INDEL_LENGTH_TSV="${OUTPUT_DIR}/indel_lengths_${REGION}.tsv"
    REF_INDEL_LENGTH_TSV="${OUTPUT_DIR}/indel_lengths_ref_${REGION}.tsv"
    NONREF_INDEL_LENGTH_TSV="${OUTPUT_DIR}/indel_lengths_nonref_${REGION}.tsv"
    
    # Where do we put the plots?
    LENGTH_PLOT_FILE="${OUTPUT_DIR}/plot_indel_lengths_${REGION}.${PLOT_FILETYPE}"
    HISTOGRAM_PLOT_FILE="${OUTPUT_DIR}/plot_indel_histogram_${REGION}.${PLOT_FILETYPE}"
    LONG_HISTOGRAM_PLOT_FILE="${OUTPUT_DIR}/plot_indel_histogram_long_${REGION}.${PLOT_FILETYPE}"
    REF_HISTOGRAM_PLOT_FILE="${OUTPUT_DIR}/plot_indel_histogram_ref_${REGION}.${PLOT_FILETYPE}"
    NONREF_HISTOGRAM_PLOT_FILE="${OUTPUT_DIR}/plot_indel_histogram_nonref_${REGION}.${PLOT_FILETYPE}"
    
    # Empty the plotting data files
    rm -f "${INDEL_LENGTH_TSV}"
    touch "${INDEL_LENGTH_TSV}"
    rm -f "${REF_INDEL_LENGTH_TSV}"
    touch "${REF_INDEL_LENGTH_TSV}"
    rm -f "${NONREF_INDEL_LENGTH_TSV}"
    touch "${NONREF_INDEL_LENGTH_TSV}"

    for GRAPH in `ls ${INPUT_DIR}/${REGION} | xargs -n 1 basename`
    do
        # For every graph we ran on this region
    
        for SAMPLE_FILE in `ls ${INPUT_DIR}/${REGION}/${GRAPH}/*_sample_preprocessed.vcf`
        do
    
            # For all the un-normalized sample call VCFs for the graph
            # Don't use the normalized ones because they are duplicating and weirdly left-shifting deletions.
    
            # Get all the indel lengths
            scripts/indelLengths.py --in_file ${SAMPLE_FILE} --indels_only --distinguish > "${OUTPUT_DIR}/temp.tsv"
    
            # Get all the indel lengths, tack on the graph name, and save them to the file
            cat "${OUTPUT_DIR}/temp.tsv" | awk "{print \"${GRAPH}\t\", \$1}" >> "${INDEL_LENGTH_TSV}"
            
            # Do it for ref only
            cat ${SAMPLE_FILE} | grep "XREF" | scripts/indelLengths.py --indels_only --distinguish > "${OUTPUT_DIR}/temp.tsv"
            cat "${OUTPUT_DIR}/temp.tsv" | awk "{print \"${GRAPH}\t\", \$1}" >> "${REF_INDEL_LENGTH_TSV}"
            
            # And for nonref only
            cat ${SAMPLE_FILE} | grep -v "XREF" | scripts/indelLengths.py --indels_only --distinguish > "${OUTPUT_DIR}/temp.tsv"
            cat "${OUTPUT_DIR}/temp.tsv" | awk "{print \"${GRAPH}\t\", \$1}" >> "${NONREF_INDEL_LENGTH_TSV}"
            
            
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
        --range --log_y --min 0.1 --no_n --font_size 20 \
        "${PLOT_PARAMS[@]}"
        
    # Do a length histogram
    ./scripts/histogram.py "${INDEL_LENGTH_TSV}" \
        --title "$(printf "All indel lengths in ${HR_REGION}")" \
        --x_label "Length (bp)" --y_label "Indel count" --save "${HISTOGRAM_PLOT_FILE}" \
        --line --no_n \
        --bins 30 --no_zero_ends --fake_zero --split_at_zero --log_counts --x_min -150 --x_max 150 --y_min 0.9 \
        --legend_overlay "upper right" \
        --style "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" \
        --width 6 --height 2 \
        "${PLOT_PARAMS[@]}"
    
    # And one that's over a wider range
    ./scripts/histogram.py "${INDEL_LENGTH_TSV}" \
        --title "$(printf "All indel lengths in ${HR_REGION}")" \
        --x_label "Length (bp)" --y_label "Indel count" --save "${LONG_HISTOGRAM_PLOT_FILE}" \
        --line --no_n \
        --bins 50 --no_zero_ends --fake_zero --split_at_zero --log_counts --x_min -500 --x_max 500 --y_min 0.9 \
        --legend_overlay "upper right" \
        --ticks -500 -150 0 150 500 \
        --style "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" \
        --width 6 --height 2 \
        "${PLOT_PARAMS[@]}"
        
    # And histograms for ref only and nonref only
    ./scripts/histogram.py "${REF_INDEL_LENGTH_TSV}" \
        --title "$(printf "Reference indel lengths in ${HR_REGION}")" \
        --x_label "Length (bp)" --y_label "Indel count" --save "${REF_HISTOGRAM_PLOT_FILE}" \
        --line --no_n \
        --bins 30 --no_zero_ends --fake_zero --split_at_zero --log_counts --x_min -150 --x_max 150 --y_min 0.9 \
        --legend_overlay "upper right" \
        --style "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" \
        --width 6 --height 2 \
        "${PLOT_PARAMS[@]}"
        
    ./scripts/histogram.py "${NONREF_INDEL_LENGTH_TSV}" \
        --title "$(printf "Nonreference indel lengths in ${HR_REGION}")" \
        --x_label "Length (bp)" --y_label "Indel count" --save "${NONREF_HISTOGRAM_PLOT_FILE}" \
        --line --no_n \
        --bins 30 --no_zero_ends --fake_zero --split_at_zero --log_counts --x_min -150 --x_max 150 --y_min 0.9 \
        --legend_overlay "upper right" \
        --style "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" \
        --width 6 --height 2 \
        "${PLOT_PARAMS[@]}"
        
done

