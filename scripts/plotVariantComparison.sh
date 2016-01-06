#!/usr/bin/env bash
# Run after evaluateVariantCalls.py.
# Makes plots comparing the variant call concordances across graphs.

set -e

# Grab the input directory to look in
INPUT_DIR=${1}

if [[ ! -d "${INPUT_DIR}" ]]
then
    echo "Specify input directory!"
    exit 1
fi

# Set up the plot parameters
# Include both versions of the 1kg SNPs graph name
PLOT_PARAMS=(
    --categories
    snp1kg
    snp1000g
    haplo1kg
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
    "1KG Haplo"
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
    Control
    --colors
    "#fb9a99"
    "#fb9a99"
    "#fdbf6f"
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

# Where are the stats files
STATS_DIR="${INPUT_DIR}/stats"

# Where do we put the plots?
PLOTS_ROOT_DIR="${INPUT_DIR}/plots"
mkdir -p "${PLOTS_ROOT_DIR}"

# Where are the GATK output files by region?
GATK_DIR="${STATS_DIR}/gatk"

for REGION in `ls ${GATK_DIR}`
do
    # For every region we ran
    
    # Work out its directory
    REGION_DIR="${GATK_DIR}/${REGION}"
    
    # Plot the bases dropped
    DROPPED_TSV="${PLOTS_ROOT_DIR}/dropped-${REGION}.tsv"
    DROPPED_PNG="${PLOTS_ROOT_DIR}/plot-dropped-${REGION}.png"
    cat "${STATS_DIR}/bases_dropped.tsv" | grep "${REGION}" | cut -f2,4 > "${DROPPED_TSV}"
    scripts/barchart.py "${DROPPED_TSV}" --save "${DROPPED_PNG}" --title "Bases Dropped in ${REGION^^}" --x_label "Graph" --y_label "Bases" --x_sideways --log_y "${PLOT_PARAMS[@]}"
    
    # Plot the concordances somehow
    
    # Where will the concordance data for comparable sites go?
    CONCORDANCE_TSV="${PLOTS_ROOT_DIR}/concordance-${REGION}.tsv"
    CONCORDANCE_PNG="${PLOTS_ROOT_DIR}/plot-concordance-${REGION}.png"
    
    # And the fraction of true variants that were able to be used for comparison
    COMPARABLE_TSV="${PLOTS_ROOT_DIR}/comparable-${REGION}.tsv"
    COMPARABLE_PNG="${PLOTS_ROOT_DIR}/plot-comparable-${REGION}.png"
    
    # And precision and recall for there being a variant?
    PRECISION_RECALL_TSV="${PLOTS_ROOT_DIR}/precisionrecall-${REGION}.tsv"
    PRECISION_RECALL_PNG="${PLOTS_ROOT_DIR}/plot-precisionrecall-${REGION}.png"
    
    # Clear up any old values
    truncate -s 0 "${CONCORDANCE_TSV}"
    truncate -s 0 "${COMPARABLE_TSV}"
    truncate -s 0 "${PRECISION_RECALL_TSV}"
    
    for GRAPH in `ls ${REGION_DIR}`
    do
        # For every graph
        GRAPH_DIR="${REGION_DIR}/${GRAPH}"
    
    
        for SAMPLE in `ls ${GRAPH_DIR}`
        do
            # For each grp statistics file for that graph type
            
            # Grab the sample name
            SAMPLE_NAME="${SAMPLE%.grp}"
            
            # Grab the file name
            SAMPLE_FILE="${GRAPH_DIR}/${SAMPLE}"
            
            echo "Processing ${SAMPLE_FILE} concordance..."
            
            # Add graph and concordance value for comparable sites to the concordance TSV
            printf "${GRAPH}\t" >> "${CONCORDANCE_TSV}"
            cat ${SAMPLE_FILE} | tail -n 7 | head -n 1 | awk '{print $4}' >> "${CONCORDANCE_TSV}"
            
            echo "Processing ${SAMPLE_FILE} comparable sites..."
            
            # Add graph and portion of sites comparable value to the comparable TSV
            # Calculate sites with matching alleles / sites with any variant in the truth set
            printf "${GRAPH}\t" >> "${COMPARABLE_TSV}"
            cat ${SAMPLE_FILE} | tail -n 2 | head -n 1 | awk '{print ($1 / ($1 + $2 + $3 + $4 + $6))}' >> "${COMPARABLE_TSV}"
            
            echo "Processing ${SAMPLE_FILE} precision and recall..."
            
            # Calculate precision and recall (to handle zero denominators)
            PR_LINE=`cat ${SAMPLE_FILE} | tail -n 2 | head -n 1`
            
            PR_NUM=`echo "${PR_LINE}" | awk '{print ($1 + $2 + $3 + $4)}'`
            PRECISION_DENOM=`echo "${PR_LINE}" | awk '{print ($1 + $2 + $3 + $4 + $5)}'`
            RECALL_DENOM=`echo "${PR_LINE}" | awk '{print ($1 + $2 + $3 + $4 + $6)}'`
            
            if [ "${PRECISION_DENOM}" == "0" ]
            then
                # Limit of precision is 1
                PRECISION=1
            else
                PRECISION=`echo "${PR_NUM} / ${PRECISION_DENOM}" | bc -l`
            fi
            
            if [ "${RECALL_DENOM}" == "0" ]
            then
                # Limit of recall is 1
                RECALL=1
            else
                RECALL=`echo "${PR_NUM} / ${RECALL_DENOM}" | bc -l`
            fi
            
            
            # Do precision (portion of sites with any variant in the truth when there's a variant in the sample)
            printf "${GRAPH}\t${RECALL}\t${PRECISION}\n" >> "${PRECISION_RECALL_TSV}"
            
            
        done
        
    done
    
    # Now plot the plots
    scripts/boxplot.py "${CONCORDANCE_TSV}" --save "${CONCORDANCE_PNG}" --title "Genotype Concordance in ${REGION^^}" --x_label "Graph" --y_label "$(printf "Portion of identical variants\nwith concordant genotypes")" --x_sideways "${PLOT_PARAMS[@]}"
    
    scripts/boxplot.py "${COMPARABLE_TSV}" --save "${COMPARABLE_PNG}" --title "Identical Variants in ${REGION^^}" --x_label "Graph" --y_label "$(printf "Portion of true variants\nwith correct alleles identified")" --x_sideways "${PLOT_PARAMS[@]}"
    
    scripts/scatter.py "${PRECISION_RECALL_TSV}" --save "${PRECISION_RECALL_PNG}" --title "Variant Existence Precision vs. Recall in ${REGION^^}" --x_label "Recall" --y_label "Precision" "${PLOT_PARAMS[@]}"
    
done











