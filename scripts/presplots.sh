#!/usr/bin/env bash
# Make one-off plots for presentation

set -e

# Grab the input directory to look in
INPUT_DIR=${1}

if [[ ! -d "${INPUT_DIR}" ]]
then
    echo "Specify input directory!"
    exit 1
fi

# Set up the plot parameters
# USe a ColorBrewer scheme, except drop the super light color
PLOT_PARAMS=(
    --categories "cactus" "camel" "vg" "curoverse" "simons" "snp1000g" "prg" "debruijn-k31" "debruijn-k63"
    "refonly" "trivial" "level1" "level2" "level3"
    --category_labels "Cactus" "Camel" "VG"  "Curoverse" "Simons SNPs" "1000 GSNPs" "PRG" "k=31" "k=63"
    "RefOnly" "Trivial" "Level1" "Level2" "Level3"
    --colors "#7fc97f" "#fdc086" "k" "k" "#386cb0" "#beaed4" "k" "#8000FF" "#f0027f" "#bf5b17" "#666666" "k" "k" "k"
)

# Where are the stats files
STATS_DIR="${INPUT_DIR}/stats"

# Where do we put the plots?
PLOTS_DIR="${INPUT_DIR}/presplots"
mkdir -p "${PLOTS_DIR}"

# We need overall files for mapped and multimapped
OVERALL_MAPPING_FILE="${PLOTS_DIR}/mapping.tsv"
OVERALL_MAPPING_PLOT="${PLOTS_DIR}/mapping.png"
OVERALL_MULTIMAPPING_FILE="${PLOTS_DIR}/multimapping.tsv"
OVERALL_MULTIMAPPING_PLOT="${PLOTS_DIR}/multimapping.png"

for REGION in `ls ${STATS_DIR}`
do
    # For every region we ran
    
    # We have intermediate data files for plotting from
    MAPPING_FILE="${PLOTS_DIR}/mapping.${REGION}.tsv"
    MAPPING_PLOT="${PLOTS_DIR}/mapping.${REGION}.png"
    MULTIMAPPING_FILE="${PLOTS_DIR}/multimapping.${REGION}.tsv"
    MULTIMAPPING_PLOT="${PLOTS_DIR}/multimapping.${REGION}.png"
    RUNTIME_FILE="${PLOTS_DIR}/runtime.${REGION}.tsv"
    RUNTIME_PLOT="${PLOTS_DIR}/runtime.${REGION}.png"
    
    # Make them empty
    :>"${MAPPING_FILE}"
    :>"${MULTIMAPPING_FILE}"
    :>"${RUNTIME_FILE}"
    
    
    echo "Plotting ${REGION^^}..."
    
    for GRAPH_NAME in `ls ${STATS_DIR}/${REGION}`
    do
        # For every graph we ran for it
        
        for STATS_FILE in `ls ${STATS_DIR}/${REGION}/${GRAPH_NAME}`
        do
            # For each sample run, parse its JSON and add a point to the region
            # TSV for the appropriate graph.
            
            # First build the path to the JSON file to look at
            JSON_FILE="${STATS_DIR}/${REGION}/${GRAPH_NAME}/${STATS_FILE}"
            
            # We need to account for the well-mapped/well-multimapped identity thresholds
            
            # First: portion mapped with <=2 mismatches out of 100 expected length
            printf "${GRAPH_NAME}\t" >> "${MAPPING_FILE}"
            # We need the 0 + in case there are no sufficiently good mappings
            cat "${JSON_FILE}" | jq -r '(0 + .primary_mismatches."0" + .primary_mismatches."1" + .primary_mismatches."2") / .total_reads' >> "${MAPPING_FILE}"
            
            # Next: portion NOT multimapped with <=2 mismatches out of 100 expected length
            printf "${GRAPH_NAME}\t" >> "${MULTIMAPPING_FILE}"
            # We need the 0 + in case there are no sufficiently good mappings
            cat "${JSON_FILE}" | jq -r '1 - ((0 + .secondary_mismatches."0" + .secondary_mismatches."1" + .secondary_mismatches."2") / .total_reads)' >> "${MULTIMAPPING_FILE}"
            
            # Next: runtime in seconds
            printf "${GRAPH_NAME}\t" >> "${RUNTIME_FILE}"
            cat "${JSON_FILE}" | jq -r '.run_time' >> "${RUNTIME_FILE}"
        done
    done
    
    ./boxplot.py "${MAPPING_FILE}" \
        --title "" \
        --x_label "Graph" --y_label "Portion well-mapped" --save "${MAPPING_PLOT}" \
        --x_sideways \
        "${PLOT_PARAMS[@]}" \
        --font_size 20 --dpi 300 --line_width 3
        
    ./boxplot.py "${MULTIMAPPING_FILE}" \
        --title "" \
        --x_label "Graph" --y_label "Portion not-well-multimapped" --save "${MULTIMAPPING_PLOT}" \
        --x_sideways \
        "${PLOT_PARAMS[@]}" \
        --font_size 20 --dpi 300 --line_width 3
        
    ./boxplot.py "${RUNTIME_FILE}" \
        --title "" \
        --x_label "Graph" --y_label "Runtime per sample (seconds)" --save "${RUNTIME_PLOT}" \
        --x_sideways \
        "${PLOT_PARAMS[@]}" \
        --font_size 20 --dpi 300 --line_width 3
    
done

# Aggregate the overall files
cat "${PLOTS_DIR}"/mapping.*.tsv > "${OVERALL_MAPPING_FILE}"
cat "${PLOTS_DIR}"/multimapping.*.tsv > "${OVERALL_MULTIMAPPING_FILE}"

# Make the overall plots
./boxplot.py "${OVERALL_MAPPING_FILE}" \
    --title "" \
    --x_label "Graph" --y_label "Portion well-mapped" --save "${OVERALL_MAPPING_PLOT}" \
    --x_sideways \
    "${PLOT_PARAMS[@]}" \
    --font_size 20 --dpi 300 --line_width 3
        
./boxplot.py "${OVERALL_MULTIMAPPING_FILE}" \
    --title "" \
    --x_label "Graph" --y_label "Portion not-well-multimapped" --save "${OVERALL_MULTIMAPPING_PLOT}" \
    --x_sideways \
    "${PLOT_PARAMS[@]}" \
    --font_size 20 --dpi 300 --line_width 3

