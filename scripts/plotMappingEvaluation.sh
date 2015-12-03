#!/usr/bin/env bash
# Run after mapping_evaluations.sh or parallelMappingEvaluation.py.
# Makes plots comparing the graphs in each.

set -e

# Grab the input directory to look in
INPUT_DIR=${1}

if [[ ! -d "${INPUT_DIR}" ]]
then
    echo "Specify input directory!"
    exit 1
fi

# Set up the plot parameters
PLOT_PARAMS=(
    --categories
    snp1kg
    haplo1kg
    sbg
    cactus
    camel
    curoverse
    Debruijn-k31
    Debruijn-k63
    level1
    level2
    level3
    prg
    refonly
    simons
    trivial
    vglr
    --category_labels 
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
    Reference
    SGDP
    Trivial
    VGLR
    --colors
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
    --font_size 20 --dpi 90 --no_n
)

# Where are the stats files
STATS_DIR="${INPUT_DIR}/stats"

# Where do we put the plots?
PLOTS_DIR="${INPUT_DIR}/plots"
mkdir -p "${PLOTS_DIR}"

# We need overall files for mapped and multimapped
OVERALL_MAPPING_FILE="${PLOTS_DIR}/mapping.tsv"
OVERALL_MAPPING_PLOT="${PLOTS_DIR}/mapping.ALL.png"
OVERALL_PERFECT_FILE="${PLOTS_DIR}/perfect.tsv"
OVERALL_PERFECT_PLOT="${PLOTS_DIR}/perfect.ALL.png"
OVERALL_ONE_ERROR_FILE="${PLOTS_DIR}/oneerror.tsv"
OVERALL_ONE_ERROR_PLOT="${PLOTS_DIR}/oneerror.ALL.png"
OVERALL_SINGLE_MAPPING_FILE="${PLOTS_DIR}/singlemapping.tsv"
OVERALL_SINGLE_MAPPING_PLOT="${PLOTS_DIR}/singlemapping.ALL.png"

for REGION in `ls ${PLOTS_DIR}/mapping.*.tsv | xargs -n 1 basename | sed 's/mapping.\(.*\).tsv/\1/'`
do
    # For every region we ran
    
    # We have intermediate data files for plotting from
    MAPPING_FILE="${PLOTS_DIR}/mapping.${REGION}.tsv"
    MAPPING_PLOT="${PLOTS_DIR}/mapping.${REGION}.png"
    PERFECT_FILE="${PLOTS_DIR}/perfect.${REGION}.tsv"
    PERFECT_PLOT="${PLOTS_DIR}/perfect.${REGION}.png"
    ONE_ERROR_FILE="${PLOTS_DIR}/oneerror.${REGION}.tsv"
    ONE_ERROR_PLOT="${PLOTS_DIR}/oneerror.${REGION}.png"
    SINGLE_MAPPING_FILE="${PLOTS_DIR}/singlemapping.${REGION}.tsv"
    SINGLE_MAPPING_PLOT="${PLOTS_DIR}/singlemapping.${REGION}.png"
    RUNTIME_FILE="${PLOTS_DIR}/runtime.${REGION}.tsv"
    RUNTIME_PLOT="${PLOTS_DIR}/runtime.${REGION}.png"
    
    echo "Plotting ${REGION^^}..."
    
    # Remove underscores from region names to make them human readable
    HR_REGION=`echo ${REGION^^} | sed 's/_/ /g'`
    
    # TODO: you need to run collateStatistics.py to build the per-region-and-
    # graph stats files. We expect them to exist and only concatenate the final
    # overall files and make the plots.
    
    ./scripts/boxplot.py "${MAPPING_FILE}" \
        --title "$(printf "Mapped (<=2 mismatches)\nreads in ${HR_REGION}")" \
        --x_label "Graph" --y_label "Portion mapped" --save "${MAPPING_PLOT}" \
        --x_sideways --hline_median refonly \
        "${PLOT_PARAMS[@]}"
        
    ./scripts/boxplot.py "${PERFECT_FILE}" \
        --title "$(printf "Perfectly mapped\nreads in ${HR_REGION}")" \
        --x_label "Graph" --y_label "Portion perfectly mapped" --save "${PERFECT_PLOT}" \
        --x_sideways --hline_median refonly \
        "${PLOT_PARAMS[@]}"
        
    ./scripts/boxplot.py "${ONE_ERROR_FILE}" \
        --title "$(printf "One-error (<=1 mismatch)\nreads in ${HR_REGION}")" \
        --x_label "Graph" --y_label "Portion" --save "${ONE_ERROR_PLOT}" \
        --x_sideways --hline_median refonly \
        "${PLOT_PARAMS[@]}"
        
    ./scripts/boxplot.py "${SINGLE_MAPPING_FILE}" \
        --title "$(printf "Uniquely mapped (<=2 mismatches)\nreads in ${HR_REGION}")" \
        --x_label "Graph" --y_label "Portion uniquely mapped" --save "${SINGLE_MAPPING_PLOT}" \
        --x_sideways --hline_median refonly \
        "${PLOT_PARAMS[@]}"
        
    ./scripts/boxplot.py "${RUNTIME_FILE}" \
        --title "$(printf "Per-read runtime\n in ${HR_REGION}")" \
        --x_label "Graph" --y_label "Runtime per read (seconds)" --save "${RUNTIME_PLOT}" \
        --x_sideways --max_max 0.006 \
        "${PLOT_PARAMS[@]}"
    
done

# Aggregate the overall files
cat "${PLOTS_DIR}"/mapping.*.tsv > "${OVERALL_MAPPING_FILE}"
cat "${PLOTS_DIR}"/perfect.*.tsv > "${OVERALL_PERFECT_FILE}"
cat "${PLOTS_DIR}"/oneerror.*.tsv > "${OVERALL_ONE_ERROR_FILE}"
cat "${PLOTS_DIR}"/singlemapping.*.tsv > "${OVERALL_SINGLE_MAPPING_FILE}"

# Make the overall plots
./scripts/boxplot.py "${OVERALL_MAPPING_FILE}" \
    --title "$(printf "Mapped (<=2 mismatches)\nreads")" \
    --x_label "Graph" --y_label "Portion mapped" --save "${OVERALL_MAPPING_PLOT}" \
    --x_sideways  --hline_median trivial \
    "${PLOT_PARAMS[@]}"
    
./scripts/boxplot.py "${OVERALL_PERFECT_FILE}" \
    --title "$(printf "Perfectly mapped\nreads")" \
    --x_label "Graph" --y_label "Portion perfectly mapped" --save "${OVERALL_PERFECT_PLOT}" \
    --x_sideways --hline_median trivial \
    "${PLOT_PARAMS[@]}"
    
./scripts/boxplot.py "${OVERALL_ONE_ERROR_FILE}" \
    --title "$(printf "One-error (<=1 mismatch)\nreads")" \
    --x_label "Graph" --y_label "Portion" --save "${OVERALL_ONE_ERROR_PLOT}" \
    --x_sideways --hline_median trivial \
    "${PLOT_PARAMS[@]}"

./scripts/boxplot.py "${OVERALL_SINGLE_MAPPING_FILE}" \
    --title "$(printf "Uniquely mapped (<=2 mismatches)\nreads")" \
    --x_label "Graph" --y_label "Portion uniquely mapped" --save "${OVERALL_SINGLE_MAPPING_PLOT}" \
    --x_sideways --hline_median refonly \
    "${PLOT_PARAMS[@]}"

