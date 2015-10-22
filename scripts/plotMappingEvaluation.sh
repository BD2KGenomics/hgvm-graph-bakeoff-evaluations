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
    --categories "cactus" "camel" "vg" "curoverse" "simons" "snp1000g" "prg" "debruijn-k31" "debruijn-k63"
    "sbg" "refonly" "trivial" "level1" "level2" "level3"
    --category_labels "Cactus" "Camel" "VG"  "Curoverse" "Simons" "1000 Genomes" "PRG" "k=31" "k=63"
    "7 Bridges" "RefOnly" "Trivial" "Level1" "Level2" "Level3"
    --colors "#5C755E" "#C19A6B" "#000099" "#31184A" "#384DA0" "k" "#353C47" "r" "m"
    "#93AC2B" "c" "b" "c" "m" "y"
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
    
    # TODO: you need to run collateStatistics.py to build the per-region-and-
    # graph stats files. We expect them to exist and only concatenate the final
    # overall files and make the plots.
    
    ./boxplot.py "${MAPPING_FILE}" \
        --title "$(printf "Mapped (<=2 mismatches)\nreads in ${REGION^^}")" \
        --x_label "Graph" --y_label "Portion mapped" --save "${MAPPING_PLOT}" \
        --x_sideways --hline_median trivial \
        "${PLOT_PARAMS[@]}"
        
    ./boxplot.py "${PERFECT_FILE}" \
        --title "$(printf "Perfectly mapped\nreads in ${REGION^^}")" \
        --x_label "Graph" --y_label "Portion perfectly mapped" --save "${PERFECT_PLOT}" \
        --x_sideways --hline_median trivial \
        "${PLOT_PARAMS[@]}"
        
    ./boxplot.py "${ONE_ERROR_FILE}" \
        --title "$(printf "One-error (<=1 mismatch)\nreads in ${REGION^^}")" \
        --x_label "Graph" --y_label "Portion" --save "${ONE_ERROR_PLOT}" \
        --x_sideways --hline_median trivial \
        "${PLOT_PARAMS[@]}"
        
    ./boxplot.py "${SINGLE_MAPPING_FILE}" \
        --title "$(printf "Uniquely mapped (<=2 mismatches)\nreads in ${REGION^^}")" \
        --x_label "Graph" --y_label "Portion uniquely mapped" --save "${SINGLE_MAPPING_PLOT}" \
        --x_sideways --hline_median refonly \
        "${PLOT_PARAMS[@]}"
        
    ./boxplot.py "${RUNTIME_FILE}" \
        --title "$(printf "Per-read runtime\n in ${REGION^^}")" \
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
./boxplot.py "${OVERALL_MAPPING_FILE}" \
    --title "$(printf "Mapped (<=2 mismatches)\nreads")" \
    --x_label "Graph" --y_label "Portion mapped" --save "${OVERALL_MAPPING_PLOT}" \
    --x_sideways  --hline_median trivial \
    "${PLOT_PARAMS[@]}"
    
./boxplot.py "${OVERALL_PERFECT_FILE}" \
    --title "$(printf "Perfectly mapped\nreads")" \
    --x_label "Graph" --y_label "Portion perfectly mapped" --save "${OVERALL_PERFECT_PLOT}" \
    --x_sideways --hline_median trivial \
    "${PLOT_PARAMS[@]}"
    
./boxplot.py "${OVERALL_ONE_ERROR_FILE}" \
    --title "$(printf "One-error (<=1 mismatch)\nreads")" \
    --x_label "Graph" --y_label "Portion" --save "${OVERALL_ONE_ERROR_PLOT}" \
    --x_sideways --hline_median trivial \
    "${PLOT_PARAMS[@]}"

./boxplot.py "${OVERALL_SINGLE_MAPPING_FILE}" \
    --title "$(printf "Uniquely mapped (<=2 mismatches)\nreads")" \
    --x_label "Graph" --y_label "Portion uniquely mapped" --save "${OVERALL_SINGLE_MAPPING_PLOT}" \
    --x_sideways --hline_median refonly \
    "${PLOT_PARAMS[@]}"

