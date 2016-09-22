#!/usr/bin/env bash
# Run after biasDetector.py.
# Makes plots comparing the populations for each graph.

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
    --font_size 20 --dpi 90 --no_n
)

# Set up proper graph names
declare -A HR_NAMES

HR_NAMES["snp1kg"]="1KG"
HR_NAMES["snp1000g"]="1KG"
HR_NAMES["haplo1kg30"]="1KG Haplo 30"
HR_NAMES["haplo1kg50"]="1KG Haplo"
HR_NAMES["sbg"]="7BG"
HR_NAMES["cactus"]="Cactus"
HR_NAMES["camel"]="Camel"
HR_NAMES["curoverse"]="Curoverse"
HR_NAMES["debruijn-k31"]="De Bruijn 31"
HR_NAMES["debruijn-k63"]="De Bruijn 63"
HR_NAMES["level1"]="Level1"
HR_NAMES["level2"]="Level2"
HR_NAMES["level3"]="Level3"
HR_NAMES["prg"]="PRG"
HR_NAMES["refonly"]="Reference"
HR_NAMES["simons"]="SGDP"
HR_NAMES["trivial"]="Trivial"
HR_NAMES["vglr"]="VGLR"
HR_NAMES["shifted1kg"]="Scrambled"
HR_NAMES["snp1kg_kp"]="1KG KP"

for MODE in normalized_distributions distributions
do

    # Are we absolute or normalized?
    NORMALIZED="absolute"

    # Do we want portion (default) or absolute deviations
    DEVIATIONS=""

    PORTION="Portion"
    if [ "${MODE}" == "normalized_distributions" ]
    then
        PORTION="Difference in portion"
        NORMALIZED="normalized"
        # No sense normalizing twice
        DEVIATIONS="--absolute_deviation"
    fi

    # Where are the bias files
    DISTRIBUTION_DIR="${INPUT_DIR}/bias/${MODE}"

    for REGION in `ls "${DISTRIBUTION_DIR}" | grep -v '\.png$' | grep -v '\.svg$'`
    do
        # For every region, plot every graph
        
        # Remove underscores from region names to make them human readable
        HR_REGION=`echo ${REGION^^} | sed 's/_/ /g'`
        
        REGION_DIR="${DISTRIBUTION_DIR}/${REGION}"
        
        for GRAPH_TSV in `ls "${REGION_DIR}" | grep '\.tsv$'`
        do
            # Every TSV becomes a boxplot
            
            # Pull out the graph name
            GRAPH=`basename ${GRAPH_TSV} | sed 's/\(.*\).tsv/\1/'` 
            
            # Ge tthe actual path to the graph TSV
            GRAPH_TSV_PATH="${REGION_DIR}/${GRAPH_TSV}"
            
            # Where should the plot go?
            PLOT_PATH="${DISTRIBUTION_DIR}/${NORMALIZED}_bias_${REGION}_${GRAPH}.svg"
            
            # Get the human readable graph name
            HR_GRAPH=${HR_NAMES["${GRAPH}"]}
            
            echo "Plotting ${HR_REGION} ${HR_GRAPH} (${GRAPH})"
            
            ./scripts/boxplot.py "${GRAPH_TSV_PATH}" \
                --title "$(printf "Perfectly mapped\nreads in ${HR_REGION} ${HR_GRAPH}")" \
                --x_label "Population" --y_label "${PORTION} mapped" --save "${PLOT_PATH}" \
                --x_sideways --hline_median EUR \
                ${DEVIATIONS} \
                "${PLOT_PARAMS[@]}"
        
        done

    done
done


        
