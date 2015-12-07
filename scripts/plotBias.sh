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
HR_NAMES["haplo1kg"]="1KG Haplo"
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

# Where are the bias files
DISTRIBUTION_DIR="${INPUT_DIR}/bias/normalized_distributions"

for REGION in `ls "${DISTRIBUTION_DIR}"`
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
        PLOT_PATH="${REGION_DIR}/${GRAPH}.png"
        
        # Get the human readable graph name
        HR_GRAPH=${HR_NAMES["${GRAPH}"]}
        
        echo "Plotting ${HR_REGION} ${HR_GRAPH} (${GRAPH})"
        
        ./scripts/boxplot.py "${GRAPH_TSV_PATH}" \
            --title "$(printf "Mapped (<=2 mismatches)\nreads in ${HR_REGION} ${HR_GRAPH}")" \
            --x_label "Population" --y_label "Relative portion mapped" --save "${PLOT_PATH}" \
            --x_sideways \
            "${PLOT_PARAMS[@]}"
    
    done

done


        
