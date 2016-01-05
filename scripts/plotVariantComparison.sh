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
    scripts/barchart.py "${DROPPED_TSV}" --save "${DROPPED_PNG}" --title "Bases Dropped in ${REGION^^}" --x_label "Graph" --y_label "Bases" --x_sideways "${PLOT_PARAMS[@]}"
    
done
