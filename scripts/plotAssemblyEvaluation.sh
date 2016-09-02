#!/usr/bin/env bash
# Makes plots comparing the assembly realignment statistics in each region.

set -ex

# What plot filetype should we produce?
PLOT_FILETYPE="png"

# Grab the input directory to look in
INPUT_DIR=${1}

if [[ ! -d "${INPUT_DIR}" ]]
then
    echo "Specify input directory!"
    exit 1
fi

# Set to "call" or "genotype" for comparison experiment
PARAM_SET="defray"

# What evaluation are we?
EVAL="assembly_sd"

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
    platypus
    freebayes
    samtools
    empty
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
    Platypus
    Freebayes
    Samtools
    Empty
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
    "#9E7C72"
    "#FF00FF"
    "#2F4F4F"
    "#FFFF00"
    --dpi 90 --font_size 14 --no_n
)

# Set up combined stats vars. They have to be associative arrays because we need
# to stick stuff in for different graphs.
declare -A ALL_INSERT_BASES
declare -A ALL_DELETE_BASES
declare -A ALL_SUBSTITUTE_BASES
declare -A ALL_UNVISITED_BASES

declare -A ALL_INSERT_COUNT
declare -A ALL_DELETE_COUNT
declare -A ALL_SUBSTITUTE_COUNT
declare -A ALL_UNVISITED_COUNT

for GRAPH in empty snp1kg refonly shifted1kg freebayes platypus samtools; do
    # Fill global stats arrays with 0
    ALL_INSERT_BASES[$GRAPH]=0
    ALL_DELETE_BASES[$GRAPH]=0
    ALL_SUBSTITUTE_BASES[$GRAPH]=0
    ALL_UNVISITED_BASES[$GRAPH]=0
    
    ALL_INSERT_COUNT[$GRAPH]=0
    ALL_DELETE_COUNT[$GRAPH]=0
    ALL_SUBSTITUTE_COUNT[$GRAPH]=0
    ALL_UNVISITED_COUNT[$GRAPH]=0
done

# Clear out old stats
true > "${INPUT_DIR}/evals/${EVAL}/plots/bp/stats/ALL-insertions-${PARAM_SET}.tsv"
true > "${INPUT_DIR}/evals/${EVAL}/plots/bp/stats/ALL-deletions-${PARAM_SET}.tsv"
true > "${INPUT_DIR}/evals/${EVAL}/plots/bp/stats/ALL-substitutions-${PARAM_SET}.tsv"
true > "${INPUT_DIR}/evals/${EVAL}/plots/bp/stats/ALL-unvisited-${PARAM_SET}.tsv"

true > "${INPUT_DIR}/evals/${EVAL}/plots/count/stats/ALL-insertions-${PARAM_SET}.tsv"
true > "${INPUT_DIR}/evals/${EVAL}/plots/count/stats/ALL-deletions-${PARAM_SET}.tsv"
true > "${INPUT_DIR}/evals/${EVAL}/plots/count/stats/ALL-substitutions-${PARAM_SET}.tsv"
true > "${INPUT_DIR}/evals/${EVAL}/plots/count/stats/ALL-unvisited-${PARAM_SET}.tsv"


for REGION in brca1 brca2 lrc_kir mhc; do
      
    mkdir -p "${INPUT_DIR}/evals/${EVAL}/plots/bp"  
    mkdir -p "${INPUT_DIR}/evals/${EVAL}/plots/count"
    
    # Keep stats files in their own directories
    mkdir -p "${INPUT_DIR}/evals/${EVAL}/plots/bp/stats"  
    mkdir -p "${INPUT_DIR}/evals/${EVAL}/plots/count/stats"
    
    
    # Clear out old stats
    true > "${INPUT_DIR}/evals/${EVAL}/plots/bp/stats/${REGION}-insertions-${PARAM_SET}.tsv"
    true > "${INPUT_DIR}/evals/${EVAL}/plots/bp/stats/${REGION}-deletions-${PARAM_SET}.tsv"
    true > "${INPUT_DIR}/evals/${EVAL}/plots/bp/stats/${REGION}-substitutions-${PARAM_SET}.tsv"
    true > "${INPUT_DIR}/evals/${EVAL}/plots/bp/stats/${REGION}-unvisited-${PARAM_SET}.tsv"
    
    true > "${INPUT_DIR}/evals/${EVAL}/plots/count/stats/${REGION}-insertions-${PARAM_SET}.tsv"
    true > "${INPUT_DIR}/evals/${EVAL}/plots/count/stats/${REGION}-deletions-${PARAM_SET}.tsv"
    true > "${INPUT_DIR}/evals/${EVAL}/plots/count/stats/${REGION}-substitutions-${PARAM_SET}.tsv"
    true > "${INPUT_DIR}/evals/${EVAL}/plots/count/stats/${REGION}-unvisited-${PARAM_SET}.tsv"


    for GRAPH in empty snp1kg refonly shifted1kg freebayes platypus samtools; do
    
        # Parse all the counts from the stats files
        
        GRAPH_STATS="${INPUT_DIR}/evals/${EVAL}/stats/${REGION}/${GRAPH}-${PARAM_SET}.txt"
    
        if [[ -e ${GRAPH_STATS} ]]; then
        
            # Do plot input files by base count
            INSERT_BASES=$(cat "${GRAPH_STATS}" | grep "Insertions" | sed -E 's/.* ([0-9]+) bp.*/\1/g')
            DELETE_BASES=$(cat "${GRAPH_STATS}" | grep "Deletions" | sed -E 's/.* ([0-9]+) bp.*/\1/g')
            SUBSTITUTE_BASES=$(cat "${GRAPH_STATS}" | grep "Substitutions" | sed -E 's/.* ([0-9]+) bp.*/\1/g')
            UNVISITED_BASES=$(cat "${GRAPH_STATS}" | grep "Unvisited" | sed -E 's/.* \(([0-9]+) bp\).*/\1/g')
            
            ALL_INSERT_BASES[$GRAPH]=$((ALL_INSERT_BASES[$GRAPH] + INSERT_BASES))
            ALL_DELETE_BASES[$GRAPH]=$((ALL_DELETE_BASES[$GRAPH] + DELETE_BASES))
            ALL_SUBSTITUTE_BASES[$GRAPH]=$((ALL_SUBSTITUTE_BASES[$GRAPH] + SUBSTITUTE_BASES))
            ALL_UNVISITED_BASES[$GRAPH]=$((ALL_UNVISITED_BASES[$GRAPH] + UNVISITED_BASES))
            
            printf "${GRAPH}\t${INSERT_BASES}\n" \
                >> "${INPUT_DIR}/evals/${EVAL}/plots/bp/stats/${REGION}-insertions-${PARAM_SET}.tsv"
            printf "${GRAPH}\t${DELETE_BASES}\n" \
                >> "${INPUT_DIR}/evals/${EVAL}/plots/bp/stats/${REGION}-deletions-${PARAM_SET}.tsv"
            printf "${GRAPH}\t${SUBSTITUTE_BASES}\n" \
                >> "${INPUT_DIR}/evals/${EVAL}/plots/bp/stats/${REGION}-substitutions-${PARAM_SET}.tsv"
            printf "${GRAPH}\t${UNVISITED_BASES}\n" \
                >> "${INPUT_DIR}/evals/${EVAL}/plots/bp/stats/${REGION}-unvisited-${PARAM_SET}.tsv"
            
            # Now do by event count
            INSERT_COUNT=$(cat "${GRAPH_STATS}" | grep "Insertions" | sed -E 's/.* ([0-9]+) read events.*/\1/g')
            DELETE_COUNT=$(cat "${GRAPH_STATS}" | grep "Deletions" | sed -E 's/.* ([0-9]+) read events.*/\1/g')
            SUBSTITUTE_COUNT=$(cat "${GRAPH_STATS}" | grep "Substitutions" | sed -E 's/.* ([0-9]+) read events.*/\1/g')
            UNVISITED_COUNT=$(cat "${GRAPH_STATS}" | grep "Unvisited" | sed -E 's/.* ([0-9]+)\/.*/\1/g')
            
            ALL_INSERT_COUNT[$GRAPH]=$((ALL_INSERT_COUNT[$GRAPH] + INSERT_COUNT))
            ALL_DELETE_COUNT[$GRAPH]=$((ALL_DELETE_COUNT[$GRAPH] + DELETE_COUNT))
            ALL_SUBSTITUTE_COUNT[$GRAPH]=$((ALL_SUBSTITUTE_COUNT[$GRAPH] + SUBSTITUTE_COUNT))
            ALL_UNVISITED_COUNT[$GRAPH]=$((ALL_UNVISITED_COUNT[$GRAPH] + UNVISITED_COUNT))
            
            printf "${GRAPH}\t${INSERT_COUNT}\n" \
                >> "${INPUT_DIR}/evals/${EVAL}/plots/count/stats/${REGION}-insertions-${PARAM_SET}.tsv"
            printf "${GRAPH}\t${DELETE_COUNT}\n" \
                >> "${INPUT_DIR}/evals/${EVAL}/plots/count/stats/${REGION}-deletions-${PARAM_SET}.tsv"
            printf "${GRAPH}\t${SUBSTITUTE_COUNT}\n" \
                >> "${INPUT_DIR}/evals/${EVAL}/plots/count/stats/${REGION}-substitutions-${PARAM_SET}.tsv"
            printf "${GRAPH}\t${UNVISITED_COUNT}\n" \
                >> "${INPUT_DIR}/evals/${EVAL}/plots/count/stats/${REGION}-unvisited-${PARAM_SET}.tsv"
                
        fi

    done
    
    # Plot by bp
    ./scripts/barchart.py "${INPUT_DIR}/evals/${EVAL}/plots/bp/stats/${REGION}-insertions-${PARAM_SET}.tsv" \
        --title "Inserted bases relative to sample ${REGION^^} (${PARAM_SET})" \
        --x_label "Graph type" \
        --y_label "Inserted bases in assembly re-alignment" \
        --save "${INPUT_DIR}/evals/${EVAL}/plots/bp/${REGION}-insertions-${PARAM_SET}.${PLOT_FILETYPE}" \
        "${PLOT_PARAMS[@]}"
        
    ./scripts/barchart.py "${INPUT_DIR}/evals/${EVAL}/plots/bp/stats/${REGION}-deletions-${PARAM_SET}.tsv" \
        --title "Deleted bases relative to sample ${REGION^^} (${PARAM_SET})" \
        --x_label "Graph type" \
        --y_label "Deleted bases in assembly re-alignment" \
        --save "${INPUT_DIR}/evals/${EVAL}/plots/bp/${REGION}-deletions-${PARAM_SET}.${PLOT_FILETYPE}" \
        "${PLOT_PARAMS[@]}"
        
    ./scripts/barchart.py "${INPUT_DIR}/evals/${EVAL}/plots/bp/stats/${REGION}-substitutions-${PARAM_SET}.tsv" \
        --title "Substituted bases relative to sample ${REGION^^} (${PARAM_SET})" \
        --x_label "Graph type" \
        --y_label "Substituted bases in assembly re-alignment" \
        --save "${INPUT_DIR}/evals/${EVAL}/plots/bp/${REGION}-substitutions-${PARAM_SET}.${PLOT_FILETYPE}" \
        "${PLOT_PARAMS[@]}"
        
    ./scripts/barchart.py "${INPUT_DIR}/evals/${EVAL}/plots/bp/stats/${REGION}-unvisited-${PARAM_SET}.tsv" \
        --title "Unvisited node length in sample ${REGION^^} (${PARAM_SET})" \
        --x_label "Graph type" \
        --y_label "Length of unvisited called nodes" \
        --save "${INPUT_DIR}/evals/${EVAL}/plots/bp/${REGION}-unvisited-${PARAM_SET}.${PLOT_FILETYPE}" \
        "${PLOT_PARAMS[@]}"
        
    # Plot by event count
    ./scripts/barchart.py "${INPUT_DIR}/evals/${EVAL}/plots/count/stats/${REGION}-insertions-${PARAM_SET}.tsv" \
        --title "Insertion events relative to sample ${REGION^^} (${PARAM_SET})" \
        --x_label "Graph type" \
        --y_label "Insertions in assembly re-alignment" \
        --save "${INPUT_DIR}/evals/${EVAL}/plots/count/${REGION}-insertions-${PARAM_SET}.${PLOT_FILETYPE}" \
        "${PLOT_PARAMS[@]}"
        
    ./scripts/barchart.py "${INPUT_DIR}/evals/${EVAL}/plots/count/stats/${REGION}-deletions-${PARAM_SET}.tsv" \
        --title "Deletion events relative to sample ${REGION^^} (${PARAM_SET})" \
        --x_label "Graph type" \
        --y_label "Deletions in assembly re-alignment" \
        --save "${INPUT_DIR}/evals/${EVAL}/plots/count/${REGION}-deletions-${PARAM_SET}.${PLOT_FILETYPE}" \
        "${PLOT_PARAMS[@]}"
        
    ./scripts/barchart.py "${INPUT_DIR}/evals/${EVAL}/plots/count/stats/${REGION}-substitutions-${PARAM_SET}.tsv" \
        --title "Substitution events relative to sample ${REGION^^} (${PARAM_SET})" \
        --x_label "Graph type" \
        --y_label "Substitutions in assembly re-alignment" \
        --save "${INPUT_DIR}/evals/${EVAL}/plots/count/${REGION}-substitutions-${PARAM_SET}.${PLOT_FILETYPE}" \
        "${PLOT_PARAMS[@]}"
        
    ./scripts/barchart.py "${INPUT_DIR}/evals/${EVAL}/plots/count/stats/${REGION}-unvisited-${PARAM_SET}.tsv" \
        --title "Unvisited nodes in sample ${REGION^^} (${PARAM_SET})" \
        --x_label "Graph type" \
        --y_label "Number of unvisited called nodes" \
        --save "${INPUT_DIR}/evals/${EVAL}/plots/count/${REGION}-unvisited-${PARAM_SET}.${PLOT_FILETYPE}" \
        "${PLOT_PARAMS[@]}"
    
done

for GRAPH in "${!ALL_INSERT_BASES[@]}"; do

    # Make overall stats files
    printf "${GRAPH}\t${ALL_INSERT_BASES[$GRAPH]}\n" \
        >> "${INPUT_DIR}/evals/${EVAL}/plots/bp/stats/ALL-insertions-${PARAM_SET}.tsv"
    printf "${GRAPH}\t${ALL_DELETE_BASES[$GRAPH]}\n" \
        >> "${INPUT_DIR}/evals/${EVAL}/plots/bp/stats/ALL-deletions-${PARAM_SET}.tsv"
    printf "${GRAPH}\t${ALL_SUBSTITUTE_BASES[$GRAPH]}\n" \
        >> "${INPUT_DIR}/evals/${EVAL}/plots/bp/stats/ALL-substitutions-${PARAM_SET}.tsv"
    printf "${GRAPH}\t${ALL_UNVISITED_BASES[$GRAPH]}\n" \
        >> "${INPUT_DIR}/evals/${EVAL}/plots/bp/stats/ALL-unvisited-${PARAM_SET}.tsv"
        
    printf "${GRAPH}\t${ALL_INSERT_COUNT[$GRAPH]}\n" \
        >> "${INPUT_DIR}/evals/${EVAL}/plots/count/stats/ALL-insertions-${PARAM_SET}.tsv"
    printf "${GRAPH}\t${ALL_DELETE_COUNT[$GRAPH]}\n" \
        >> "${INPUT_DIR}/evals/${EVAL}/plots/count/stats/ALL-deletions-${PARAM_SET}.tsv"
    printf "${GRAPH}\t${ALL_SUBSTITUTE_COUNT[$GRAPH]}\n" \
        >> "${INPUT_DIR}/evals/${EVAL}/plots/count/stats/ALL-substitutions-${PARAM_SET}.tsv"
    printf "${GRAPH}\t${ALL_UNVISITED_COUNT[$GRAPH]}\n" \
        >> "${INPUT_DIR}/evals/${EVAL}/plots/count/stats/ALL-unvisited-${PARAM_SET}.tsv"

done

# Plot totals by bp
./scripts/barchart.py "${INPUT_DIR}/evals/${EVAL}/plots/bp/stats/ALL-insertions-${PARAM_SET}.tsv" \
    --title "Inserted bases relative to sample overall (${PARAM_SET})" \
    --x_label "Graph type" \
    --y_label "Inserted bases in assembly re-alignment" \
    --save "${INPUT_DIR}/evals/${EVAL}/plots/bp/ALL-insertions-${PARAM_SET}.${PLOT_FILETYPE}" \
    "${PLOT_PARAMS[@]}"
    
./scripts/barchart.py "${INPUT_DIR}/evals/${EVAL}/plots/bp/stats/ALL-deletions-${PARAM_SET}.tsv" \
    --title "Deleted bases relative to sample overall (${PARAM_SET})" \
    --x_label "Graph type" \
    --y_label "Deleted bases in assembly re-alignment" \
    --save "${INPUT_DIR}/evals/${EVAL}/plots/bp/ALL-deletions-${PARAM_SET}.${PLOT_FILETYPE}" \
    "${PLOT_PARAMS[@]}"
    
./scripts/barchart.py "${INPUT_DIR}/evals/${EVAL}/plots/bp/stats/ALL-substitutions-${PARAM_SET}.tsv" \
    --title "Substituted bases relative to sample overall (${PARAM_SET})" \
    --x_label "Graph type" \
    --y_label "Substituted bases in assembly re-alignment" \
    --save "${INPUT_DIR}/evals/${EVAL}/plots/bp/ALL-substitutions-${PARAM_SET}.${PLOT_FILETYPE}" \
    "${PLOT_PARAMS[@]}"
    
./scripts/barchart.py "${INPUT_DIR}/evals/${EVAL}/plots/bp/stats/ALL-unvisited-${PARAM_SET}.tsv" \
    --title "Unvisited node length in sample overall (${PARAM_SET})" \
    --x_label "Graph type" \
    --y_label "Length of unvisited called nodes" \
    --save "${INPUT_DIR}/evals/${EVAL}/plots/bp/ALL-unvisited-${PARAM_SET}.${PLOT_FILETYPE}" \
    "${PLOT_PARAMS[@]}"
    
# Plot totals by event count
./scripts/barchart.py "${INPUT_DIR}/evals/${EVAL}/plots/count/stats/ALL-insertions-${PARAM_SET}.tsv" \
    --title "Insertion events relative to sample overall (${PARAM_SET})" \
    --x_label "Graph type" \
    --y_label "Insertions in assembly re-alignment" \
    --save "${INPUT_DIR}/evals/${EVAL}/plots/count/ALL-insertions-${PARAM_SET}.${PLOT_FILETYPE}" \
    "${PLOT_PARAMS[@]}"
    
./scripts/barchart.py "${INPUT_DIR}/evals/${EVAL}/plots/count/stats/ALL-deletions-${PARAM_SET}.tsv" \
    --title "Deletion events relative to sample overall (${PARAM_SET})" \
    --x_label "Graph type" \
    --y_label "Deletions in assembly re-alignment" \
    --save "${INPUT_DIR}/evals/${EVAL}/plots/count/ALL-deletions-${PARAM_SET}.${PLOT_FILETYPE}" \
    "${PLOT_PARAMS[@]}"
    
./scripts/barchart.py "${INPUT_DIR}/evals/${EVAL}/plots/count/stats/ALL-substitutions-${PARAM_SET}.tsv" \
    --title "Substitution events relative to sample overall (${PARAM_SET})" \
    --x_label "Graph type" \
    --y_label "Substitutions in assembly re-alignment" \
    --save "${INPUT_DIR}/evals/${EVAL}/plots/count/ALL-substitutions-${PARAM_SET}.${PLOT_FILETYPE}" \
    "${PLOT_PARAMS[@]}"
    
./scripts/barchart.py "${INPUT_DIR}/evals/${EVAL}/plots/count/stats/ALL-unvisited-${PARAM_SET}.tsv" \
    --title "Unvisited nodes in sample overall (${PARAM_SET})" \
    --x_label "Graph type" \
    --y_label "Number of unvisited called nodes" \
    --save "${INPUT_DIR}/evals/${EVAL}/plots/count/ALL-unvisited-${PARAM_SET}.${PLOT_FILETYPE}" \
    "${PLOT_PARAMS[@]}"
    
        
