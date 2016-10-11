#!/usr/bin/env bash
# Makes plots comparing the assembly realignment statistics in each region.

set -eux

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
    "GRCh38"
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
    "#edd91b"
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

declare -A ALL_FALSE_POSITIVES
declare -A ALL_FALSE_NEGATIVES
declare -A ALL_TRUE_POSITIVES

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
    
    ALL_FALSE_POSITIVES[$GRAPH]=0
    ALL_FALSE_NEGATIVES[$GRAPH]=0
    ALL_TRUE_POSITIVES[$GRAPH]=0
done

# We also need to know the sizes of the regions in bp
declare -A REGION_SIZES
REGION_SIZES["brca1"]=81189
REGION_SIZES["brca2"]=84989
REGION_SIZES["lrc_kir"]=1058685
REGION_SIZES["mhc"]=4970458
REGION_SIZES["sma"]=2397625

# Make plot directories
mkdir -p "${INPUT_DIR}/evals/${EVAL}/plots/bp"  
mkdir -p "${INPUT_DIR}/evals/${EVAL}/plots/count"

# Keep stats files in their own directories
mkdir -p "${INPUT_DIR}/evals/${EVAL}/plots/bp/stats"  
mkdir -p "${INPUT_DIR}/evals/${EVAL}/plots/count/stats"

# Clear out old stats
true > "${INPUT_DIR}/evals/${EVAL}/plots/bp/stats/ALL-insertions-${PARAM_SET}.tsv"
true > "${INPUT_DIR}/evals/${EVAL}/plots/bp/stats/ALL-deletions-${PARAM_SET}.tsv"
true > "${INPUT_DIR}/evals/${EVAL}/plots/bp/stats/ALL-substitutions-${PARAM_SET}.tsv"
true > "${INPUT_DIR}/evals/${EVAL}/plots/bp/stats/ALL-unvisited-${PARAM_SET}.tsv"
true > "${INPUT_DIR}/evals/${EVAL}/plots/bp/stats/ALL-indels-${PARAM_SET}.tsv"

true > "${INPUT_DIR}/evals/${EVAL}/plots/bp/stats/PR-ALL-${PARAM_SET}.tsv"
true > "${INPUT_DIR}/evals/${EVAL}/plots/bp/stats/F1-ALL-${PARAM_SET}.tsv"

true > "${INPUT_DIR}/evals/${EVAL}/plots/count/stats/ALL-insertions-${PARAM_SET}.tsv"
true > "${INPUT_DIR}/evals/${EVAL}/plots/count/stats/ALL-deletions-${PARAM_SET}.tsv"
true > "${INPUT_DIR}/evals/${EVAL}/plots/count/stats/ALL-substitutions-${PARAM_SET}.tsv"
true > "${INPUT_DIR}/evals/${EVAL}/plots/count/stats/ALL-unvisited-${PARAM_SET}.tsv"
true > "${INPUT_DIR}/evals/${EVAL}/plots/count/stats/ALL-indels-${PARAM_SET}.tsv"

# Count up otoal length of all regions
TOTAL_LENGTH=0

# Regions with good assembly coverage:
# brca1 brca2 lrc_kir mhc
for REGION in  brca1 brca2 mhc lrc_kir; do
      
    # Clear out old stats
    true > "${INPUT_DIR}/evals/${EVAL}/plots/bp/stats/${REGION}-insertions-${PARAM_SET}.tsv"
    true > "${INPUT_DIR}/evals/${EVAL}/plots/bp/stats/${REGION}-deletions-${PARAM_SET}.tsv"
    true > "${INPUT_DIR}/evals/${EVAL}/plots/bp/stats/${REGION}-substitutions-${PARAM_SET}.tsv"
    true > "${INPUT_DIR}/evals/${EVAL}/plots/bp/stats/${REGION}-unvisited-${PARAM_SET}.tsv"
    true > "${INPUT_DIR}/evals/${EVAL}/plots/bp/stats/${REGION}-indels-${PARAM_SET}.tsv"
    
    true > "${INPUT_DIR}/evals/${EVAL}/plots/bp/stats/PR-${REGION}-${PARAM_SET}.tsv"
    true > "${INPUT_DIR}/evals/${EVAL}/plots/bp/stats/F1-${REGION}-${PARAM_SET}.tsv"
    
    true > "${INPUT_DIR}/evals/${EVAL}/plots/count/stats/${REGION}-insertions-${PARAM_SET}.tsv"
    true > "${INPUT_DIR}/evals/${EVAL}/plots/count/stats/${REGION}-deletions-${PARAM_SET}.tsv"
    true > "${INPUT_DIR}/evals/${EVAL}/plots/count/stats/${REGION}-substitutions-${PARAM_SET}.tsv"
    true > "${INPUT_DIR}/evals/${EVAL}/plots/count/stats/${REGION}-unvisited-${PARAM_SET}.tsv"
    true > "${INPUT_DIR}/evals/${EVAL}/plots/count/stats/${REGION}-indels-${PARAM_SET}.tsv"

    # Add region length in to total
    TOTAL_LENGTH=$((TOTAL_LENGTH + REGION_SIZES[$REGION]))

    for GRAPH in empty snp1kg refonly shifted1kg cactus prg freebayes platypus samtools; do
    
        # Parse all the counts from the stats files
        
        GRAPH_STATS="${INPUT_DIR}/evals/${EVAL}/stats/${REGION}/${GRAPH}-${PARAM_SET}.txt"
        
        # We need to rerun the stats collection, in case we want new stats.
        RTEMP="${INPUT_DIR}/evals/${EVAL}/temp/${PARAM_SET}/${REGION}"
        TEMP="${RTEMP}/graph/${GRAPH}"
        vg stats "${TEMP}/sample.vg" -v -a "${TEMP}/assembly_aligned.gam" > "${GRAPH_STATS}" 2>&1
    
        if [[ -e ${GRAPH_STATS} ]]; then
        
            # Do plot input files by base count
            INSERT_BASES=$(cat "${GRAPH_STATS}" | grep "Insertions" | sed -E 's/.* ([0-9]+) bp.*/\1/g')
            DELETE_BASES=$(cat "${GRAPH_STATS}" | grep "Deletions" | sed -E 's/.* ([0-9]+) bp.*/\1/g')
            SUBSTITUTE_BASES=$(cat "${GRAPH_STATS}" | grep "Substitutions" | sed -E 's/.* ([0-9]+) bp.*/\1/g')
            UNVISITED_BASES=$(cat "${GRAPH_STATS}" | grep "Unvisited" | sed -E 's/.* \(([0-9]+) bp\).*/\1/g')
            SINGLE_VISITED_BASES=$(cat "${GRAPH_STATS}" | grep "Unvisited" | sed -E 's/.* \(([0-9]+) bp\).*/\1/g')
            
            ALL_INSERT_BASES[$GRAPH]=$((ALL_INSERT_BASES[$GRAPH] + INSERT_BASES))
            ALL_DELETE_BASES[$GRAPH]=$((ALL_DELETE_BASES[$GRAPH] + DELETE_BASES))
            ALL_SUBSTITUTE_BASES[$GRAPH]=$((ALL_SUBSTITUTE_BASES[$GRAPH] + SUBSTITUTE_BASES))
            ALL_UNVISITED_BASES[$GRAPH]=$((ALL_UNVISITED_BASES[$GRAPH] + UNVISITED_BASES))
            
            # Normalize base counts by region size for plotting
            INSERT_PORTION=$(echo "${INSERT_BASES} / ${REGION_SIZES[$REGION]}" | bc -l)
            DELETE_PORTION=$(echo "${DELETE_BASES} / ${REGION_SIZES[$REGION]}" | bc -l)
            SUBSTITUTE_PORTION=$(echo "${SUBSTITUTE_BASES} / ${REGION_SIZES[$REGION]}" | bc -l)
            UNVISITED_PORTION=$(echo "${UNVISITED_BASES} / ${REGION_SIZES[$REGION]}" | bc -l)
            
            # Combine indels
            INDEL_PORTION=$(echo "${INSERT_PORTION} + ${DELETE_PORTION}" | bc -l)
            
            printf "${GRAPH}\t${INSERT_PORTION}\n" \
                >> "${INPUT_DIR}/evals/${EVAL}/plots/bp/stats/${REGION}-insertions-${PARAM_SET}.tsv"
            printf "${GRAPH}\t${DELETE_PORTION}\n" \
                >> "${INPUT_DIR}/evals/${EVAL}/plots/bp/stats/${REGION}-deletions-${PARAM_SET}.tsv"
            printf "${GRAPH}\t${SUBSTITUTE_PORTION}\n" \
                >> "${INPUT_DIR}/evals/${EVAL}/plots/bp/stats/${REGION}-substitutions-${PARAM_SET}.tsv"
            printf "${GRAPH}\t${UNVISITED_PORTION}\n" \
                >> "${INPUT_DIR}/evals/${EVAL}/plots/bp/stats/${REGION}-unvisited-${PARAM_SET}.tsv"
                
            printf "${GRAPH}\t${INDEL_PORTION}\n" \
                >> "${INPUT_DIR}/evals/${EVAL}/plots/bp/stats/${REGION}-indels-${PARAM_SET}.tsv"
            
            # Now do by event count
            INSERT_COUNT=$(cat "${GRAPH_STATS}" | grep "Insertions" | sed -E 's/.* ([0-9]+) read events.*/\1/g')
            DELETE_COUNT=$(cat "${GRAPH_STATS}" | grep "Deletions" | sed -E 's/.* ([0-9]+) read events.*/\1/g')
            SUBSTITUTE_COUNT=$(cat "${GRAPH_STATS}" | grep "Substitutions" | sed -E 's/.* ([0-9]+) read events.*/\1/g')
            UNVISITED_COUNT=$(cat "${GRAPH_STATS}" | grep "Unvisited" | sed -E 's/.* ([0-9]+)\/.*/\1/g')
            
            ALL_INSERT_COUNT[$GRAPH]=$((ALL_INSERT_COUNT[$GRAPH] + INSERT_COUNT))
            ALL_DELETE_COUNT[$GRAPH]=$((ALL_DELETE_COUNT[$GRAPH] + DELETE_COUNT))
            ALL_SUBSTITUTE_COUNT[$GRAPH]=$((ALL_SUBSTITUTE_COUNT[$GRAPH] + SUBSTITUTE_COUNT))
            ALL_UNVISITED_COUNT[$GRAPH]=$((ALL_UNVISITED_COUNT[$GRAPH] + UNVISITED_COUNT))
            
            # Normalize event counts by region size for plotting
            INSERT_FREQUENCY=$(echo "${INSERT_COUNT} / ${REGION_SIZES[$REGION]}" | bc -l)
            DELETE_FREQUENCY=$(echo "${DELETE_COUNT} / ${REGION_SIZES[$REGION]}" | bc -l)
            SUBSTITUTE_FREQUENCY=$(echo "${SUBSTITUTE_COUNT} / ${REGION_SIZES[$REGION]}" | bc -l)
            UNVISITED_FREQUENCY=$(echo "${UNVISITED_COUNT} / ${REGION_SIZES[$REGION]}" | bc -l)
            
            # Combine indels
            INDEL_FREQUENCY=$(echo "${INSERT_FREQUENCY} + ${DELETE_FREQUENCY}" | bc -l)
            
            printf "${GRAPH}\t${INSERT_FREQUENCY}\n" \
                >> "${INPUT_DIR}/evals/${EVAL}/plots/count/stats/${REGION}-insertions-${PARAM_SET}.tsv"
            printf "${GRAPH}\t${DELETE_FREQUENCY}\n" \
                >> "${INPUT_DIR}/evals/${EVAL}/plots/count/stats/${REGION}-deletions-${PARAM_SET}.tsv"
            printf "${GRAPH}\t${SUBSTITUTE_FREQUENCY}\n" \
                >> "${INPUT_DIR}/evals/${EVAL}/plots/count/stats/${REGION}-substitutions-${PARAM_SET}.tsv"
            printf "${GRAPH}\t${UNVISITED_FREQUENCY}\n" \
                >> "${INPUT_DIR}/evals/${EVAL}/plots/count/stats/${REGION}-unvisited-${PARAM_SET}.tsv"
                
            printf "${GRAPH}\t${INDEL_FREQUENCY}\n" \
                >> "${INPUT_DIR}/evals/${EVAL}/plots/count/stats/${REGION}-indels-${PARAM_SET}.tsv"
                
            # Get the length of the sample graph. TODO: use this as the
            # denominator instead for portions and rates of events?
            GRAPH_BASES=$(vg stats -l "${TEMP}/sample.vg" | cut -f2)
            
            # How many bases are in nodes visited twice?
            DOUBLE_VISITED_NODE_BASES=$((GRAPH_BASES - UNVISITED_BASES - SINGLE_VISITED_BASES))
            # How many visitings of nodes are there, weighted by node bases?
            NODE_BASE_VISITS=$((2 * DOUBLE_VISITED_NODE_BASES + SINGLE_VISITED_BASES))
            # Of that, how many actual vitist to bases are there? We need to
            # subtract out deleted and substituted bases, which detract from
            # bases actually touched by the assemblies.
            BASE_VISITS=$((NODE_BASE_VISITS - DELETE_BASES - SUBSTITUTE_BASES))
            
            # That gives us our true positives
            TRUE_POSITIVES=${BASE_VISITS}
            # False positives is unvisited bases, plus those deleted (because we
            # had asserted a copy incorrectly per deletion), plus those
            # substituted (because we had asserted a copy incorrectly per
            # substitution). We assume unvisited bases are only worth 1, because
            # if there wasn't a way around them in the graph (i.e. if they were
            # called homozygous) they would be deletions/substitutions and not
            # just unvisited.
            FALSE_POSITIVES=$((UNVISITED_BASES + DELETE_BASES + SUBSTITUTE_BASES))
            # False negatives are insertions and substitutions: things we didn't
            # call as existing.
            FALSE_NEGATIVES=$((INSERT_BASES + SUBSTITUTE_BASES))
            
            # Add to global arrays
            ALL_FALSE_POSITIVES[$GRAPH]=$((ALL_FALSE_POSITIVES[$GRAPH] + FALSE_POSITIVES))
            ALL_FALSE_NEGATIVES[$GRAPH]=$((ALL_FALSE_NEGATIVES[$GRAPH] + FALSE_NEGATIVES))
            ALL_TRUE_POSITIVES[$GRAPH]=$((ALL_TRUE_POSITIVES[$GRAPH] + TRUE_POSITIVES))
            
            # Calculate precision and recall
            PRECISION=$(echo "${TRUE_POSITIVES} / ( ${TRUE_POSITIVES} + ${FALSE_POSITIVES} )" | bc -l)
            RECALL=$(echo "${TRUE_POSITIVES} / ( ${TRUE_POSITIVES} + ${FALSE_NEGATIVES} )" | bc -l)
            
            # Prepare the plot file for a precision vs. recall plot (where precision is Y)
            printf "${GRAPH}\t${RECALL}\t${PRECISION}\n" >> \
                "${INPUT_DIR}/evals/${EVAL}/plots/bp/stats/PR-${REGION}-${PARAM_SET}.tsv"
                
            # Calculate an inverse F1 score
            F1_SCORE_INV=$(echo "1 - 2 * ( ${PRECISION} * ${RECALL} ) / ( ${PRECISION} + ${RECALL} )" | bc -l)
            
            # Prepare the plot file for an F1 bar chart
            printf "${GRAPH}\t${F1_SCORE_INV}\n" >> \
                "${INPUT_DIR}/evals/${EVAL}/plots/bp/stats/F1-${REGION}-${PARAM_SET}.tsv"
            
        fi

    done
    
    # Plot by bp
    ./scripts/barchart.py "${INPUT_DIR}/evals/${EVAL}/plots/bp/stats/${REGION}-insertions-${PARAM_SET}.tsv" \
        --title "Inserted bases relative to sample ${REGION^^}" \
        --x_label "Graph type" \
        --y_label "Inserted bases per region base" \
        --x_sideways \
        --save "${INPUT_DIR}/evals/${EVAL}/plots/bp/${REGION}-insertions-${PARAM_SET}.${PLOT_FILETYPE}" \
        "${PLOT_PARAMS[@]}"
        
    ./scripts/barchart.py "${INPUT_DIR}/evals/${EVAL}/plots/bp/stats/${REGION}-deletions-${PARAM_SET}.tsv" \
        --title "Deleted bases relative to sample ${REGION^^}" \
        --x_label "Graph type" \
        --y_label "Deleted bases per region base" \
        --x_sideways \
        --save "${INPUT_DIR}/evals/${EVAL}/plots/bp/${REGION}-deletions-${PARAM_SET}.${PLOT_FILETYPE}" \
        "${PLOT_PARAMS[@]}"
        
    ./scripts/barchart.py "${INPUT_DIR}/evals/${EVAL}/plots/bp/stats/${REGION}-substitutions-${PARAM_SET}.tsv" \
        --title "Substituted bases relative to sample ${REGION^^}" \
        --x_label "Graph type" \
        --y_label "Substituted bases per region base" \
        --x_sideways \
        --save "${INPUT_DIR}/evals/${EVAL}/plots/bp/${REGION}-substitutions-${PARAM_SET}.${PLOT_FILETYPE}" \
        "${PLOT_PARAMS[@]}"
        
    ./scripts/barchart.py "${INPUT_DIR}/evals/${EVAL}/plots/bp/stats/${REGION}-unvisited-${PARAM_SET}.tsv" \
        --title "Unvisited bases relative to sample ${REGION^^}" \
        --x_label "Graph type" \
        --y_label "Unvisited bases per region base" \
        --x_sideways \
        --save "${INPUT_DIR}/evals/${EVAL}/plots/bp/${REGION}-unvisited-${PARAM_SET}.${PLOT_FILETYPE}" \
        "${PLOT_PARAMS[@]}"
        
    ./scripts/barchart.py "${INPUT_DIR}/evals/${EVAL}/plots/bp/stats/${REGION}-indels-${PARAM_SET}.tsv" \
        --title "Indel bases relative to sample ${REGION^^}" \
        --x_label "Graph type" \
        --y_label "Indel bases per region base" \
        --x_sideways \
        --save "${INPUT_DIR}/evals/${EVAL}/plots/bp/${REGION}-indels-${PARAM_SET}.${PLOT_FILETYPE}" \
        "${PLOT_PARAMS[@]}"
        
    if [[ "${REGION^^}" == "BRCA1" || "${REGION^^}" == "BRCA2" ]]; then
        PR_BOUNDS="--max_x 1 --min_x 0.998 --max_y 1 --min_y 0.998"
    else
        PR_BOUNDS="--max_x 1 --min_x 0.992 --max_y 1 --min_y 0.995"
    fi
        
    ./scripts/scatter.py "${INPUT_DIR}/evals/${EVAL}/plots/bp/stats/PR-${REGION}-${PARAM_SET}.tsv" \
        --title "Precision vs. Recall in ${REGION^^}" \
        --x_label "Recall" \
        --y_label "Precision" \
        --width 6 --height 6 \
        ${PR_BOUNDS} \
        --no_legend --annotate \
        --save "${INPUT_DIR}/evals/${EVAL}/plots/bp/PR-${REGION}-${PARAM_SET}.${PLOT_FILETYPE}" \
        "${PLOT_PARAMS[@]}"
        
    ./scripts/barchart.py "${INPUT_DIR}/evals/${EVAL}/plots/bp/stats/F1-${REGION}-${PARAM_SET}.tsv" \
        --title "F1 score in ${REGION^^}" \
        --x_label "Graph" \
        --y_label "1 - F1" \
        --x_sideways \
        --save "${INPUT_DIR}/evals/${EVAL}/plots/bp/F1-${REGION}-${PARAM_SET}.${PLOT_FILETYPE}" \
        "${PLOT_PARAMS[@]}"
        
    # Plot by event count
    ./scripts/barchart.py "${INPUT_DIR}/evals/${EVAL}/plots/count/stats/${REGION}-insertions-${PARAM_SET}.tsv" \
        --title "Insertion events relative to sample ${REGION^^}" \
        --x_label "Graph type" \
        --y_label "Insertions per region base" \
        --x_sideways \
        --save "${INPUT_DIR}/evals/${EVAL}/plots/count/${REGION}-insertions-${PARAM_SET}.${PLOT_FILETYPE}" \
        "${PLOT_PARAMS[@]}"
        
    ./scripts/barchart.py "${INPUT_DIR}/evals/${EVAL}/plots/count/stats/${REGION}-deletions-${PARAM_SET}.tsv" \
        --title "Deletion events relative to sample ${REGION^^}" \
        --x_label "Graph type" \
        --y_label "Deletions per region base" \
        --x_sideways \
        --save "${INPUT_DIR}/evals/${EVAL}/plots/count/${REGION}-deletions-${PARAM_SET}.${PLOT_FILETYPE}" \
        "${PLOT_PARAMS[@]}"
        
    ./scripts/barchart.py "${INPUT_DIR}/evals/${EVAL}/plots/count/stats/${REGION}-substitutions-${PARAM_SET}.tsv" \
        --title "Substitution events relative to sample ${REGION^^}" \
        --x_label "Graph type" \
        --y_label "Substitutions per region base" \
        --x_sideways \
        --save "${INPUT_DIR}/evals/${EVAL}/plots/count/${REGION}-substitutions-${PARAM_SET}.${PLOT_FILETYPE}" \
        "${PLOT_PARAMS[@]}"
        
    ./scripts/barchart.py "${INPUT_DIR}/evals/${EVAL}/plots/count/stats/${REGION}-unvisited-${PARAM_SET}.tsv" \
        --title "Unvisited nodes in sample ${REGION^^}" \
        --x_label "Graph type" \
        --y_label "Unvisited nodes per region base" \
        --x_sideways \
        --save "${INPUT_DIR}/evals/${EVAL}/plots/count/${REGION}-unvisited-${PARAM_SET}.${PLOT_FILETYPE}" \
        "${PLOT_PARAMS[@]}"
        
    ./scripts/barchart.py "${INPUT_DIR}/evals/${EVAL}/plots/count/stats/${REGION}-indels-${PARAM_SET}.tsv" \
        --title "Indel events relative to sample ${REGION^^}" \
        --x_label "Graph type" \
        --y_label "Indels per region base" \
        --x_sideways \
        --save "${INPUT_DIR}/evals/${EVAL}/plots/count/${REGION}-indels-${PARAM_SET}.${PLOT_FILETYPE}" \
        "${PLOT_PARAMS[@]}"
    
done

for GRAPH in "${!ALL_INSERT_BASES[@]}"; do

    # Normalize base counts by total size for plotting
    INSERT_PORTION=$(echo "${ALL_INSERT_BASES[$GRAPH]} / ${TOTAL_LENGTH}" | bc -l)
    DELETE_PORTION=$(echo "${ALL_DELETE_BASES[$GRAPH]} / ${TOTAL_LENGTH}" | bc -l)
    SUBSTITUTE_PORTION=$(echo "${ALL_SUBSTITUTE_BASES[$GRAPH]} / ${TOTAL_LENGTH}" | bc -l)
    UNVISITED_PORTION=$(echo "${ALL_UNVISITED_BASES[$GRAPH]} / ${TOTAL_LENGTH}" | bc -l)
    
    # Combine indels
    INDEL_PORTION=$(echo "${INSERT_PORTION} + ${DELETE_PORTION}" | bc -l)

    # Make overall stats files
    printf "${GRAPH}\t${INSERT_PORTION}\n" \
        >> "${INPUT_DIR}/evals/${EVAL}/plots/bp/stats/ALL-insertions-${PARAM_SET}.tsv"
    printf "${GRAPH}\t${DELETE_PORTION}\n" \
        >> "${INPUT_DIR}/evals/${EVAL}/plots/bp/stats/ALL-deletions-${PARAM_SET}.tsv"
    printf "${GRAPH}\t${SUBSTITUTE_PORTION}\n" \
        >> "${INPUT_DIR}/evals/${EVAL}/plots/bp/stats/ALL-substitutions-${PARAM_SET}.tsv"
    printf "${GRAPH}\t${UNVISITED_PORTION}\n" \
        >> "${INPUT_DIR}/evals/${EVAL}/plots/bp/stats/ALL-unvisited-${PARAM_SET}.tsv"
        
    printf "${GRAPH}\t${INDEL_PORTION}\n" \
        >> "${INPUT_DIR}/evals/${EVAL}/plots/bp/stats/ALL-indels-${PARAM_SET}.tsv"
        
    # Normalize base counts by total size for plotting
    INSERT_FREQUENCY=$(echo "${ALL_INSERT_COUNT[$GRAPH]} / ${TOTAL_LENGTH}" | bc -l)
    DELETE_FREQUENCY=$(echo "${ALL_DELETE_COUNT[$GRAPH]} / ${TOTAL_LENGTH}" | bc -l)
    SUBSTITUTE_FREQUENCY=$(echo "${ALL_SUBSTITUTE_COUNT[$GRAPH]} / ${TOTAL_LENGTH}" | bc -l)
    UNVISITED_FREQUENCY=$(echo "${ALL_UNVISITED_COUNT[$GRAPH]} / ${TOTAL_LENGTH}" | bc -l)
    
    # Combine indels
    INDEL_FREQUENCY=$(echo "${INSERT_FREQUENCY} + ${DELETE_FREQUENCY}" | bc -l)
    
    printf "${GRAPH}\t${INSERT_FREQUENCY}\n" \
        >> "${INPUT_DIR}/evals/${EVAL}/plots/count/stats/ALL-insertions-${PARAM_SET}.tsv"
    printf "${GRAPH}\t${DELETE_FREQUENCY}\n" \
        >> "${INPUT_DIR}/evals/${EVAL}/plots/count/stats/ALL-deletions-${PARAM_SET}.tsv"
    printf "${GRAPH}\t${SUBSTITUTE_FREQUENCY}\n" \
        >> "${INPUT_DIR}/evals/${EVAL}/plots/count/stats/ALL-substitutions-${PARAM_SET}.tsv"
    printf "${GRAPH}\t${UNVISITED_FREQUENCY}\n" \
        >> "${INPUT_DIR}/evals/${EVAL}/plots/count/stats/ALL-unvisited-${PARAM_SET}.tsv"
        
    printf "${GRAPH}\t${INDEL_FREQUENCY}\n" \
        >> "${INPUT_DIR}/evals/${EVAL}/plots/count/stats/ALL-indels-${PARAM_SET}.tsv"
        
    # Calculate precision and recall
    PRECISION=$(echo "${ALL_TRUE_POSITIVES[$GRAPH]} / ( ${ALL_TRUE_POSITIVES[$GRAPH]} + ${ALL_FALSE_POSITIVES[$GRAPH]} )" | bc -l)
    RECALL=$(echo "${ALL_TRUE_POSITIVES[$GRAPH]} / ( ${ALL_TRUE_POSITIVES[$GRAPH]} + ${ALL_FALSE_NEGATIVES[$GRAPH]} )" | bc -l)
    
    # Prepare the plot file for a precision vs. recall plot (where precision is Y)
    printf "${GRAPH}\t${RECALL}\t${PRECISION}\n" >> \
        "${INPUT_DIR}/evals/${EVAL}/plots/bp/stats/PR-ALL-${PARAM_SET}.tsv"
        
    # Calculate an inverse F1 score
    F1_SCORE_INV=$(echo "1 - 2 * ( ${PRECISION} * ${RECALL} ) / ( ${PRECISION} + ${RECALL} )" | bc -l)
    
    # Prepare the plot file for an F1 bar chart
    printf "${GRAPH}\t${F1_SCORE_INV}\n" >> \
        "${INPUT_DIR}/evals/${EVAL}/plots/bp/stats/F1-ALL-${PARAM_SET}.tsv"

done

# Plot totals by bp
./scripts/barchart.py "${INPUT_DIR}/evals/${EVAL}/plots/bp/stats/ALL-insertions-${PARAM_SET}.tsv" \
    --title "Inserted bases relative to sample overall" \
    --x_label "Graph type" \
    --y_label "Inserted bases per region base" \
    --x_sideways \
    --save "${INPUT_DIR}/evals/${EVAL}/plots/bp/ALL-insertions-${PARAM_SET}.${PLOT_FILETYPE}" \
    "${PLOT_PARAMS[@]}"
    
./scripts/barchart.py "${INPUT_DIR}/evals/${EVAL}/plots/bp/stats/ALL-deletions-${PARAM_SET}.tsv" \
    --title "Deleted bases relative to sample overall" \
    --x_label "Graph type" \
    --y_label "Deleted bases per region base" \
    --x_sideways \
    --save "${INPUT_DIR}/evals/${EVAL}/plots/bp/ALL-deletions-${PARAM_SET}.${PLOT_FILETYPE}" \
    "${PLOT_PARAMS[@]}"
    
./scripts/barchart.py "${INPUT_DIR}/evals/${EVAL}/plots/bp/stats/ALL-substitutions-${PARAM_SET}.tsv" \
    --title "Substituted bases relative to sample overall" \
    --x_label "Graph type" \
    --y_label "Substituted bases per region base" \
    --x_sideways \
    --save "${INPUT_DIR}/evals/${EVAL}/plots/bp/ALL-substitutions-${PARAM_SET}.${PLOT_FILETYPE}" \
    "${PLOT_PARAMS[@]}"
    
./scripts/barchart.py "${INPUT_DIR}/evals/${EVAL}/plots/bp/stats/ALL-unvisited-${PARAM_SET}.tsv" \
    --title "Unvisited node length in sample overall" \
    --x_label "Graph type" \
    --y_label "Unvisited bases per region base" \
    --x_sideways \
    --save "${INPUT_DIR}/evals/${EVAL}/plots/bp/ALL-unvisited-${PARAM_SET}.${PLOT_FILETYPE}" \
    "${PLOT_PARAMS[@]}"
    
./scripts/barchart.py "${INPUT_DIR}/evals/${EVAL}/plots/bp/stats/ALL-indels-${PARAM_SET}.tsv" \
    --title "Indle bases relative to sample overall" \
    --x_label "Graph type" \
    --y_label "Indel bases per region base" \
    --x_sideways \
    --save "${INPUT_DIR}/evals/${EVAL}/plots/bp/ALL-indels-${PARAM_SET}.${PLOT_FILETYPE}" \
    "${PLOT_PARAMS[@]}"
    
# Plot total precision/recall
./scripts/scatter.py "${INPUT_DIR}/evals/${EVAL}/plots/bp/stats/PR-ALL-${PARAM_SET}.tsv" \
    --title "Accuracy against PacBio assembly" \
    --x_label "Recall" \
    --y_label "Precision" \
    --width 6.5 --height 6 \
    --max_x 1 --min_x 0.992 --max_y 1 --min_y 0.995 \
    --no_legend --annotate \
    --save "${INPUT_DIR}/evals/${EVAL}/plots/bp/PR-ALL-${PARAM_SET}.${PLOT_FILETYPE}" \
    "${PLOT_PARAMS[@]}"

# Plot total F1
./scripts/barchart.py "${INPUT_DIR}/evals/${EVAL}/plots/bp/stats/F1-ALL-${PARAM_SET}.tsv" \
    --title "F1 score overall" \
    --x_label "Graph" \
    --y_label "1 - F1" \
    --x_sideways \
    --save "${INPUT_DIR}/evals/${EVAL}/plots/bp/F1-ALL-${PARAM_SET}.${PLOT_FILETYPE}" \
    "${PLOT_PARAMS[@]}"
    
# Plot totals by event count
./scripts/barchart.py "${INPUT_DIR}/evals/${EVAL}/plots/count/stats/ALL-insertions-${PARAM_SET}.tsv" \
    --title "Insertion events relative to sample overall" \
    --x_label "Graph type" \
    --y_label "Insertions per region base" \
    --x_sideways \
    --save "${INPUT_DIR}/evals/${EVAL}/plots/count/ALL-insertions-${PARAM_SET}.${PLOT_FILETYPE}" \
    "${PLOT_PARAMS[@]}"
    
./scripts/barchart.py "${INPUT_DIR}/evals/${EVAL}/plots/count/stats/ALL-deletions-${PARAM_SET}.tsv" \
    --title "Deletion events relative to sample overall" \
    --x_label "Graph type" \
    --y_label "Deletions per region base" \
    --x_sideways \
    --save "${INPUT_DIR}/evals/${EVAL}/plots/count/ALL-deletions-${PARAM_SET}.${PLOT_FILETYPE}" \
    "${PLOT_PARAMS[@]}"
    
./scripts/barchart.py "${INPUT_DIR}/evals/${EVAL}/plots/count/stats/ALL-substitutions-${PARAM_SET}.tsv" \
    --title "Substitution events relative to sample overall" \
    --x_label "Graph type" \
    --y_label "Substitutions per region base" \
    --x_sideways \
    --save "${INPUT_DIR}/evals/${EVAL}/plots/count/ALL-substitutions-${PARAM_SET}.${PLOT_FILETYPE}" \
    "${PLOT_PARAMS[@]}"
    
./scripts/barchart.py "${INPUT_DIR}/evals/${EVAL}/plots/count/stats/ALL-unvisited-${PARAM_SET}.tsv" \
    --title "Unvisited nodes in sample overall" \
    --x_label "Graph type" \
    --y_label "Unvisited nodes per region base" \
    --x_sideways \
    --save "${INPUT_DIR}/evals/${EVAL}/plots/count/ALL-unvisited-${PARAM_SET}.${PLOT_FILETYPE}" \
    "${PLOT_PARAMS[@]}"
    
./scripts/barchart.py "${INPUT_DIR}/evals/${EVAL}/plots/count/stats/ALL-indels-${PARAM_SET}.tsv" \
    --title "Indel events relative to sample overall" \
    --x_label "Graph type" \
    --y_label "Indels per region base" \
    --x_sideways \
    --save "${INPUT_DIR}/evals/${EVAL}/plots/count/ALL-indels-${PARAM_SET}.${PLOT_FILETYPE}" \
    "${PLOT_PARAMS[@]}"
    
        
