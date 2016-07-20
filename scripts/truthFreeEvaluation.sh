#!/usr/bin/env bash
# Evaluates variant calls without recourse to a truth set
# Run after parallelMappingEvaluation.sh
# Usage: ./scripts/truthFreeEvaluation.sh <alignment input directory> <read source directory>

set -ex

# What plot filetype should we produce?
PLOT_FILETYPE="svg"

# Grab the input directory to look in
INPUT_DIR=${1}

if [[ ! -d "${INPUT_DIR}" ]]
then
    echo "Specify input directory!"
    exit 1
fi

# Grab the read directory to look in
READ_DIR=${2}

if [[ ! -d "${READ_DIR}" ]]
then
    echo "Specify read directory!"
    exit 1
fi



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
    _truth_
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
    "'Truth'"
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
    "#fdbf6f"
    --dpi 90 --no_n
)

# Extract the graphs
EXTRACTED_DIR="${INPUT_DIR}/evals/truthfree/extracted"
if [[ ! -d "${EXTRACTED_DIR}" ]]; then
    # Extract all the graphs
    mkdir -p "${EXTRACTED_DIR}"
    ./scripts/extractGraphs.py "${INPUT_DIR}"/indexes/*/*/*.tar.gz "${EXTRACTED_DIR}"
fi

# Make a place to index reads
READ_INDEX_DIR="${INPUT_DIR}/evals/truthfree/read_indexes"

# Where do we put the plots?
PLOTS_ROOT_DIR="${INPUT_DIR}/evals/truthfree/plots"
mkdir -p "${PLOTS_ROOT_DIR}"

# And the stats files to make the plots?
STATS_ROOT_DIR="${INPUT_DIR}/evals/truthfree/stats"
mkdir -p "${STATS_ROOT_DIR}"

for REGION in brca2; do
    # Go through all the regions (TODO: hardcoded for now)
    
    # Make a file to put per-graph stats in
    mkdir -p "${STATS_ROOT_DIR}/${REGION}"
    REGION_STATS_FILE="${STATS_ROOT_DIR}/${REGION}/stats.tsv"
    # Don't clear out the stats because we only want to run each graph once.
    #true > ${REGION_STATS_FILE}
    
    
    for SAMPLE in NA12878; do
        # Go through all the samples (TODO: hardcoded for now)
    
        for GRAPH in snp1kg refonly cactus prg sbg camel _truth_; do
            # Go through all the graphs (TODO: hardcoded for now)
    
            # TODO: shifted1kg fails due to not having any alt paths in its graph somehow.
        
            # We need a directory for the .vcf to live in
            VCF_DIR="${INPUT_DIR}/evals/truthfree/stats/${REGION}/${GRAPH}/${SAMPLE}"
            
            if [[ -d "${VCF_DIR}" ]]; then
                # Assume we already ran this one and put it in the stats TSV.
                # TODO: add a rerun flag somehow?
                continue
            fi
            
            mkdir -p "${VCF_DIR}"
        
            if [[ "${GRAPH}" == "_truth_" ]]; then
                # Special graph: grab from the truth set
                
                cat "data/platinum/${SAMPLE}/${REGION^^}.vcf" | grep "#" > "${VCF_DIR}/unsorted.vcf"
                echo "##contig=<ID=ref,length=10000000>" >> "${VCF_DIR}/unsorted.vcf"
                
                # Subtract out a different base offset to convert to region-relative coordinates for every region.
                if [[ "${REGION}" == "brca2" ]]; then
                    OFFSET=32314860
                elif [[ "${REGION}" == "mhc" ]]; then
                    OFFSET=28510119
                elif [[ "${REGION}" == "lrc_kir" ]]; then
                    OFFSET=54025633
                else
                    echo "Need base offset for region ${REGION}!"
                    exit
                fi
                
                # Adjust positions and only take passing variants                
                cat "data/platinum/${SAMPLE}/${REGION^^}.vcf" | grep -v "^#" | awk -v offset="${OFFSET}" -F $'\t' 'BEGIN {OFS = FS} { $1 = "ref"; $2 -= offset; print $0 }' | grep "PASS"  >> "${VCF_DIR}/unsorted.vcf"
                
                cat "${VCF_DIR}/unsorted.vcf" | sort -n -k2 > "${VCF_DIR}/sorted.vcf"
                
                rm -f "${VCF_DIR}/sorted.vcf.gz"
                bgzip "${VCF_DIR}/sorted.vcf" -c > "${VCF_DIR}/sorted.vcf.gz"
                tabix -f -p vcf "${VCF_DIR}/sorted.vcf.gz"
                
            else
                # Use the genotypes we get from this graph
            
            
                # Where do we index these reads?
                INDEX_DIR="${READ_INDEX_DIR}/${REGION}/${GRAPH}/${SAMPLE}/reads.index"
                
                if [[ ! -d "${INDEX_DIR}" ]]; then
                    # Index the reads
                    mkdir -p "${INDEX_DIR}"
                    FILTERED_FILE=`mktemp`
                    # Get just the primary alignments
                    vg view -aj "${INPUT_DIR}/alignments/${REGION}/${GRAPH}/${SAMPLE}.gam" | jq -c 'select(.is_secondary | not)' | vg view -JGa - > "${FILTERED_FILE}"
                    # Index them into the index
                    vg index -d "${INDEX_DIR}" -N "${FILTERED_FILE}"
                    rm "${FILTERED_FILE}"
                fi
                
                # Do the genotyping
                # We re-do this every time because the genotype changes pretty fast
                vg genotype "${EXTRACTED_DIR}/${GRAPH}-${REGION}.vg" "${INDEX_DIR}" -v -q -i -C > "${VCF_DIR}/called.vcf"
                
                # Sort, compress, and index the VCF
                cat "${VCF_DIR}/called.vcf" | sort -n -k2 > "${VCF_DIR}/sorted.vcf"
                rm -f "${VCF_DIR}/sorted.vcf.gz"
                bgzip "${VCF_DIR}/sorted.vcf" -c > "${VCF_DIR}/sorted.vcf.gz"
                tabix -f -p vcf "${VCF_DIR}/sorted.vcf.gz"
                
            fi
            
            # Grab the reference FASTA
            ln -f "data/altRegions/${REGION^^}/ref.fa" "${VCF_DIR}/ref.fa"
            
            # Make a sample graph
            vg construct -r "${VCF_DIR}/ref.fa" -v "${VCF_DIR}/sorted.vcf.gz" -a -f > "${VCF_DIR}/reconstructed.vg"
            vg mod -v "${VCF_DIR}/sorted.vcf.gz" "${VCF_DIR}/reconstructed.vg" > "${VCF_DIR}/sample.vg"
            
            # Remap reads
            vg map -V "${VCF_DIR}/sample.vg" -f "${READ_DIR}/${REGION^^}/${SAMPLE}/${SAMPLE}.bam.fq" -k 16 -i > "${VCF_DIR}/realigned.gam"
            
            # Collect scores
            vg view -aj "${VCF_DIR}/realigned.gam" | jq '.score' | sed 's/null/0/' > "${VCF_DIR}/scores.txt"
            
            # Collect portion of het sites that are significantly biased
            # Sed escaped parens are the regex capture ones for some reason...
            PORTION_BIASED=$(vg stats "${VCF_DIR}/sample.vg" -a "${VCF_DIR}/realigned.gam" | tail -n1 | sed 's/.*(\(.*\)%).*/\1/g')
            
            # Find the mean score
            MEAN_SCORE=$(cat "${VCF_DIR}/scores.txt" | ./scripts/mean.sh)
            
            # Put them in as a point for this graph
            printf "${GRAPH}\t${PORTION_BIASED}\t${MEAN_SCORE}\n" >> "${REGION_STATS_FILE}"
 
        done
    done
    
    # Now make the plot for the region
    ./scripts/scatter.py "${STATS_ROOT_DIR}/${REGION}/stats.tsv" \
        --save "${PLOTS_ROOT_DIR}/${REGION}.png" \
        --x_label "Significantly Allele-Biased Hets (%)" \
        --y_label "Mean Read Alignment Score (points)" \
        --title "$(printf "Mean read score vs.\nsite bias in ${REGION}")" \
        --legend_overlay best \
        "${PLOT_PARAMS[@]}"
    
done

