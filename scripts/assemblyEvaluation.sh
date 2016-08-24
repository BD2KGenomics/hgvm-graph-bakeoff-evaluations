#!/usr/bin/env bash
# Realigns mole assembly to sample graphs to evaluate them.

set -ex

# Grab the input directory to look in
INPUT_DIR=${1}

if [[ ! -d "${INPUT_DIR}" ]]
then
    echo "Specify input directory!"
    exit 1
fi

# Get the absloute path so we can use the extract script, which apparently gets
# upset if used with arbitrary input paths
ABSOLUTE_INPUT_DIR="$(cd ${INPUT_DIR} && pwd)"
LOCAL_INPUT_DIR="$(basename ${INPUT_DIR})"

rm "${INPUT_DIR}"/extracted/*.vg || true
mkdir -p "${INPUT_DIR}/extracted"
# Make sure the glob actually activates.
./scripts/extractGraphs.py "${INPUT_DIR}"/indexes/*/*/*.tar.gz "${INPUT_DIR}/extracted"

# What evaluation are we?
EVAL="assembly_sd"

# What samples do we use for our synthetic diploid
SAMPLES=(
    "CHM1"
    "CHM13"
)

for PARAM_SET in call defray; do

    for REGION in brca1 brca2 sma lrc_kir mhc; do
        REF_FASTA="data/altRegions/${REGION^^}/ref.fa"
        
        
            
        # This holds all the region-level data that is constant across graphs
        RTEMP="${INPUT_DIR}/evals/${EVAL}/temp/${PARAM_SET}/${REGION}"
        mkdir -p "${RTEMP}"
            
        # We need to grab the relevant parts of the assembly from the FASTA
        true > "${RTEMP}/reads.txt"
        for SAMPLE in "${SAMPLES[@]}"; do
            # Get the relevant assembly regions from each sample.
            
            # Look up assembly by sample
            if [[ "${SAMPLE}" == "CHM1" ]]; then
                ASSEMBLY_FASTA="mole_assembly/GCA_001297185.1_PacBioCHM1_r2_GenBank_08312015_genomic.fna"
            elif [[ "${SAMPLE}" == "CHM13" ]]; then
                ASSEMBLY_FASTA="mole_assembly/GCA_000983455.2_CHM13_Draft_Assembly_genomic.fna"
            else
                echo "Unknown sample assembly ${SAMPLE}"
                exit 1
            fi
            
            bedtools getfasta -fi "${ASSEMBLY_FASTA}" -bed "mole_regions/${SAMPLE}/${REGION^^}.bed" -fo "${RTEMP}/relevant.fa"
            cat "${RTEMP}/relevant.fa" | ./scripts/fasta2reads.py --uppercase >> "${RTEMP}/reads.txt"
        done

        for GRAPH in empty snp1kg refonly shifted1kg; do

            # This holds all the temporary files for this graph for this region.
            TEMP="${RTEMP}/graph/${GRAPH}"
            mkdir -p "${TEMP}"

            if [[ "${GRAPH}" == "freebayes" ]]; then
                # Grab the Freebayes calls for this region
                # TODO: replace these with pooled CHM1/CHM13 freebayes calls
                
                cp "/cluster/home/hickey/ga4gh/hgvm-graph-bakeoff-evalutations/platinum_classic3/freebayes/${SAMPLE}/${REGION^^}.vcf.ref" "${TEMP}/on_ref_sorted.vcf"
                rm -f "${TEMP}/on_ref.vcf.gz"
                bgzip "${TEMP}/on_ref_sorted.vcf" -c > "${TEMP}/on_ref.vcf.gz"
                tabix -f -p vcf "${TEMP}/on_ref.vcf.gz"
                vg construct -r "${REF_FASTA}" -v "${TEMP}/on_ref.vcf.gz" -a -f > "${TEMP}/reconstructed.vg"
                vg validate "${TEMP}/reconstructed.vg"
                vg mod -v "${TEMP}/on_ref.vcf.gz" "${TEMP}/reconstructed.vg" > "${TEMP}/sample.vg"
            
            elif [[ "${GRAPH}" == "empty" ]]; then
                # Don't use a vcf, just vg construct
                # We can't actually tabix index an empty VCF I think
                
                vg construct -r "${REF_FASTA}" -a -f > "${TEMP}/sample.vg"
                
            else
                # Do the genotyping since this is a real graph

                VGFILE="${INPUT_DIR}/extracted/${GRAPH}-${REGION}.vg"
                XGFILE="${INPUT_DIR}/extracted/${GRAPH}-${REGION}.xg"
                
                # Make a combined GAM
                GAM="${TEMP}/combined.gam"
                true > "${GAM}"
                
                for SAMPLE in "${SAMPLES[@]}"; do
                    SAMPLE_GAM="${INPUT_DIR}/alignments/${REGION}/${GRAPH}/${SAMPLE}.gam"
                    # Concatenate all the sample GAMs together.
                    # TODO: balance read counts or something
                    cat "${SAMPLE_GAM}" >> "${GAM}"
                done
                
                if [[ "${PARAM_SET}" == "call" ]]; then
                    
                    vg filter -r 0.90 -d 0.00 -e 0.00 -afu -s 10000 -q 15 -o 0 -E 4 "${GAM}" > "${TEMP}/filtered.gam"
                    vg pileup -w 40 -m 10 -q 10 -a "${VGFILE}" "${TEMP}/filtered.gam" > "${TEMP}/pileup.vgpu"
                    
                    # Guess ref path, because if we ask for output on "ref" and input isn't on "ref" we get in trouble
                    REF_PATH="$(vg view -j ${VGFILE} | jq -r '.path[].name' | grep -v 'GI' | head -n 1)"
                    
                    vg call -b 5 -s 1 -d 1 -f 0 -D 10 -H 3 -n 1 -F 0.2 -B 250 -R 4 --contig ref -r "${REF_PATH}" "${VGFILE}" "${TEMP}/pileup.vgpu" > "${TEMP}/calls.vcf"
                    cat "${TEMP}/calls.vcf" | sort -n -k2 | uniq | ./scripts/vcfFilterQuality.py - 5 --ad > "${TEMP}/on_ref_sorted.vcf"
                
                elif [[ "${PARAM_SET}" == "defray" ]]; then
                    # Just like call but with an xg index and --defray-ends on the filter step.
                    
                    if [[ (! -e "${XGFILE}")  || ("${VGFILE}" -nt "${XGFILE}") ]]; then
                        # No XG file, or VG file is newer than it. Index.
                        vg index -x "${XGFILE}" "${VGFILE}"
                    fi
                    
                    vg filter -x "${XGFILE}" -r 0.90 -d 0.00 -e 0.00 -afu -s 10000 -q 15 -o 0 -E 4 --defray-ends 40 "${GAM}" > "${TEMP}/filtered.gam"
                    vg pileup -w 40 -m 10 -q 10 -a "${VGFILE}" "${TEMP}/filtered.gam" > "${TEMP}/pileup.vgpu"
                    
                    # Guess ref path, because if we ask for output on "ref" and input isn't on "ref" we get in trouble
                    REF_PATH="$(vg view -j ${VGFILE} | jq -r '.path[].name' | grep -v 'GI' | head -n 1)"
                    
                    vg call -b 5 -s 1 -d 1 -f 0 -D 10 -H 3 -n 1 -F 0.2 -B 250 -R 4 --contig ref -r "${REF_PATH}" "${VGFILE}" "${TEMP}/pileup.vgpu" > "${TEMP}/calls.vcf"
                    cat "${TEMP}/calls.vcf" | sort -n -k2 | uniq | ./scripts/vcfFilterQuality.py - 5 --ad > "${TEMP}/on_ref_sorted.vcf"
                    
                elif [[ "${PARAM_SET}" == "genotype" ]]; then
                
                    # Use vg genotype
                    
                    # Make sure to drop all secondaries
                    # TODO: can't vg filter just do that?
                    vg filter -r 0.9 -d 0.00 -e 0.00 -afu -s 10000 -q 15 -o 0 "${GAM}" | vg view -aj - | jq 'select(.is_secondary | not)' | vg view -JGa - > "${TEMP}/filtered.gam"
                
                    rm -Rf "${TEMP}/reads.index"
                    vg index -d "${TEMP}/reads.index" -N "${TEMP}/filtered.gam"
                    vg genotype "${VGFILE}" "${TEMP}/reads.index" -C -q -i -v --contig ref --min_per_strand 1 --het_prior_denom 10 > "${TEMP}/calls.vcf" 2>"${TEMP}/log.txt"
                    cat "${TEMP}/calls.vcf" | sort -n -k2 | uniq | ./scripts/vcfFilterQuality.py - 5 --ad > "${TEMP}/on_ref_sorted.vcf"
                
                else
                
                    echo "Unknown parameter set ${PARAM_SET}"
                    exit 1
                
                fi
                
                
            
                rm -f "${TEMP}/on_ref.vcf.gz"
                bgzip "${TEMP}/on_ref_sorted.vcf" -c > "${TEMP}/on_ref.vcf.gz"
                tabix -f -p vcf "${TEMP}/on_ref.vcf.gz"
                vg construct -r "${REF_FASTA}" -v "${TEMP}/on_ref.vcf.gz" -a -f > "${TEMP}/reconstructed.vg"
                vg validate "${TEMP}/reconstructed.vg"
                vg mod -v "${TEMP}/on_ref.vcf.gz" "${TEMP}/reconstructed.vg" > "${TEMP}/sample.vg"
            
            fi
            
            
            vg validate "${TEMP}/sample.vg"
            
            # Index the sample graph and align the assembly contigs
            vg index -x "${TEMP}/sample.xg" -g "${TEMP}/sample.gcsa" -k 16 "${TEMP}/sample.vg"
            vg map -x "${TEMP}/sample.xg" -g "${TEMP}/sample.gcsa" -r "${RTEMP}/reads.txt" > "${TEMP}/assembly_aligned.gam"
            
            mkdir -p "${INPUT_DIR}/evals/${EVAL}/stats/${REGION}"
            
            vg stats "${TEMP}/sample.vg" -v -a "${TEMP}/assembly_aligned.gam" > "${INPUT_DIR}/evals/${EVAL}/stats/${REGION}/${GRAPH}-${PARAM_SET}.txt" 2>&1
            
        done
        
    done
    
done
