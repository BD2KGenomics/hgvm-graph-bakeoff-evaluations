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

rm -f "${INPUT_DIR}"/extracted/*.vg || true
mkdir -p "${INPUT_DIR}/extracted"
# Make sure the glob actually activates.
./scripts/extractGraphs.py "${INPUT_DIR}"/indexes/*/*/*.tar.gz "${INPUT_DIR}/extracted"

# What evaluation are we?
EVAL="assembly_sd"

# What sample is our pre-merged balanced synthetic diploid?
SAMPLE="SYNDIP"

# What assemblies are its two halves
ASSEMBLIES=(
    "CHM1"
    "CHM13"
)

for PARAM_SET in defray; do

    for REGION in lrc_kir brca1 brca2 mhc; do
        REF_FASTA="data/altRegions/${REGION^^}/ref.fa"
        
        # How many bases should we knock off the start of the region? Some of
        # our regions don't have good assembly coverage all the way to the start
        # in both assemblies?
        TRIM_START=0
        if [[ "${REGION}" == "lrc_kir" ]]; then
            # LRC_KIR has the second assembly come in about here, according to
            # my BED files I made from my alignment.
            TRIM_START=87796
        fi
            
        # This holds all the region-level data that is constant across graphs
        RTEMP="${INPUT_DIR}/evals/${EVAL}/temp/${PARAM_SET}/${REGION}"
        mkdir -p "${RTEMP}"
            
        # We need to grab the relevant parts of the assembly from the FASTA
        true > "${RTEMP}/reads.txt"
        for ASSEMBLY in "${ASSEMBLIES[@]}"; do
            # Get the relevant assembly regions from each sample.
            
            # Look up assembly by sample
            if [[ "${ASSEMBLY}" == "CHM1" ]]; then
                ASSEMBLY_FASTA="mole_assembly/GCA_001297185.1_PacBioCHM1_r2_GenBank_08312015_genomic.fna"
            elif [[ "${ASSEMBLY}" == "CHM13" ]]; then
                ASSEMBLY_FASTA="mole_assembly/GCA_000983455.2_CHM13_Draft_Assembly_genomic.fna"
            else
                echo "Unknown assembly ${ASSEMBLY}"
                exit 1
            fi
            
            bedtools getfasta -fi "${ASSEMBLY_FASTA}" -bed "mole_regions/${ASSEMBLY}/${REGION^^}.bed" -fo "${RTEMP}/relevant.fa"
            cat "${RTEMP}/relevant.fa" | ./scripts/fasta2reads.py --uppercase >> "${RTEMP}/reads.txt"
        done

        # Graphs we like:
        # empty snp1kg refonly shifted1kg cactus prg freebayes platypus samtools
        for GRAPH in empty snp1kg refonly shifted1kg cactus prg freebayes platypus samtools; do

            # This holds all the temporary files for this graph for this region.
            TEMP="${RTEMP}/graph/${GRAPH}"
            mkdir -p "${TEMP}"

            # Trim the FASTA. I could write a FASTA trimmer script, or
            # depend on a FASTA trimmer script, but using inline Python is
            # so much easier.
            cat "${REF_FASTA}" | python -c 'import sys; import Bio.SeqIO; Bio.SeqIO.write(Bio.SeqIO.read(sys.stdin, "fasta")[int(sys.argv[1]):], sys.stdout, "fasta")' ${TRIM_START} > "${TEMP}/ref.fa"

            if [[ "${GRAPH}" == "empty" ]]; then
                # Don't use a vcf, just vg construct
                # We can't actually tabix index an empty VCF I think
                
                vg construct -r "${TEMP}/ref.fa" -a -f > "${TEMP}/sample.vg"
            else
                # We will use a vcf. Each if branch has to output an on_ref_sorted.vcf.

                if [[ "${GRAPH}" == "freebayes" ]]; then
                    # Grab the Freebayes calls for this region
                    
                    cp "/cluster/home/hickey/ga4gh/hgvm-graph-bakeoff-evalutations/syndip_classic/freebayes/${SAMPLE}/${REGION^^}.vcf.ref" "${TEMP}/on_ref_sorted.vcf"
                
                elif [[ "${GRAPH}" == "platypus" ]]; then
                    # Grab the Platypus calls for this region
                    
                    cp "/cluster/home/hickey/ga4gh/hgvm-graph-bakeoff-evalutations/syndip_classic/platypus/${SAMPLE}/${REGION^^}.vcf.ref" "${TEMP}/on_ref_sorted.vcf"
                
                elif [[ "${GRAPH}" == "samtools" ]]; then
                    # Grab the Samtools calls for this region
                    
                    cp "/cluster/home/hickey/ga4gh/hgvm-graph-bakeoff-evalutations/syndip_classic/samtools/${SAMPLE}/${REGION^^}.vcf.ref" "${TEMP}/on_ref_sorted.vcf"
               
                else
                    # Do the genotyping since this is a real graph

                    VGFILE="${INPUT_DIR}/extracted/${GRAPH}-${REGION}.vg"
                    XGFILE="${INPUT_DIR}/extracted/${GRAPH}-${REGION}.xg"
                    
                    # Find the sample's GAM
                    GAM="${INPUT_DIR}/alignments/${REGION}/${GRAPH}/${SAMPLE}.gam"
                    
                    # Now we need to pick a calling method. Each of these will produce an unsorted calls.vcf.
                    
                    if [[ "${PARAM_SET}" == "call" ]]; then
                        
                        vg filter -r 0.90 -d 0.00 -e 0.00 -afu -s 10000 -q 15 -o 0 -E 4 "${GAM}" > "${TEMP}/filtered.gam"
                        vg pileup -w 40 -m 10 -q 10 -a "${VGFILE}" "${TEMP}/filtered.gam" > "${TEMP}/pileup.vgpu"
                        
                        # Guess ref path, because if we ask for output on "ref" and input isn't on "ref" we get in trouble
                        REF_PATH="$(vg view -j ${VGFILE} | jq -r '.path[].name' | grep -v 'GI' | head -n 1)"
                        
                        vg call -b 5 -s 1 -d 1 -f 0 -D 10 -H 3 -n 1 -F 0.2 -B 250 -R 4 --contig ref -r "${REF_PATH}" "${VGFILE}" "${TEMP}/pileup.vgpu" > "${TEMP}/calls.vcf"
                    
                    elif [[ "${PARAM_SET}" == "defray" ]]; then
                        # Just like call but with an xg index and --defray-ends on the filter step.
                        
                        if [[ (! -e "${XGFILE}")  || ("${VGFILE}" -nt "${XGFILE}") ]]; then
                            # No XG file, or VG file is newer than it. Index.
                            vg index -x "${XGFILE}" "${VGFILE}"
                        fi
                        
                        vg filter -x "${XGFILE}" -r 0.90 -d 0.00 -e 0.00 -afu -s 10000 -q 15 -o 0 -E 4 --defray-ends 999 "${GAM}" > "${TEMP}/filtered.gam"
                        vg pileup -w 40 -m 10 -q 10 -a "${VGFILE}" "${TEMP}/filtered.gam" > "${TEMP}/pileup.vgpu"
                        
                        # Guess ref path, because if we ask for output on "ref" and input isn't on "ref" we get in trouble
                        REF_PATH="$(vg view -j ${VGFILE} | jq -r '.path[].name' | grep -v 'GI' | head -n 1)"
                        
                        # Use default vg call parameters that Glenn set to be his favorite
                        vg call --contig ref -r "${REF_PATH}" "${VGFILE}" "${TEMP}/pileup.vgpu" > "${TEMP}/calls.vcf"
                        
                    elif [[ "${PARAM_SET}" == "genotype" ]]; then
                    
                        # Use vg genotype
                        
                        # Make sure to drop all secondaries
                        # TODO: can't vg filter just do that?
                        vg filter -r 0.9 -d 0.00 -e 0.00 -afu -s 10000 -q 15 -o 0 "${GAM}" | vg view -aj - | jq 'select(.is_secondary | not)' | vg view -JGa - > "${TEMP}/filtered.gam"
                    
                        rm -Rf "${TEMP}/reads.index"
                        vg index -d "${TEMP}/reads.index" -N "${TEMP}/filtered.gam"
                        vg genotype "${VGFILE}" "${TEMP}/reads.index" -C -q -i -v --contig ref --min_per_strand 1 --het_prior_denom 10 > "${TEMP}/calls.vcf" 2>"${TEMP}/log.txt"
                    
                    else
                    
                        echo "Unknown parameter set ${PARAM_SET}"
                        exit 1
                    
                    fi
                    
                    # Sort and filter on quality
                    cat "${TEMP}/calls.vcf" | sort -n -k2 | uniq | ./scripts/vcfFilterQuality.py - 5 --ad > "${TEMP}/on_ref_sorted.vcf"
                    
                fi
                    
                # Now drop failing records, trim the VCF and shift all the variants down
                cat "${TEMP}/on_ref_sorted.vcf" | \
                    awk -F $'\t' 'BEGIN {OFS = FS} { if($1 !~ "#") { if($7 == "PASS" || $7 == ".") { print } } else print; }' | \
                    awk -v offset=${TRIM_START} -F $'\t' 'BEGIN {OFS = FS} {if($1 !~ "#") { $2 -= offset; if($2 >= 0) print} else {print}}' > "${TEMP}/on_ref_trimmed.vcf"
                
                rm -f "${TEMP}/on_ref.vcf.gz"
                bgzip "${TEMP}/on_ref_trimmed.vcf" -c > "${TEMP}/on_ref.vcf.gz"
                tabix -f -p vcf "${TEMP}/on_ref.vcf.gz"
                vg construct -r "${TEMP}/ref.fa" -v "${TEMP}/on_ref.vcf.gz" -a -f > "${TEMP}/reconstructed.vg"
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
