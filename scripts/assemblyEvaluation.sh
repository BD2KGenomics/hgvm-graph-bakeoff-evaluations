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

rm -R "${INPUT_DIR}/extracted"
./scripts/extractGraphs.py "${INPUT_DIR}/indexes/*/*/*.tar.gz" "${INPUT_DIR}/extracted"


for REGION in brca1 brca2; do
    REF_FASTA="data/altRegions/${REGION^^}/ref.fa"
    ASSEMBLY_FASTA="mole_assembly/GCA_001297185.1_PacBioCHM1_r2_GenBank_08312015_genomic.fna"
        
    RTEMP=`mktemp -d`
        
    # We need to grab the relevant parts of the assembly from the FASTA
    bedtools getfasta -fi "${ASSEMBLY_FASTA}" -bed "mole_regions/CHM1/${REGION^^}.bed" -fo "${RTEMP}/relevant.fa"
    cat "${RTEMP}/relevant.fa" | ./scripts/fasta2reads.py --uppercase > "${RTEMP}/reads.txt"

    for GRAPH in snp1kg refonly shifted1kg; do

        TEMP=`mktemp -d`

        VGFILE="mole_graphs_extracted/${GRAPH}-${REGION}.vg"
        GAM="mole_alignments_updated/alignments/${REGION}/${GRAPH}/CHM1.gam"
        
        
        vg filter -r 0.9 -d 0.00 -e 0.00 -afu -s 10000 -q 15 -o 0 "${GAM}" > "${TEMP}/filtered.gam"
        vg pileup -w 40 -m 10 -q 10 -a "${VGFILE}" "${TEMP}/filtered.gam" > "${TEMP}/pileup.vgpu"
        vg call -r 0.0001 -b 5 -s 1 -d 1 -f 0 -l "${VGFILE}" "${TEMP}/pileup.vgpu" --calls "${TEMP}/calls.tsv" -l > "${TEMP}/augmented.vg"
        glenn2vcf --depth 10 --max_het_bias 3 --min_count 1 --min_fraction 0.2 --contig ref "${TEMP}/augmented.vg" "${TEMP}/calls.tsv" > "${TEMP}/calls.vcf"
        
        cat "${TEMP}/calls.vcf" | sort -n -k2 | uniq > "${TEMP}/on_ref_sorted.vcf"
        rm -f "${TEMP}/on_ref.vcf.gz"
        bgzip "${TEMP}/on_ref_sorted.vcf" -c > "${TEMP}/on_ref.vcf.gz"
        tabix -f -p vcf "${TEMP}/on_ref.vcf.gz"
        vg construct -r "${REF_FASTA}" -v "${TEMP}/on_ref.vcf.gz" -a -f > "${TEMP}/reconstructed.vg"
        vg validate "${TEMP}/reconstructed.vg"
        vg mod -v "${TEMP}/on_ref.vcf.gz" "${TEMP}/reconstructed.vg" > "${TEMP}/sample.vg"
        vg validate "${TEMP}/sample.vg"
        
        vg index -x "${TEMP}/sample.xg" -g "${TEMP}/sample.gcsa" -k 16 "${TEMP}/sample.vg"
        vg map -x "${TEMP}/sample.xg" -g "${TEMP}/sample.gcsa" -r "${RTEMP}/reads.txt" > "${TEMP}/assembly_aligned.gam"
        
        mkdir -p "mole_alignments_updated/evals/assembly/stats/${REGION}"
        
        vg stats "${TEMP}/sample.vg" -v -a "${TEMP}/assembly_aligned.gam" > "mole_alignments_updated/evals/assembly/stats/${REGION}/${GRAPH}.txt" 2>&1
        
        rm -R "${TEMP}"
    
    done
    
    rm -R "${RTEMP}"

done
