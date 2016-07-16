#!/bin/bash

# run some other variant callers on linear (non-graph) alignmetns as points of comparison
# use mostly default options.  

# expects there to be a rocksdb index for graphs in <graph_name>.index

if [ "$#" -ne 4 ]; then
	 echo "Syntax $0 <graphs_dir> <alginments_dir> <reads_dir> <out_dir>"
	 exit 1
fi

GRAPHS_DIR=$1
ALIGNMENT_DIR=$2
READS_DIR=$3
OUT_DIR=$4


# must be absolute
FA_DIR="/cluster/home/hickey/ga4gh/hgvm-graph-bakeoff-evalutations/data/altRegions"
PLATYPUS_CMD="python /cluster/home/hickey/tools/Platypus/bin/Platypus.py"

REGIONS=( "BRCA2" "MHC" "BRCA1" "SMA" "LRC_KIR" )
#SAMPLES=( "NA12878" "NA12877" "NA12879" )

#REGIONS=( "LRC_KIR" )
SAMPLES=( "NA12878" )
SAMPLES=( "CHM1" )

PLATYPUS_OPTS=" --mergeClusteredVariants=1"
FREEBAYES_OPTS=" --strict-vcf"
MPILEUP_OPTS=""
BCFTOOLS_CALL_OPTS=""
BWA_OPTS=" -t 20 -p"
#BWA_OPTS=" -t 20"
SURJECT_OPTS=" -p ref -b -t 30"
OVERWRITE=1

function primary_graph {
	 local REGION=$1

	 echo ${GRAPHS_DIR}/refonly-${REGION,,}.vg
}

function primary_gam {
	 local REGION=$1
	 local SAMPLE=$2

	 echo ${ALIGNMENT_DIR}/${REGION,,}/refonly/${SAMPLE}.gam
}

# get the vcf -R coordinates for a region
function get_chrom {

    local REGION=$1
    # get the bed path
    local BED=data/g1kvcf/${REGION}.bed
    # get contig
    local CONTIG=`cat ${BED} | awk '{print $1}'`
	 echo $CONTIG
}

function get_offset {

	 local REGION=$1
	 local BED=data/g1kvcf/${REGION}.bed
    #get coordinates (convert BED into 1-based, inclusive)
    local START=`cat ${BED} | awk '{print $2}'`
	 echo $START
}

# fix vcf coordinates
function change_coordinates {
	 local IN_VCF=$1
	 local REGION=$2
	 local SAMPLE=$3

	 local CHROM=`get_chrom $REGION`
	 local OFFSET=`get_offset $REGION`

	 cp $IN_VCF $IN_VCF.ref
	 echo "##fileformat=VCFv4.1" > $IN_VCF.fix
	 echo "##contig=<ID=5,length=181538259>" >> $IN_VCF.fix
	 echo "##contig=<ID=6,length=170805979>" >> $IN_VCF.fix
	 echo "##contig=<ID=13,length=114364328>" >> $IN_VCF.fix
	 echo "##contig=<ID=17,length=83257441>" >> $IN_VCF.fix
	 echo "##contig=<ID=19,length=58617616>" >> $IN_VCF.fix
	 bcftools view -h $IN_VCF | grep "##" | grep -v contig | grep -v fileformat >> $IN_VCF.fix
	 bcftools view -h $IN_VCF | grep -v "##" | sed -e s/unknown/${SAMPLE}/g >> $IN_VCF.fix

	 bcftools view -H $IN_VCF | awk -v OFS='\t' -v c=$CHROM -v o=$OFFSET '{$1=c; $2=$2+o; print $0}' >> $IN_VCF.fix
	 mv $IN_VCF.fix $IN_VCF
}

# zap quality scores
function set_mapq {
	 local IN_BAM=$1
	 local QUAL=$2
	 local OUT_BAM=$3

	 samtools view -H $IN_BAM > $OUT_BAM
	 samtools view $IN_BAM | awk -v OFS='\t' -v q=$QUAL '{$5=q; print $0}' >> $OUT_BAM
	 samtools view -b $OUT_BAM > $OUT_BAM.temp
	 mv $OUT_BAM.temp $OUT_BAM

	 samtools sort $OUT_BAM -o ${OUT_BAM}.sort
	 mv ${OUT_BAM}.sort $OUT_BAM
	 samtools index -b $OUT_BAM
}

# zap pairing
function remove_pairing {
	 local IN_BAM=$1
	 local OUT_BAM=$2

	 samtools view -H $IN_BAM > $OUT_BAM
	 samtools view $IN_BAM | awk -v OFS='\t' '{$8=1; $9=0; print $0}' >> $OUT_BAM
	 samtools view -b $OUT_BAM > $OUT_BAM.temp
	 mv $OUT_BAM.temp $OUT_BAM

	 samtools sort $OUT_BAM -o ${OUT_BAM}.sort
	 mv ${OUT_BAM}.sort $OUT_BAM
	 samtools index -b $OUT_BAM
}


# do the bwa index, never overwrite
function bwa_index {
	 local FA=$1
	 local WORK_DIR=$2

	 local IDX_FILE=${WORK_DIR}/$(basename $FA).bwt
	 
	 if [ ! -f $IDX_FILE ]
	 then
		  ln -fs $FA $WORK_DIR
		  bwa index ${WORK_DIR}/$(basename $FA)
	 fi
}

# do bwa on all the reads
function bwa_mem {
	 local FA=$1
	 local READS=$2
	 local OUT_BAM=$3

	 if [[ ! -f $OUT_BAM ]] || [[ $OVERWRITE -eq 1 ]]
	 then
		  echo "bwa mem $BWA_OPTS $FA $READS | samtools view -1 - > $OUT_BAM"
		  bwa mem $BWA_OPTS $FA $READS | samtools view -1 - > $OUT_BAM
		  samtools sort $OUT_BAM -o ${OUT_BAM}.sort
		  mv ${OUT_BAM}.sort $OUT_BAM
		  samtools index -b $OUT_BAM
	 fi
}

# get bam from surject
function surject_gam {
	 local INDEX=$1.index
	 local GAM=$2
	 local OUT_BAM=$3

	 if [[ ! -f $OUT_BAM ]] || [[ $OVERWRITE -eq 1 ]]
	 then
		  vg surject $SURJECT_OPTS -d $INDEX $GAM > $OUT_BAM
		  samtools sort $OUT_BAM -o ${OUT_BAM}.sort
		  mv ${OUT_BAM}.sort $OUT_BAM
		  samtools index -b $OUT_BAM
	 fi	 
}

# run platypus
function platypus {
	 local IN_BAM=$1
	 local IN_REF=$2
	 local REGION=$3
	 local SAMPLE=$4
	 local OUT_VCF=$5

	 if [[ ! -f $OUT_VCF ]] || [[ $OVERWRITE -eq 1 ]]
	 then
		  samtools faidx $IN_REF
		  $PLATYPUS_CMD callVariants --refFile $IN_REF --bamFiles $IN_BAM $PLATYPUS_OPTS -o $OUT_VCF
		  change_coordinates $OUT_VCF $REGION $SAMPLE
	 fi	 
}

# run freebayes
function free_bayes {
	 local IN_BAM=$1
	 local IN_REF=$2
	 local REGION=$3
	 local SAMPLE=$4
	 local OUT_VCF=$5

	 if [[ ! -f $OUT_VCF ]] || [[ $OVERWRITE -eq 1 ]]
	 then
		  samtools faidx $IN_REF
		  freebayes -f  $IN_REF $IN_BAM $FREEBAYES_OPTS > $OUT_VCF
		  change_coordinates $OUT_VCF $REGION $SAMPLE
	 fi	 
}

# run samtools
function sam_tools {
	 local IN_BAM=$1
	 local IN_REF=$2
	 local REGION=$3
	 local SAMPLE=$4
	 local OUT_VCF=$5

	 if [[ ! -f $OUT_VCF ]] || [[ $OVERWRITE -eq 1 ]]
	 then
		  samtools faidx $IN_REF
		  samtools mpileup -uf $IN_REF $IN_BAM $MPILEUP_OPTS | bcftools call -mv -Ov $BCFTOOLS_CALL_OPTS > $OUT_VCF
		  change_coordinates $OUT_VCF $REGION $SAMPLE
	 fi	 
}


# do bwa on all the read sets
function align_all {

	 mkdir $OUT_DIR 2> /dev/null
	 mkdir ${OUT_DIR}/platypus 2> /dev/null
	 mkdir ${OUT_DIR}/freebayes 2> /dev/null
	 mkdir ${OUT_DIR}/samtools 2> /dev/null
	 mkdir ${OUT_DIR}/platypus_surject 2> /dev/null
	 mkdir ${OUT_DIR}/freebayes_surject 2> /dev/null
	 mkdir ${OUT_DIR}/samtools_surject 2> /dev/null
	 mkdir ${OUT_DIR}/platypus_mapq 2> /dev/null
	 mkdir ${OUT_DIR}/freebayes_mapq 2> /dev/null
	 mkdir ${OUT_DIR}/samtools_mapq 2> /dev/null
	 mkdir ${OUT_DIR}/platypus_pair 2> /dev/null
	 mkdir ${OUT_DIR}/freebayes_pair 2> /dev/null
	 mkdir ${OUT_DIR}/samtools_pair 2> /dev/null
	 mkdir ${OUT_DIR}/platypus_pairo 2> /dev/null
	 mkdir ${OUT_DIR}/freebayes_pairo 2> /dev/null
	 mkdir ${OUT_DIR}/samtools_pairo 2> /dev/null
	 
	 for REGION in "${REGIONS[@]}"
	 do
		  local BASE_DIR=${OUT_DIR}/${REGION}
		  mkdir ${BASE_DIR} 2> /dev/null

		  local FASTA=${FA_DIR}/${REGION}/ref.fa
		  bwa_index $FASTA $BASE_DIR

		  for SAMPLE in "${SAMPLES[@]}"
		  do
				local FASTQ=${READS_DIR}/${REGION}/${SAMPLE}/${SAMPLE}.bam.fq
				local INPUT_REF=${BASE_DIR}/$(basename $FASTA)
				local OUTPUT_BAM=${BASE_DIR}/${SAMPLE}.bam

				# bwa-mem pipelines
				bwa_mem $INPUT_REF $FASTQ $OUTPUT_BAM

                echo "Platypus for ${SAMPLE} ${REGION}"
				mkdir ${OUT_DIR}/platypus/${SAMPLE} 2> /dev/null
				local PLAT_OUTPUT=${OUT_DIR}/platypus/${SAMPLE}/${REGION}.vcf
				platypus $OUTPUT_BAM $INPUT_REF $REGION $SAMPLE $PLAT_OUTPUT

                echo "Freebayes for ${SAMPLE} ${REGION}"
				mkdir ${OUT_DIR}/freebayes/${SAMPLE} 2> /dev/null
				local FREEBAYES_OUTPUT=${OUT_DIR}/freebayes/${SAMPLE}/${REGION}.vcf
				free_bayes $OUTPUT_BAM $INPUT_REF $REGION $SAMPLE $FREEBAYES_OUTPUT

				mkdir ${OUT_DIR}/samtools/${SAMPLE} 2> /dev/null
				local SAMTOOLS_OUTPUT=${OUT_DIR}/samtools/${SAMPLE}/${REGION}.vcf
				sam_tools $OUTPUT_BAM $INPUT_REF $REGION $SAMPLE $SAMTOOLS_OUTPUT

				# as above, but zap map qualities to 60
				local MAPQ_BAM=${BASE_DIR}/${SAMPLE}_mapq.bam
				set_mapq $OUTPUT_BAM 60 $MAPQ_BAM
				
				echo "Platypus no MAPQ for ${SAMPLE} ${REGION}"
				mkdir ${OUT_DIR}/platypus_mapq/${SAMPLE} 2> /dev/null
				local PLAT_OUTPUT_MAPQ=${OUT_DIR}/platypus_mapq/${SAMPLE}/${REGION}.vcf
				platypus $MAPQ_BAM $INPUT_REF $REGION $SAMPLE $PLAT_OUTPUT_MAPQ

				mkdir ${OUT_DIR}/freebayes_mapq/${SAMPLE} 2> /dev/null
				local FREEBAYES_OUTPUT_MAPQ=${OUT_DIR}/freebayes_mapq/${SAMPLE}/${REGION}.vcf
				free_bayes $MAPQ_BAM $INPUT_REF $REGION $SAMPLE $FREEBAYES_OUTPUT_MAPQ

				mkdir ${OUT_DIR}/samtools_mapq/${SAMPLE} 2> /dev/null
				local SAMTOOLS_OUTPUT_MAPQ=${OUT_DIR}/samtools_mapq/${SAMPLE}/${REGION}.vcf
				sam_tools $MAPQ_BAM $INPUT_REF $REGION $SAMPLE $SAMTOOLS_OUTPUT_MAPQ

				# remove pairing information
				local PAIRO_BAM=${BASE_DIR}/${SAMPLE}_pairo.bam
				remove_pairing $OUTPUT_BAM $PAIRO_BAM
				
				mkdir ${OUT_DIR}/platypus_pairo/${SAMPLE} 2> /dev/null
				local PLAT_OUTPUT_PAIRO=${OUT_DIR}/platypus_pairo/${SAMPLE}/${REGION}.vcf
				platypus $PAIRO_BAM $INPUT_REF $REGION $SAMPLE $PLAT_OUTPUT_PAIRO

				mkdir ${OUT_DIR}/freebayes_pairo/${SAMPLE} 2> /dev/null
				local FREEBAYES_OUTPUT_PAIRO=${OUT_DIR}/freebayes_pairo/${SAMPLE}/${REGION}.vcf
				free_bayes $PAIRO_BAM $INPUT_REF $REGION $SAMPLE $FREEBAYES_OUTPUT_PAIRO

				mkdir ${OUT_DIR}/samtools_pairo/${SAMPLE} 2> /dev/null
				local SAMTOOLS_OUTPUT_PAIRO=${OUT_DIR}/samtools_pairo/${SAMPLE}/${REGION}.vcf
				sam_tools $PAIRO_BAM $INPUT_REF $REGION $SAMPLE $SAMTOOLS_OUTPUT_PAIRO

				# zap map qualities to 60 and remove pairing information
				local PAIR_BAM=${BASE_DIR}/${SAMPLE}_pair.bam
				remove_pairing $MAPQ_BAM $PAIR_BAM
				
				mkdir ${OUT_DIR}/platypus_pair/${SAMPLE} 2> /dev/null
				local PLAT_OUTPUT_PAIR=${OUT_DIR}/platypus_pair/${SAMPLE}/${REGION}.vcf
				platypus $PAIR_BAM $INPUT_REF $REGION $SAMPLE $PLAT_OUTPUT_PAIR

				mkdir ${OUT_DIR}/freebayes_pair/${SAMPLE} 2> /dev/null
				local FREEBAYES_OUTPUT_PAIR=${OUT_DIR}/freebayes_pair/${SAMPLE}/${REGION}.vcf
				free_bayes $PAIR_BAM $INPUT_REF $REGION $SAMPLE $FREEBAYES_OUTPUT_PAIR

				mkdir ${OUT_DIR}/samtools_pair/${SAMPLE} 2> /dev/null
				local SAMTOOLS_OUTPUT_PAIR=${OUT_DIR}/samtools_pair/${SAMPLE}/${REGION}.vcf
				sam_tools $PAIR_BAM $INPUT_REF $REGION $SAMPLE $SAMTOOLS_OUTPUT_PAIR

				# surject bams and run on those
				local SURJECT_BAM=${BASE_DIR}/${SAMPLE}_surject.bam
				local GRAPH_PATH=`primary_graph $REGION`
				local GAM_PATH=`primary_gam $REGION $SAMPLE`
				surject_gam $GRAPH_PATH $GAM_PATH $SURJECT_BAM.surject
				set_mapq $SURJECT_BAM.surject 60 $SURJECT_BAM

				mkdir ${OUT_DIR}/platypus_surject/${SAMPLE} 2> /dev/null
				local PLAT_OUTPUT_SURJECT=${OUT_DIR}/platypus_surject/${SAMPLE}/${REGION}.vcf
				platypus $SURJECT_BAM $INPUT_REF $REGION $SAMPLE $PLAT_OUTPUT_SURJECT

				mkdir ${OUT_DIR}/freebayes_surject/${SAMPLE} 2> /dev/null
				local FREEBAYES_OUTPUT_SURJECT=${OUT_DIR}/freebayes_surject/${SAMPLE}/${REGION}.vcf
				free_bayes $SURJECT_BAM $INPUT_REF $REGION $SAMPLE $FREEBAYES_OUTPUT_SURJECT

				mkdir ${OUT_DIR}/samtools_surject/${SAMPLE} 2> /dev/null
				local SAMTOOLS_OUTPUT_SURJECT=${OUT_DIR}/samtools_surject/${SAMPLE}/${REGION}.vcf
				sam_tools $SURJECT_BAM $INPUT_REF $REGION $SAMPLE $SAMTOOLS_OUTPUT_SURJECT

		  done
	 done				
}

align_all

#temp
exit 0

