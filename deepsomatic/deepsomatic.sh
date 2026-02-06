#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../mamba/bin:${PATH}

HG=${BASEDIR}/../genome/hg38.fa

BIN_VERSION="1.9.0"
BASE=/tmp

# Set up input and output directory data
INPUT_DIR="${BASE}/input/data"
OUTPUT_DIR="${BASE}/output"

## DeepSomatic
docker pull google/deepsomatic:"${BIN_VERSION}"
for BAM in ${BASEDIR}/../alignment/*.bam
do
    ID=`echo ${BAM} | sed 's/^.*\///' | sed 's/.bam$//'`
    echo ${ID}
    if [ ! -f ${ID}.bcf ]
    then
	## Create local directory structure
	mkdir -p "${INPUT_DIR}"
	mkdir -p "${OUTPUT_DIR}"
	if [ ! -f ${INPUT_DIR}/hg38.fa ]
	then
	    cp ${HG} ${INPUT_DIR}/
	    cp ${HG}.fai ${INPUT_DIR}/
	fi
	if [ ! -f ${INPUT_DIR}/${ID}.bam ]
	then
	    cp ${BAM} ${INPUT_DIR}/
	    cp ${BAM}.bai ${INPUT_DIR}/
	fi
	
	echo "Run docker"
	for CHR in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY
	do
	    docker run -v ${INPUT_DIR}:${INPUT_DIR} -v ${OUTPUT_DIR}:${OUTPUT_DIR} google/deepsomatic:"${BIN_VERSION}" run_deepsomatic --model_type=ONT_TUMOR_ONLY --ref=${INPUT_DIR}/hg38.fa --reads_tumor=${INPUT_DIR}/${ID}.bam --output_vcf=${OUTPUT_DIR}/${ID}.${CHR}.vcf.gz --sample_name_tumor="${ID}" --num_shards=8 --logging_dir=${OUTPUT_DIR}/logs --intermediate_results_dir=${OUTPUT_DIR}/intermediate_results_dir --use_default_pon_filtering=true --regions=${CHR}
	done
	bcftools concat -a -O b -o ${ID}.bcf ${OUTPUT_DIR}/${ID}.chr*.vcf.gz
	bcftools index ${ID}.bcf
	rm -rf ${INPUT_DIR}/${ID}.bam ${INPUT_DIR}/${ID}.bam.bai
    fi
done
