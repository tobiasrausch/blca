#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")

# Nextflow parameters
export NXF_JVM_ARGS="-Xms32g -Xmx64g"
export NXF_SINGULARITY_CACHEDIR=`pwd`/singularity
export NXF_CACHE_DIR=`pwd`/nxfcache

if [ ! -f nextflow ]
then
    curl -s https://get.nextflow.io | bash
fi

## Collect samples
if [ ! -f samples.csv ]
then
    echo "patient,sex,status,sample,lane,fastq_1,fastq_2" > samples.csv
    for BAM in ${BASEDIR}/../control/*.bam
    do
	NAME=`echo ${BAM} | sed 's/^.*\///' |  sed 's/.bam$//'`
	echo ${NAME}
	FQ1=${NAME}.R1.fastq.gz
	FQ2=${NAME}.R2.fastq.gz
	samtools collate -uO ${BASEDIR}/../control/*.bam | samtools fastq -1 ${FQ1} -2 ${FQ2} -0 /dev/null -s /dev/null -
	if [ -f ${FQ1} ]
	then
	    if [ -f ${FQ2} ]
	    then
		SEX="NA"
		ID=`echo ${FQ1} | sed 's/^.*\///' |  sed 's/.R1.fastq.gz$//' | sed 's/benigntissue.*$//' | sed 's/germline.*$//'`
		TYPE=`echo ${FQ1} | sed 's/^.*\///' |  sed 's/.R1.fastq.gz$//' | sed "s/^${ID}//"`
		FQ2=`echo ${FQ1} | sed 's/.R1.fastq.gz$/.R2.fastq.gz/'`

		echo ${ID}
		if [ `cat sex.tsv | grep -w ${ID} | cut -f 2` == "F" ]; then SEX="XX"; fi
		if [ `cat sex.tsv | grep -w ${ID} | cut -f 2` == "M" ]; then SEX="XY"; fi
		STATUS="NA"
		if [ ${TYPE} == "germline" ]; then STATUS=0; fi
		if [ ${TYPE} == "benigntissue" ]; then STATUS=1; fi
		if [ ${TYPE} == "tumor" ]; then STATUS=1; fi
		echo ${ID} ${TYPE} ${STATUS} ${SEX}
		if [ ! -d ${ID} ]
		then
		    echo "${ID},${SEX},${STATUS},${ID}${TYPE},lane_1,${FQ1},${FQ2}" >> samples.csv
		fi
	    fi
	fi
    done
fi

## nf-core/sarek
./nextflow run nf-core/sarek -resume -profile singularity -revision 3.5.0 --input samples.csv --genome GATK.GRCh38 --outdir sarek/ --tools 'manta,mutect2,strelka,vep'
