#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")

# Nextflow parameters
export NXF_JVM_ARGS="-Xms50g -Xmx140g"
export NXF_SINGULARITY_CACHEDIR=`pwd`/singularity
export NXF_CACHE_DIR=`pwd`/nxfcache

# Temporary directories
mkdir -p `pwd`/tmp
export SINGULARITYENV_TMPDIR=`pwd`/tmp/
export APPTAINERENV_TMPDIR=`pwd`/tmp/
export SINGULARITYENV_NXF_TASK_WORKDIR=`pwd`/tmp/
export APPTAINERENV_NXF_TASK_WORKDIR=`pwd`/tmp/
export SINGULARITYENV_NXF_DEBUG=`pwd`/tmp/
export APPTAINERENV_NXF_DEBUG=`pwd`/tmp/

if [ ! -f nextflow ]
then
    curl -s https://get.nextflow.io | bash
fi

if [ ! -f Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz ]
then
    wget https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
fi

if [ ! -f Homo_sapiens.GRCh38.113.gtf.gz ]
then
    wget https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz
fi

## Collect samples
if [ ! -f samples.csv ]
then
    echo "sample,fastq_1,fastq_2,strandedness" > samples.csv
    for FQ1 in ${BASEDIR}/../rna/*.1.fq.gz
    do
	ID=`echo ${FQ1} | sed 's/^.*\///' |  sed 's/.1.fq.gz$//' | sed 's/$/tumor/'`
	FQ2=`echo ${FQ1} | sed 's/.1.fq.gz/.2.fq.gz/'`
	echo ${ID}
	if [ -f ${FQ1} ]
	then
	    if [ -f ${FQ2} ]
	    then
		echo "${ID},${FQ1},${FQ2},auto" >> samples.csv
	    fi
	fi
    done
fi

./nextflow run nf-core/rnaseq -resume -profile singularity -revision 3.18.0 --input samples.csv --fasta Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz --gtf Homo_sapiens.GRCh38.113.gtf.gz --outdir rnaseq/
