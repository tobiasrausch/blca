#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")

export PATH=${BASEDIR}/mamba/bin:${PATH}

if [ ! -d ${BASEDIR}/L1EM ]
then
    git clone https://github.com/FenyoLab/L1EM    
fi

if [ ! -f hg38.fa ]
then
    ## Build GRCh38 reference
    wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
    zcat hg38.fa.gz > hg38.fa
    samtools faidx hg38.fa
    bwa index hg38.fa
    rm hg38.fa.gz

    ## Build L1EM reference
    cd ${BASEDIR}/L1EM/ && bash generate_L1EM_fasta_and_index.sh ${BASEDIR}/hg38.fa && cd ../
fi

## Run L1EM
for FQ1 in ${BASEDIR}/../rna/*.1.fq.gz
do
    FQ2=`echo ${FQ1} | sed 's/.1.fq.gz$/.2.fq.gz/'`
    if [ -f ${FQ2} ]
    then
	ID=`echo ${FQ1} | sed 's/^.*\///' | sed 's/.1.fq.gz$/tumor/'`
	echo ${ID}
	if [ ! -d ${ID} ]
	then
	    if [ ! -f ${ID}.bam ]
	    then
		bwa mem -R "@RG\tID:${ID}\tSM:${ID}" -t 24 ${BASEDIR}/hg38.fa ${FQ1} ${FQ2} | samtools fixmate -m - - | samtools sort -o ${ID}.bam -
		samtools index ${ID}.bam
	    fi
	    bash -e L1EM/run_L1EM.sh ${BASEDIR}/${ID}.bam ${BASEDIR}/L1EM/ ${BASEDIR}/hg38.fa
	    mkdir -p ${ID}
	    mv X_final.pkl names_final.pkl filter_L1HS_FPM.txt l1hs_transcript_counts.txt full_counts.txt baminfo.txt ${ID}/
	    rm -rf G_of_R/ split_fqs/ idL1reads/
	    rm -rf ${ID}.bam ${ID}.bam.bai
	fi
    fi
done
