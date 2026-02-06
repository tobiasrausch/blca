#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/mamba/bin:${PATH}

GENOME=${BASEDIR}/../../genome/hg38.fa
for F in ../*.consensus.bcf
do
    ID=`echo ${F} | sed 's/^.*\///' | sed 's/.consensus.bcf//'`
    GERM=`echo ${ID} | sed 's/tumor/germline/' | sed 's/B060/B60/'`
    if [ -f ${BASEDIR}/../../control/${GERM}.bam ]
    then
	echo ${ID}
	FILES=""
	for CHR in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY
	do
	    bcftools view ${F} ${CHR} | bgzip > ${ID}.${CHR}.consensus.vcf.gz
	    tabix ${ID}.${CHR}.consensus.vcf.gz
	    freebayes -r ${CHR} -l -@ ${ID}.${CHR}.consensus.vcf.gz --no-partial-observations --min-repeat-entropy 1 --report-genotype-likelihood-max --min-alternate-fraction 0.01 --fasta-reference ${GENOME} --genotype-qualities ${BASEDIR}/../../control/${GERM}.bam -v ${ID}.${CHR}.vcf
	    rm ${ID}.${CHR}.consensus.vcf.gz ${ID}.${CHR}.consensus.vcf.gz.tbi
	    bgzip ${ID}.${CHR}.vcf
	    tabix ${ID}.${CHR}.vcf.gz
	    FILES=${FILES}" "${ID}.${CHR}.vcf.gz
	done
	bcftools concat -a -O b -o ${ID}.bcf ${FILES}
	bcftools index ${ID}.bcf
	rm ${ID}.chr*.vcf.gz ${ID}.chr*.vcf.gz.tbi
    fi
done
