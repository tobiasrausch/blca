#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/mamba/bin:${PATH}

HG=${BASEDIR}/../genome/hg38.fa

## Merge MEI annotations and intersect with somatic SVs
for ID in `ls *.ins.vcf.gz *.del.vcf.gz | sed 's/.ins.vcf.gz//' | sed 's/.del.vcf.gz//' | sort | uniq`
do
    echo ${ID}
    if [ ! -f ${ID}.vcf.gz ]
    then
	bcftools concat -a ${ID}.del.vcf.gz ${ID}.ins.vcf.gz | grep "^#" > ${ID}.vcf
	zcat ${ID}.del.vcf.gz ${ID}.ins.vcf.gz | grep -v "^#" | grep -w -Ff <(bcftools view ${BASEDIR}/../consensusSV/${ID}.somatic.bcf | cut -f 1,2,3,4,5) >> ${ID}.vcf
	bcftools sort ${ID}.vcf | bgzip > ${ID}.vcf.gz
	tabix ${ID}.vcf.gz
	rm ${ID}.vcf
    fi
done

## Statistics
for BCF in *tumor.vcf.gz
do
    bcftools query -f "%ID\t%ITYPE_N\t%FAM_N\n" ${BCF} | grep "^INS" | cut -f 2- | sort | uniq -c
done
