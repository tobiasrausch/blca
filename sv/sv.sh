#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../mamba/bin:${PATH}

HG=${BASEDIR}/../genome/hg38.fa

if [ ! -f human_GRCh38_no_alt_analysis_set.trf.bed ]
then
    wget https://raw.githubusercontent.com/KolmogorovLab/Severus/refs/heads/main/vntrs/human_GRCh38_no_alt_analysis_set.trf.bed
fi

if [ ! -f PoN_1000G_hg38.tsv.gz ]
then
    wget https://github.com/KolmogorovLab/Severus/raw/refs/heads/main/pon/PoN_1000G_hg38.tsv.gz
fi

for BAM in ${BASEDIR}/../alignment/*.bam
do
    if [ -f ${BAM} ]
    then
	ID=`echo ${BAM} | sed 's/^.*\///' | sed 's/.bam$//'`
	echo ${ID}

	## Severus
	severus --target-bam ${BAM} --out-dir ${ID} -t 32 --vntr-bed human_GRCh38_no_alt_analysis_set.trf.bed --PON PoN_1000G_hg38.tsv.gz
	bgzip ${ID}/somatic_SVs/severus_somatic.vcf
	tabix ${ID}/somatic_SVs/severus_somatic.vcf.gz

	## Delly
	delly lr -g ${HG} -o ${ID}.delly.bcf ${BAM}
    fi
done
