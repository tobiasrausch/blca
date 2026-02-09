#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../mamba/bin:${PATH}

GENOME=${BASEDIR}/../genome/hg38.fa
EXCL=human.hg38.excl.tsv
MAP=Homo_sapiens.GRCh38.dna.primary_assembly.fa.r101.s501.blacklist.gz
WIN=10000

if [ ! -f ${MAP} ]
then
    wget https://gear-genomics.embl.de/data/delly/Homo_sapiens.GRCh38.dna.primary_assembly.fa.r101.s501.blacklist.gz.gzi
    wget https://gear-genomics.embl.de/data/delly/Homo_sapiens.GRCh38.dna.primary_assembly.fa.r101.s501.blacklist.gz.fai
    wget https://gear-genomics.embl.de/data/delly/Homo_sapiens.GRCh38.dna.primary_assembly.fa.r101.s501.blacklist.gz
fi

if [ ! -f rd.R ]
then
    wget https://raw.githubusercontent.com/dellytools/delly/refs/heads/main/R/rd.R
fi

if [ ! -f ${EXCL} ]
then
    wget https://raw.githubusercontent.com/dellytools/delly/refs/heads/main/excludeTemplates/human.hg38.excl.tsv
fi

for F in ${BASEDIR}/../control/*.bam
do
    ID=`echo ${F} | sed 's/^.*\///' |  sed 's/.bam$//'`
    ## SV calling
    if [ ! -f ${ID}.delly.bcf ]
    then
	delly call -g ${GENOME} -x ${EXCL} -o ${ID}.delly.bcf ${F}
    fi
    ## Coverage
    if [ ! -f ${ID}.cov.gz ]
    then
	delly cnv -i ${WIN} -j ${WIN} -w ${WIN} -g ${GENOME} -m ${MAP} -o ${ID}.cov.bcf -c ${ID}.cov.gz ${F}
    fi
    ## Read-depth profile
    if [ ! -d ${ID} ]
    then
	echo ${ID}
	rm -rf ${ID}/
	mkdir -p ${ID}
	Rscript rd.R ${ID}.cov.gz
	mv plot.*.png ${ID}/
    fi
done
