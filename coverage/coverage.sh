#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../mamba/bin:${PATH}

GENOME=${BASEDIR}/../genome/hg38.fa
MAP=Homo_sapiens.GRCh38.dna.primary_assembly.fa.r101.s501.blacklist.gz
WIN=25000

if [ ! -f ${MAP} ]
then
    wget https://gear-genomics.embl.de/data/delly/Homo_sapiens.GRCh38.dna.primary_assembly.fa.r101.s501.blacklist.gz.gzi
    wget https://gear-genomics.embl.de/data/delly/Homo_sapiens.GRCh38.dna.primary_assembly.fa.r101.s501.blacklist.gz.fai
    wget https://gear-genomics.embl.de/data/delly/Homo_sapiens.GRCh38.dna.primary_assembly.fa.r101.s501.blacklist.gz
fi

for F in ${BASEDIR}/../alignment/*.bam
do
    ID=`echo ${F} | sed 's/^.*\///' |  sed 's/.bam$//'`
    if [ ! -f ${ID}.cov.gz ]
    then
	delly cnv -i ${WIN} -j ${WIN} -w ${WIN} -g ${GENOME} -m ${MAP} -o ${ID}.cov.bcf -c ${ID}.cov.gz ${F}
    fi
    if [ ! -d ${ID} ]
    then
	echo ${ID}
	rm -rf ${ID}/
	mkdir -p ${ID}
	Rscript rd.R ${ID}.cov.gz
	mv plot.*.png ${ID}/
    fi
done
