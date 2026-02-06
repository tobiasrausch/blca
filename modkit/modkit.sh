#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../mamba/bin:${PATH}

if [ ! -f ./dist_modkit_v0.6.1_481e3c9/modkit ]
then
    wget https://github.com/nanoporetech/modkit/releases/download/v0.6.1/modkit_v0.6.1_u16_x86_64.tar.gz
    tar -xzf modkit_v0.6.1_u16_x86_64.tar.gz 
fi

HG=${BASEDIR}/../genome/hg38.fa
for BAM in ${BASEDIR}/../alignment/*.bam
do
    ID=`echo ${BAM} | sed 's/^.*\///' | sed 's/.bam$//'`
    echo ${ID}
    if [ ! -f ${ID}.pileup.bed.gz ]
    then
	./dist_modkit_v0.6.1_481e3c9/modkit pileup --cpg --reference ${HG} --modified-bases 5mC --combine-strands --threads 8 --log-filepath ${ID}.log ${BAM} ${ID}.pileup.bed
	bgzip ${ID}.pileup.bed
	tabix -p bed ${ID}.pileup.bed.gz
    fi
done
