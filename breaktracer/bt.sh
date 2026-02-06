#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../mamba/bin:${PATH}

if [ ! -f breaktracer_v0.0.5_linux_x86_64bit ]
then
    wget https://github.com/tobiasrausch/breaktracer/releases/download/v0.0.5/breaktracer_v0.0.5_linux_x86_64bit
    chmod a+x breaktracer_v0.0.5_linux_x86_64bit
fi

REF=${BASEDIR}/../genome/hg38.fa
for BAM in ${BASEDIR}/../alignment/*.bam
do
    ID=`echo ${BAM} | sed 's/^.*\///' | sed 's/.bam$//'`
    if [ ! -f ${ID}.telseq.tsv ]
    then
	echo ${ID}
	./breaktracer_v0.0.5_linux_x86_64bit find -o ${ID}.l1.tsv -g ${REF} ${BAM} 
    fi
done
