#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/mamba/bin:${PATH}

source activate clairs-to

## Call
REF=${BASEDIR}/../genome/hg38.fa
for F in ${BASEDIR}/../alignment/*.bam
do
    ID=`echo ${F} | sed 's/^.*\///' | sed 's/.bam$//'`
    echo ${ID}
    if [ ! -d ${ID} ]
    then
	./ClairS-TO/run_clairs_to -s ${ID} -T ${F} -R ${REF} -o ${ID} -t 48 -p 'ont_r10_dorado_sup_4khz' > ${ID}.log 2> ${ID}.err
	rm -rf ${ID}/logs/ ${ID}/tmp/
    fi
done
