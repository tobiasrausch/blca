#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../mamba/bin:${PATH}

for DELLY in ${BASEDIR}/../sv/*.delly.bcf
do
    ID=`echo ${DELLY} | sed 's/^.*\///' | sed 's/.delly.bcf$//'`
    PID=`echo ${ID} | sed 's/tumor//' | sed 's/B060/B60/'`
    if [ -f ${DELLY} ]
    then
	SEV=`ls ../sv/${ID}/somatic_SVs/severus_somatic.vcf.gz`
	if [ -f ${SEV} ]
	then
	    echo ${ID}
	    ## Consensus
	    sansa compvcf --nosvt -a ${SEV} -e 0 -m 0 -n 250000000 ${DELLY}
	    cat out.sv.classification | awk '$2=="TP"' | cut -f 1 | sort | uniq > delly.tp
	    bcftools view ${DELLY} | grep -v "^#" | cut -f 3 | sort | uniq > all.sv
	    cat all.sv | grep -v -w -Ff delly.tp > remove.sv
	    bcftools view ${DELLY} chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY | grep -v -w -Ff remove.sv | bcftools view -O b -o ${ID}.somatic.bcf -
	    rm delly.tp all.sv remove.sv out.sv.classification out.tsv
	    bcftools index -f ${ID}.somatic.bcf

	    ## Implant and genotype SVs using short-reads
	    if [ -f genotypingSR/${PID}.depth.regions.bed.gz ]
	    then
		COV=`grep -P "^total\t" genotypingSR/${PID}.depth.mosdepth.summary.txt | cut -f 4`
		zcat genotypingSR/${PID}.depth.regions.bed.gz | sed "s/$/\t${COV}/" | awk '$5/$6 > 0.5' | cut -f 4 | sort | uniq > remove.sv
		bcftools view ${ID}.somatic.bcf | grep -v -w -Ff remove.sv | bcftools view -O b -o ${ID}.somatic.bcf.tmp -
		rm remove.sv
		mv ${ID}.somatic.bcf.tmp ${ID}.somatic.bcf
		bcftools index -f ${ID}.somatic.bcf
	    fi
	fi
    fi
done

## Statistics
for VCF in *.somatic.bcf
do
    ID=`echo ${VCF} | sed 's/.somatic.bcf//'`
    echo ${ID} `bcftools view ${VCF} | grep -v "^#" | cut -f 3 | cut -c 1-3 | sort | uniq -c`
    bcftools view ${VCF} | grep -v "^#" | wc -l
done
