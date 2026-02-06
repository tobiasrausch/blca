#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../mamba/bin:${PATH}

for F in ../deepsomatic/*tumor.bcf
do
    ID=`echo ${F} | sed 's/^.*\///' | sed 's/tumor.bcf//'`
    if [ `echo ${ID} | grep "_" | wc -l` -eq 1 ]; then continue; fi
    echo ${ID}
    DS=../deepsomatic/${ID}tumor.bcf
    SC=../clairS/${ID}tumor/snv.vcf.gz
    SI=../clairS/${ID}tumor/indel.vcf.gz
    if [ -f ${DS} ]
    then
	if [ -f ${SC} ]
	then
	    if [ -f ${SI} ]
	    then
		if [ ! -f ${ID}tumor.consensus.bcf ]
		then
		    ## Concat SNVs and InDels for clairS
		    bcftools concat -a -O b -o ${ID}tumor.raw.clair3.bcf ${SC} ${SI}
		    bcftools index ${ID}tumor.raw.clair3.bcf
		    
		    ## clairS: min. 5x coverage
		    bcftools view -f 'PASS,.' -v snps,indels -m 2 -M 2 -i 'sum(AD) > 5' -O b -o ${ID}tumor.filtered.clair3.bcf ${ID}tumor.raw.clair3.bcf chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX
		    bcftools index ${ID}tumor.filtered.clair3.bcf
		    rm ${ID}tumor.raw.clair3.bcf ${ID}tumor.raw.clair3.bcf.csi
		    
		    ## DeepSomatic: SNVs on autosomes, min. 5x coverage
		    bcftools view -f 'PASS,.' -v snps,indels -m 2 -M 2 -i 'sum(AD) > 5' -O b -o ${ID}tumor.filtered.ds.bcf ${DS} chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX
		    bcftools index ${ID}tumor.filtered.ds.bcf
		    
		    ## Consensus
		    bcftools isec -O b -o ${ID}tumor.consensus.bcf -n=2 -w1 ${ID}tumor.filtered.clair3.bcf ${ID}tumor.filtered.ds.bcf
		    bcftools index ${ID}tumor.consensus.bcf
		    rm ${ID}tumor.filtered.ds.bcf ${ID}tumor.filtered.ds.bcf.csi
		    rm ${ID}tumor.filtered.clair3.bcf ${ID}tumor.filtered.clair3.bcf.csi
		fi

		## Filter against short-read controls
		if [ ! -f ${ID}tumor.filtered.bcf ]
		then
		    if [ -f genotypeSR/${ID}tumor.bcf ]
		    then
			bcftools view -i 'AC>0' -O b -o ${ID}.blood.bcf genotypeSR/${ID}tumor.bcf
			bcftools index ${ID}.blood.bcf

			## Remove blood calls
			bcftools isec -O b -o ${ID}tumor.filtered.bcf ${ID}tumor.consensus.bcf ${ID}.blood.bcf -n -1 -w 1
			bcftools index ${ID}tumor.filtered.bcf
			rm ${ID}.blood.bcf ${ID}.blood.bcf.csi
		    fi
		fi
	    fi
	fi
    fi
done
