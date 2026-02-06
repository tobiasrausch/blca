#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/mamba/bin:${PATH}

HG=${BASEDIR}/../../genome/hg38.fa
for F in ../*.somatic.bcf
do
    ID=`echo ${F} | sed 's/^.*\///' | sed 's/.somatic.bcf//'`
    PID=`echo ${ID} | sed 's/tumor//'`
    GERM=`echo ${ID} | sed 's/tumor/germline/' | sed 's/B060/B60/'`
    if [ -f ${BASEDIR}/../../control/${GERM}.bam ]
    then
	echo ${ID}

	## Implant SVs
	if [ ! -f ${ID}.new_ref.bed ]
	then
	    bcftools view -e 'ALT~"<.*>"' ${F} | bcftools view -i 'SVTYPE=="INS"' - | bcftools norm --check-ref s -m- -f ${HG} - | bgzip > ${ID}.norm.vcf.gz
	    tabix ${ID}.norm.vcf.gz
	    cat ${HG} | vcf-consensus ${ID}.norm.vcf.gz > ${ID}.new_ref.fa
	
	    ## Translate coordinates
	    bcftools query -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\n" ${ID}.norm.vcf.gz > ${ID}.coord.tsv
	    OLDCHR="NA"
	    CUMSUM=0
	    rm -f ${ID}.new_ref.bed
	    while read CHR POS SVID REF ALT
	    do
		if [ ${OLDCHR} != ${CHR} ]
		then
		    CUMSUM=0
		    OLDCHR=${CHR}
		fi
		ILEN=`echo -e "${REF}\t${ALT}" | awk '{print length($2)-length($1);}'`
		NEWS=`expr $POS + $CUMSUM + 1`
		NEWE=`expr $NEWS + $ILEN`
		CUMSUM=`expr ${CUMSUM} + ${ILEN}`
		echo -e "${CHR}\t${NEWS}\t${NEWE}\t${SVID}" >> ${ID}.new_ref.bed
	    done < ${ID}.coord.tsv
	fi
	## Remap to new reference
	if [ -f ${ID}.new_ref.fa ]
	then
	    ## Convert to FASTQ
	    samtools collate -uO ${BASEDIR}/../../control/${GERM}.bam | samtools fastq -1 ${ID}.1.fq -2 ${ID}.2.fq -0 /dev/null -s /dev/null -n -
	    FQ1=${ID}.1.fq
	    FQ2=${ID}.2.fq
	    
	    ## Remap to new reference
	    if [ -f ${FQ1} ]
	    then
		if [ -f ${FQ2} ]
		then
		    if [ ! -f ${ID}.new_ref.fa.pac ]
		    then
			bwa index ${ID}.new_ref.fa
		    fi
		    
		    if [ ! -f ${PID}.bam ]
		    then
			if [ ! -f ${PID}.srt.bam ]
			then
			    bwa mem -R "@RG\tID:${PID}\tSM:${PID}" -t 16 ${ID}.new_ref.fa ${FQ1} ${FQ2} | samtools fixmate -m - - | samtools sort -o ${PID}.srt.bam -
			    samtools index ${PID}.srt.bam
			fi
			
			## Mark duplicates
			samtools markdup ${PID}.srt.bam ${PID}.bam
			samtools index ${PID}.bam
		    fi
		    
		    ## Depth
		    mosdepth -m -q 1 -n --by ${ID}.new_ref.bed ${PID}.depth ${PID}.bam
		fi
	    fi
	fi
    fi
done
