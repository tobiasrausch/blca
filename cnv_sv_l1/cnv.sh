#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../mamba/bin:${PATH}

HG=${BASEDIR}/../genome/hg38.fa

## PCAWG somatic L1s
if [ ! -f SupplT2.csv ]
then
    if [ ! -f 41588_2019_562_MOESM3_ESM.xlsx ]
    then
	wget https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-019-0562-0/MediaObjects/41588_2019_562_MOESM3_ESM.xlsx
    fi
    echo "Please export the suppl_table2 sheet in CSV format from 41588_2019_562_MOESM3_ESM.xlsx as SupplT2.csv and rerun the script."
    exit;
fi

## PCAWG somatic SV and CNV calls
if [[ -z "${PCAWG_PATH}" ]]; then
    echo "Please set the PCAWG_PATH environment variable to the access-controlled somatic SV and CNV calls."
    exit;
fi

## Get mobile elements
if [ ! -f mei_per_sample.tsv ]
then
    cat SupplT2.csv  | cut -f 1-3 -d ',' | grep -v '"' | sed 's/,/\t/g' | sort | uniq -c  | awk '$1>=1 {print $2"\t"$3"\t"$4"\t"$1;}' > mei_per_sample.tsv
    cat mei_per_sample.tsv  | cut -f 3 | sort | uniq -c | awk '$1>=20 {print $2;}' > selected_hist
    echo "Tumor histologies with at least 20 samples that have at least 1 somatic MEI: " `cat selected_hist | wc -l`
    cat mei_per_sample.tsv | grep -w -Ff selected_hist > mei_per_sample.tsv.tmp
    mv mei_per_sample.tsv.tmp mei_per_sample.tsv
    rm selected_hist
fi

## Process each sample
echo -e "key\tsample\thistology\tsvdist\tl1dist\tl1count\tl1sv\tl1loose\tenrichment" > summary.tsv
while read icgc_donor_id tumor_wgs_icgc_sample_id histology l1count
do
    ID=`zcat release_may2016.v1.2.tsv.gz | cut -f 27,28 | grep -w "^${tumor_wgs_icgc_sample_id}" | cut -f 2 | head -n 1`
    if [ `echo ${ID} | awk '{print length($1);}'` -gt 10 ]
    then
	SVFILE=${PCAWG_PATH}/pcawg_consensus_1.6.161116.somatic_svs/${ID}.pcawg_consensus_1.6.161116.somatic.sv.bedpe.gz
	CNVFILE=${PCAWG_PATH}/consensus.20170119.somatic.cna/${ID}.consensus.20170119.somatic.cna.txt
	if [ -f ${SVFILE} ]
	then
	    if [ -f ${CNVFILE} ]
	    then
		if [ `cat ${CNVFILE} | wc -l` -gt 0 ]
		then
		    ## Get the L1 insertion positions
		    cat SupplT2.csv  | cut -f 1-5 -d ',' | grep -w "${icgc_donor_id}" | grep -w "${tumor_wgs_icgc_sample_id}" | grep -w "${histology}" | cut -f 4,5 -d ',' | sed 's/,/\t/' | sort -k1,1V -k2,2n | awk '{print "chr"$1"\t"($2)"\t"($2+1);}' | grep -P "^chr[0-9X]*\t" | sed 's/^chr//' > ${ID}.l1.bed
		    echo ${l1count} ${ID} `wc -l ${ID}.l1.bed`
		    ./analyze_l1_cna_sv.py --cna ${CNVFILE} --sv ${SVFILE} --l1 ${ID}.l1.bed | sed "s/.l1.bed/\t${histology}/" >> summary.tsv
		    rm ${ID}.l1.bed
		fi
	    fi
	fi
    fi
done < mei_per_sample.tsv

## Plot results
cat summary.tsv  | egrep "^key|Stats" > summary.table
Rscript cnv_sv_l1.R
