#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=`pwd`/mamba/bin:${PATH}

HG=${BASEDIR}/../genome/hg38.fa

if [ ! -d SVAN/ ]
then
    git clone https://github.com/REPBIO-LAB/SVAN.git 
fi

if [ ! -f trf409.linux64 ]
then
    wget https://github.com/Benson-Genomics-Lab/TRF/releases/download/v4.09.1/trf409.linux64
    chmod a+x trf409.linux64
fi
    
for BCF in ${BASEDIR}/../sv/*.delly.bcf
do
    ID=`echo ${BCF} | sed 's/^.*\///' | sed 's/.delly.bcf//'`
    echo ${ID}

    ## Insertions
    mkdir -p ${ID}
    bcftools view ${BCF} | head -n 1 > ${ID}.ins.vcf
    awk '{print "##contig=<ID=" $1 ",length=" $2 ">"}' ${HG}.fai >> ${ID}.ins.vcf
    bcftools view -i '(QUAL>300) && (STRLEN(ALT)>(STRLEN(REF)+50)) && (SVTYPE=="INS")' ${BCF} chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY | grep -v "^##contig" | grep -v "^##fileformat" >> ${ID}.ins.vcf
    PYTHONPATH=`pwd`/SVAN/GAPI:`pwd`/SVAN python3 SVAN/scripts/ins2fasta.py ${ID}.ins.vcf ${ID}
    ./trf409.linux64 ${ID}/insertions_seq.fa 2 7 7 80 10 10 500 -h -d -ngs 1 > ${ID}/ins_trf.out
    PYTHONPATH=`pwd`/SVAN/GAPI:`pwd`/SVAN python3 SVAN/SVAN-INS.py ${ID}.ins.vcf ${ID}/ins_trf.out hg38/VNTR_hg38.bed hg38/EXONS_hg38.bed hg38/REPEATS_hg38.bed hg38/CONSENSUS.fa ${HG} ${ID}
    mv ${ID}.vcf ${ID}.ins.vcf
    bgzip ${ID}.ins.vcf
    tabix ${ID}.ins.vcf.gz
    rm -rf ${ID}/ tmp/

    ## Deletions
    mkdir -p ${ID}
    bcftools view ${BCF} | head -n 1 > ${ID}.del.vcf
    awk '{print "##contig=<ID=" $1 ",length=" $2 ">"}' ${HG}.fai >> ${ID}.del.vcf
    bcftools view -i '(QUAL>300) && (STRLEN(REF)>(STRLEN(ALT)+50)) && (SVTYPE=="DEL")' ${BCF} chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY | grep -v "^##contig" | grep -v "^##fileformat" >> ${ID}.del.vcf
    PYTHONPATH=`pwd`/SVAN/GAPI:`pwd`/SVAN python3 SVAN/scripts/del2fasta.py ${ID}.del.vcf ${ID}
    ./trf409.linux64 ${ID}/deletions_seq.fa 2 7 7 80 10 10 500 -h -d -ngs 1 > ${ID}/del_trf.out
    PYTHONPATH=`pwd`/SVAN/GAPI:`pwd`/SVAN python3 SVAN/SVAN-DEL.py ${ID}.del.vcf ${ID}/del_trf.out hg38/VNTR_hg38.bed hg38/EXONS_hg38.bed hg38/REPEATS_hg38.bed hg38/CONSENSUS.fa ${HG} ${ID}
    mv ${ID}.vcf ${ID}.del.vcf
    bgzip ${ID}.del.vcf
    tabix ${ID}.del.vcf.gz
    rm -rf ${ID}/ tmp/
done
