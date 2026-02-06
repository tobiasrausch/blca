#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/mamba/bin:${PATH}

HG=${BASEDIR}/../genome/hg38.fa

REFPANEL=CoRAL/reference/hg38full_ref_5k.cnn
REF=hg38
GAINTH=4.0
MINBPS=2

if [ ! -d CoRAL ]
then
    git clone https://github.com/AmpliconSuite/CoRAL
    cd CoRAL/
    git checkout a7f8ece28a1452b8682bb69e8f68c548b54e025a
    cd ../
fi

FD=outcnv
mkdir -p ${FD}
for BAM in ${BASEDIR}/../alignment/*.bam
do
    ID=`echo ${BAM} | sed 's/^.*\///' | sed 's/.bam$//'`
    echo ${ID}
    if [ ! -f ${FD}/${ID}/${ID}.cns ]
    then
	cnvkit.py batch $BAM --seq-method wgs --drop-low-coverage --reference ${REFPANEL} --scatter --diagram -d ${FD}/${ID}
    fi

    SEG=${FD}/${ID}/${ID}.cns
    ls ${SEG}
    if [ -f ${SEG} ]
    then
	if [ ! -f ${ID}_reconstruct.log ]
	then
	    ## Seed
	    python3 ./CoRAL/src/CoRAL.py seed --cn_seg ${SEG} --out ${ID}_CNV_SEEDS.bed --gain ${GAINTH} --min_seed_size 100000 --max_seg_gap 300000
	    ## Reconstruct
	    python3 ./CoRAL/src/CoRAL.py reconstruct --lr_bam ${BAM} --cnv_seed ${ID}_CNV_SEEDS.bed --output_prefix ${ID} --cn_seg ${SEG} --log_fn ${ID}_reconstruct.log --min_bp_support ${MINBPS} --cycle_decomp_time_limit 15000 --cycle_decomp_threads 16
	fi
    fi

    ## Plotting & convert to BED
    for G in ${ID}_*_graph.txt
    do
	if [ -f ${G} ]
	then
	    OUTP=`echo ${G} | sed 's/_graph.txt$//'`
	    echo ${OUTP}
	    python3 ./CoRAL/src/CoRAL.py plot --ref ${REF} --bam ${BAM} --plot_graph --graph ${G} --output_prefix ${OUTP}
	fi
    done
    for C in ${ID}_*_cycles.txt
    do
	if [ -f ${C} ]
	then
	    OUTP=`echo ${C} | sed 's/_cycles.txt$//'`
	    echo ${OUTP}
	    python3 ./CoRAL/src/CoRAL.py plot --ref ${REF} --bam ${BAM} --plot_cycles --cycles ${C} --output_prefix ${OUTP}
	fi
    done
done
