#!/bin/bash

# GENOME="phaCin"
# DATA_PATH="/data/tusers/boxu/lrft/rawdata/nanopore/Koala"
# OUT_PATH="/data/tusers/boxu/lrft/result/nanopore/Koala"
# ANNO_PATH="/data/tusers/boxu/annotation"
# TE_REF=${ANNO_PATH}/${GENOME}/phaCin.transposon.fa
# SAMPLES="Koala2_20191231_Adelaide_rep1"
# CPU=10

# GENOME="hs37d5"
# DATA_PATH="/data/tusers/boxu/lrft/rawdata/nanopore/human"
# OUT_PATH="/data/tusers/boxu/lrft/result/nanopore/human"
# ANNO_PATH="/data/tusers/boxu/annotation"
# TE_REF=${ANNO_PATH}/${GENOME}/ALUL1SVA.2.fa
# SAMPLES="SRR11669560"
# CPU=10

# GENOME="dm6_clean"
# DATA_PATH="/data/tusers/boxu/lrft/rawdata/nanopore/fly_gut"
# OUT_PATH="/data/tusers/boxu/lrft/result/nanopore/fly_gut"
# ANNO_PATH="/data/tusers/boxu/annotation"
# TE_REF=${ANNO_PATH}/${GENOME}/dm6_clean.transposon.2.fa
# SAMPLES="ERR4920933,ERR4920934,ERR4920935,ERR4920936"
# CPU=10


# GENOME="hs37d5"
# BIN_PATH="/data/tusers/boxu/lrft/scripts"
# DATA_PATH="/data/tusers/boxu/lrft/rawdata/nanopore"
# # DATA_PATH="/data/tusers/boxu/lrft/result/nanopore/human/simulation/nanosim/simulate3"
# OUT_PATH="/data/tusers/boxu/lrft/result/nanopore"
# # OUT_PATH="/data/tusers/boxu/lrft/result/nanopore/human/simulation/nanosim/simulate3"
# ANNO_PATH="/data/tusers/boxu/annotation"
# CPU=10
# # Koala_sample="Koala2_20191231_Adelaide_rep1"
# # human_sample="SRR11669560"

# # human_sample="SRR11669560"
# human_sample="HG002_ONT-UL_GIAB_20200122"



while getopts ":g:i:o:n:t:s:c:" OPTION;
do
	case $OPTION in
		g)	GENOME=${OPTARG};;
		i)  DATA_PATH=${OPTARG};;
		o)	OUT_PATH=${OPTARG};;
		n)	ANNO_PATH=${OPTARG};;
		t)  TE_REF=${OPTARG};;
		s)	SAMPLES=${OPTARG};;
		c)	CPU=${OPTARG};;
	esac
done
# lrft_map.sh -g hs37d5 -i data/tusers/boxu/lrft/rawdata/nanopore/human -o /data/tusers/boxu/lrft/result/nanopore/human -n /data/tusers/boxu/annotation -t /data/tusers/boxu/annotation/hs37d5/ALUL1SVA.2.fa -s SRR11669560 -c 10
# lrft_map.sh -g phaCin -i /data/tusers/boxu/lrft/rawdata/nanopore/Koala -o /data/tusers/boxu/lrft/result/nanopore/Koala -n /data/tusers/boxu/annotation -t /data/tusers/boxu/annotation/hs37d5/phaCin.transposon.fa -s Koala2_20191231_Adelaide_rep1 -c 10

# 
if_uniq=
if_remap=
if_genome_first=
# lrft_temp.sh -g hs37d5 -i /data/tusers/boxu/lrft/rawdata/nanopore -o /data/tusers/boxu/lrft/result/nanopore -s SRR11669560


TE_temp=${TE_REF##*/}
TE=${TE_temp%.*}
GENOME_REF=${ANNO_PATH}/${GENOME}/${GENOME}.fa
SAMPLES=(${SAMPLES//,/ })

BIN_PATH=$(cd "$(dirname "$0")";pwd)

SAMPLE="HG002_ONT-UL_GIAB_20200122"

${BIN_PATH}/lrft.map.sh -g ${GENOME} -i ${DATA_PATH} -o ${OUT_PATH}/lrft -n ${ANNO_PATH} -t ${ANNO_PATH}/${GENOME}/ALUL1SVA.2.fa -s ${SAMPLE} -c ${CPU}
python ${BIN_PATH}/lrft_insertion.v8.py ${OUT_PATH}/lrft ${SAMPLE}.mapped.sorted.tcg.bam ${ANNO_PATH}/${GENOME}/${GENOME}.fa ${SAMPLE} 2>${OUT_PATH}/lrft.log


# tldr -b /data/tusers/boxu/lrft/result/nanopore/human/map/${human_sample}_0.mapped.sorted.tcg.bam -e /data/tusers/boxu/annotation/hs37d5/ALUL1SVA.2.fa -r /data/tusers/boxu/annotation/hs37d5/hs37d5.fa --color_consensus -o /data/tusers/boxu/lrft/result/nanopore/human/insertion/${human_sample}_0 > /data/tusers/boxu/lrft/log/${human_sample}_simulate.tldr.log 2>&1





