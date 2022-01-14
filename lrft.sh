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
BIN_PATH="/data/tusers/boxu/lrft/scripts"
DATA_PATH="data/tusers/boxu/lrft/rawdata/nanopore"
OUT_PATH="/data/tusers/boxu/lrft/result/nanopore"
ANNO_PATH="/data/tusers/boxu/annotation"
CPU=10
# Koala_sample="Koala2_20191231_Adelaide_rep1"
# human_sample="SRR11669560"

human_sample="SRR11669560"




# map
${BIN_PATH}/lrft.map.sh -g hs37d5 -i ${DATA_PATH}/human -o ${OUT_PATH}/human -n ${ANNO_PATH} -t ${ANNO_PATH}/hs37d5/ALUL1SVA.2.fa -s ${human_sample} -c ${CPU}

# cluster reads
# python ${BIN_PATH}/lrft_insertion.py ${OUT_PATH}/human ${human_sample}.mapped.sorted.tcg.uniq.bam ${human_sample}.uniq 
# python ${BIN_PATH}/lrft.cluster.py

# # merge insertion found by reads
# ${BIN_PATH}/lrft.merge.sh

# # get consensus sequence
# python ${BIN_PATH}/lrft.consensus.py

python ${BIN_PATH}/lrft_insertion.py /data/tusers/boxu/lrft/result/nanopore/human SRR11669560_0.mapped.sorted.tcg.bam /data/tusers/boxu/annotation/hs37d5/hs37d5.fa SRR11669560_0

