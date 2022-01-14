#!/bin/bash

version='v6'
version='q0'
DATA_PATH="/data/tusers/boxu/lrft/result/nanopore/human/figure"
RESULT_PATH="/data/tusers/boxu/lrft/result/nanopore/human/figure/${version}"
insertion_PATH="/data/tusers/boxu/lrft/result/nanopore/human/insertion"
# echo "" > ${DATA_PATH}/chr1.intersect.txt 
# echo "" > ${DATA_PATH}/all.intersect.txt 


[ -f ${RESULT_PATH}/all.intersect.txt ] && rm ${RESULT_PATH}/all.intersect.txt

# awk 'BEGIN{FS=OFS="\t"}{split($5,a,"|");split(a[3],b,"_");if($3-$2>=10){type="germline"}else{type="novel"}; print $1,$2,$3,b[1],$4,"lrft:"type}' ${insertion_PATH}/SRR11669560.${version}.lrft.table.txt | sed 1d  > ${RESULT_PATH}/lrft.all.bed
awk 'BEGIN{FS=OFS="\t"}{split($5,a,"|");split(a[3],b,"_"); print $1,$2,$3,b[1],$4,"lrft:"$9}' ${insertion_PATH}/SRR11669560_0.lrft.table.txt | sed 1d > ${RESULT_PATH}/lrft.all.bed

awk 'BEGIN{FS=OFS="\t"}{split($4,a,":");if($5 >= 0.1){print $1,$2,$3,$6, a[1], "TEMP2:"$7}}'  ${DATA_PATH}/TEMP2.bed  | sed 1d | grep 1p1 > ${RESULT_PATH}/TEMP2.all.bed
awk 'BEGIN{FS=OFS="\t"}{print $2,$3,$4,$5,$6,"tldr"}' ${DATA_PATH}/tldr.bed   | sed 1d  > ${RESULT_PATH}/tldr.all.bed


cat ${RESULT_PATH}/lrft.all.bed ${RESULT_PATH}/TEMP2.all.bed ${RESULT_PATH}/tldr.all.bed > ${RESULT_PATH}/ttl.all.bed
LC_COLLATE=C sort -k1,1 -k2,2n ${RESULT_PATH}/ttl.all.bed  | bedtools merge -i - -c 4,5,6 -d 20 -o collapse -delim ";"  > ${RESULT_PATH}/merge.all.bed


for meth in tldr TEMP2 lrft
do
    # grep ${meth} -n ${DATA_PATH}/merge.bed | awk -v me=${meth} 'BEGIN{FS=OFS="\t"}{split($1,a,":");print me,a[1]}' >> ${RESULT_PATH}/chr1.intersect.txt 
    grep ${meth} -n ${RESULT_PATH}/merge.all.bed | awk -v me=${meth} 'BEGIN{FS=OFS="\t"}{split($1,a,":");print me,a[1]}' >> ${RESULT_PATH}/all.intersect.txt 
done


Rscript ./intersection.R

