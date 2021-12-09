#!/bin/bash

PROJECT_PATH="/data/tusers/boxu/lrft/result/nanopore/human"
BED_PATH="${PROJECT_PATH}/insertion/tmp"

for i in ` ls ${BED_PATH} `
do
    LC_COLLATE=C sort -k1,1 -k2,2n ${i} | bedtools merge -i - -c 6 -o collapse -delim ";" | head 
done




