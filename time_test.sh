#!/bin/bash



GENOME="hs37d5"
DATA_PATH="/data/tusers/boxu/lrft/result/nanopore/human"
OUT_PATH="/data/tusers/boxu/lrft/result/nanopore/human/time_test"
ANNO_PATH="/data/tusers/boxu/annotation"
TE_REF=${ANNO_PATH}/${GENOME}/ALUL1SVA.2.fa
SAMPLE="SRR11669560"
CPU=10
TE_temp=${TE_REF##*/}
TE=${TE_temp%.*}

filter=${DATA_PATH}/filter
PREFIX=${OUT_PATH}/${SAMPLE}


[ -f time_test_result.txt ] && rm time_test_result.txt


# start=$(date +%s)
# # map to genome
# minimap2 -a ${ANNO_PATH}/${GENOME}/${GENOME}.mmi ${filter}/${SAMPLE}.clean.fastq > ${PREFIX}.alignment.genome.sam
# samtools sort -@ 8 -O bam -o ${PREFIX}.alignment.genome.bam ${PREFIX}.alignment.genome.sam
# samtools index -@ 8 ${PREFIX}.alignment.genome.bam
# end=$(date +%s)
# time1=$(( $end - $start ))
# echo -e "\n-----------\n" >> time_test_result.txt
# echo "map to genome:" >> time_test_result.txt
# echo ${time1}'s' >> time_test_result.txt



# start=$(date +%s)
# # map to TE
# minimap2 -a ${ANNO_PATH}/${GENOME}/${TE}.mmi ${filter}/${SAMPLE}.clean.fastq > ${PREFIX}.alignment.TE.sam 
# samtools sort -@ ${CPU} -O bam -o ${PREFIX}.alignment.TE.sorted.bam ${PREFIX}.alignment.TE.sam
# samtools index -@ ${CPU} ${PREFIX}.alignment.TE.sorted.bam
# BAM=${PREFIX}.alignment.TE.sorted.bam

# samtools view -bhSF 4 ${BAM} | samtools view -bhSq 0 > ${PREFIX}.mapped.TE.bam
# samtools sort -@ ${CPU} -o ${PREFIX}.mapped.sorted.TE.bam ${PREFIX}.mapped.TE.bam
# samtools index -@ ${CPU} ${PREFIX}.mapped.sorted.TE.bam
# BAM=${PREFIX}.mapped.sorted.TE.bam

# python lrft_clip_reads.py ${BAM} ${PREFIX}.mapped.sorted.TE.clip.fq ${PREFIX}.mapped.sorted.TE.clip.bed ${PREFIX}.clip.reads.pkl

# minimap2 -a ${ANNO_PATH}/${GENOME}/${GENOME}.mmi ${PREFIX}.mapped.sorted.TE.clip.fq > ${PREFIX}.mapped.TE.clip.genome.sam 	
# samtools sort -@ 8 -O bam -o ${PREFIX}.mapped.TE.clip.genome.sorted.bam ${PREFIX}.mapped.TE.clip.genome.sam
# samtools index -@ 8 ${PREFIX}.mapped.TE.clip.genome.sorted.bam
# BAM=${PREFIX}.mapped.TE.clip.genome.sorted.bam

# samtools view -bhSF 4 ${BAM} > ${PREFIX}.mapped.tcg.bam
# samtools sort -@ ${CPU} -o ${PREFIX}.mapped.sorted.tcg.bam ${PREFIX}.mapped.tcg.bam
# samtools index -@ ${CPU} ${PREFIX}.mapped.sorted.tcg.bam

# end=$(date +%s)
# time2=$(( $end - $start ))
# echo -e "\n-----------\n" >> time_test_result.txt
# echo "map to TE first and get tcg.bam" >> time_test_result.txt
# echo ${time2}'s' >> time_test_result.txt


# BAM=${PREFIX}.alignment.genome.bam
# samtools view -bhSF 4 ${BAM} > ${PREFIX}.mapped.genome.bam
# samtools sort -@ ${CPU} -o ${PREFIX}.mapped.genome.sorted.bam ${PREFIX}.mapped.genome.bam
# samtools index -@ ${CPU} ${PREFIX}.mapped.genome.sorted.bam


start=$(date +%s)
tldr -b ${PREFIX}.mapped.sorted.tcg.bam -e /data/tusers/boxu/annotation/hs37d5/ALUL1SVA.2.fa -r /data/tusers/boxu/annotation/hs37d5/hs37d5.fa --color_consensus -o ${PREFIX} > ${PREFIX}.tldr.log 2>&1
end=$(date +%s)
time2=$(( $end - $start ))
echo -e "\n-----------\n" >> time_test_result.txt
echo "tldr" >> time_test_result.txt
echo ${time2}'s' >> time_test_result.txt


echo "DONE"
