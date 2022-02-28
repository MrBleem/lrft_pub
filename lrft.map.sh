#!/bin/bash

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
# if_remap=
if_genome_first=1
# lrft_temp.sh -g hs37d5 -i /data/tusers/boxu/lrft/rawdata/nanopore -o /data/tusers/boxu/lrft/result/nanopore -s SRR11669560


TE_temp=${TE_REF##*/}
TE=${TE_temp%.*}
GENOME_REF=${ANNO_PATH}/${GENOME}/${GENOME}.fa
SAMPLES=(${SAMPLES//,/ })

BIN_PATH=$(cd "$(dirname "$0")";pwd)

TE_QC="0"


# GENOME="dm6_clean"
# DATA_PATH="/data/tusers/boxu/lrft/rawdata/gold"
# OUT_PATH="/data/tusers/boxu/lrft/result/gold"
# ANNO_PATH="/data/tusers/boxu/annotation"
# TE_REF=${ANNO_PATH}/${GENOME}/dm6_clean.transposon.2.fa
# GENOME_REF=${ANNO_PATH}/${GENOME}/${GENOME}.fa
# SAMPLES=('ERR4920934' 'ERR4920935' 'ERR4920936')


# checking dic
QC=${OUT_PATH}/QC
nanoplot=${OUT_PATH}/nanoplot
summary=${nanoplot}/summary
filter=${OUT_PATH}/filter
map=${OUT_PATH}/map
insertion=${OUT_PATH}/insertion

[ ! -d ${OUT_PATH} ] && mkdir ${OUT_PATH}
[ ! -d ${QC} ] && mkdir ${QC}
[ ! -d ${nanoplot} ] && mkdir ${nanoplot}
[ ! -d ${summary} ] && mkdir ${summary}
[ ! -d ${filter} ] && mkdir ${filter}
[ ! -d ${map} ] && mkdir ${map}
[ ! -d ${insertion} ] && mkdir ${insertion}

# checking index
if [ ! -f ${ANNO_PATH}/${GENOME}/${GENOME}.mmi ];then
	minimap2 -d ${ANNO_PATH}/${GENOME}/${GENOME}.mmi ${GENOME_REF}
fi
	
if [ ! -f ${ANNO_PATH}/${GENOME}/${TE}.mmi ];then
	minimap2 -d ${ANNO_PATH}/${GENOME}/${TE}.mmi ${TE_REF}
fi



for SAMPLE in ${SAMPLES[*]}
do
# 	echo "### ${SAMPLE}"
	echo -e "\033[40;32;1m### ${SAMPLE}\033[0m"
	PREFIX=${OUT_PATH}/map/${SAMPLE}_${TE_QC}
	
	echo ">>>>>>>>>>>-----step1 : QC-----<<<<<<<<<<"
	[ ! -d ${QC}/${SAMPLE} ] && mkdir ${QC}/${SAMPLE}
	if [ ! -f ${QC}/${SAMPLE}/nanoQC.html ];then
		nanoQC ${DATA_PATH}/${SAMPLE}*.f*q -o ${QC}/${SAMPLE}
		echo "skip"
	else
		echo "already QC"
	fi
 	
 	# step3--nanofilt
 	echo ">>>>>>>>>>>-----step2 : NanoFilt filtering reads-----<<<<<<<<<<"
 	if [ ! -f ${filter}/${SAMPLE}.clean.fastq ];then
 		NanoFilt -q 7 -l 1000 --headcrop 50 --tailcrop 50 ${DATA_PATH}/${SAMPLE}*.f*q  > ${filter}/${SAMPLE}.clean.fastq
 	else
 		echo "already clean"
 	fi
 	

	
	echo ">>>>>>>>>>>-----step3 : mapping-----<<<<<<<<<<"
	if [ ! -f ${PREFIX}.alignment.TE.sorted.bam ];then
		echo "map to transposon"
		minimap2 -a ${ANNO_PATH}/${GENOME}/${TE}.mmi ${filter}/${SAMPLE}.clean.fastq > ${PREFIX}.alignment.TE.sam 
		# minimap2 -a ${ANNO_PATH}/${GENOME}/${TE}.mmi ${DATA_PATH}/${SAMPLE}*.f*q > ${PREFIX}.alignment.TE.sam 
		samtools sort -@ ${CPU} -O bam -o ${PREFIX}.alignment.TE.sorted.bam ${PREFIX}.alignment.TE.sam
		samtools index -@ ${CPU} ${PREFIX}.alignment.TE.sorted.bam
	else
		echo "already map to transposon"
	fi
	BAM=${PREFIX}.alignment.TE.sorted.bam
	
	
	echo ">>>>>>>>>>>-----step4 : filtering--TE -----<<<<<<<<<<"
	if [ ! -f ${PREFIX}.mapped.sorted.TE.bam ];then
		echo "filter aligned reads"
		# samtools view -bhSf 4 ${BAM} > ${PREFIX}.unmapped.bam # 提取未比对到参考序列上的比对结果
		# 提取比对到参考序列上的比对结果
		# 比对到 TE 上的结果可以稍微放宽一点
		# samtools view -bhSF 4 ${BAM} | samtools view -bhSq 60 > ${PREFIX}.mapped.TE.bam
		samtools view -bhSF 4 ${BAM} | samtools view -bhSq ${TE_QC} > ${PREFIX}.mapped.TE.bam
		
		samtools sort -@ ${CPU} -o ${PREFIX}.mapped.sorted.TE.bam ${PREFIX}.mapped.TE.bam
		samtools index -@ ${CPU} ${PREFIX}.mapped.sorted.TE.bam
	else
		echo "already filter"
	fi
	BAM=${PREFIX}.mapped.sorted.TE.bam
	
	
	
	############################### un uniq  ####################################################v##########################
	BAM=${PREFIX}.mapped.sorted.TE.bam
	echo ">>>>>>>>>>>-----step? : extract clip reads -----<<<<<<<<<<"
	if [ ! -f ${PREFIX}.mapped.sorted.TE.clip.fq ];then
		python ${BIN_PATH}/lrft_clip_reads.py ${BAM} ${PREFIX}.mapped.sorted.TE.clip.fq ${PREFIX}.mapped.sorted.TE.clip.bed ${PREFIX}.clip.reads.pkl
		# python ${BIN_PATH}/bed2dic.py ${PREFIX}.mapped.sorted.TE.clip.bed ${PREFIX}.clip.reads.pkl
	else
		echo "already clip"
	fi


 	echo ">>>>>>>>>>>-----step? : map mapped reads to genome -----<<<<<<<<<<"
 	if [ ! -f ${PREFIX}.mapped.TE.clip.genome.sorted.bam ];then
 		minimap2 -a ${ANNO_PATH}/${GENOME}/${GENOME}.mmi ${PREFIX}.mapped.sorted.TE.clip.fq > ${PREFIX}.mapped.TE.clip.genome.sam 	
 		samtools sort -@ 8 -O bam -o ${PREFIX}.mapped.TE.clip.genome.sorted.bam ${PREFIX}.mapped.TE.clip.genome.sam
		samtools index -@ 8 ${PREFIX}.mapped.TE.clip.genome.sorted.bam
	else
		echo "already mapped"
	fi
	BAM=${PREFIX}.mapped.TE.clip.genome.sorted.bam
	
	
	echo ">>>>>>>>>>>-----step? : filtering--genome -----<<<<<<<<<<"
	if [ ! -f ${PREFIX}.mapped.sorted.tcg.bam ];then
		echo "filter aligned reads"
		samtools view -bhSF 4 ${BAM} > ${PREFIX}.mapped.tcg.bam
		samtools sort -@ ${CPU} -o ${PREFIX}.mapped.sorted.tcg.bam ${PREFIX}.mapped.tcg.bam
		samtools index -@ ${CPU} ${PREFIX}.mapped.sorted.tcg.bam
	else
		echo "already filter"
	fi
	BAM=${PREFIX}.mapped.sorted.tcg.bam
	##################################################################################################################################
	




	
	############################### uniq  ####################################################v##########################
	
	if [ ${if_uniq} ];then
		BAM=${PREFIX}.mapped.sorted.TE.bam
		echo ">>>>>>>>>>>-----step-uniq : filtering--TE -----<<<<<<<<<<"
		if [ ! -f ${PREFIX}.mapped.sorted.TE.uniq.bam ];then
			echo "filter aligned reads"
# 			samtools view -h -q 60 -F 4 ${BAM} | grep -v XA:Z |grep -v SA:Z | samtools view -b - > ${PREFIX}.mapped.TE.uniq.bam 
			samtools view -bSh -q 60 -F 4 ${BAM} > ${PREFIX}.mapped.TE.uniq.bam 

			samtools sort -@ ${CPU} -o ${PREFIX}.mapped.sorted.TE.uniq.bam ${PREFIX}.mapped.TE.uniq.bam 
			samtools index ${PREFIX}.mapped.sorted.TE.uniq.bam 
		else
			echo "already filter"
		fi
		BAM=${PREFIX}.mapped.sorted.TE.uniq.bam 
	
		echo ">>>>>>>>>>>-----step-uniq : extract clip reads -----<<<<<<<<<<"
		if [ ! -f ${PREFIX}.mapped.sorted.TE.uniq.clip.fq ];then
			python ${BIN_PATH}/lrft_clip_reads.py ${BAM} ${PREFIX}.mapped.sorted.TE.uniq.clip.fq ${PREFIX}.mapped.sorted.TE.uniq.clip.bed ${PREFIX}.uniq.clip.reads.pkl
			# python ${BIN_PATH}/bed2dic.py ${PREFIX}.mapped.sorted.TE.uniq.clip.bed ${PREFIX}.uniq.clip.reads.pkl
		else
			echo "already clip"
		fi
	
		echo ">>>>>>>>>>>-----step-uniq : map mapped reads to genome -----<<<<<<<<<<"
 		if [ ! -f ${PREFIX}.mapped.TE.uniq.clip.genome.sorted.bam ];then
 			minimap2 -a ${ANNO_PATH}/${GENOME}/${GENOME}.mmi ${PREFIX}.mapped.sorted.TE.uniq.clip.fq > ${PREFIX}.mapped.TE.uniq.clip.genome.sam 	
 			samtools sort -@ 8 -O bam -o ${PREFIX}.mapped.TE.uniq.clip.genome.sorted.bam ${PREFIX}.mapped.TE.uniq.clip.genome.sam
			samtools index -@ 8 ${PREFIX}.mapped.TE.uniq.clip.genome.sorted.bam
		else
			echo "already mapped"
		fi
		BAM=${PREFIX}.mapped.TE.uniq.clip.genome.sorted.bam
	
		echo ">>>>>>>>>>>-----step-uniq : filtering--genome -----<<<<<<<<<<"
		if [ ! -f ${PREFIX}.mapped.sorted.tcg.uniq.bam ];then
			echo "filter aligned reads"
			samtools view -bhSF 4 -q 60 ${BAM} > ${PREFIX}.mapped.tcg.uniq.bam
			samtools sort -@ ${CPU} -o ${PREFIX}.mapped.sorted.tcg.uniq.bam ${PREFIX}.mapped.tcg.uniq.bam
			samtools index -@ ${CPU} ${PREFIX}.mapped.sorted.tcg.uniq.bam
		else
			echo "already filter"
		fi
		BAM=${PREFIX}.mapped.sorted.tcg.uniq.bam
	fi
	##################################################################################################################################
	
	
	
	
	
	
	############################### remap  ####################################################v##########################
	if [ ${if_remap} ];then
		echo ">>>>>>>>>>>-----step? : remap mapped reads to TE -----<<<<<<<<<<"
 		if [ ! -f ${PREFIX}.realignment.TE.clip.bam ];then
 			minimap2 -a ${ANNO_PATH}/${GENOME}/${TE}.mmi ${PREFIX}.mapped.sorted.TE.clip.fq > ${PREFIX}.realignment.TE.clip.sam 	
 			samtools sort -@ 8 -O bam -o ${PREFIX}.realignment.TE.clip.bam ${PREFIX}.realignment.TE.clip.sam
			samtools index -@ 8 ${PREFIX}.realignment.TE.clip.bam
		else
			echo "already mapped"
		fi
		BAM=${PREFIX}.realignment.TE.clip.bam
	
		echo ">>>>>>>>>>>-----step? : filtering--re_TE -----<<<<<<<<<<"
		if [ ! -f ${PREFIX}.remapped.TE.clip.sorted.bam ];then
			echo "filter aligned reads"
			samtools view -bhSF 4 ${BAM} > ${PREFIX}.remapped.TE.clip.bam
			samtools sort -@ ${CPU} -o ${PREFIX}.remapped.TE.clip.sorted.bam ${PREFIX}.remapped.TE.clip.bam
			samtools index -@ ${CPU} ${PREFIX}.remapped.TE.clip.sorted.bam
		else
			echo "already filter"
		fi
		BAM=${PREFIX}.remapped.TE.clip.sorted.bam
	
	
		echo ">>>>>>>>>>>-----step? : re extract clip reads 2-----<<<<<<<<<<"
		if [ ! -f ${PREFIX}.remapped.TE.clip.fq ];then
			python ${BIN_PATH}/lrft_clip_reads.py ${BAM} ${PREFIX}.remapped.TE.clip.fq ${PREFIX}.remapped.TE.clip.bed ${PREFIX}.remapped.clip.reads.pkl
			# python ${BIN_PATH}/bed2dic.py ${PREFIX}.remapped.TE.clip.bed ${PREFIX}.remapped.clip.reads.pkl
		else
			echo "already clip"
		fi
	
		echo ">>>>>>>>>>>-----step? : map remapped reads to genome -----<<<<<<<<<<"
 		if [ ! -f ${PREFIX}.remapped.TE.clip.genome.sorted.bam ];then
 			minimap2 -a ${ANNO_PATH}/${GENOME}/${GENOME}.mmi ${PREFIX}.remapped.TE.clip.fq > ${PREFIX}.remapped.TE.clip.genome.sam 	
 			samtools sort -@ 8 -O bam -o ${PREFIX}.remapped.TE.clip.genome.sorted.bam ${PREFIX}.remapped.TE.clip.genome.sam
			samtools index -@ 8 ${PREFIX}.remapped.TE.clip.genome.sorted.bam
		else
			echo "already mapped"
		fi
		BAM=${PREFIX}.remapped.TE.clip.genome.sorted.bam
	
	
		echo ">>>>>>>>>>>-----step? : refiltering--genome -----<<<<<<<<<<"
		if [ ! -f ${PREFIX}.remapped.sorted.tcg.bam ];then
			echo "filter aligned reads"
			samtools view -bhSF 4 ${BAM} > ${PREFIX}.remapped.tcg.bam
			samtools sort -@ ${CPU} -o ${PREFIX}.remapped.sorted.tcg.bam ${PREFIX}.remapped.tcg.bam
			samtools index -@ ${CPU} ${PREFIX}.remapped.sorted.tcg.bam
		else
			echo "already filter"
		fi
		BAM=${PREFIX}.remapped.sorted.tcg.bam
	fi
	##################################################################################################################################
	
	
	
	
	############################### map to genome first  ####################################################v##########################
	if [ ${if_genome_first} ];then
		echo ">>>>>>>>>>>-----step? : map mapped reads to genome -----<<<<<<<<<<"
 		if [ ! -f ${PREFIX}.alignment.genome.bam ];then
 			minimap2 -a ${ANNO_PATH}/${GENOME}/${GENOME}.mmi ${filter}/${SAMPLE}.clean.fastq > ${PREFIX}.alignment.genome.sam
 			samtools sort -@ 8 -O bam -o ${PREFIX}.alignment.genome.bam ${PREFIX}.alignment.genome.sam
			samtools index -@ 8 ${PREFIX}.alignment.genome.bam
		else
			echo "already mapped"
		fi
		BAM=${PREFIX}.alignment.genome.bam
	
	
		echo ">>>>>>>>>>>-----step? : filtering--genome -----<<<<<<<<<<"
		if [ ! -f ${PREFIX}.mapped.genome.sorted.bam ];then
			echo "filter aligned reads"
			samtools view -bhSF 4 ${BAM} > ${PREFIX}.mapped.genome.bam
			samtools sort -@ ${CPU} -o ${PREFIX}.mapped.genome.sorted.bam ${PREFIX}.mapped.genome.bam
			samtools index -@ ${CPU} ${PREFIX}.mapped.genome.sorted.bam
		else
			echo "already filter"
		fi
		BAM=${PREFIX}.mapped.genome.sorted.bam
	
		echo ">>>>>>>>>>>-----step? : test tldr -----<<<<<<<<<<"
		if [ ! -f ${insertion}/${SAMPLE}.table.txt ];then
			tldr -b ${BAM} -e ${TE_REF} -r ${GENOME_REF} --color_consensus -o ${insertion}/${SAMPLE}
		else
			echo "alrady detect insertion"
		fi
	fi
	##################################################################################################################################
	
	
done

echo "find insertion"

echo -e "\033[40;32;1m>?< DONE >_<\033[0m"

