#! /usr/lib/R/bin/Rscript --vanilla


# Args <- commandArgs()
# TE=Args[7]
# data_path='/home/boxu/Re-analyzing/ChIPseq/chipseq_pipeline/map_to_transposon'
# out_path='/home/boxu/Re-analyzing/ChIPseq/chipseq_pipeline/figure'

# pdf(paste(out_path,'/chip_',TE,'.pdf',sep=''), width = 10, height = 5)


library('ggplot2')
library(VennDiagram)

merge.all <- read.table("/data/tusers/boxu/lrft/result/nanopore/human/figure/all.intersect.txt")

merge.all.tldr <- c( merge.all[merge.all$V1=="tldr", ]$V2)
merge.all.TEMP2 <- c(merge.all[merge.all$V1=="TEMP2", ]$V2)
merge.all.lrft <- c(merge.all[merge.all$V1=="lrft", ]$V2)

venn_data <- list(tldr=merge.all.tldr, TEMP2=merge.all.TEMP2, lrft=merge.all.lrft)


venn.plot <- venn.diagram(venn_data,
  filename = "/data/tusers/boxu/lrft/result/nanopore/human/figure/merge.all.png",##韦恩图的名字
  lty = 1,
  lwd = 0,
  col = "black",  ##圈的颜色
  fill = c("#1262EB", "#ff3535", "#FFCD00"),##对应每个圈的颜色，有几个数据集，就需要有相应数量的颜色
  cat.col = "black",##此处设置每个数据集的名称颜色，也可以使用c（）函数输入三种颜色
  cat.cex = 1.5,
  cat.fontface = "bold",
  margin = 0.07,
  cex = 1.5
)


