#!/home/jhy/anaconda2/envs/scanpy3/bin/Rscript
  
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

args = commandArgs(trailingOnly = TRUE)
inputdir <- args[1]
gene_list <- args[2]

data <- readRDS(file = paste0(inputdir,"/preprocessing.rds"))
gene_list <- read.table(gene_list,sep='\t',header=F)
gene_list <- as.character(gene_list)
gene1 <- gene_list[1:9]
gene2 <- gene_list[10:18]
gene3 <- gene_list[19:27]
gene4 <- gene_list[28:length(gene_list)]

pdf(paste0(inputdir,'Dotplot_221019_1.pdf'))
#SpatialFeaturePlot(data, features=c("TP53","CD36","WFDC2","MUC16","CLDN4","COL1A1","VIM","FN1","POSTN"))
SpatialFeaturePlot(data, features=gene1)
dev.off()

pdf(paste0(inputdir,'Dotplot_221019_2.pdf'))
SpatialFeaturePlot(data, features=gene2)
dev.off()

pdf(paste0(inputdir,'Dotplot_221019_3.pdf'))
SpatialFeaturePlot(data, features=gene3)
dev.off()

pdf(paste0(inputdir,'Dotplot_221019_4.pdf'))
SpatialFeaturePlot(data, features=gene4)
dev.off()


