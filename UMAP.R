#! /usr/bin/env Rscript
  
## /home/jhy/spatial_transcriptomics/220905/code/SpatialInferCNV/UMAP.R /home/jhy/spatial_transcriptomics/220905/result/seurat/CPS_OV1RtOV3/ CPS_OV1RtOV3

library(Seurat)
library(dplyr)

args = commandArgs(trailingOnly = TRUE)
sample = args[2]

pbmc <- readRDS(paste0(args[1],'/dimensional_reduction_dim.rds'))
#data <- read.table(paste0(args[1],'/SpatialInferCNV_231213/InferCNV_outputs_k3/infercnv.observation_groupings.txt'),sep=' ',header=T)
data <- read.table(paste0(args[1],'/SpatialInferCNV_231213/InferCNV_outputs_k5/infercnv.observation_groupings.txt'),sep=' ',header=T)
#data <- read.table(paste0(args[1],'/SpatialInferCNV_231122/InferCNV_outputs_test3/infercnv.observation_groupings.txt'),sep=' ',header=T)
dim(pbmc)
dim(data)

#rownames(data) <- gsub("^CPS_OV34RtOV1_","",rownames(data))
#rownames(data) <- gsub("^CPS_OV10RTOV4-1_","",rownames(data))
rownames(data) <- gsub(paste0("^",sample,"_"),"",rownames(data))
rownames(data) <- sub(".1$","-1",rownames(data))
p <- colnames(pbmc)
test <- data
test <- test[,'Dendrogram.Group', drop=FALSE]
test <- test %>% slice(match(p,rownames(test)))
missing_rows <- setdiff(p,rownames(test))
missing_data <- data.frame(value = rep(NA, length(missing_rows)), row.names = missing_rows)
rownames(missing_data) <- missing_rows
colnames(missing_data) <- c('Dendrogram.Group')
test <- rbind(test,missing_data)
test <- test%>% slice(match(p, rownames(test)))
pbmc$cnv <- test$Dendrogram.Group
Idents(pbmc) <- 'cnv'
head(pbmc)

#pdf(paste0(args[1],'/SpatialInferCNV_231213/InferCNV_outputs_k3/SpatialDimplot_k3_noborder_240223.pdf'))
pdf(paste0(args[1],'/SpatialInferCNV_231213/InferCNV_outputs_k3/SpatialDimplot_k5_noborder_240223.pdf'))
#pdf(paste0(args[1],'/SpatialInferCNV_231122/InferCNV_outputs_test3/SpatialDimplot_test3.pdf'))
SpatialDimPlot(pbmc,stroke=NA)
dev.off()

