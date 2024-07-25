#! /usr/bin/env Rscript

## Rscript /home/jhy/spatial_transcriptomics/220905/code/image/edgeR.R /home/jhy/spatial_transcriptomics/220905/result/seurat/combine/ /home/jhy/spatial_transcriptomics/220905/result/seurat/combine//image/bilaterality_ovary/ /home/jhy/spatial_transcriptomics/220905/result/seurat/combine//image/bilaterality_ovary/bilaterality_ovary.csv bilaterality_ovary 

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

args = commandArgs(trailingOnly = TRUE)
seuratdir <- args[1]
inputdir <- args[2]
samples <- args[3]
image <- args[4]

data10 <- readRDS(file = paste0(seuratdir,"/dimensional_reduction_dim.rds"))

data10
head(data10)

samples <- read.table(samples,header=T,sep=',')
samples$ids[samples$ids == "CPS_OV1"] <- "CPS_OV1RtOV3"
samples$ids[samples$ids == "CPS_OV5"] <- "CPS_OV5LtOV4"
samples$ids[samples$ids == "CPS_OV10"] <- "CPS_OV10RTOV4-1"
samples$ids[samples$ids == "CPS_OV19"] <- "CPS_OV19_LtOV1"
samples$ids[samples$ids == "CPS_OV20"] <- "CPS_OV20RtOV4"
samples$ids[samples$ids == "CPS_OV24"] <- "CPS_OV24RTOV4"
samples$ids[samples$ids == "CPS_OV34"] <- "CPS_OV34RtOV1"
samples$ids[samples$ids == "CPS_OV71"] <- "CPS_OV71_1"
Case <- subset(samples,samples$type == "Case")
Case
Case$ids
Case <- subset(data10, subset = sample %in% Case$ids)
Case
Case$group <- 'Case'
head(Case)
Control <- subset(samples,samples$type == "Control")
Control
Control <- subset(data10, subset = sample %in% Control$ids)
Control
Control$group <- 'Control'
head(Control)

Combined <- merge(Case, Control)
Combined
head(Combined)
tail(Combined)

Idents(Combined) <- Combined$group

head(Combined$SCT@counts)



