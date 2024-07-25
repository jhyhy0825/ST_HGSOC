#!/usr/bin/env Rscript

## /home/jhy/spatial_transcriptomics/TCIA/OV/TCGA-OV/code/heatmap.R /home/jhy/spatial_transcriptomics/TCIA/OV/TCGA-OV/DEG/binary1_age/

library(made4)
library(ggplot2)
library(dplyr)

args = commandArgs(trailingOnly = TRUE)
inputdir <- args[1]

data <- read.table(paste0(inputdir,"DEG_q0.1_logFC2.txt"), sep = '\t', header = TRUE)
data <- data[,-1]
data <- data %>% select(gene_name, everything())
data <- data %>% select(1:41)
head(data)
rownames(data) <- data[,1]
data <- data[,-1]
head(data)

pdf(paste0(inputdir,"heatplot_q0.1_logFC2.pdf"),width=14,height=16)
heatplot(data, dend =c("both", "row", "column", "none"), margins=c(10,20),cols.default = TRUE, lowcol = "green", highcol = "red", scale="none",  classvec=NULL, classvecCol=NULL, classvec2=NULL, distfun=NULL, returnSampleTree=FALSE,method="ave", dualScale=TRUE, zlim=c(-3,3),scaleKey=TRUE)
dev.off()

data <- read.table(paste0(inputdir,"DEG_q0.05_log2FC1.txt"), sep = '\t', header = TRUE)
data <- data[,-1]
data <- data %>% select(gene_name, everything())
data <- data %>% select(1:41)
head(data)
rownames(data) <- data[,1]
data <- data[,-1]
head(data)

pdf(paste0(inputdir,"heatplot_q0.05_log2FC1.pdf"),width=12,height=25)
heatplot(data, dend = c("both", "row", "column", "none"), margins=c(15,20),cols.default = TRUE, lowcol = "green", highcol = "red", scale="none",  classvec=NULL, classvecCol=NULL, classvec2=NULL, distfun=NULL, returnSampleTree=FALSE,method="ave", dualScale=TRUE, zlim=c(-3,3),scaleKey=TRUE)
dev.off()


