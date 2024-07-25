#! /usr/bin/env Rscript
##/home/jhy/anaconda2/envs/SpatialInferCNV/bin/Rscript

## /home/jhy/spatial_transcriptomics/220905/code/SpatialInferCNV/SpatialInferCNV.R /home/jhy/spatial_transcriptomics/220905/result/spaceranger/CPS_OV1RtOV3/CPS_OV1RtOV3/outs/filtered_feature_bc_matrix.h5 /home/jhy/spatial_transcriptomics/220905/result/seurat/CPS_OV1RtOV3/dimensional_reduction_dim.rds CPS_OV1RtOV3 /home/jhy/test/230104_scCNVenv/CPS_OV1RtOV3/

## /home/jhy/spatial_transcriptomics/220905/code/SpatialInferCNV/SpatialInferCNV.R /home/jhy/spatial_transcriptomics/220905/result/spaceranger/CPS_OV34RtOV1/CPS_OV34RtOV1/outs/filtered_feature_bc_matrix.h5 /home/jhy/spatial_transcriptomics/220905/result/seurat/CPS_OV34RtOV1/dimensional_reduction_dim.rds CPS_OV34RtOV1 /home/jhy/test/230104_scCNVenv/CPS_OV34RtOV1/

## /home/jhy/spatial_transcriptomics/220905/code/SpatialInferCNV/SpatialInferCNV.R /home/jhy/spatial_transcriptomics/220905/result/spaceranger/CPS_OV34RtOV1/CPS_OV34RtOV1/outs/filtered_feature_bc_matrix.h5 /home/jhy/spatial_transcriptomics/220905/result/seurat/CPS_OV34RtOV1/dimensional_reduction_dim.rds CPS_OV34RtOV1 /home/jhy/spatial_transcriptomics/220905/result/seurat/CPS_OV34RtOV1/SpatialInferCNV/

library(Seurat)
library(SpatialInferCNV)
library(tidyverse)

args = commandArgs(trailingOnly = TRUE)
spaceranger <- args[1]
seurat <- args[2]
sample <- args[3]
outputdir <- args[4]

count <- ImportCountData(sample, spaceranger)

data <- readRDS(seurat)
cluster <- data.frame(data$seurat_clusters)
cluster$Barcode <- rownames(cluster)
colnames(cluster) <- c("Histology","Barcode")
cluster <- cluster[,c(2,1)]
write.csv(cluster,paste0(outputdir,"Histology_",sample,".csv"),row.names=F,quote=F)
cluster <- paste0(outputdir,"Histology_",sample,".csv")
cluster <- ImportHistologicalAnnotations(sample, cluster)

#########
#Because an error, I used the original code.
#joined_count <- Userguide_10xBreast_Joined_Counts <- MergingCountAndAnnotationData('CPS_OV1RtOV3',count,histology)
SectionName <- sample
InputAnnotationFile <- cluster
InputCountFile <- count
formerge <- select(InputAnnotationFile, -Histology)
MergedAnnotationsandCounts <- inner_join(formerge, InputCountFile)
MergedAnnotationsandCounts <- remove_rownames(MergedAnnotationsandCounts)
MergedAnnotationsandCounts <- column_to_rownames(MergedAnnotationsandCounts, "Barcode")
MergedAnnotationsandCounts$Total <- rowSums(MergedAnnotationsandCounts)
MergedAnnotationsandCounts <- MergedAnnotationsandCounts %>% filter(Total >= 500)
MergedAnnotationsandCounts <- select(MergedAnnotationsandCounts, -Total)
MergedAnnotationsandCounts <- as.data.frame(t(MergedAnnotationsandCounts))
MergedAnnotationsandCounts <- MergedAnnotationsandCounts[,colSums(is.na(MergedAnnotationsandCounts))<nrow(MergedAnnotationsandCounts)]
MergedAnnotationsandCounts <- tibble::rownames_to_column(MergedAnnotationsandCounts, "Genes")
joined_count <- MergedAnnotationsandCounts
#########
head(joined_count)
joined_count <- joined_count  %>% column_to_rownames(var = "Genes")
head(joined_count)
dim(joined_count)


FinalAnnotationsForExport <- FinalAnnotations(cluster, joined_count)
head(FinalAnnotationsForExport)
dim(FinalAnnotationsForExport)

write.table(joined_count,paste0(outputdir,"Joined_count.tsv"),sep='\t')
write.table(FinalAnnotationsForExport,paste0(outputdir,"FinalAnnotationsForExport.tsv"),sep='\t',quote = FALSE,col.names = FALSE,row.names = FALSE)

infCNV <- infercnv::CreateInfercnvObject(raw_counts_matrix=paste0(outputdir,"Joined_count.tsv"), gene_order_file=paste0("/home/jhy/test/230104_scCNVenv/gene_pos.txt"), annotations_file=paste0(outputdir,"FinalAnnotationsForExport.tsv"), delim='\t', ref_group_names=NULL)

infCNV = infercnv::run(infCNV,cutoff=0.1,out_dir=paste0(outputdir,"/InferCNV_outputs/"),cluster_by_groups=FALSE,HMM = FALSE,denoise=TRUE)
#infCNV = infercnv::run(infCNV,cutoff=0.1,out_dir=paste0(outputdir,"/InferCNV_outputs/"),cluster_by_groups=FALSE,HMM=TRUE,denoise=TRUE)


