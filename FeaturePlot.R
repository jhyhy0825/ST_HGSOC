#! /usr/bin/env Rscript

## /home/jhy/spatial_transcriptomics/220905/code/Giotto/FeaturePlot.R

library(Seurat)

enrichment_data <- read.table('/home/jhy/spatial_transcriptomics/220905/result/seurat/CPS_OV10RTOV4-1/Giotto_GSE180661_GSE211956_bulk_GSE146026/enrichment_result2.txt',header=T,sep=' ')

sp <- readRDS('/home/jhy/spatial_transcriptomics/220905/result/seurat/CPS_OV10RTOV4-1/dimensional_reduction_dim.rds')

Cancer = as.vector(enrichment_data$Ovarian.cancer.cell)
Fibroblast = as.vector(enrichment_data$Fibroblast)
Endothelial = as.vector(enrichment_data$Endothelial.cell)
Tcell = as.vector(enrichment_data$T.cell)
B_Plasma = as.vector(enrichment_data$B.Plasma.cell)
Dendritic = as.vector(enrichment_data$Dendritic.cell)
Myeloid = as.vector(enrichment_data$Myeloid.cell)
Macrophage = as.vector(enrichment_data$Macrophage)

sp <- AddMetaData(sp,metadata=Cancer,col.name="Ovarian_Cancer")
sp <- AddMetaData(sp,metadata=Fibroblast,col.name="Fibroblast")
sp <- AddMetaData(sp,metadata=Endothelial,col.name="Endothelial_cell")
sp <- AddMetaData(sp,metadata=Tcell,col.name="T_cell")
sp <- AddMetaData(sp,metadata=B_Plasma,col.name="B_Plasma_cell")
sp <- AddMetaData(sp,metadata=Dendritic,col.name="Dendritic_cell")
sp <- AddMetaData(sp,metadata=Myeloid,col.name="Myeloid_cell")
sp <- AddMetaData(sp,metadata=Macrophage,col.name="Macrophage")

pdf('/home/jhy/spatial_transcriptomics/220905/result/seurat/CPS_OV10RTOV4-1/Giotto_GSE180661_GSE211956_bulk_GSE146026/CellTypeEnrichment_230704.pdf')
SpatialFeaturePlot(sp,features=c("Ovarian_Cancer","Fibroblast","Endothelial_cell","T_cell","B_Plasma_cell","Dendritic_cell","Myeloid_cell","Macrophage"))
dev.off()

