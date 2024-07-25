#!/home/jhy/anaconda2/envs/scanpy3/bin/Rscript

library(Seurat)
#library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)

args = commandArgs(trailingOnly = TRUE)
inputdir <- args[1]
outputdir <- args[2]
project <- args[3]

## Import_input_data
data <- Load10X_Spatial(data.dir=inputdir, filename="filtered_feature_bc_matrix.h5", assay="Spatial", filter.matrix=TRUE)
data$orig.ident <- project
Idents(data) <- data$orig.ident
data
head(data)

## Because_of_10X_input
data@images[["slice1"]]@coordinates[["tissue"]] <- as.integer(data@images[["slice1"]]@coordinates[["tissue"]])
data@images[["slice1"]]@coordinates[["row"]] <- as.integer(data@images[["slice1"]]@coordinates[["row"]])
data@images[["slice1"]]@coordinates[["col"]] <- as.integer(data@images[["slice1"]]@coordinates[["col"]])
data@images[["slice1"]]@coordinates[["imagerow"]] <- as.integer(data@images[["slice1"]]@coordinates[["imagerow"]])
data@images[["slice1"]]@coordinates[["imagecol"]] <- as.integer(data@images[["slice1"]]@coordinates[["imagecol"]])


## nCount_plot
pdf(paste0(outputdir,'Violin_and_Feature.pdf'),width=10,height=5)
plot1 <- VlnPlot(data, features="nCount_Spatial", pt.size=0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(data, features="nCount_Spatial") + theme(legend.position="right")
wrap_plots(plot1,plot2)
dev.off()


## SCTransform which builds regularized negative binomial models of gene expression in order to account for technical artifacts while preserving biological variance
data <- SCTransform(data, assay="Spatial", verbose = FALSE)
data
head(data)

pdf(paste0(outputdir,'Violin_and_Feature_after_sct.pdf'),width=10,height=5)
plot1 <- VlnPlot(data, features="nCount_Spatial", pt.size=0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(data, features="nCount_Spatial") + theme(legend.position="right")
wrap_plots(plot1,plot2)
dev.off()

pdf(paste0(outputdir,'Feature_TP53.pdf'),width=10,height=5)
p1 <- SpatialFeaturePlot(data, features="TP53")
p2 <- SpatialFeaturePlot(data, features="TP53", pt.size.factor=1)
p3 <- SpatialFeaturePlot(data, features="TP53", alpha=c(0.1,1))
p1+p2+p3
dev.off()

saveRDS(data, file = paste0(outputdir,"/preprocessing.rds"))

