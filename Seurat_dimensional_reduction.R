#!/home/jhy/anaconda2/envs/scanpy3/bin/Rscript
  
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

args = commandArgs(trailingOnly = TRUE)
inputdir <- args[1]
dim <- str(args[2])

data <- readRDS(file = paste0(inputdir,"/preprocessing.rds"))

data <- RunPCA(data, assay = "SCT", verbose = FALSE)


### Dim30 ###
data <- FindNeighbors(data, reduction = "pca", dims = 1:24)
data <- FindClusters(data, verbose = FALSE)
data <- RunUMAP(data, reduction = "pca", dims = 1:24)
pdf(paste0(inputdir,"Dimplot_and_SpatialDimplot_dim",dim,".pdf"),width=10,height=5)
p1 <- DimPlot(data, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(data, label = TRUE, label.size = 3)
p1 + p2
dev.off()

### Dim20 ###
#data20 <- FindNeighbors(data, reduction = "pca", dims = 1:20)
#data20 <- FindClusters(data20, verbose = FALSE)
#data20 <- RunUMAP(data20, reduction = "pca", dims = 1:20)
#pdf(paste0(inputdir,"Dimplot_and_SpatialDimplot_dim20.pdf"),width=10,height=5)
#p1 <- DimPlot(data20, reduction = "umap", label = TRUE)
#p2 <- SpatialDimPlot(data20, label = TRUE, label.size = 3)
#p1 + p2
#dev.off()

### Dim10 ###
#data10 <- FindNeighbors(data, reduction = "pca", dims = 1:10)
#data10 <- FindClusters(data10, verbose = FALSE)
#data10 <- RunUMAP(data10, reduction = "pca", dims = 1:10)
#pdf(paste0(inputdir,"Dimplot_and_SpatialDimplot_dim10.pdf"),width=10,height=5)
#p1 <- DimPlot(data10, reduction = "umap", label = TRUE)
#p2 <- SpatialDimPlot(data10, label = TRUE, label.size = 3)
#p1 + p2
#dev.off()


### Identification of Spatially Variable Features ###
de_markers <- FindAllMarkers(data, only.pos=TRUE, min.pct = 0.25, logfc.threshold = 0.25)
de_markers <- subset(de_markers, subset = p_val_adj <= 0.05)
top100 <- de_markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)
write.table(data.frame(top100), paste0(inputdir,'top100_markers_log2FC0.25_adjP0.05_dim',dim,'.txt'))

top10 <- de_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top10
pdf(paste0(inputdir,"top10_genes_heatmap_dim",dim,".pdf"), width=10, height=10)
DoHeatmap(data, features = top10$gene)
dev.off()


#de_markers20 <- FindAllMarkers(data20, only.pos=TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#de_markers20 <- subset(de_markers20, subset = p_val_adj <= 0.05)
#top100 <- de_markers20 %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)
#write.table(data.frame(top100), paste0(inputdir,'top100_markers_log2FC0.25_adjP0.05_dim20.txt'))

#top10 <- de_markers20 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
#top10
#pdf(paste0(inputdir,"top10_genes_heatmap_dim20.pdf"), width=10, height=10)
#DoHeatmap(data20, features = top10$gene)
#dev.off()


#de_markers10 <- FindAllMarkers(data10, only.pos=TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#de_markers10 <- subset(de_markers10, subset = p_val_adj <= 0.05)
#top100 <- de_markers10 %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)
#write.table(data.frame(top100), paste0(inputdir,'top100_markers_log2FC0.25_adjP0.05_dim10.txt'))

#top10 <- de_markers10 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
#top10
#pdf(paste0(inputdir,"top10_genes_heatmap_dim10.pdf"), width=10, height=10)
#DoHeatmap(data10, features = top10$gene)
#dev.off()


markers <- FindAllMarkers(data, only.pos=TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers <- subset(markers, subset = p_val_adj <= 0.05)
write.table(data.frame(markers), paste0(inputdir,'markers_log2FC0.25_adjP0.05_dim',dim,'.txt'))


pdf(file=paste0(inputdir,"FeaturePlot_dim",dim,".pdf"), height=6, width=8)
FeaturePlot(data, features=c("TP53","CD36","WFDC2","MUC16","CLDN4","COL1A1","VIM","FN1","POSTN"))
dev.off()

pdf(file=paste0(inputdir,"SpatialFeaturePlot_dim",dim,".pdf"), height=6, width=8)
SpatialFeaturePlot(data, features=c("TP53","CD36","WFDC2","MUC16","CLDN4","COL1A1","VIM","FN1","POSTN"))
dev.off()


new <- subset(markers,markers$gene == 'TP53' | markers$gene == 'CD36' | markers$gene == 'WFDC2' | markers$gene == 'MUC16' | markers$gene == 'CLDN4' | markers$gene == 'COL1A1' | markers$gene == 'VIM' | markers$gene == 'FN1' | markers$gene == 'POSTN')
write.table(data.frame(new), paste0(inputdir,'tumor_stroma_markers_dim',dim,'.txt'))


saveRDS(data, file = paste0(inputdir,"/dimensional_reduction_dim",dim,".rds"))
#saveRDS(data20, file = paste0(inputdir,"/dimensional_reduction_dim20.rds"))
#saveRDS(data10, file = paste0(inputdir,"/dimensional_reduction_dim10.rds"))



