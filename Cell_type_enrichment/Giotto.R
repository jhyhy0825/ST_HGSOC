#! /usr/bin/env Rscript

## /home/jhy/spatial_transcriptomics/220905/code/Giotto/Giotto4.R /home/jhy/single_cell_lung_mouse/HumanValidation/rawdata/GSE180661/Adnexa_count_matrix.txt.gz /home/jhy/single_cell_lung_mouse/HumanValidation/rawdata/GSE180661/Adnexa_cell_type.txt /home/jhy/spatial_transcriptomics/220905/result/seurat/CPS_OV10RTOV4-1/Giotto_GSE180661_matrix/ /home/jhy/spatial_transcriptomics/220905/result/seurat/CPS_OV10RTOV4-1/

library(Seurat)
library(Giotto)

args = commandArgs(trailingOnly = TRUE)

GEX = args[1]
ann = args[2]
outputdir = args[3]
rawdir = args[4]


### SPdata ###
data <- readRDS(file = paste0(rawdir,"/dimensional_reduction_dim.rds"))
count <- data@assays$Spatial@counts
count <- data.matrix(count)
head(count)
dim(count)
location <- data@images$slice1@coordinates
location <- location[,4:5]
location <- data.matrix(location)
head(location)
dim(location)

my_giotto_object <- createGiottoObject(raw_exprs=count,spatial_locs=location)

seurat_clusters <- data$seurat_clusters
my_giotto_object <- addCellMetadata(my_giotto_object,seurat_clusters)
colnames(pDataDT(my_giotto_object))
my_giotto_object <- normalizeGiotto(gobject = my_giotto_object)


### SCdata ###
#single_cell_DT = read.table(GEX,header=T,sep='\t')
library(data.table)
single_cell_DT <- fread(GEX, sep='\t', header=T)
dim(single_cell_DT)
head(colnames(single_cell_DT))
head(rownames(single_cell_DT))
#head(single_cell_DT)
single_cell_matrix = Giotto:::dt_to_matrix(single_cell_DT)
single_cell_matrix[1:4,1:4]
dim(single_cell_matrix)

cell_annotations = read.table(file = ann)
head(cell_annotations)
dim(cell_annotations)
cell_annotations = as.vector(cell_annotations$V1)
head(cell_annotations)
tail(cell_annotations)
length(cell_annotations)


### test ###
#ex <- c("Other","T.cell","Endothelial.cell","B.cell","Mast.cell","Dendritic.cell","Myeloid.cell","Plasma.cell")
#idx <- which(unlist(cell_annotations) %in% ex)
#cell_annotations <- cell_annotations[!cell_annotations %in% ex]
#unique(cell_annotations)
#length(cell_annotations)
#single_cell_matrix <- single_cell_matrix[,-idx]
#dim(single_cell_matrix)


### EnrichTest ###
rank_matrix = makeSignMatrixRank(sc_matrix = single_cell_matrix, sc_cluster_ids = cell_annotations)
my_giotto_object = runRankEnrich(gobject = my_giotto_object, sign_matrix = rank_matrix)
colnames(pDataDT(my_giotto_object))


### CellType ###
cell_types_subset = unique(cell_annotations)
cell_types_subset


### Score ###
Giotto <- my_giotto_object
enrichment_results = Giotto@spatial_enrichment$rank
enrich_cell_types = colnames(enrichment_results)
enrich_cell_types = enrich_cell_types[enrich_cell_types != 'cell_ID']
head(enrichment_results)
enrichment_results2 = as.data.frame(enrichment_results)
rownames(enrichment_results2) = enrichment_results2[,1]
enrichment_results2 = enrichment_results2[,-1]
head(enrichment_results2)
write.table(enrichment_results2,paste0(outputdir,'enrichment_score.txt'))

data1 <- data$seurat_clusters
data2 <- as.data.frame(data1)
new <- merge(enrichment_results2,data2,by='row.names')
head(new)
write.table(new,paste0(outputdir,'enrichment_score_with_seurat_clusters.txt'))


### SpatPlot ###
pdf(paste0(outputdir,'spatplot.pdf'), width=8, height=8)
spatCellPlot(gobject = my_giotto_object, spat_enr_names = 'rank', cell_annotation_values = cell_types_subset, cow_n_col = 3, coord_fix_ratio = NULL, point_size = 0.05, point_shape='no_border', cell_color_gradient = c("blue", "white", "red"))
dev.off()

pdf(paste0(outputdir,'spatplot2.pdf'), width=8, height=8)
spatCellPlot(gobject = my_giotto_object, spat_enr_names = 'rank', cell_annotation_values = cell_types_subset, cow_n_col = 3, coord_fix_ratio = NULL, point_size = 0.05, point_shape='no_border', cell_color_gradient = c('gray90', 'gray90', 'red'))
dev.off()


### SpatPlot_using_Seurat ###
sp <- data
enrichment_data <- new

Cancer = as.vector(enrichment_data$Ovarian.cancer.cell)
Fibroblast = as.vector(enrichment_data$Fibroblast)
Endothelial = as.vector(enrichment_data$Endothelial.cell)
Tcell = as.vector(enrichment_data$T.cell)
Bcell = as.vector(enrichment_data$B.cell)
Dendritic = as.vector(enrichment_data$Dendritic.cell)
Myeloid = as.vector(enrichment_data$Myeloid.cell)
Mast = as.vector(enrichment_data$Mast.cell)
Plasma = as.vector(enrichment_data$Plasma.cell)
Other = as.vector(enrichment_data$Other)

sp <- AddMetaData(sp,metadata=Cancer,col.name="Ovarian_Cancer")
sp <- AddMetaData(sp,metadata=Fibroblast,col.name="Fibroblast")
sp <- AddMetaData(sp,metadata=Endothelial,col.name="Endothelial_cell")
sp <- AddMetaData(sp,metadata=Tcell,col.name="T_cell")
sp <- AddMetaData(sp,metadata=Bcell,col.name="B_cell")
sp <- AddMetaData(sp,metadata=Dendritic,col.name="Dendritic_cell")
sp <- AddMetaData(sp,metadata=Myeloid,col.name="Myeloid_cell")
sp <- AddMetaData(sp,metadata=Mast,col.name="Mast_cell")
sp <- AddMetaData(sp,metadata=Plasma,col.name="Plasma_cell")
sp <- AddMetaData(sp,metadata=Other,col.name="Other")

pdf(paste0(outputdir,'spatplot3.pdf'), width=8, height=8)
SpatialFeaturePlot(sp,features=c("Ovarian_Cancer","Fibroblast","Endothelial_cell","T_cell","B_cell","Plasma_cell","Myeloid_cell","Dendritic_cell","Mast_cell"))
dev.off()


### Heatmap ###
pdf(paste0(outputdir,'heatmap_by_seurat_clusters_zscores.pdf'))
plotMetaDataCellsHeatmap(gobject = my_giotto_object, metadata_cols = 'seurat_clusters', value_cols = cell_types_subset, spat_enr_names = 'rank', x_text_size = 8, y_text_size = 8, show_values = 'zscores')
dev.off()

pdf(paste0(outputdir,'heatmap_by_seurat_clusters_zscores_rescaled.pdf'))
plotMetaDataCellsHeatmap(gobject = my_giotto_object, metadata_cols = 'seurat_clusters', value_cols = cell_types_subset, spat_enr_names = 'rank', x_text_size = 8, y_text_size = 8, show_values = 'zscores_rescaled')
dev.off()

pdf(paste0(outputdir,'heatmap_by_seurat_clusters_original.pdf'))
plotMetaDataCellsHeatmap(gobject = my_giotto_object, metadata_cols = 'seurat_clusters', value_cols = cell_types_subset, spat_enr_names = 'rank', x_text_size = 8, y_text_size = 8, show_values = 'original')
dev.off()


### Score_by_clusters ###
table <- calculateMetaTableCells(my_giotto_object, metadata_cols = 'seurat_clusters', value_cols = cell_types_subset, spat_enr_names = 'rank')
write.table(table,paste0(outputdir,'enrichment_score_by_seurat_clusters_original.txt'),sep='\t',quote=F)

table2 <- table[, zscores := scale(value), by = c('variable')]
write.table(table2,paste0(outputdir,'enrichment_score_by_seurat_clusters_zscores.txt'),sep='\t',quote=F)

table3 <- table2[, zscores_rescaled_per_gene := scales::rescale(zscores, to = c(-1,1)), by = c('variable')]
write.table(table3,paste0(outputdir,'enrichment_score_by_seurat_clusters_zscores_rescaled..txt'),sep='\t',quote=F)

