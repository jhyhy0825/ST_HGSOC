#! /usr/bin/env Rscript
  
## /home/jhy/spatial_transcriptomics/220905/code/SpatialInferCNV/extract.R CPS_OV1RtOV3 stroma

library(Seurat)
library(dplyr)

args = commandArgs(trailingOnly = TRUE)
sample <- args[1]
region <- args[2]

data <- read.table(paste0('/home/jhy/spatial_transcriptomics/220905/result/seurat/',sample,'/SpatialInferCNV_231122/InferCNV_outputs_group/infercnv.observations.txt'),header=T,sep=' ')
your_data_frame <- as.data.frame(data)
head(your_data_frame)
#max_value_location <- which(your_data_frame == max(your_data_frame), arr.ind = TRUE)
max_value_location <- which(your_data_frame > 1, arr.ind = TRUE)
row_with_max_value <- rownames(your_data_frame)[max_value_location[, "row"]]
col_with_max_value <- colnames(your_data_frame)[max_value_location[, "col"]]
result <- list()
for(i in 1:length(col_with_max_value)){result[[i]] <- paste0(row_with_max_value[i],',',col_with_max_value[i],',',your_data_frame[row_with_max_value[i],col_with_max_value[i]])} 
#####result <- lapply(seq_along(col_with_max_value), function(i) {paste0(row_with_max_value[i], ',', col_with_max_value[i], ',', your_data_frame[row_with_max_value[i], col_with_max_value[i]])})
#####result <- mapply(function(row, col, value) paste0(row, ',', col, ',', value),row_with_max_value, col_with_max_value,your_data_frame[row_with_max_value, col_with_max_value], SIMPLIFY = FALSE)
#####result <- mapply(function(row, col, value) paste0(row, ',', col, ',', value), row_with_max_value, col_with_max_value, your_data_frame[row_with_max_value, col_with_max_value], SIMPLIFY = FALSE)
head(result)
obj <- readRDS(paste0('/home/jhy/spatial_transcriptomics/220905/result/seurat/',sample,'/SpatialInferCNV_231122/new_object.rds'))
head(obj)
stroma <- subset(obj, info==region)
cells <- colnames(stroma)
result_df <- do.call(rbind, lapply(result, function(x) unlist(strsplit(x, ","))))
head(result_df)
#new <- sub("^","CPS_OV1RtOV3_",cells)
new <- sub("^",paste0(sample,'_'),cells)
new <- sub("-1$",".1",new)
test <- subset(result_df,result_df[,2] %in% new)
head(test)
counts <- table(test[,1])
gtf <- read.table('/home/jhy/spatial_transcriptomics/220905/spaceranger/ref/refdata-gex-GRCh38-2020-A/genes/genes.gtf',sep='\t',header=F)
gtf <- subset(gtf, gtf$V3 == 'gene')
gene_ids <- sub(".*gene_id\\s(\\S+);.*", "\\1", gtf$V9)
gene_info <- gsub(".*gene_id\\s(\\S+);.*gene_name\\s(\\S+);.*", "\\1 \\2", gtf$V9)
gene_info <- strsplit(gene_info, " ")
result_df <- data.frame()
for (i in seq_along(gene_ids)) {
	gene_id <- gene_ids[i]
	if (gene_id %in% rownames(counts)) {
		position_info <- gtf[i, c("V1", "V4", "V5")]
		#gene_counts <- counts[counts$Var1 == gene_id, "Freq"]
		gene_counts <- counts[gene_id]
		gene_name <- gene_info[i][[1]][2]
		result_df <- rbind(result_df, cbind(gene_id = gene_id, gene_name = gene_name, position_info, Freq = gene_counts))
	}
}
head(result_df)
write.table(result_df,paste0('/home/jhy/spatial_transcriptomics/220905/result/seurat/',sample,'/SpatialInferCNV_231122/InferCNV_outputs_group/',sample,'_max_',region,'_counts_pos_test.txt'),sep=',',quote=F,row.names=F)


#data <- read.table(paste0('/home/jhy/spatial_transcriptomics/220905/result/seurat/',sample,'/SpatialInferCNV_231122/InferCNV_outputs_group/infercnv.observations.txt'),header=T,sep=' ')
#your_data_frame <- as.data.frame(data)
#max_value_location <- which(your_data_frame == min(your_data_frame), arr.ind = TRUE)
max_value_location <- which(your_data_frame < 1, arr.ind = TRUE)
row_with_max_value <- rownames(your_data_frame)[max_value_location[, "row"]]
col_with_max_value <- colnames(your_data_frame)[max_value_location[, "col"]]
result <- list()
for(i in 1:length(col_with_max_value)){result[[i]] <- paste0(row_with_max_value[i],',',col_with_max_value[i],',',your_data_frame[row_with_max_value[i],col_with_max_value[i]])}
#obj <- readRDS(paste0('/home/jhy/spatial_transcriptomics/220905/result/seurat/',sample,'/SpatialInferCNV_231122/new_object.rds'))
#stroma <- subset(obj, info==region)
#cells <- colnames(stroma)
result_df <- do.call(rbind, lapply(result, function(x) unlist(strsplit(x, ","))))
#new <- sub("^","CPS_OV1RtOV3_",cells)
new <- sub("^",paste0(sample,'_'),cells)
new <- sub("-1$",".1",new)
test <- subset(result_df,result_df[,2] %in% new)
counts <- table(test[,1])
head(counts)
#gtf <- read.table('/home/jhy/spatial_transcriptomics/220905/spaceranger/ref/refdata-gex-GRCh38-2020-A/genes/genes.gtf',sep='\t',header=F)
#gtf <- subset(gtf, gtf$V3 == 'gene')
#gene_ids <- sub(".*gene_id\\s(\\S+);.*", "\\1", gtf$V9)
#gene_info <- gsub(".*gene_id\\s(\\S+);.*gene_name\\s(\\S+);.*", "\\1 \\2", gtf$V9)
#gene_info <- strsplit(gene_info, " ")
result_df <- data.frame()
for (i in seq_along(gene_ids)) {
	gene_id <- gene_ids[i]
	if (gene_id %in% rownames(counts)) {
		position_info <- gtf[i, c("V1", "V4", "V5")]
		gene_counts <- counts[gene_id]
		gene_name <- gene_info[i][[1]][2]
		result_df <- rbind(result_df, cbind(gene_id = gene_id, gene_name = gene_name, position_info, Freq = gene_counts))}}
head(result_df)
write.table(result_df,paste0('/home/jhy/spatial_transcriptomics/220905/result/seurat/',sample,'/SpatialInferCNV_231122/InferCNV_outputs_group/',sample,'_min_',region,'_counts_pos_test.txt'),sep=',',quote=F,row.names=F)


