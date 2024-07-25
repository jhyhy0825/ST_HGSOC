#data <- read.table('/home/jhy/spatial_transcriptomics/220905/result/seurat/CellPhoneDB/output_non_0.2/non0.2_sigpvals.txt',sep='\t',header=T)
data <- read.table('/home/jhy/spatial_transcriptomics/220905/result/seurat/CellPhoneDB/output_recur_0.2/recur0.2_sigpvals.txt',sep='\t',header=T)
#cells <- colnames(data)[14:77]
cells <- colnames(data)[14:94]
cells <- lapply(cells, function(x) gsub('cell.','cell - ',x))
cells <- lapply(cells, function(x) gsub('Fibroblast.','Fibroblast - ',x))
cells <- lapply(cells, function(x) gsub('\\.', ' ', x))
colnames(data) <- c(colnames(data)[1:13],cells)
create_combinations <- function(list_a) {
	result <- character()
	for (i in seq_along(list_a)) {
		for (j in seq_along(list_a)) {
			result <- c(result, paste0(list_a[i], " - ", list_a[j]))
		}
	}
	return(result)
}
#a <- c('Ovarian cancer cell','Fibroblast','T cell','B cell','Plasma cell','Mast cell','Dendritic cell','Myeloid cell')
a <- c('Ovarian cancer cell','Fibroblast','Endothelial cell','T cell','B cell','Plasma cell','Mast cell','Dendritic cell','Myeloid cell')       
b <- create_combinations(a)
colname <- c(colnames(data)[1:13],b)
new <- data[,colname]

library(reshape2)
#test <- t(new[,14:77])
test <- t(new[,14:94])
colnames(test) <- new$interacting_pair
melted_df <- melt(test, id.vars = rownames(test))
colnames(melted_df) <- c("cell_cell", "gene_gene", "pval")
head(melted_df)

#data <- read.table('/home/jhy/spatial_transcriptomics/220905/result/seurat/CellPhoneDB/output_non_0.2/non0.2_sigmeans2.txt',header=T,sep='\t')
data <- read.table('/home/jhy/spatial_transcriptomics/220905/result/seurat/CellPhoneDB/output_recur_0.2/recur0.2_sigmeans2.txt',header=T,sep='\t')
#cells <- colnames(data)[14:77]
cells <- colnames(data)[14:94]
cells <- lapply(cells, function(x) gsub('cell.','cell - ',x))
cells <- lapply(cells, function(x) gsub('Fibroblast.','Fibroblast - ',x))
cells <- lapply(cells, function(x) gsub('\\.', ' ', x)) 
colnames(data) <- c(colnames(data)[1:13],cells)
create_combinations <- function(list_a) {
	result <- character()
	for (i in seq_along(list_a)) {
		for (j in seq_along(list_a)) {
			result <- c(result, paste0(list_a[i], " - ", list_a[j]))
		}
	}
	return(result)
}
#a <- c('Ovarian cancer cell','Fibroblast','T cell','B cell','Plasma cell','Mast cell','Dendritic cell','Myeloid cell')
a <- c('Ovarian cancer cell','Fibroblast','Endothelial cell','T cell','B cell','Plasma cell','Mast cell','Dendritic cell','Myeloid cell')
b <- create_combinations(a)
colname <- c(colnames(data)[1:13],b)
new <- data[,colname]
#test <- t(new[,14:77])
test <- t(new[,14:94])
colnames(test) <- new$interacting_pair
melted_df2 <- melt(test, id.vars = rownames(test))
colnames(melted_df2) <- c("cell_Cell", "gene_gene", "mean")
head(melted_df2)

melted_df$mean <- melted_df2$mean
head(melted_df)

data <- melted_df

split_cell_cell <- strsplit(as.character(data$cell_cell), " - ")
data$cell1 <- sapply(split_cell_cell, `[`, 1)
data$cell2 <- sapply(split_cell_cell, `[`, 2)
order_sequence <- c('Ovarian cancer cell', 'Fibroblast', 'Endothelial cell', 'T cell', 'B cell', 'Plasma cell', 'Mast cell', 'Dendritic cell', 'Myeloid cell')
#order_sequence <- c('Ovarian cancer cell','Fibroblast','T cell','B cell','Plasma cell','Mast cell','Dendritic cell','Myeloid cell')
data$cell3 <- as.numeric(factor(data$cell1, levels = order_sequence))
data$cell_cell <- paste(data$cell2, data$cell3, sep='   ')

unique_values <- unique(data$cell_cell)
data$cell_cell <- factor(data$cell_cell, levels = unique_values)

neg_log_pval <- -log10(data$pval)
data$neglogp <- ifelse(is.infinite(neg_log_pval), 3, neg_log_pval)

#pdf('/home/jhy/spatial_transcriptomics/220905/result/seurat/CellPhoneDB/output_non_0.2/method2_dotplot.pdf',width=12,height=10)
pdf('/home/jhy/spatial_transcriptomics/220905/result/seurat/CellPhoneDB/output_recur_0.2/method2_dotplot.pdf',width=15,height=12)
ggplot(data=data,aes(x=cell_cell, y=gene_gene, color=log2(mean), size=neglogp)) + geom_point() + scale_color_viridis_c(name = 'log2 (mean)') +cowplot::theme_cowplot() + theme(axis.line  = element_blank()) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab('') + theme(axis.ticks = element_blank())
#ggplot(data=data,aes(x=cell_cell, y=gene_gene, color=log2(mean), size=-log10(pval))) + geom_point() + scale_color_viridis_c(name = 'log2 (mean)') +cowplot::theme_cowplot() + theme(axis.line  = element_blank()) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab('') + theme(axis.ticks = element_blank())
dev.off()


