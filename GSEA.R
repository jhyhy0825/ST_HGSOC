## Rscript /home/jhy/spatial_transcriptomics/220905/code/image/GSEA.R /home/jhy/spatial_transcriptomics/TCIA/OV/TCGA-OV/DEG/bilaterality_ovary/ /home/jhy/spatial_transcriptomics/220905/result/seurat/image/bilaterality_ovary/ bilaterality_ovary /home/jhy/spatial_transcriptomics/220905/result/seurat/image/Stroma.rds stroma

## Rscript /home/jhy/spatial_transcriptomics/220905/code/image/GSEA.R /home/jhy/spatial_transcriptomics/TCIA/OV/TCGA-OV/DEG/binary1_longest_diameter_ovary/ /home/jhy/spatial_transcriptomics/220905/result/seurat/image/binary1_longest_diameter_ovary/ binary1_longest_diameter_ovary /home/jhy/spatial_transcriptomics/220905/result/seurat/image/Stroma.rds stroma

## Rscript /home/jhy/spatial_transcriptomics/220905/code/image/GSEA.R /home/jhy/spatial_transcriptomics/TCIA/OV/TCGA-OV/DEG/binary2_longest_diameter_ovary/ /home/jhy/spatial_transcriptomics/220905/result/seurat/image/binary2_longest_diameter_ovary/ binary2_longest_diameter_ovary /home/jhy/spatial_transcriptomics/220905/result/seurat/image/Stroma.rds stroma

## Rscript /home/jhy/spatial_transcriptomics/220905/code/image/GSEA.R /home/jhy/spatial_transcriptomics/TCIA/OV/TCGA-OV/DEG/seeding_pattern/ /home/jhy/spatial_transcriptomics/220905/result/seurat/image/seeding_pattern/ seeding_pattern /home/jhy/spatial_transcriptomics/220905/result/seurat/image/Stroma.rds stroma

## Rscript /home/jhy/spatial_transcriptomics/220905/code/image/GSEA.R /home/jhy/spatial_transcriptomics/TCIA/OV/TCGA-OV/DEG/seeding_paracolic_gutter/ /home/jhy/spatial_transcriptomics/220905/result/seurat/image/seeding_paracolic_gutter/ seeding_paracolic_gutter /home/jhy/spatial_transcriptomics/220905/result/seurat/image/Stroma.rds stroma

## Rscript /home/jhy/spatial_transcriptomics/220905/code/image/GSEA.R /home/jhy/spatial_transcriptomics/TCIA/OV/TCGA-OV/DEG/binary1_number_seeding_location/ /home/jhy/spatial_transcriptomics/220905/result/seurat/image/binary1_number_seeding_location/ binary1_number_seeding_location /home/jhy/spatial_transcriptomics/220905/result/seurat/image/Stroma.rds stroma

## Rscript /home/jhy/spatial_transcriptomics/220905/code/image/GSEA.R /home/jhy/spatial_transcriptomics/TCIA/OV/TCGA-OV/DEG/peritoneal_disease/ /home/jhy/spatial_transcriptomics/220905/result/seurat/image/peritoneal_disease/ peritoneal_disease /home/jhy/spatial_transcriptomics/220905/result/seurat/image/Stroma.rds stroma

## Rscript /home/jhy/spatial_transcriptomics/220905/code/image/GSEA.R /home/jhy/spatial_transcriptomics/TCIA/OV/TCGA-OV/DEG/pleural_effusion/ /home/jhy/spatial_transcriptomics/220905/result/seurat/image/pleural_effusion/ pleural_effusion /home/jhy/spatial_transcriptomics/220905/result/seurat/image/Stroma.rds stroma

## Rscript /home/jhy/spatial_transcriptomics/220905/code/image/GSEA.R /home/jhy/spatial_transcriptomics/220905/result/seurat/image/Recur/ /home/jhy/spatial_transcriptomics/220905/result/seurat/image/Recur/ Recur /home/jhy/spatial_transcriptomics/220905/result/seurat/image/Cancer.rds tumor

library(Seurat)
library(presto)
library(msigdbr)
library(fgsea)
library(dplyr)
library(ggplot2)
library(tibble)

args = commandArgs(trailingOnly = TRUE)
inputdir1 <- args[1]
inputdir2 <- args[2]
factor <- args[3]
rds <- args[4]
region <- args[5]

#data1 <- read.table(paste0(inputdir1,"edgeR_result_with_gene_name.txt"), sep = '\t', header = TRUE)
#head(data1)
data1 <- read.table(paste0(inputdir2,factor,'_case_vs_control_markers_log2FC0.25_adjP0.05_dim_',region,'region_230919.txt'),sep=' ',header=T)
data1$gene_name <- rownames(data1)

m_df <- msigdbr(species = "Homo sapiens", category = "H")
head(m_df)
gs <- m_df$gs_name
gs <- unique(gs)

fgsea_sets <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

#genes1 <- data1 %>%
#	arrange(desc(logFC)) %>%
#	dplyr::select(gene_name, logFC)
#colnames(genes1) <- c('feature','logFC')
#head(genes1)
genes1 <- data1 %>%
	arrange(desc(avg_log2FC)) %>%
	dplyr::select(gene_name, avg_log2FC)
colnames(genes1) <- c('feature','logFC')
head(genes1)

data2 <- readRDS(rds)
samples <- read.table(paste0(inputdir2,factor,'.csv'),header=T,sep=',')
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
Case <- subset(data2, subset = orig.ident %in% Case$ids)
Case
Case$group <- 'Case'
head(Case)
Control <- subset(samples,samples$type == "Control")
Control
Control <- subset(data2, subset = orig.ident %in% Control$ids)
Control
Control$group <- 'Control'
head(Control)

Combined <- merge(Case, Control)
Combined
head(Combined)
tail(Combined)
Idents(Combined) <- Combined$group
Combined <- PrepSCTFindMarkers(Combined, verbose = FALSE)
pbmc.genes <- wilcoxauc(Combined, 'group', seurat_assay='SCT')
pbmc.genes <- subset(pbmc.genes, group == 'Case')
head(pbmc.genes)
dplyr::count(pbmc.genes, group)

genes2 <- pbmc.genes %>%
	arrange(desc(logFC)) %>%
	dplyr::select(feature, logFC)

#genes <- rbind(genes1,genes2)
#genes <- genes1

ranks1 <- deframe(genes1)
head(ranks1)
length(ranks1)

fgseaRes1 <- fgsea(fgsea_sets, stats = ranks1, nperm = 1000)

fgseaResTidy1 <- fgseaRes1 %>%
	as_tibble() %>%
	arrange(desc(NES))

fgseaResTidy1 %>%
	dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>%
	arrange(padj) %>%
	head()
tail(fgseaResTidy1)

result1 <- fgseaResTidy1 %>%  mutate(leadingEdge = sapply(leadingEdge, function(x) paste(x, collapse = ",")))
result1 <- bind_rows(result1)
write.table(result1,paste0(inputdir2,'/Hall_pathways_StMary_',region,'_usinglog2FC0.25adjP0.05.txt'),sep='\t',quote=F,row.names=F)

pdf(paste0(inputdir2,'/Hall_pathways_p0.05_StMary_',region,'_usinglog2FC0.25adjP0.05.pdf'), width=8, height=8)
ggplot(fgseaResTidy1 %>% head(n= 100), aes(reorder(pathway, NES), NES)) +
	geom_col(aes(fill= padj < 0.05)) +
	coord_flip() +
	labs(x="Pathway", y="Normalized Enrichment Score", title=paste0("Hallmark pathway ",factor,' ',region)) +
	theme_minimal()
dev.off()


ranks2 <- deframe(genes2)
head(ranks2)
length(ranks2)

fgseaRes2 <- fgsea(fgsea_sets, stats = ranks2, nperm = 1000)

fgseaResTidy2 <- fgseaRes2 %>%
	as_tibble() %>%
	arrange(desc(NES))

fgseaResTidy2 %>%
	dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>%
	arrange(padj) %>%
	head()
tail(fgseaResTidy2)

result2 <- fgseaResTidy2 %>%  mutate(leadingEdge = sapply(leadingEdge, function(x) paste(x, collapse = ",")))
result2 <- bind_rows(result2)
write.table(result2,paste0(inputdir2,'/Hall_pathways_StMary_',region,'.txt'),sep='\t',quote=F,row.names=F)
pdf(paste0(inputdir2,'/Hall_pathways_p0.05_StMary_',region,'.pdf'), width=8, height=8)
ggplot(fgseaResTidy2 %>% head(n= 100), aes(reorder(pathway, NES), NES)) +
	geom_col(aes(fill= padj < 0.05)) +
	coord_flip() +
	labs(x="Pathway", y="Normalized Enrichment Score", title=paste0("Hallmark pathway ",factor,' ',region)) +
	theme_minimal()
dev.off()


