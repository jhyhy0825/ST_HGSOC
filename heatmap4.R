#!/home/jhy/anaconda2/envs/scanpy3/bin/Rscript

## Rscript /home/jhy/spatial_transcriptomics/220905/code/image/heatmap4.R /home/jhy/spatial_transcriptomics/220905/result/seurat/image/ /home/jhy/spatial_transcriptomics/220905/result/seurat/image/binary1_longest_diameter_ovary/ /home/jhy/spatial_transcriptomics/220905/result/seurat/image/binary1_longest_diameter_ovary/binary1_longest_diameter_ovary.csv binary1_longest_diameter_ovary tumor Cancer

## Rscript /home/jhy/spatial_transcriptomics/220905/code/image/heatmap4.R /home/jhy/spatial_transcriptomics/220905/result/seurat/image/ /home/jhy/spatial_transcriptomics/220905/result/seurat/image/bilaterality_ovary/ /home/jhy/spatial_transcriptomics/220905/result/seurat/image/bilaterality_ovary/bilaterality_ovary.csv bilaterality_ovary tumor Cancer

## Rscript /home/jhy/spatial_transcriptomics/220905/code/image/heatmap4.R /home/jhy/spatial_transcriptomics/220905/result/seurat/image/ /home/jhy/spatial_transcriptomics/220905/result/seurat/image/peritoneal_disease/ /home/jhy/spatial_transcriptomics/220905/result/seurat/image/peritoneal_disease/peritoneal_disease.csv peritoneal_disease tumor Cancer

## Rscript /home/jhy/spatial_transcriptomics/220905/code/image/heatmap4.R /home/jhy/spatial_transcriptomics/220905/result/seurat/image/ /home/jhy/spatial_transcriptomics/220905/result/seurat/image/binary1_longest_diameter_ovary/ /home/jhy/spatial_transcriptomics/220905/result/seurat/image/binary1_longest_diameter_ovary/binary1_longest_diameter_ovary.csv binary1_longest_diameter_ovary tumor Cancer

## Rscript /home/jhy/spatial_transcriptomics/220905/code/image/heatmap4.R /home/jhy/spatial_transcriptomics/220905/result/seurat/image/ /home/jhy/spatial_transcriptomics/220905/result/seurat/image/binary2_longest_diameter_ovary/ /home/jhy/spatial_transcriptomics/220905/result/seurat/image/binary2_longest_diameter_ovary/binary2_longest_diameter_ovary.csv binary2_longest_diameter_ovary tumor Cancer

## Rscript /home/jhy/spatial_transcriptomics/220905/code/image/heatmap4.R /home/jhy/spatial_transcriptomics/220905/result/seurat/image/ /home/jhy/spatial_transcriptomics/220905/result/seurat/image/seeding_pattern/ /home/jhy/spatial_transcriptomics/220905/result/seurat/image/seeding_pattern/seeding_pattern.csv seeding_pattern tumor Cancer

## Rscript /home/jhy/spatial_transcriptomics/220905/code/image/heatmap4.R /home/jhy/spatial_transcriptomics/220905/result/seurat/image/ /home/jhy/spatial_transcriptomics/220905/result/seurat/image/seeding_paracolic_gutter/ /home/jhy/spatial_transcriptomics/220905/result/seurat/image/seeding_paracolic_gutter/seeding_paracolic_gutter.csv seeding_paracolic_gutter tumor Cancer

## Rscript /home/jhy/spatial_transcriptomics/220905/code/image/heatmap4.R /home/jhy/spatial_transcriptomics/220905/result/seurat/image/ /home/jhy/spatial_transcriptomics/220905/result/seurat/image/binary1_number_seeding_location/ /home/jhy/spatial_transcriptomics/220905/result/seurat/image/binary1_number_seeding_location/binary1_number_seeding_location.csv binary1_number_seeding_location tumor Cacner

## Rscript /home/jhy/spatial_transcriptomics/220905/code/image/heatmap4.R /home/jhy/spatial_transcriptomics/220905/result/seurat/image/ /home/jhy/spatial_transcriptomics/220905/result/seurat/image/pleural_effusion/ /home/jhy/spatial_transcriptomics/220905/result/seurat/image/pleural_effusion/pleural_effusion.csv pleural_effusion tumor Cancer

## Rscript /home/jhy/spatial_transcriptomics/220905/code/image/heatmap4.R /home/jhy/spatial_transcriptomics/220905/result/seurat/image/ /home/jhy/spatial_transcriptomics/220905/result/seurat/image/Recur/ /home/jhy/spatial_transcriptomics/220905/result/seurat/image/Recur/Recur.csv Recur tumor Cancer

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

args = commandArgs(trailingOnly = TRUE)
seuratdir <- args[1]
inputdir <- args[2]
samples <- args[3]
image <- args[4]
region <-args[5]
region2 <- args[6]

#data10 <- readRDS(file = paste0(seuratdir,"/Cancer.rds"))
data10 <- readRDS(file = paste0(seuratdir,"/",region2,".rds"))
#data20 <- readRDS(file = paste0(inputdir,"/dimensional_reduction_dim20.rds"))
#data30 <- readRDS(file = paste0(inputdir,"/dimensional_reduction_dim30.rds"))

#Idents(data10)
#Idents(data10) <- data10$orig.ident
#Idents(data10)
data10
head(data10)
#Case <- subset(data10, subset = sample=='CPS_OV10RTOV4-1' | sample=='CPS_OV5LtOV4' | sample=='CPS_OV20RtOV4'| sample=='CPS_OV24RTOV4' | sample=='CPS_OV18CDS1')
#Control <- subset(data10, subset = sample=='CPS_OV19_LtOV1' | sample=='CPS_OV1RtOV3' | sample=='CPS_OV34RtOV1' | sample=='CPS_OV71_1')
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
#Case <- subset(data10, subset = sample %in% Case$ids)
Case <- subset(data10, subset = orig.ident %in% Case$ids)
Case
Case$group <- 'Case'
head(Case)
Control <- subset(samples,samples$type == "Control")
Control
#Control <- subset(data10, subset = sample %in% Control$ids)
Control <- subset(data10, subset = orig.ident %in% Control$ids)
Control
Control$group <- 'Control' 
head(Control)

Combined <- merge(Case, Control)
Combined
head(Combined)
tail(Combined)

Idents(Combined) <- Combined$group

#Combined <- SCTransform(Combined, assay='Spatial')
#Combined <- PrepSCTFindMarkers(Combined, verbose = FALSE)

#markers <- FindMarkers(Combined, ident.1='Case', ident.2='Control', min.pct = 0.25, logfc.threshold = 0.25)
#markers <- subset(markers, subset = p_val_adj <= 0.05)
#write.table(data.frame(markers), paste0(inputdir,'/',image,'_case_vs_control_markers_log2FC0.25_adjP0.05_dim.txt'))

#markers2 <- FindMarkers(Combined, ident.1='Case', ident.2='Control', min.pct = 0.25, logfc.threshold = 2)
#markers2 <- subset(markers2, subset = p_val_adj <= 0.05)
#write.table(data.frame(markers2), paste0(inputdir,'/',image,'_case_vs_control_markers_log2FC2_adjP0.05_dim.txt'))

#markers3 <- FindMarkers(Combined, ident.1='Case', ident.2='Control', min.pct = 0.25, logfc.threshold = 1)
#markers3 <- subset(markers3, subset = p_val_adj <= 0.05)
#write.table(data.frame(markers3), paste0(inputdir,'/',image,'_case_vs_control_markers_log2FC1_adjP0.05_dim.txt'))

#pbmc.markers <- FindAllMarkers(object = Combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use=method)


###
#pbmc.markers <- FindAllMarkers(object = Combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#head(pbmc.markers)
#p <- subset(pbmc.markers, subset = p_val_adj <= 0.05)
#head(p)
#markers <- FindMarkers(Combined, ident.1='Case', ident.2='Control', min.pct = 0.25, logfc.threshold = 0.25)
#p1 <- subset(markers, subset = p_val_adj <= 0.05)
#p2 <- subset(markers, subset = avg_log2FC >= 1)
#p3 <- p1 %>% group_by(group) %>% top_n(n = 50, wt = avg_log2FC)
#head(p3)
#tail(p3)
Combined = ScaleData(Combined)
###

#genes1 <- read.table(paste0(inputdir,image,'_case_vs_control_markers_log2FC0.25_adjP0.05_dim_',region,'region_230919.txt'),header=T,sep=' ')
#genes1 <- genes1$gene
#head(genes1)
#genes1 <- genes1 %>% arrange(desc(avg_log2FC))

#genes2 <- read.table(paste0(inputdir,image,'_case_vs_control_markers_log2FC1_adjP0.05_dim_',region,'region_230919.txt'),header=T,sep=' ')
#genes2 <- genes2$gene
#genes2 <- genes2 %>% arrange(desc(avg_log2FC))

genes3 <- read.table(paste0(inputdir,image,'_case_vs_control_markers_log2FC0.25_adjP0.05_dim_',region,'region_231115.txt'),header=T,sep=' ')
genes3 <- genes3 %>% top_n(n=100, wt=abs(avg_log2FC))
genes3<- genes3 %>% arrange(desc(avg_log2FC))

#pbmc.markers <- FindMarkers(object = Combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use=method)
#top50 <- markers %>% group_by(group) %>% top_n(n = 50, wt = avg_log2FC)
#top50 <- markers %>% top_n(n = 50, wt = avg_log2FC)

#head(genes1$gene)
#head(rownames(genes2))
#pdf(file=paste0(inputdir,'/',image,"_case_vs_control_",region,"_heatmap_log2FC0.25_adjP0.05.pdf"), height=20, width=6)
#DoHeatmap(Combined, features=rownames(genes1), group.by='group')
#DoHeatmap(Combined, features=p1$gene)
#dev.off()

#pdf(file=paste0(inputdir,'/',image,"_case_vs_control_",region,"_heatmap_log2FC1_adjP0.05.pdf"), height=12, width=6)
#DoHeatmap(Combined, features=p2$gene)
#DoHeatmap(Combined, features=rownames(genes2),group.by='group')
#dev.off()

pdf(file=paste0(inputdir,'/',image,"_case_vs_control_",region,"_heatmap_log2FC0.25_adjP0.05_top100.pdf"), height=15, width=6)
#DoHeatmap(Combined, features=p3$gene)
DoHeatmap(Combined, features=rownames(genes3),group.by='group')
dev.off()

#Idents(data20)
#Idents(data20) <- data20$orig.ident
#Idents(data20)
#data20
#head(data20)
#markers20 <- FindMarkers(data20, ident.1='recur', ident.2='non', min.pct = 0.25, logfc.threshold = 0.25)
#markers20 <- subset(markers20, subset = p_val_adj <= 0.05)
#write.table(data.frame(markers20), paste0(inputdir,'markers_log2FC0.25_adjP0.05_dim20.txt'))

#Idents(data30)
#Idents(data30) <- data30$orig.ident
#Idents(data30)
#data30
#head(data30)
#markers30 <- FindMarkers(data10, ident.1='recur', ident.2='non', min.pct = 0.25, logfc.threshold = 0.25)
#markers30 <- subset(markers30, subset = p_val_adj <= 0.05)
#write.table(data.frame(markers30), paste0(inputdir,'markers_log2FC0.25_adjP0.05_dim30.txt'))


