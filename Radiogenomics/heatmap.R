#!/home/jhy/anaconda2/envs/scanpy3/bin/Rscript

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


Combined = ScaleData(Combined)
###

genes3 <- read.table(paste0(inputdir,image,'_case_vs_control_markers_log2FC0.25_adjP0.05_dim_',region,'region_231115.txt'),header=T,sep=' ')
genes3 <- genes3 %>% top_n(n=100, wt=abs(avg_log2FC))
genes3<- genes3 %>% arrange(desc(avg_log2FC))

#pbmc.markers <- FindMarkers(object = Combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use=method)
#top50 <- markers %>% group_by(group) %>% top_n(n = 50, wt = avg_log2FC)
#top50 <- markers %>% top_n(n = 50, wt = avg_log2FC)

#pdf(file=paste0(inputdir,'/',image,"_case_vs_control_",region,"_heatmap_log2FC1_adjP0.05.pdf"), height=12, width=6)
#DoHeatmap(Combined, features=p2$gene)
#DoHeatmap(Combined, features=rownames(genes2),group.by='group')
#dev.off()

pdf(file=paste0(inputdir,'/',image,"_case_vs_control_",region,"_heatmap_log2FC0.25_adjP0.05_top100.pdf"), height=15, width=6)
#DoHeatmap(Combined, features=p3$gene)
DoHeatmap(Combined, features=rownames(genes3),group.by='group')
dev.off()
