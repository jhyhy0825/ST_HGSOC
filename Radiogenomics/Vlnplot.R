## Rscript /home/jhy/spatial_transcriptomics/220905/code/image/Vlnplot.R /home/jhy/spatial_transcriptomics/220905/result/seurat/image/ /home/jhy/spatial_transcriptomics/220905/result/seurat/image/binary1_longest_diameter_ovary/ /home/jhy/spatial_transcriptomics/220905/result/seurat/image/binary1_longest_diameter_ovary/binary1_longest_diameter_ovary.csv binary1_longest_diameter_ovary tumor Cancer

## Rscript /home/jhy/spatial_transcriptomics/220905/code/image/Vlnplot.R /home/jhy/spatial_transcriptomics/220905/result/seurat/image/ /home/jhy/spatial_transcriptomics/220905/result/seurat/image/bilaterality_ovary/ /home/jhy/spatial_transcriptomics/220905/result/seurat/image/bilaterality_ovary/bilaterality_ovary.csv bilaterality_ovary tumor Cancer

## Rscript /home/jhy/spatial_transcriptomics/220905/code/image/Vlnplot.R /home/jhy/spatial_transcriptomics/220905/result/seurat/image/ /home/jhy/spatial_transcriptomics/220905/result/seurat/image/pleural_effusion/ /home/jhy/spatial_transcriptomics/220905/result/seurat/image/pleural_effusion/pleural_effusion.csv pleural_effusion tumor Cancer

## Rscript /home/jhy/spatial_transcriptomics/220905/code/image/Vlnplot.R /home/jhy/spatial_transcriptomics/220905/result/seurat/image/ /home/jhy/spatial_transcriptomics/220905/result/seurat/image/Recur/ /home/jhy/spatial_transcriptomics/220905/result/seurat/image/Recur/Recur.csv Recur tumor Cancer

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

data10 <- readRDS(file = paste0(seuratdir,"/",region2,".rds"))

data10
head(data10)

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

#Combined = ScaleData(Combined)

library(patchwork)
library(ggpubr)

#genes <- c('ABCA1','EFNB2','EPHB3','APLP2','APOA1','C3','EFNA1','EFNB1','EPHA1','EPHA2','CXCL14','CXCR4','NTN4')
genes <- c('MGST1','PPP1R14A','SH3KBP1','FOXM1','TCEAL2','SHKBP1','MEG3','PTGDS','ALPL')
vp_case1 <- function(gene_signature, file_name, test_sign){
	plot_case1 <- function(signature, y_max = 3.5){
		VlnPlot(Combined, features = signature, pt.size = 0, group.by = "group", y.max = y_max) + stat_compare_means(comparisons = test_sign, method = "wilcox.test", label = "p.signif")}
plot_list <- list()
for (gene in gene_signature) {
	plot_list[[gene]] <- plot_case1(gene)
}
cowplot::plot_grid(plotlist = plot_list)
file_name <- paste0(file_name, ".pdf")
ggsave(file_name, width = 12, height = 12)
}
gene_sig <- genes
gene_sig
comparisons1 <- list(c("Case","Control"))
vp_case1(gene_signature = gene_sig, file_name = paste0(inputdir,'/',image,"_case_vs_control_",region,"_vlnplot_pt0_max3.5_240219"), test_sign = comparisons1)


#pdf(file=paste0(inputdir,'/',image,"_case_vs_control_",region,"_vlnplotpt0_with_stat.pdf"), height=4, width=6)
#VlnPlot(Combined, features =c('S100A4','HOXD9','MGST1'), pt.size = 0,y.max=6) + stat_compare_means(comparisons = list(c("Case","Control")), method = "wilcox.test")
#dev.off()




