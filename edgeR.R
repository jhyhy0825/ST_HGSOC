#! /usr/bin/env Rscript
## Rscript /home/jhy/test/211027_edgeR/code/step1_edgeR_multi.R /home/jhy/test/211027_edgeR/result/ GSE155192_TPM_Ndufs4KO.txt

## Rscript /home/jhy/spatial_transcriptomics/TCIA/OV/TCGA-OV/code/edgeR.R /home/jhy/spatial_transcriptomics/TCIA/OV/TCGA-OV/DEG/presence_carcinomatosis_peritonei/ presence_carcinomatosis_peritonei_tpm.txt

library(edgeR)

args = commandArgs(trailingOnly = TRUE)
inputdir <- args[1]
file <- paste0(inputdir,args[2])
#files = list.files(path=inputdir, full.names=T, recursive=F, include.dirs=F)
#files = files[grep(sample, files)]

data <- read.table(file,header=F, nrow=2)
#head(data)
groups <- data[2,-1]
groups
#length(groups)
samples <- data[1,-1]
samples
#length(samples)
col <- data[2,]
data <- read.table(file,header=F, skip=2)
#head(data)
#tail(data)
#dim(data)
colnames(data)<-col
row.names(data)<-data$GeneID
data = data[,-1]
#head(data)
#tail(data)
#dim(data)

### DEG ###
dge <- DGEList(counts=data, group=groups)
#dge
#dim(dge)

#dge$samples$group <- factor(dge$samples$group, levels = c("Case","Control"))
dge$samples$group <- relevel(dge$samples$group, ref="Control")
#dge$samples$group

design <- model.matrix(~dge$samples$group)
#design
rownames(design) <- rownames(dge$samples)
#design
dge_dispersion <- estimateGLMCommonDisp(dge,design)
#dge_dispersion
dge_tagwise <- estimateGLMTagwiseDisp(dge_dispersion, design)
#dge_tagwise
fit <- glmFit(dge_tagwise,design,dispersion=dge_tagwise$tagwise.dispersion)
#fit
#lrt <- glmLRT(fit, coef=c(2,3,4))
lrt <- glmLRT(fit)
#lrt
norm <- round(cpm(dge_tagwise,normalized.lib.sizes = T))
## top10000
#normExp <- merge(norm, topTags(lrt,n=10000), by.x="row.names", by.y="row.names")
## all
normExp <- merge(norm, topTags(lrt, n=nrow(data)), by.x="row.names", by.y="row.names")
names(normExp)[1]
normExp <- normExp[order(normExp$PValue),]

outputfile <- paste0(inputdir,"edgeR_result_using_raw.txt")
write.table(normExp, outputfile, row.names=FALSE, sep='\t')
