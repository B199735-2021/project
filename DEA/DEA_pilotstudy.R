setwd("/localdisk/home/s2172876/project/DEA")
library(recount3)
library(stringr)
library(DESeq2)
library(edgeR)
library(tibble)

study_id <- read.csv('../annotation/pre_re_annotation.csv',sep = '\t')[,c(1,2)]
experiment <- read.csv('../annotation/group.txt',sep = '\t')
abstract <- read.csv('/localdisk/home/s2172876/project/DEA/recount3_selection_2022-05-26 17_28_49.csv')[,c(1,3)]
study_id <-  merge(study_id,abstract,by.x = 'study',by.y = 'project')
human_gene_annotation <- read.table('human_gene_annotation.txt',sep = '\t',header = TRUE)
human_gene_annotation <- human_gene_annotation[which(human_gene_annotation$Gene.type=='protein_coding'),]
mouse_gene_annotation <- read.table('mouse_gene_annotation.txt',sep = ',',header = TRUE)
mouse_gene_annotation <- mouse_gene_annotation[which(mouse_gene_annotation$Gene.type=='protein_coding'),]

#########################------------------------------------------------------############################
# differentially expressed gene analysis for GSE79661
# wildtype VS knockout
# sensitivity and specificity calculation
#########################------------------------------------------------------############################
groups <- experiment[which(experiment[,1]=='SRP072457'),]
dif <- vector()
for(j in 1:8){
    dif <- c(dif,length(table(groups[,j])))
}
coldata<-data.frame()
# control row number
cn <- which(grepl('wildtype', groups[,which(dif!=1)]))
# other rows
otn <- which(!grepl('wildtype', groups[,which(dif!=1)]))
k=1
g1 <- str_split(groups[cn,10],',')[[1]]
g2 <- str_split(groups[k,10],',')[[1]]
coldata <- data.frame(row.names=c(str_split(groups[cn,10],',')[[1]],
                                                    str_split(groups[k,10],',')[[1]]), 
                                    factor = c(rep(groups[cn,which(dif!=1)],
                                                    length(str_split(groups[cn,10],',')[[1]])),
                                                rep(groups[k,which(dif!=1)],
                                                    length(str_split(groups[k,10],',')[[1]]))))
coldata$factor <- gsub(' ', "_", coldata$factor)
colnames(coldata)[1] <- colnames(groups)[which(dif!=1)]
Expr <- read.table('/localdisk/home/s2172876/project/Count/SRP072457.gz',stringsAsFactors = FALSE)
# subset of the expression matrix
Expr <- Expr[,c(g1,g2)]
#  filter protein coding gene
Expr <- Expr[which(sapply(str_split(row.names(Expr), '\\.'),"[[",1) %in% mouse_gene_annotation$Gene.stable.ID),]

dds <- DESeqDataSetFromMatrix(countData = Expr,colData = coldata,design = as.formula(paste('~',colnames(coldata),sep = ' ')))
# Apply filters to remove genes that do not pass quality filters
# first filter based upon counts
keep <- rowSums(counts(dds)) > 1
dds <- dds[keep,]
nrow(dds)

# second filter-require 25% in >= 2 samples
keep <- rowSums(cpm(dds) >= quantile(cpm(dds),obs=0.25)[2]) >= min(length(g1),length(g2))/2
dds <- dds[keep,]
nrow(dds)

# size factor estimation for normalisation
dds <- estimateSizeFactors(dds)
sizeFactors(dds)

# Defferential gene expression
deg <- DESeq(dds,parallel=TRUE)
res <- results(deg, contrast = c(colnames(coldata), levels(factor(coldata[,1]))[1], levels(factor(coldata[,1]))[2]))
deg <- data.frame(res)
# deg filter
deg_res <- deg[which((abs(deg$log2FoldChange)>=1) & (deg$padj<=0.05)),]
dim(deg_res)
# compare with supplementary file
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6261680/bin/NIHMS1507725-supplement-Sup_Table_3.xlsx
GSE79661 <- read.csv('GSE79661.txt',sep='\t',header=TRUE)[,c(3,10,13)]
GSE79661 <- GSE79661[which(abs(GSE79661$log2.fold_change.)>1),]
GSE79661 <- merge(GSE79661,mouse_gene_annotation,by.x='gene',by.y='Gene.name')
# sensitivity and specificity
TP <- length(intersect(sapply(str_split(row.names(deg_res), '\\.'),"[[",1),GSE79661[,4]))
FN <- length(setdiff(GSE79661[,4],sapply(str_split(row.names(deg_res), '\\.'),"[[",1)))
FP <- length(setdiff(sapply(str_split(row.names(deg_res), '\\.'),"[[",1),GSE79661[,4]))
TN <- nrow(dds) - length(union(sapply(str_split(row.names(deg_res), '\\.'),"[[",1), GSE79661[,4]))
sensitivity <- TP/(TP+FN)
specificity <- TN/(TN+FP)
sensitivity
specificity

# save.image(file = "/localdisk/home/s2172876/project/DEA/pilotstudy.RData")
