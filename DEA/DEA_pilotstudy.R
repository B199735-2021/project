setwd("/localdisk/home/s2172876/project/DEA")
library(recount3)
library(stringr)
library(DESeq2)
library(edgeR)
library(BiocParallel)
library(tibble)
library(doParallel) 
library(BiocParallel) 
register(MulticoreParam(10))
no_cores <- 10
study_id <- read.csv('../annotation/pre_re_annotation.csv',sep = '\t')[,c(1,2)]
experiment <- read.csv('../annotation/group.txt',sep = '\t')
abstract <- read.csv('/localdisk/home/s2172876/project/DEA/recount3_selection_2022-05-26 17_28_49.csv')[,c(1,3)]
study_id <-  merge(study_id,abstract,by.x = 'study',by.y = 'project')
human_gene_annotation <- read.table('human_gene_annotation.txt',sep = '\t',header = TRUE)
human_gene_annotation <- human_gene_annotation[which(human_gene_annotation$Gene.type=='protein_coding'),]
mouse_gene_annotation <- read.table('mouse_gene_annotation.txt',sep = ',',header = TRUE)
mouse_gene_annotation <- mouse_gene_annotation[which(mouse_gene_annotation$Gene.type=='protein_coding'),]


projects <- rbind(
    recount3::available_projects("human"),

    recount3::available_projects("mouse")
)

# if single factor, two level
DEA <- function(g1, g2, cdata, project){
    
}
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

###############################################################################################################
for(i in 'SRP029367'){

}
groups <- experiment[which(experiment[,1]=='SRP072457'),]
dif <- vector()
for(j in 1:8){
    dif <- c(dif,length(table(groups[,j])))
}
coldata<-data.frame()
# single factor
if(length(which(dif!=1))==1){
    if(dif[dif!=1]==2){ # two levels
        print('two level')
    }else{ # multiple levels
        if(which(grepl('wildtype', groups[,which(dif!=1)]))) { # include control group
            # control row number
            cn <- which(grepl('wildtype', groups[,which(dif!=1)]))
            # other rows
            otn <- which(!grepl('wildtype', groups[,which(dif!=1)]))
            # for循环，每个和它比
            for(k in 1){
                coldata <- data.frame(row.names=c(str_split(groups[cn,10],',')[[1]],
                                                    str_split(groups[k,10],',')[[1]]), 
                                    factor = c(rep(groups[cn,which(dif!=1)],
                                                    length(str_split(groups[cn,10],',')[[1]])),
                                                rep(groups[k,which(dif!=1)],
                                                    length(str_split(groups[k,10],',')[[1]]))))
                coldata$factor <- gsub(' ', "_", coldata$factor)
                colnames(coldata)[1] <- colnames(groups)[which(dif!=1)]
            }
    }else{ # all combination
        print('all combination')
    }

}
}

save.image(file = "/localdisk/home/s2172876/project/DEA/pilotstudy.RData")


k1=2
g1 <- str_split(groups[cn,10],',')[[1]]
g2 <- str_split(groups[k1,10],',')[[1]]
k2=3
g3 <- str_split(groups[k2,10],',')[[1]]
coldata <- data.frame(row.names=c(str_split(groups[cn,10],',')[[1]],
                                                    str_split(groups[k1,10],',')[[1]],
                                                    str_split(groups[k2,10],',')[[1]]), 
                                    factor = c(rep(groups[cn,which(dif!=1)],
                                                    length(str_split(groups[cn,10],',')[[1]])),
                                                rep(groups[k1,which(dif!=1)],
                                                    length(str_split(groups[k1,10],',')[[1]])),
                                                rep(groups[k2,which(dif!=1)],
                                                    length(str_split(groups[k2,10],',')[[1]]))))
coldata$factor <- gsub(' ', "_", coldata$factor)
colnames(coldata)[1] <- colnames(groups)[which(dif!=1)]
Expr <- read.table('/localdisk/home/s2172876/project/Count/SRP029367.gz',stringsAsFactors = FALSE)
# subset of the expression matrix
Expr <- Expr[,c(g1,g2,g3)]
#  filter protein coding gene
Expr <- Expr[which(sapply(str_split(row.names(Expr), '\\.'),"[[",1) %in% human_gene_annotation$Gene.stable.ID),]
dim(Expr)
dds <- DESeqDataSetFromMatrix(countData = Expr,colData = coldata,design = as.formula(paste('~',colnames(coldata),sep = ' ')))
# Apply filters to remove genes that do not pass quality filters
# first filter based upon counts
keep <- rowSums(counts(dds)) > 0
dds <- dds[keep,]
nrow(dds)

# second filter-require >5 in >= 2 samples
keep <- rowSums(cpm(dds) >= round(quantile(cpm(dds),obs=0.25)[2],0)) >= min(length(g1),length(g2))/2
dds <- dds[keep,]
nrow(dds)

# size factor estimation for normalisation
dds <- estimateSizeFactors(dds)
sizeFactors(dds)

# Defferential gene expression
deg <- DESeq(dds,parallel=TRUE)
resultsNames(deg)
res <- results(deg, contrast = c(0,0.5,0.5))
deg <- data.frame(res)
res
