setwd("/localdisk/home/s2172876/project/DEA")
library(stringr)
library(DESeq2)
library(edgeR)
library(BiocParallel)
library(tibble)
# import data
study_id <- read.csv('../annotation/pre_re_annotation.csv',sep = '\t')[,c(1,2)]
experiment <- read.csv('../annotation/group.txt',sep = '\t')
abstract <- read.csv('/localdisk/home/s2172876/project/DEA/recount3_selection_2022-05-26 17_28_49.csv')[,c(1,3)]
study_id <-  merge(study_id,abstract,by.x = 'study',by.y = 'project')
single_factor_multi_levels <- read.csv('single_factor_multi_levels.txt',header = FALSE)

# gene annotations
human_gene_annotation <- read.table('human_gene_annotation.txt',sep = '\t',header = TRUE)
human_gene_annotation <- human_gene_annotation[which(human_gene_annotation$Gene.type=='protein_coding'),]
mouse_gene_annotation <- read.table('mouse_gene_annotation.txt',sep = ',',header = TRUE)
mouse_gene_annotation <- mouse_gene_annotation[which(mouse_gene_annotation$Gene.type=='protein_coding'),]

write.table(matrix(c('ExpID','ProjectID','Comparison','DEG','Up','Down','non_DEG','Group1','Group2'), nrow=1), 
            file="/localdisk/home/s2172876/project/DEA_Result/summary_1.txt", append=TRUE, sep = '\t',quote = FALSE,row.names = FALSE,col.names = FALSE)
            
k = 1 # ExpID counter
# multi_level <- function(h){
for(i in single_factor_multi_levels[,1]){
    try({
        # i = single_factor_multi_levels[h,1]
        # print(i)
        groups <- experiment[which(experiment[,1]==i),]
        # trying to find which annotation is different
        dif <- vector()
        for(j in 1:9){
            dif <- c(dif,length(table(groups[,j])))
        }
        # construct design matrix
        coldata<-data.frame()
        rn <- vector()
        factor <- vector()
        for(j in 1:nrow(groups)){
            rn <- c(rn,c(str_split(groups[j,10],',')[[1]]))
            factor <- c(factor,rep(groups[j,which(dif!=1)],length(str_split(groups[j,10],',')[[1]])))
        }
        coldata <- data.frame(row.names=rn, factor)
        coldata[1] <- gsub('\\+','(pos)',coldata[,1])
        coldata[1] <- gsub('-','(neg)',coldata[,1])
        colnames(coldata)[1] <- colnames(groups)[which(dif!=1)]
        # print(coldata)
        # import expression matrix
        Expr <- read.table(paste0('/localdisk/home/s2172876/project/Count/',i,'.gz'),stringsAsFactors = FALSE)
        # subset of the expression matrix
        Expr <- Expr[,rn]
        if(ncol(Expr)<50){
            #  filter protein coding gene
            if(abstract[which(abstract[,2]==i),1] == 'mouse'){
                Expr <- Expr[which(sapply(str_split(row.names(Expr), '\\.'),"[[",1) %in% mouse_gene_annotation$Gene.stable.ID),]
            }else{
                Expr <- Expr[which(sapply(str_split(row.names(Expr), '\\.'),"[[",1) %in% human_gene_annotation$Gene.stable.ID),]
            }
            # construct DEA object
            dds <- DESeqDataSetFromMatrix(countData = Expr,colData = coldata,design = as.formula(paste('~0+',colnames(coldata),sep = ' ')))
            # Apply filters to remove genes that do not pass quality filters
            # first filter based upon counts
            keep <- rowSums(counts(dds)) > 1
            dds <- dds[keep,]
            nrow(dds)
            # second filter-require 25% in >= 2 samples
            keep <- rowSums(cpm(dds) >= quantile(cpm(dds),obs=0.25)[2]) >= min(table(coldata))/2
            dds <- dds[keep,]
            nrow(dds)
            
            # Defferential gene expression
            deg <- DESeq(dds)
            # resultsNames(deg)
            contrast_vector <- rep(-1/(length(unique(coldata[,1]))-1),length(unique(coldata[,1])))
            coldata[1] <- gsub('\\(pos\\)','\\+',coldata[,1])
            coldata[1] <- gsub('\\(neg\\)','-',coldata[,1])
            for(j in 1:length(unique(coldata[,1]))){
                # contrast information
                contrast_description <- (paste0(unique(coldata[,1])[j],' vs ','(',paste(unique(coldata[,1])[-j],collapse=" + "),')/',length(unique(coldata[,1]))-1))
                contrast_vector[j] <- 1
                res <- results(deg, contrast = contrast_vector)
                res <- data.frame(res)
                # remove gene id version
                row.names(res) <- sapply(str_split(row.names(res), '\\.'),"[[",1)
                # round
                res[,1:4] <- round(res[,1:4],2)
                res[,5:6] <- signif(res[,5:6],2)
                # DEG filter
                deg_res <- res[which((abs(res$log2FoldChange)>=1) & (res$padj<=0.05)),]
                up <- res[which((res$log2FoldChange>=1) & (res$padj<=0.05)),]
                contrast_vector <- rep(-1/(length(unique(coldata[,1]))-1),length(unique(coldata[,1])))
                exp_info <- matrix(c(paste0('Exp',k),i,contrast_description,
                                        dim(deg_res)[1],dim(up)[1],dim(deg_res)[1]-dim(up)[1],
                                        nrow(dds)-dim(deg_res)[1],groups[which(groups[,which(dif!=1)]==unique(coldata[,1])[j]),'external_id'],
                                        paste(apply(groups[which(groups[,which(dif!=1)]!=unique(coldata[,1])[j]),c(which(dif!=1),10)] , 1 , paste , collapse = ":" ),collapse = ';')), nrow=1)
                write.table(exp_info, file="/localdisk/home/s2172876/project/DEA_Result/summary_1.txt", append=TRUE, sep = '\t',quote = FALSE,row.names = FALSE,col.names = FALSE)
                write.table(res,file = gzfile(paste0('/localdisk/home/s2172876/project/DEA_Result/singlefactor_multilevel/',paste0('Exp',k),'.gz')),sep = '\t',quote = FALSE)
                k = k+1
            }
        }
    })
}

