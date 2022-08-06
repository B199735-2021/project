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
single_factor_two_levels <- read.csv('single_factor_two_levels.txt',header = FALSE)

# gene annotations
human_gene_annotation <- read.table('human_gene_annotation.txt',sep = '\t',header = TRUE)
human_gene_annotation <- human_gene_annotation[which(human_gene_annotation$Gene.type=='protein_coding'),]
mouse_gene_annotation <- read.table('mouse_gene_annotation.txt',sep = ',',header = TRUE)
mouse_gene_annotation <- mouse_gene_annotation[which(mouse_gene_annotation$Gene.type=='protein_coding'),]


exp_info <- data.frame(matrix(ncol = 9))
colnames(exp_info) <- c('ExpID','ProjectID','Comparison','DEG','Up','Down','non_DEG','Group1','Group2')
write.table(matrix(c('ExpID','ProjectID','Comparison','DEG','Up','Down','non_DEG','Group1','Group2'), nrow=1), 
            file="/localdisk/home/s2172876/project/DEA_Result/summary.txt", append=TRUE, sep = '\t',quote = FALSE,row.names = FALSE,col.names = FALSE)
t1 <- Sys.time()
k=1

for(i in single_factor_two_levels[,1]){ 
    try({
        exp_info[k,1] <- paste0('Exp',k)
        exp_info[k,2] <- i
        groups <- experiment[which(experiment[,1]==i),]
        # trying to find which annotation is different
        dif <- vector()
        for(j in 1:9){
            dif <- c(dif,length(table(groups[,j])))
        }
        # construct design matrix
        coldata <- data.frame()
        coldata <- data.frame(row.names=c(str_split(groups[1,10],',')[[1]],
                                                            str_split(groups[2,10],',')[[1]]), 
                                            factor = c(rep(groups[1,which(dif!=1)],
                                                            length(str_split(groups[1,10],',')[[1]])),
                                                        rep(groups[2,which(dif!=1)],
                                                            length(str_split(groups[2,10],',')[[1]]))))
        coldata[1] <- gsub('\\+','(pos)',coldata[,1])
        coldata[1] <- gsub('-','(neg)',coldata[,1])
        colnames(coldata)[1] <- colnames(groups)[which(dif!=1)]
        # import expression matrix
        Expr <- read.table(paste0('/localdisk/home/s2172876/project_backup/Count/',i,'.gz'),stringsAsFactors = FALSE)
        g1 <- str_split(groups[1,10],',')[[1]]
        g2 <- str_split(groups[2,10],',')[[1]]
        # subset of the expression matrix
        Expr <- Expr[,c(g1,g2)]
        if(ncol(Expr)<1000){
            #  filter protein coding gene
            if(abstract[which(abstract[,2]==i),1] == 'mouse'){
                Expr <- Expr[which(sapply(str_split(row.names(Expr), '\\.'),"[[",1) %in% mouse_gene_annotation$Gene.stable.ID),]
            }else{
                Expr <- Expr[which(sapply(str_split(row.names(Expr), '\\.'),"[[",1) %in% human_gene_annotation$Gene.stable.ID),]
            }
            # construct DEA object
            dds <- DESeqDataSetFromMatrix(countData = Expr,colData = coldata,design = as.formula(paste('~',colnames(coldata),sep = ' ')))
            # Apply filters to remove genes that do not pass quality filters
            # first filter based upon counts
            keep <- rowSums(counts(dds)) > 1
            dds <- dds[keep,]
            nrow(dds)
            # second filter-require 25% in >= 2 samples
            library(edgeR)
            keep <- rowSums(cpm(dds) >= quantile(cpm(dds),obs=0.25)[2]) >= min(length(g1),length(g2))/2
            dds <- dds[keep,]
            exp_info[k,7] <- nrow(dds)

            # size factor estimation for normalisation
            dds <- estimateSizeFactors(dds)

            # Defferential gene expression
            deg <- DESeq(dds)
            res <- results(deg, contrast = c(colnames(coldata), levels(factor(coldata[,1]))[1], levels(factor(coldata[,1]))[2]))
            deg <- data.frame(res)
            # remove gene id version
            row.names(deg) <- sapply(str_split(row.names(deg), '\\.'),"[[",1)
            # round
            deg[,1:4] <- round(deg[,1:4],2)
            deg[,5:6] <- signif(deg[,5:6],2)
            # deg filter
            exp_info[k,4] <- dim(deg[which((abs(deg$log2FoldChange)>=1) & (deg$padj<=0.05)),])[1]
            exp_info[k,5] <- dim(deg[which((deg$log2FoldChange)>=1 & (deg$padj<=0.05)),])[1]
            exp_info[k,6] <- dim(deg[which((deg$log2FoldChange)<=(-1) & (deg$padj<=0.05)),])[1]
            exp_info[k,7] <- exp_info[k,7]-exp_info[k,4]
            info <- mcols(res, use.names=TRUE)
            exp_info[k,3] <- word(str_match(info[2,2],"\\(MLE\\)\\: (.*?)$")[2],2,-1)
            exp_info[k,3] <- gsub('\\(pos\\)','\\+',str_match(info[2,2],"\\(MLE\\)\\: (.*?)$")[2])
            exp_info[k,3] <- gsub('\\(neg\\)','-',exp_info[k,3])
            coldata[1] <- gsub('\\(pos\\)','\\+',coldata[,1])
            coldata[1] <- gsub('\\(neg\\)','-',coldata[,1])
            exp_info[k,8] <- groups[which(groups[,which(dif!=1)]==levels(factor(coldata[,1]))[1]),'external_id']
            exp_info[k,9] <- groups[which(groups[,which(dif!=1)]==levels(factor(coldata[,1]))[2]),'external_id']
            write.table(deg,file = gzfile(paste0('/localdisk/home/s2172876/project/DEA_Result/singlefactor_twolevel/',exp_info[k,1],'.gz')),sep = '\t',quote = FALSE)
            write.table(exp_info[k,], file="/localdisk/home/s2172876/project/DEA_Result/summary.txt", append=TRUE, sep = '\t',quote = FALSE,row.names = FALSE,col.names = FALSE)
            k = k+1
         }
        }
       
    )
}
t2 <- Sys.time()
t2-t1
