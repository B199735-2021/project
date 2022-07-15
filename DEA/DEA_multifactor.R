setwd("/localdisk/home/s2172876/project/DEA")
library(stringr)
library(DESeq2)
library(edgeR)
library(tibble)

# import data
study_id <- read.csv('../annotation/pre_re_annotation.csv',sep = '\t')[,c(1,2)]
experiment <- read.csv('../annotation/group.txt',sep = '\t')
abstract <- read.csv('/localdisk/home/s2172876/project/DEA/recount3_selection_2022-05-26 17_28_49.csv')[,c(1,3)]
study_id <-  merge(study_id,abstract,by.x = 'study',by.y = 'project')
two_factor <- read.csv('two_factor.txt',header = FALSE)
three_factor <- read.csv('three_factor.txt',header = FALSE)
multiple_factor <- rbind(two_factor,three_factor)
# gene annotations
human_gene_annotation <- read.table('human_gene_annotation.txt',sep = '\t',header = TRUE)
human_gene_annotation <- human_gene_annotation[which(human_gene_annotation$Gene.type=='protein_coding'),]
mouse_gene_annotation <- read.table('mouse_gene_annotation.txt',sep = ',',header = TRUE)
mouse_gene_annotation <- mouse_gene_annotation[which(mouse_gene_annotation$Gene.type=='protein_coding'),]

write.table(matrix(c('ExpID','ProjectID','Comparison','DEG','Up','Down','non_DEG','Group1','Group2'), nrow=1), 
            file="/localdisk/home/s2172876/project/DEA_Result/summary_2.txt", append=TRUE, sep = '\t',quote = FALSE,row.names = FALSE,col.names = FALSE)
   
k = 1
for(i in multiple_factor[,1]){
    try({
        groups <- experiment[which(experiment[,1]==i),]
        # trying to find which annotation is different
        dif <- vector()
        for(j in 1:9){
            dif <- c(dif,length(table(groups[,j])))
        }
        groups[,11] <- apply( groups[,which(dif!=1)] , 1 , paste , collapse = "_" )
        colnames(groups)[11] <- paste(colnames(groups)[which(dif!=1)],collapse='_')
        groups <- groups[,-which(dif!=1)]
        # construct design matrix
        coldata<-data.frame()
        rn <- vector()
        factor <- vector()
        for(j in 1:nrow(groups)){
            rn <- c(rn,rep(str_split(groups$external_id[j],',')[[1]]))
            factor <- c(factor,rep(groups[j,ncol(groups)],length(str_split(groups$external_id[j],',')[[1]])))
        }
        coldata <- data.frame(row.names=rn, factor)
        coldata[1] <- gsub('\\+','(pos)',coldata[,1])
        coldata[1] <- gsub('-','(neg)',coldata[,1])
        colnames(coldata)[1] <- colnames(groups)[ncol(groups)]
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

            # size factor estimation for normalisation
            dds <- estimateSizeFactors(dds)

            # Defferential gene expression
            deg <- DESeq(dds)
            if(length(table(coldata[,1]))==2){
                res <- results(deg, contrast = c(colnames(coldata), levels(factor(coldata[,1]))[1], levels(factor(coldata[,1]))[2]))
                info <- mcols(res, use.names=TRUE)
                res <- data.frame(res)
                # remove gene id version
                row.names(res) <- sapply(str_split(row.names(res), '\\.'),"[[",1)
                # round
                res[,1:4] <- round(res[,1:4],2)
                res[,5:6] <- signif(res[,5:6],2)
                # remove gene id version
                row.names(res) <- sapply(str_split(row.names(res), '\\.'),"[[",1)
                # round
                res[,1:4] <- round(res[,1:4],2)
                res[,5:6] <- signif(res[,5:6],2)
                # DEG filter
                deg_res <- res[which((abs(res$log2FoldChange)>=1) & (res$padj<=0.05)),]
                up <- res[which((res$log2FoldChange>=1) & (res$padj<=0.05)),]
                contrast_vector <- rep(-1/(length(unique(coldata[,1]))-1),length(unique(coldata[,1])))

                coldata[1] <- gsub('\\(pos\\)','\\+',coldata[,1])
                coldata[1] <- gsub('\\(neg\\)','-',coldata[,1])

                exp_info <- matrix(c(paste0('Exp',k),i,word(str_match(info[2,2],"\\(MLE\\)\\: (.*?)$")[2],2,-1),
                                    dim(deg_res)[1],dim(up)[1],dim(deg_res)[1]-dim(up)[1],nrow(dds)-dim(deg_res)[1],
                                    groups[which(groups[,ncol(groups)]==levels(factor(coldata[,1]))[1]),'external_id'],
                                    groups[which(groups[,ncol(groups)]==levels(factor(coldata[,1]))[2]),'external_id']), nrow=1)
                write.table( exp_info, file="/localdisk/home/s2172876/project/DEA_Result/summary_2.txt", append=TRUE, sep = '\t',quote = FALSE,row.names = FALSE,col.names = FALSE)
                write.table(res,file = gzfile(paste0('/localdisk/home/s2172876/project/DEA_Result/multifactor/',exp_info[1],'.gz')),sep = '\t',quote = FALSE)
                k = k+1   
            }
            if(length(table(coldata[,1]))>2){
                # resultsNames(deg)
                contrast_vector <- rep(-1/(length(unique(coldata[,1]))-1),length(unique(coldata[,1])))
                coldata[1] <- gsub('\\(pos\\)','\\+',coldata[,1])
                coldata[1] <- gsub('\\(neg\\)','-',coldata[,1])
                for(j in 1:length(unique(coldata[,1]))){
                    # contrast information
                    contrast_description <- (paste0(unique(coldata[,1])[j],' vs ','(',paste(unique(coldata[,1])[-j],collapse=" + "),')/',length(unique(coldata[,1]))-1))
                    contrast_vector[j] <- 1
                    print(contrast_vector)
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
                                        nrow(dds)-dim(deg_res)[1],
                                        groups[which(groups[,colnames(groups)[ncol(groups)]]==unique(coldata[,1])[j]),'external_id'],
                                        paste(apply(groups[which(groups[,colnames(groups)[ncol(groups)]]!=unique(coldata[,1])[j]),c(ncol(groups),ncol(groups)-1)] , 1 , paste , collapse = ":" ),collapse = ';')), nrow=1)
                


                    write.table( exp_info, file="/localdisk/home/s2172876/project/DEA_Result/summary_2.txt", append=TRUE, sep = '\t',quote = FALSE,row.names = FALSE,col.names = FALSE)
                    write.table(res,file = gzfile(paste0('/localdisk/home/s2172876/project/DEA_Result/multifactor/',exp_info[1],'.gz')),sep = '\t',quote = FALSE)
                    k = k+1
                }
        }
        }
    })    
}
