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
tissue_specific <- read.csv('tissue_specific.txt',header = FALSE)

# gene annotations
human_gene_annotation <- read.table('human_gene_annotation.txt',sep = '\t',header = TRUE)
human_gene_annotation <- human_gene_annotation[which(human_gene_annotation$Gene.type=='protein_coding'),]
mouse_gene_annotation <- read.table('mouse_gene_annotation.txt',sep = ',',header = TRUE)
mouse_gene_annotation <- mouse_gene_annotation[which(mouse_gene_annotation$Gene.type=='protein_coding'),]

k=1
for(i in tissue_specific[,1]){
    try({
        groups <- experiment[which(experiment[,1]==i),]
        # trying to find which annotation is different
        dif <- vector()
        for(j in 1:9){
            dif <- c(dif,length(table(groups[,j])))
        }
        groups[,11] <- apply( groups[,which(dif!=1)] , 1 , paste , collapse = "_" )
        colnames(groups)[11] <- paste(colnames(groups)[which(dif!=1)],collapse='_')
        groups <- groups[,sort(c(2,which(dif[1:9]==1),10,11))]
        # construct design matrix
        coldata<-data.frame()
        rn <- vector()
        factor <- vector()
        tissue <- vector()
        for(j in 1:nrow(groups)){
            rn <- c(rn,c(str_split(groups$external_id[j],',')[[1]]))
            factor <- c(factor,rep(groups[j,ncol(groups)],length(str_split(groups$external_id[j],',')[[1]])))
            tissue <- c(tissue,rep(groups[j,2],length(str_split(groups$external_id[j],',')[[1]])))
        }
        coldata <- data.frame(row.names=rn, factor,tissue)
        coldata[1] <- gsub('\\+','(pos)',coldata[,1])
        coldata[1] <- gsub('-','(neg)',coldata[,1])
        colnames(coldata)[1] <- colnames(groups)[ncol(groups)]
        coldata <- coldata[which(coldata[,2]!=''),]
        if(nrow(coldata)<500){
            # import expression matrix
            Expr <- read.table(paste0('/localdisk/home/s2172876/project/Count/',i,'.gz'),stringsAsFactors = FALSE)
            # subset of the expression matrix
            Expr <- Expr[,rn]
            #  filter protein coding gene
            if(abstract[which(abstract[,2]==i),1] == 'mouse'){
                Expr <- Expr[which(sapply(str_split(row.names(Expr), '\\.'),"[[",1) %in% mouse_gene_annotation$Gene.stable.ID),]
            }else{
                Expr <- Expr[which(sapply(str_split(row.names(Expr), '\\.'),"[[",1) %in% human_gene_annotation$Gene.stable.ID),]
            }
            # construct DEA object
            dds <- DESeqDataSetFromMatrix(countData = Expr,colData = coldata,design = as.formula(paste('~0+',colnames(coldata)[1],sep = ' ')))
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
            # resultsNames(deg)
            coldata[1] <- gsub('\\(pos\\)','\\+',coldata[,1])
            coldata[1] <- gsub('\\(neg\\)','-',coldata[,1])
            # two tissue        
            if(length(table(coldata[,2]))==2){ 
                contrast_vector <- rep(1,length(unique(coldata[,1])))
                # group1
                g1 <- unique(coldata[which(coldata[,2]==names(table(coldata[,2])[1])),1])
                g2 <- unique(coldata[which(coldata[,2]==names(table(coldata[,2])[2])),1])
                contrast_vector[unique(coldata[,1]) %in% g1] <- rep(1/length(g1),length(g1))
                contrast_vector[unique(coldata[,1]) %in% g2] <- rep(-1/length(g2),length(g2))
                res <- results(deg, contrast = contrast_vector)
                res <- data.frame(res)
                # remove gene id version
                row.names(res) <- sapply(str_split(row.names(res), '\\.'),"[[",1)
                # round
                res[,1:4] <- round(res[,1:4],2)
                res[,5:6] <- signif(res[,5:6],2)
                deg_res <- res[which((abs(res$log2FoldChange)>=1) & (res$padj<=0.05)),]
                up <- res[which((res$log2FoldChange>=1) & (res$padj<=0.05)),]
                # description
                contrast_description <- paste0('(',paste(g1,collapse=" + "),")/",length(g1)," vs ",'(',paste(g2,collapse=" + "),")/",length(g2))
                contrast_description <- gsub('\\(pos\\)','\\+',contrast_description)
                contrast_description <- gsub('\\(neg\\)','-',contrast_description)
                # tissue description
                tissue_description <- paste(names(table(coldata[,2])[1]),"vs",names(table(coldata[,2])[2]))
                if(length(g1)>1 & length(g2)>1){
                    exp_info <- matrix(c(paste0('Exp',k),i,contrast_description,tissue_description,dim(deg_res)[1],
                                        dim(up)[1],dim(deg_res)[1]-dim(up)[1],nrow(dds)-dim(deg_res)[1],
                                        paste(apply(groups[groups[,ncol(groups)] %in% g1,c(ncol(groups),ncol(groups)-1)], 1 , paste , collapse = ":" ),collapse = ';'),
                                        paste(apply(groups[groups[,ncol(groups)] %in% g2,c(ncol(groups),ncol(groups)-1)], 1 , paste , collapse = ":" ),collapse = ';')), nrow=1)
                }
                if(length(g1)==1 & length(g2)==1){
                    exp_info <- matrix(c(paste0('Exp',k),i,contrast_description,tissue_description,dim(deg_res)[1],
                                        dim(up)[1],dim(deg_res)[1]-dim(up)[1],nrow(dds)-dim(deg_res)[1],
                                        groups[groups[,ncol(groups)]==g1,ncol(groups)-1],
                                        groups[groups[,ncol(groups)]==g2,ncol(groups)-1]), nrow=1)
                }
                if(length(g1)==1 & length(g2)>1){
                    exp_info <- matrix(c(paste0('Exp',k),i,contrast_description,tissue_description,dim(deg_res)[1],
                                        dim(up)[1],dim(deg_res)[1]-dim(up)[1],nrow(dds)-dim(deg_res)[1],
                                        groups[groups[,ncol(groups)]==g1,ncol(groups)-1],
                                        paste(apply(groups[groups[,ncol(groups)] %in% g2,c(ncol(groups),ncol(groups)-1)], 1 , paste , collapse = ":" ),collapse = ';')
                                        ), nrow=1)
                }
                if(length(g1)>1 & length(g2)==1){
                    exp_info <- matrix(c(paste0('Exp',k),i,contrast_description,tissue_description,dim(deg_res)[1],
                                        dim(up)[1],dim(deg_res)[1]-dim(up)[1],nrow(dds)-dim(deg_res)[1],
                                        paste(apply(groups[groups[,ncol(groups)] %in% g1,c(ncol(groups),ncol(groups)-1)], 1 , paste , collapse = ":" ),collapse = ';'),
                                        groups[groups[,ncol(groups)]==g2,ncol(groups)-1]), nrow=1)
                }
                write.table(exp_info, file="/localdisk/home/s2172876/project/DEA_Result_version2/summary_3.txt", append=TRUE, sep = '\t',quote = FALSE,row.names = FALSE,col.names = FALSE)
                write.table(res,file = gzfile(paste0('/localdisk/home/s2172876/project/DEA_Result_version2/tissue/',exp_info[1],'.gz')),sep = '\t',quote = FALSE)
                k = k+1
            }
            # more than two tissue
            if(length(table(coldata[,2]))>2){
                for(j in unique(coldata[,2])){
                    contrast_vector <- rep(1,length(unique(coldata[,1])))
                    g1 <- unique(coldata[which(coldata[,2]==j),1])
                    g2 <- unique(coldata[which(coldata[,2]!=j),1])
                    contrast_vector[unique(coldata[,1]) %in% g1] <- rep(1/length(g1),length(g1))
                    contrast_vector[unique(coldata[,1]) %in% g2] <- rep(-1/length(g2),length(g2))
                    res <- results(deg, contrast = contrast_vector)
                    res <- data.frame(res)
                    # remove gene id version
                    row.names(res) <- sapply(str_split(row.names(res), '\\.'),"[[",1)
                    # round
                    res[,1:4] <- round(res[,1:4],2)
                    res[,5:6] <- signif(res[,5:6],2)
                    deg_res <- res[which((abs(res$log2FoldChange)>=1) & (res$padj<=0.05)),]
                    up <- res[which((res$log2FoldChange>=1) & (res$padj<=0.05)),]
                    # description
                    contrast_description <- paste0('(',paste(g1,collapse=" + "),")/",length(g1)," vs ",'(',paste(sort(g2),collapse=" + "),")/",length(g2))
                    contrast_description <- gsub('\\(pos\\)','\\+',contrast_description)
                    contrast_description <- gsub('\\(neg\\)','-',contrast_description)
                    # tissue description
                    tissue_description <- paste(names(table(coldata[,2])[1])," vs average of ",paste(names(table(coldata[,2])[-1]),collapse=","))
                    
                    exp_info <- matrix(c(paste0('Exp',k),i,contrast_description,tissue_description,dim(deg_res)[1],
                                        dim(up)[1],dim(deg_res)[1]-dim(up)[1],nrow(dds)-dim(deg_res)[1],
                                        paste(groups[which(groups[,2]==j),'external_id'],collapse = ','),
                                        paste(groups[which(groups[,2]!=j),'external_id'],collapse = ',')), nrow=1)
                    
                    write.table(exp_info, file="/localdisk/home/s2172876/project/DEA_Result/summary_3.txt", append=TRUE, sep = '\t',quote = FALSE,row.names = FALSE,col.names = FALSE)
                    write.table(res,file = gzfile(paste0('/localdisk/home/s2172876/project/DEA_Result/tissue/',exp_info[1],'.gz')),sep = '\t',quote = FALSE)
                    k = k+1
                }
            }
        }
    })
   
}
print('finish')
