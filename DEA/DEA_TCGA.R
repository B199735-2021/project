setwd('/localdisk/home/s2172876/project/DEA')
library(stringr)
library(plyr)
library(DESeq2)
library(edgeR)


TCGA_meta <- read.csv('../annotation/TCGA_meta.csv')
TCGA_meta <- TCGA_meta[,c('study','external_id','gdc_cases.project.name','cgc_sample_sample_type','gdc_cases.diagnoses.tumor_stage','gdc_cases.diagnoses.vital_status',
'gdc_cases.samples.sample_type','cgc_case_tumor_status','cgc_case_vital_status','cgc_case_pathologic_n','cgc_case_clinical_m','cgc_case_clinical_stage',
'cgc_case_pathologic_stage','cgc_case_pathologic_t','xml_primary_pathology_histological_type','gdc_cases.submitter_id')]
write.csv(TCGA_meta,'../annotation/TCGA_meta_subset.csv',quote = FALSE,row.names = FALSE)
abstract <- read.csv('/localdisk/home/s2172876/project/DEA/recount3_selection_2022-05-26 17_28_49.csv')[,c(1,3)]
# gene annotations
human_gene_annotation <- read.table('../DEA/human_gene_annotation.txt',sep = '\t',header = TRUE)
human_gene_annotation <- human_gene_annotation[which(human_gene_annotation$Gene.type=='protein_coding'),]

# grouping
# sample type
TCGA_sampletype_group <- aggregate(external_id ~ study + cgc_sample_sample_type, TCGA_meta, FUN = paste, collapse=",")
# delete groups have only one sample
TCGA_sampletype_group <- TCGA_sampletype_group[which(grepl( ',', TCGA_sampletype_group[,3], fixed = TRUE)),]
# delete projects have only one group
TCGA_sampletype_group <- TCGA_sampletype_group[-which(TCGA_sampletype_group$study %in% names(which(table(TCGA_sampletype_group[,1])==1))),]

# initiate experiment information dataframe
exp_info <- data.frame(matrix(ncol = 7))
colnames(exp_info) <- c('ExpID','ProjectID','Comparison','DEG','Up','Down','non_DEG')

k = 214

project_name <- names(sort(table(TCGA_meta[,1])))[which(names(sort(table(TCGA_meta[,1]))) %in% names(table(TCGA_sampletype_group[,1])))]
for(i in project_name[8:length(project_name)]){
    try({
    # annotation for this project
    groups <- TCGA_sampletype_group[which(TCGA_sampletype_group[,1]==i),]
    coldata<-data.frame()
    rn <- vector()
    sample_type <- vector()
    for(j in 1:nrow(groups)){
        rn <- c(rn,c(str_split(groups[j,3],',')[[1]]))
        sample_type <- c(sample_type,rep(groups[j,2],length(str_split(groups[j,3],',')[[1]])))
    }
    coldata <- data.frame(row.names=rn, sample_type)
    coldata <- merge(coldata,TCGA_meta[,c('external_id','gdc_cases.submitter_id')],by.x = 'row.names',by.y = 'external_id')
    rownames(coldata) <- coldata[,1]
    coldata <- coldata[,-1]
    coldata <- coldata[,c(2,1)]
    Expr <- read.table(gzfile(paste0('/localdisk/home/s2172876/project/Count/',i,'.gz')),check.names = FALSE,header = TRUE,stringsAsFactors = FALSE)
    comp <- combn(unique(coldata[,2]),2)
    for(j in 1:ncol(comp)){
        # construct subset of coldata
        coldata_sub <- as.data.frame(coldata[which(coldata[,2] %in% as.vector(comp[,j])),])
        rownames(coldata_sub) <- rownames(coldata)[which(coldata[,2] %in% as.vector(comp[,j]))]
        colnames(coldata_sub) <- colnames(coldata)
        # subset of the expression matrix
        Expr_sub <- Expr[,rownames(coldata_sub)]
        Expr_sub <- Expr_sub[which(sapply(str_split(row.names(Expr_sub), '\\.'),"[[",1) %in% human_gene_annotation$Gene.stable.ID),]
        # construct DEA object
        # The condition of interest should go at the end of the design formula
        dds <- DESeqDataSetFromMatrix(countData = Expr_sub,colData = coldata_sub,design = as.formula(paste('~0+',paste(colnames(coldata),collapse=' + '),sep = ' ')))
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
        sizeFactors(dds)

        # Defferential gene expression
        deg <- DESeq(dds)
        
        res <- results(deg, name = tail(resultsNames(deg[[1]]),n=1))
        # res <- results(deg, contrast = c(colnames(coldata_sub)[2], levels(factor(coldata_sub[,1]))[1], levels(factor(coldata_sub[,1]))[2]))
        res <- data.frame(res)
        # remove gene id version
        row.names(res) <- sapply(str_split(row.names(res), '\\.'),"[[",1)
        # round
        res[,1:4] <- round(res[,1:4],2)
        res[,5:6] <- signif(res[,5:6],2)
        exp_info[k,1] <- paste0('Exp',k)
        exp_info[k,2] <- i
        exp_info[k,3] <- paste(unique(coldata_sub[,2]),collapse = ' vs ')
        exp_info[k,4] <- dim(res[which((abs(res$log2FoldChange)>=1) & (res$padj<=0.05)),])[1]
        exp_info[k,5] <- dim(res[which((res$log2FoldChange)>=1 & (res$padj<=0.05)),])[1]
        exp_info[k,6] <- dim(res[which((res$log2FoldChange)<=(-1) & (res$padj<=0.05)),])[1]
        exp_info[k,7] <- nrow(dds) - exp_info[k,4]
        write.table(res,file = gzfile(paste0('/localdisk/home/s2172876/project/DEA_Result/TCGA_paired/',exp_info[k,1],'.gz')),sep = '\t',quote = FALSE)
        write.table(exp_info[k,], file="/localdisk/home/s2172876/project/DEA_Result/summary_4.txt", append=TRUE, sep = '\t',quote = FALSE,row.names = FALSE,col.names = FALSE)
        k = k+1
    }
    })
}

# # pathologic stage
# # sample type
# TCGA_pathologic_group <- aggregate(external_id ~ study + cgc_case_pathologic_stage, TCGA_meta, FUN = paste, collapse=",")
# TCGA_pathologic_group <- TCGA_pathologic_group[which(TCGA_pathologic_group$study!=''),]
# # delete groups have only one sample
# TCGA_pathologic_group <- TCGA_pathologic_group[which(grepl( ',', TCGA_pathologic_group[,3], fixed = TRUE)),]
# # delete projects have only one group
# TCGA_pathologic_group <- TCGA_pathologic_group[-which(TCGA_pathologic_group$study %in% names(which(table(TCGA_pathologic_group[,1])==1))),]
# k=1
# for(i in unique(TCGA_pathologic_group[,1])){
#         groups <- TCGA_pathologic_group[which(TCGA_pathologic_group[,1]==i),]
#         coldata<-data.frame()
#         rn <- vector()
#         sample_type <- vector()
#         for(j in 1:nrow(groups)){
#             rn <- c(rn,c(str_split(groups[j,3],',')[[1]]))
#             sample_type <- c(sample_type,rep(groups[j,2],length(str_split(groups[j,3],',')[[1]])))
#         }
#         coldata <- data.frame(row.names=rn, sample_type)
#         Expr <- read.table(gzfile(paste0('/localdisk/home/s2172876/project/Count/',i,'.gz')),check.names = FALSE,header = TRUE,stringsAsFactors = FALSE)
#         # subset of the expression matrix
#         Expr <- Expr[,rn]
#         #  filter protein coding gene
#         if(abstract[which(abstract[,2]==i),1] == 'mouse'){
#             Expr <- Expr[which(sapply(str_split(row.names(Expr), '\\.'),"[[",1) %in% mouse_gene_annotation$Gene.stable.ID),]
#         }else{
#             Expr <- Expr[which(sapply(str_split(row.names(Expr), '\\.'),"[[",1) %in% human_gene_annotation$Gene.stable.ID),]
#         }
#         # construct DEA object
#         dds <- DESeqDataSetFromMatrix(countData = Expr,colData = coldata,design = as.formula(paste('~0+',colnames(coldata),sep = ' ')))
#         # Apply filters to remove genes that do not pass quality filters
#         # first filter based upon counts
#         keep <- rowSums(counts(dds)) > 1
#         dds <- dds[keep,]
#         nrow(dds)
#         # second filter-require 25% in >= 2 samples
#         keep <- rowSums(cpm(dds) >= quantile(cpm(dds),obs=0.25)[2]) >= min(table(coldata))/2
#         dds <- dds[keep,]
#         nrow(dds)
#         # size factor estimation for normalisation
#         dds <- estimateSizeFactors(dds)
#         # Defferential gene expression
#         deg <- DESeq(dds)
#         # resultsNames(deg)
#         contrast_vector <- rep(-1/(length(unique(coldata[,1]))-1),length(unique(coldata[,1])))
#         coldata[1] <- gsub('\\(pos\\)','\\+',coldata[,1])
#         coldata[1] <- gsub('\\(neg\\)','-',coldata[,1])
        
#         for(j in 1:length(unique(coldata[,1]))){    
#             contrast_description <- (paste0(unique(coldata[,1])[j],' vs ','(',paste(unique(coldata[,1])[-j],collapse=" + "),')/',length(unique(coldata[,1]))-1))
#             contrast_vector[j] <- 1
#             res <- results(deg, contrast = contrast_vector)
#             res <- data.frame(res)
#             # remove gene id version
#             row.names(res) <- sapply(str_split(row.names(res), '\\.'),"[[",1)
#             # round
#             res[,1:4] <- round(res[,1:4],2)
#             res[,5:6] <- signif(res[,5:6],2)
#             # DEG filter
#             deg_res <- res[which((abs(res$log2FoldChange)>=1) & (res$padj<=0.05)),]
#             up <- res[which((res$log2FoldChange>=1) & (res$padj<=0.05)),]
#             contrast_vector <- rep(-1/(length(unique(coldata[,1]))-1),length(unique(coldata[,1])))
#             exp_info <- matrix(c(paste0('Exp',k),i,contrast_description,dim(deg_res)[1],dim(up)[1],
#                                         dim(deg_res)[1]-dim(up)[1],nrow(dds)-dim(deg_res)[1],
#                                         groups[j,3],
#                                         paste(apply(groups[-j,c(2,3)], 1 , paste , collapse = ":" ),collapse = ';')), nrow=1)
#             write.table(exp_info, file="/localdisk/home/s2172876/project/DEA_Result/summary_4.txt", append=TRUE, sep = '\t',quote = FALSE,row.names = FALSE,col.names = FALSE)
#             write.table(res,file = gzfile(paste0('/localdisk/home/s2172876/project/DEA_Result/TCGA/',paste0('Exp',k),'.gz')),sep = '\t',quote = FALSE)
#             # print(exp_info)
#             k = k+1
#         }
# }


# # grouping
# # histological type
# TCGA_histological_group <- aggregate(external_id ~ study + xml_primary_pathology_histological_type, TCGA_meta, FUN = paste, collapse=",")
# TCGA_histological_group <- TCGA_histological_group[which(TCGA_histological_group$study!=''),]

# # delete groups have only one sample
# TCGA_histological_group <- TCGA_histological_group[which(grepl( ',', TCGA_histological_group[,3], fixed = TRUE)),]
# # remove groups with 'yes' and 'no' as histological type
# TCGA_histological_group <- TCGA_histological_group[-which(TCGA_histological_group[,2] =='YES'|TCGA_histological_group[,2]=='NO'),]

# for(i in unique(TCGA_histological_group[,1])){
#     try({
#         groups <- TCGA_histological_group[which(TCGA_histological_group[,1]==i),]
#         coldata<-data.frame()
#         rn <- vector()
#         sample_type <- vector()
#         for(j in 1:nrow(groups)){
#             rn <- c(rn,c(str_split(groups[j,3],',')[[1]]))
#             sample_type <- c(sample_type,rep(groups[j,2],length(str_split(groups[j,3],',')[[1]])))
#         }
#         coldata <- data.frame(row.names=rn, sample_type)
#         Expr <- read.table(gzfile(paste0('/localdisk/home/s2172876/project/Count/',i,'.gz')),check.names = FALSE,header = TRUE,stringsAsFactors = FALSE)
#         # subset of the expression matrix
#         Expr <- Expr[,rn]
#         #  filter protein coding gene
#         if(abstract[which(abstract[,2]==i),1] == 'mouse'){
#             Expr <- Expr[which(sapply(str_split(row.names(Expr), '\\.'),"[[",1) %in% mouse_gene_annotation$Gene.stable.ID),]
#         }else{
#             Expr <- Expr[which(sapply(str_split(row.names(Expr), '\\.'),"[[",1) %in% human_gene_annotation$Gene.stable.ID),]
#         }
#         if(length(unique(coldata[,1]))==2){ # two level
#             # construct DEA object
#             dds <- DESeqDataSetFromMatrix(countData = Expr,colData = coldata,design = as.formula(paste('~',colnames(coldata),sep = ' ')))
#             # Apply filters to remove genes that do not pass quality filters
#             # first filter based upon counts
#             keep <- rowSums(counts(dds)) > 1
#             dds <- dds[keep,]
#             nrow(dds)
#             # second filter-require 25% in >= 2 samples
#             library(edgeR)
#             keep <- rowSums(cpm(dds) >= quantile(cpm(dds),obs=0.25)[2]) >= min(table(coldata[,1]))/2
#             dds <- dds[keep,]
#             # size factor estimation for normalisation
#             dds <- estimateSizeFactors(dds)

#             # Defferential gene expression
#             deg <- DESeq(dds)
#             res <- results(deg, contrast = c(colnames(coldata), levels(factor(coldata[,1]))[1], levels(factor(coldata[,1]))[2]))        
#             info <- mcols(res, use.names=TRUE)
#             res <- data.frame(res)
#             # remove gene id version
#             row.names(res) <- sapply(str_split(row.names(res), '\\.'),"[[",1)
#             # round
#             res[,1:4] <- round(res[,1:4],2)
#             res[,5:6] <- signif(res[,5:6],2)
#             contrast_description <- word(str_match(info[2,2],"\\(MLE\\)\\: (.*?)$")[2],2,-1)
#             # DEG filter
#             deg_res <- res[which((abs(res$log2FoldChange)>=1) & (res$padj<=0.05)),]
#             up <- res[which((res$log2FoldChange>=1) & (res$padj<=0.05)),]
#             exp_info <- matrix(c(paste0('Exp',k),i,contrast_description,dim(deg_res)[1],
#                                 dim(up)[1],dim(deg_res)[1]-dim(up)[1],nrow(dds)-dim(deg_res)[1],
#                                 groups[which(groups[,2]==str_match(contrast_description,'(.*?) vs (.*?)$')[2]),3],
#                                 groups[which(groups[,2]==str_match(contrast_description,'(.*?) vs (.*?)$')[3]),3]), nrow=1)
#             write.table(res,file = gzfile(paste0('/localdisk/home/s2172876/project/DEA_Result/TCGA/',exp_info[1],'.gz')),sep = '\t',quote = FALSE)
#             write.table(exp_info, file="/localdisk/home/s2172876/project/DEA_Result/summary_4.txt", append=TRUE, sep = '\t',quote = FALSE,row.names = FALSE,col.names = FALSE)
#             k = k+1
#         }

#         if(length(unique(coldata[,1]))>2){
#             # construct DEA object
#             dds <- DESeqDataSetFromMatrix(countData = Expr,colData = coldata,design = as.formula(paste('~0+',colnames(coldata),sep = ' ')))
#             # Apply filters to remove genes that do not pass quality filters
#             # first filter based upon counts
#             keep <- rowSums(counts(dds)) > 1
#             dds <- dds[keep,]
#             nrow(dds)
#             # second filter-require 25% in >= 2 samples
#             keep <- rowSums(cpm(dds) >= quantile(cpm(dds),obs=0.25)[2]) >= min(table(coldata))/2
#             dds <- dds[keep,]
#             nrow(dds)
            
#             # Defferential gene expression
#             deg <- DESeq(dds)
#             # resultsNames(deg)
#             contrast_vector <- rep(-1/(length(unique(coldata[,1]))-1),length(unique(coldata[,1])))
#             coldata[1] <- gsub('\\(pos\\)','\\+',coldata[,1])
#             coldata[1] <- gsub('\\(neg\\)','-',coldata[,1])
#             for(j in 1:length(unique(coldata[,1]))){    
#                 contrast_description <- (paste0(unique(coldata[,1])[j],' vs ','(',paste(unique(coldata[,1])[-j],collapse=" + "),')/',length(unique(coldata[,1]))-1))
#                 contrast_vector[j] <- 1
#                 res <- results(deg, contrast = contrast_vector)
#                 res <- data.frame(res)
#                 # remove gene id version
#                 row.names(res) <- sapply(str_split(row.names(res), '\\.'),"[[",1)
#                 # round
#                 res[,1:4] <- round(res[,1:4],2)
#                 res[,5:6] <- signif(res[,5:6],2)
#                 # DEG filter
#                 deg_res <- res[which((abs(res$log2FoldChange)>=1) & (res$padj<=0.05)),]
#                 up <- res[which((res$log2FoldChange>=1) & (res$padj<=0.05)),]
#                 contrast_vector <- rep(-1/(length(unique(coldata[,1]))-1),length(unique(coldata[,1])))
#                 exp_info <- matrix(c(paste0('Exp',k),i,contrast_description,dim(deg_res)[1],
#                                     dim(up)[1],dim(deg_res)[1]-dim(up)[1],nrow(dds)-dim(deg_res)[1],
#                                     groups[which(groups[,2]==unique(coldata[,1])[j]),3],
#                                     paste(apply(groups[which(groups[,2]==unique(coldata[,1])[j]),c(2,3)] , 1 , paste , collapse = ":" ),collapse = ';')), nrow=1)
#                 write.table(exp_info, file="/localdisk/home/s2172876/project/DEA_Result/summary_4.txt", append=TRUE, sep = '\t',quote = FALSE,row.names = FALSE,col.names = FALSE)
#                 write.table(res,file = gzfile(paste0('/localdisk/home/s2172876/project/DEA_Result/TCGA/',paste0('Exp',k),'.gz')),sep = '\t',quote = FALSE)
#                 print(exp_info)
#                 k = k+1
#             }
#         }
#     })
#     }



