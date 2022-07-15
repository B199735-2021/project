setwd('/localdisk/home/s2172876/project/DEA')
library(stringr)
library(plyr)
library(DESeq2)
library(edgeR)
library(mjcbase)
library(sva)
library(data.table)
gtex <- read.csv('../annotation/gtex_meta.csv',header = TRUE, check.names = FALSE)
gtex <- gtex[,c(2,3,4,5,7)]

# comparison_df <- data.frame(matrix(ncol = 3))
# colnames(comparison_df) <- c('tissue-specific','tissue IDs','other IDs')
# # define comparison
# k = 1
# for(i in names(sort(table(gtex$study)))){
#     comparison_df[k,] <- c(i,paste(gtex[which(gtex$study == i),2],collapse = ','),
#                             paste(sample(gtex[which(gtex$study != i),2],length(gtex[which(gtex$study == i),2])),collapse = ','))
#     k = k+1
# }
# write.table(comparison_df,'GTEx_comparison.txt',sep = '\t',quote = FALSE,row.names = FALSE)
comparison_df <- read.table('GTEx_comparison.txt',sep = '\t',header = 1,check.names = FALSE)

# counter
k = 1
# DEG analysis
sn <- scan(gzfile('../Count/GTEx.gz'), 
        what = "", nlines = 1)

for(i in 1:dim(comparison_df)[1]){
    try({
        Expr_g1 <- data.frame(fread("../Count/GTEx.gz", select = c(1,which(sn %in% str_split(comparison_df[i,2],',')[[1]])+1), check.names=FALSE), row.names=1, check.names=FALSE)
        Expr_g2 <- data.frame(fread("../Count/GTEx.gz", select = c(1,which(sn %in% str_split(comparison_df[i,3],',')[[1]])+1) , check.names=FALSE), row.names=1, check.names=FALSE)
        Expr <- merge(Expr_g1,Expr_g2,by = 'row.names')
        rownames(Expr) <- Expr[,1]
        Expr <- Expr[,-1]
        Expr <- Expr[,c(colnames(Expr_g1),colnames(Expr_g2))]
        Expr$Gene <-sapply(str_split(row.names(Expr), '\\.'),"[[",1)
        Expr <- aggregate(. ~ Gene, Expr, sum)
        row.names(Expr) <- Expr[,1]
        Expr <- Expr[,-1]
        # initialize colData
        coldata <- data.frame(row.names = colnames(Expr),
                                tissue = c(rep(comparison_df[i,1],length(str_split(comparison_df[i,2],',')[[1]])),
                                            rep('other',length(str_split(comparison_df[i,3],',')[[1]]))))
        coldata <- merge(coldata,gtex[,c(2,5)],by.x = 'row.names',by.y = 'external_id')
        rownames(coldata) <- coldata[,1]
        coldata <- coldata[,-1]
        coldata <- coldata[colnames(Expr),]
        # the following code was adapted from https://biodatascience.github.io/compbio/dist/sva.html
        t1 <- Sys.time()
        # construct DEA object
        dds <- DESeqDataSetFromMatrix(countData = Expr,colData = coldata,design = ~tissue)
        dds$batch <- factor(dds$AGE)
        table(dds$tissue, dds$batch)

        # estimate the library size correction
        dds <- estimateSizeFactors(dds)
        norm.cts <- counts(dds, normalized=TRUE)

        # the biological condition model
        mm <- model.matrix(~ tissue, colData(dds))
        # null model
        mm0 <- model.matrix(~ 1, colData(dds))
        norm.cts <- norm.cts[rowSums(norm.cts) > 0,]
        fit <- svaseq(norm.cts, mod=mm, mod0=mm0, n.sv=2)

        ddssva <- dds
        ddssva$SV1 <- fit$sv[,1]
        ddssva$SV2 <- fit$sv[,2]
        design(ddssva) <- ~ SV1 + SV2 + tissue
        deg <- DESeq(ddssva)
        res <- results(deg,parallel=TRUE)
        t2 <- Sys.time()
        # deg_time <- t2-t1
        deg <- data.frame(res)

        # round
        deg[,1:4] <- round(deg[,1:4],2)
        deg[,5:6] <- signif(deg[,5:6],2)

        contrast_description <- word(str_match(mcols(res, use.names=TRUE)[2,2],"\\(MLE\\)\\: (.*?)$")[2],2,-1)
        # DEG filter
        deg_res <- res[which((abs(res$log2FoldChange)>=1) & (res$padj<=0.05)),]
        up <- res[which((res$log2FoldChange>=1) & (res$padj<=0.05)),]
        contrast_vector <- rep(-1/(length(unique(coldata[,1]))-1),length(unique(coldata[,1])))
        if(word(contrast_description,1)=='other'){
            exp_info <- matrix(c(paste0('Exp',k),'GTEx',contrast_description,dim(deg_res)[1],dim(up)[1],
                            dim(deg_res)[1]-dim(up)[1],nrow(dds)-dim(deg_res)[1],
                            comparison_df[i,3],comparison_df[i,2]), nrow=1)
        }else{
            exp_info <- matrix(c(paste0('Exp',k),'GTEx',contrast_description,dim(deg_res)[1],dim(up)[1],
                            dim(deg_res)[1]-dim(up)[1],nrow(dds)-dim(deg_res)[1],
                            comparison_df[i,2],comparison_df[i,3]), nrow=1)
        }                
        write.table(exp_info, file="/localdisk/home/s2172876/project/DEA_Result/summary_5.txt", append=TRUE, sep = '\t',quote = FALSE,row.names = FALSE,col.names = FALSE)
        write.table(res,file = gzfile(paste0('/localdisk/home/s2172876/project/DEA_Result/GTEx/',exp_info[1],'.gz')),sep = '\t',quote = FALSE)
        k = k+1
    })
}
