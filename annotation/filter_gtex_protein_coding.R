setwd('/localdisk/home/s2172876/project/annotation')
gtex <- read.csv('gtex_meta.csv',header = TRUE, check.names = FALSE)
setwd('/localdisk/home/s2172876/project/Count')
# gene annotations
human_gene_annotation <- read.table('../DEA/human_gene_annotation.txt',sep = '\t',header = TRUE)
human_gene_annotation <- human_gene_annotation[which(human_gene_annotation$Gene.type=='protein_coding'),]

for(i in unique(gtex[,4])){
    Expr <- read.table(gzfile(paste(i,'.gz',sep = '')),check.names = FALSE)
    Expr <- Expr[which(sapply(str_split(row.names(Expr), '\\.'),"[[",1) %in% human_gene_annotation$Gene.stable.ID),]
    print(dim(Expr))
    write.table(Expr,gzfile(paste(i,'.gz',sep = '')),sep = '\t',quote = FALSE)
}

library(tibble)
i = unique(gtex[1,4])
Expr <- read.table(gzfile(paste(i,'.gz',sep = '')),check.names = FALSE)

for(i in unique(gtex[,4])[20:32]){
    temp_Expr <- read.table(gzfile(paste(i,'.gz',sep = '')),check.names = FALSE)
    Expr <- merge(Expr,temp_Expr,by = 'row.names')
    rownames(Expr) <- Expr[,1]
    Expr <- Expr[,-1]
    print(i)
}
write.table(Expr,gzfile('GTEx.gz'),sep = '\t',quote = FALSE)
