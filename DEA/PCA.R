setwd('/localdisk/home/s2172876/project/DEA')

library(stringr)
library(DESeq2)
library(edgeR)
library(factoextra)
library(ggpubr)
library(gridExtra)
# gene annotations
human_gene_annotation <- read.table('human_gene_annotation.txt',sep = '\t',header = TRUE)
human_gene_annotation <- human_gene_annotation[which(human_gene_annotation$Gene.type=='protein_coding'),]
mouse_gene_annotation <- read.table('mouse_gene_annotation.txt',sep = ',',header = TRUE)
mouse_gene_annotation <- mouse_gene_annotation[which(mouse_gene_annotation$Gene.type=='protein_coding'),]

experiment <- read.csv('../annotation/group.txt',sep = '\t')

Expr_clean <- function(i,species){
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
            
    # import expression matrix
    Expr <- read.table(paste0('/localdisk/home/s2172876/project/Count/',i,'.gz'),stringsAsFactors = FALSE)
    g1 <- str_split(groups[1,10],',')[[1]]
    g2 <- str_split(groups[2,10],',')[[1]]
    # subset of the expression matrix
    Expr <- Expr[,c(g1,g2)]
    #  filter protein coding gene
    if(species == 'mouse'){
            Expr <- Expr[which(sapply(str_split(row.names(Expr), '\\.'),"[[",1) %in% mouse_gene_annotation$Gene.stable.ID),]
        }else{
            Expr <- Expr[which(sapply(str_split(row.names(Expr), '\\.'),"[[",1) %in% human_gene_annotation$Gene.stable.ID),]
        }
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

    Expr_PCA <- counts(dds)
    return(list(Expr_PCA,coldata))
}

SRP043391 <- Expr_clean('SRP043391','human')
ERP011264 <- Expr_clean('ERP011264','human')
SRP034831 <- Expr_clean('SRP034831','human')
SRP012147 <- Expr_clean('SRP012147','mouse')

res.pca <- prcomp(t(SRP043391[[1]]), scale = TRUE)
p1 <- fviz_pca_ind(res.pca,col.ind = SRP043391[[2]][,1],
            addEllipses = TRUE, ellipse.type = "confidence",
            legend.title = "Treatment",repel = TRUE,labelsize=3) +
            labs(title ="SRP043391", x = "PC1", y = "PC2")+
            theme(legend.position = "bottom",plot.title = element_text(hjust = 0.5),
                legend.text = element_text(size = 12),
                axis.text = element_text(size = 12))
p1

res.pca <- prcomp(t(ERP011264[[1]]), scale = TRUE)
p2 <- fviz_pca_ind(res.pca,col.ind = ERP011264[[2]][,1],
            addEllipses = TRUE, ellipse.type = "confidence",
            legend.title = "Treatment",repel = TRUE,labelsize=3) +
            labs(title ="ERP011264", x = "PC1", y = "PC2")+
            theme(legend.position = "bottom",plot.title = element_text(hjust = 0.5),
                legend.text = element_text(size = 12),
                axis.text = element_text(size = 12))
p2

res.pca <- prcomp(t(SRP034831[[1]]), scale = TRUE)
p3 <- fviz_pca_ind(res.pca,col.ind = SRP034831[[2]][,1],
            addEllipses = TRUE, ellipse.type = "confidence",
            legend.title = "Treatment",repel = TRUE,labelsize=3) +
            labs(title ="SRP034831", x = "PC1", y = "PC2")+
            theme(legend.position = "bottom",plot.title = element_text(hjust = 0.5),
                legend.text = element_text(size = 12),
                axis.text = element_text(size = 12))
p3

res.pca <- prcomp(t(SRP012147[[1]]), scale = TRUE)
p4 <- fviz_pca_ind(res.pca,col.ind = SRP012147[[2]][,1],
            addEllipses = TRUE, ellipse.type = "confidence",
            legend.title = "Treatment",,repel = TRUE,labelsize=3) +
            labs(title ="SRP012147", x = "PC1", y = "PC2")+
            theme(legend.position = "bottom",plot.title = element_text(hjust = 0.5),
                legend.text = element_text(size = 12),
                axis.text = element_text(size = 12))
p4
pdf('PCA.pdf',width=13)
grid.arrange(p1, p2, p3, p4, nrow = 2, ncol =2)
dev.off()
