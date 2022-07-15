setwd('/localdisk/home/s2172876/project/DEA')
library(tximeta)
# import tcga and gtex tissue annotation
tcga_gtex_tissue <- read.csv('../annotation/TCGA_GTEX_TISSUE.txt',sep = '\t')
# filter out unpaired gtex and tcga project
tcga_gtex_tissue <- tcga_gtex_tissue[which(tcga_gtex_tissue[,3]!=''),]
# import annotation file (with sampleid)
tcga_gtex_annotation <- read.table(gzfile('../annotation/TcgaTargetGTEX_phenotype.txt.gz'),sep = '\t',header = TRUE)
tcga_gtex_annotation <- tcga_gtex_annotation[which(tcga_gtex_annotation[,2] %in% tcga_gtex_tissue[,1]),]
# import all samples expression data
tcga_gtex <- read.table(gzfile('../Count/TcgaTargetGtex_gene_expected_count.gz'),sep = '\t',check.names=TRUE,header = TRUE,row.names = 1)
# filter out unpaired gtex and tcga samples
colnames(tcga_gtex) <- gsub('\\.','-',colnames(tcga_gtex))
tcga_gtex <- tcga_gtex[,which(colnames(tcga_gtex) %in% tcga_gtex_annotation[,1])]
rownames(tcga_gtex) <- sapply(str_split(row.names(tcga_gtex), '\\.'),"[[",1)

# gene annotations
human_gene_annotation <- read.table('human_gene_annotation.txt',sep = '\t',header = TRUE)
human_gene_annotation <- human_gene_annotation[which(human_gene_annotation$Gene.type=='protein_coding'),]

# protein coding gene
tcga_gtex <- tcga_gtex[which(rownames(tcga_gtex) %in% human_gene_annotation[,3]),]
write.table(tcga_gtex,gzfile('../Count/TcgaGtex_gene_expected_count.gz'))

save.image(file = "/localdisk/home/s2172876/project/DEA/TCGA_GTEx.RData")

library(tximeta)

library(tximport)
txi.rsem <- tximport(gzfile('../Count/TcgaGtex_gene_expected_count.gz'), type = "rsem")
sampleTable <- data.frame(condition = factor(rep(c('A', 'B','C','D'), each = 3)))
dds <- DESeqDataSetFromTximport(txi.rsem, sampleTable, ~condition)


library(tximportData)
dir <- system.file("extdata", package = "tximportData")
list.files(dir)
samples <- read.table(file.path(dir, "samples.txt"), header = TRUE)
samples
files <- file.path(dir, "salmon", samples$run, "quant.sf.gz")
names(files) <- paste0("sample", 1:6)
all(file.exists(files))
files <- file.path(dir, "rsem", samples$run, paste0(samples$run, ".genes.results.gz"))
names(files) <- paste0("sample", 1:6)
txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)
head(txi.rsem$counts)
