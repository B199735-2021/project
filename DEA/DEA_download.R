library(recount3)
library(HDF5Array)
library(DESeq2)
projects <- rbind(
    recount3::available_projects("human"),
    recount3::available_projects("mouse")
)
study_id <- read.csv("/localdisk/home/s2172876/project/annotation/group.txt",sep = '\t')[1]
abstract <- read.csv('/localdisk/home/s2172876/project/DEA/recount3_selection_2022-05-26 17_28_49.csv')[,c(1,3)]
study_id <-  merge(study_id,abstract,by.x = 'study',by.y='project')
study_id <- study_id[row.names(unique(study_id[,c("study", "organism")])),]
for(i in study_id[,1]){
    try(rse_gene <- create_rse_manual(study_id[i,1],organism=study_id[i,2]))
    try(assays(rse_gene)$counts <- transform_counts(rse_gene))
    try(fn <- paste0("/localdisk/home/s2172876/project/Count/",study_id[i,1],".gz"))
    try(write.table(assays(rse_gene)$counts,gzfile(fn)))
}
for(i in unique(TCGA_meta$study)) {
    try(
    rse_gene <- recount3::create_rse_manual(
    project = i,
    project_home = "data_sources/tcga",
    organism = "human",
    annotation = "gencode_v26",
    type = "gene"))
    try(assays(rse_gene)$counts <- transform_counts(rse_gene))
    try(fn <- paste0("/localdisk/home/s2172876/project/Count/",i,".gz"))
    try(write.table(assays(rse_gene)$counts,gzfile(fn)))

}

for(i in 1:nrow(gtex_path)) {
    TCGA_meta<-str_match(gtex_path[i, 1],'gtex.gtex.(.*?).MD')[2]
    try(
    rse_gene <- recount3::create_rse_manual(
    project = gtex_meta,
    project_home = "data_sources/gtex",
    organism = "human",
    annotation = "gencode_v26",
    type = "gene"))
    try(assays(rse_gene)$counts <- transform_counts(rse_gene))
    try(fn <- paste0("/localdisk/home/s2172876/project/Count/",TCGA_meta,".gz"))
    try(write.table(assays(rse_gene)$counts,gzfile(fn)))

}
