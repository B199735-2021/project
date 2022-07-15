setwd('/localdisk/home/s2172876/project/annotation')
library(stringr)
library(plyr)
library(DESeq2)
library(edgeR)
#--------------------------------------------------------------------------------------#
# Step 1. import data
#--------------------------------------------------------------------------------------#
# load project description
pro_description <- read.csv("recount3_selection_2022-05-11 21_54_01.csv")
# tega experiment
tcga_project <- pro_description[pro_description$project_home == "data_sources/tcga", ]

# load file path
path <- read.csv("recount3_metadata_files.csv")
# geo sample file path
tcga_path <- path[grepl("/tcga/", path$project_meta, fixed = TRUE), ]

#--------------------------------------------------------------------------------------#
# Step 2. generate sample description dataframe
#--------------------------------------------------------------------------------------#
# define pre_annotation dataframe
pre_anno <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(pre_anno) <- c(
    "Project ID", "External ID",
    "Sample attributes"
)
TCGA_meta <- data.frame()
for (i in 1:nrow(tcga_path)) {
    samp_url <- tcga_path[i, 1]
    samp_id <- str_match(tcga_path[i, 3],'recount_qc.(.*?).MD')[2]
    download_command <- paste("wget -O", samp_id, samp_url)
    try(system(download_command), silent = TRUE)
    # Did it download successfully?
    if (samp_id %in% list.files()) {
        pre_meta_temp <- read.csv(gzfile(samp_id), sep = "\t")
        pre_meta_temp <- pre_meta_temp[,-which(colSums(is.na(pre_meta_temp))==nrow(pre_meta_temp))]
        # Delete file
        system(paste0("rm ", samp_id))
        # try(do.call("<-",list(str_match(tcga_path[i, 3],'recount_qc.(.*?).MD')[2], pre_meta_temp)))
        try(TCGA_meta<-rbind.fill(TCGA_meta,pre_meta_temp))
    }
}
write.csv(TCGA_meta,'TCGA_meta.csv',row.names = FALSE, quote = FALSE)
# TCGA_meta <- TCGA_meta[,c('study','tcga_barcode','gdc_cases.project.name','cgc_sample_sample_type','gdc_cases.diagnoses.tumor_stage','gdc_cases.diagnoses.vital_status',
# 'gdc_cases.samples.sample_type','cgc_case_tumor_status','cgc_case_vital_status','cgc_case_pathologic_n','cgc_case_clinical_m',
# 'cgc_case_pathologic_stage','cgc_case_pathologic_t','cgc_drug_therapy_drug_name','cgc_drug_therapy_pharmaceutical_therapy_type')]


# GTEx
#--------------------------------------------------------------------------------------#
# Step 1. import data
#--------------------------------------------------------------------------------------#
# load project description
pro_description <- read.csv("recount3_selection_2022-05-11 21_54_01.csv")
# tega experiment
gtex_project <- pro_description[pro_description$project_home == "data_sources/gtex", ]

# load file path
path <- read.csv("recount3_metadata_files.csv")
# geo sample file path
gtex_path <- path[grepl("/gtex/", path$project_meta, fixed = TRUE), ]

#--------------------------------------------------------------------------------------#
# Step 2. generate sample description dataframe
#--------------------------------------------------------------------------------------#
# define pre_annotation dataframe
pre_anno <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(pre_anno) <- c(
    "Project ID", "External ID",
    "Sample attributes"
)

setdiff(as.vector(str_match(gtex_path[,3],'recount_qc.(.*?).MD')[,2]),names(table(gtex[,3])))
which(as.vector(str_match(gtex_path[,3],'recount_qc.(.*?).MD')[,2])=='OVARY')
gtex_meta <- data.frame()
for (i in 1:nrow(gtex_path)) {
    samp_url <- gtex_path[i, 1]
    samp_id <- str_match(gtex_path[i, 3],'recount_qc.(.*?).MD')[2]
    download_command <- paste("wget -O", samp_id, samp_url)
    system(download_command)
    pre_meta_temp <- read.csv(gzfile(samp_id), sep = "\t")
    pre_meta_temp <- pre_meta_temp[,-which(colSums(is.na(pre_meta_temp))==nrow(pre_meta_temp))]
    # Delete file
    system(paste0("rm ", samp_id))
    # try(do.call("<-",list(str_match(gtex_path[i, 3],'recount_qc.(.*?).MD')[2], pre_meta_temp)))
    gtex_meta<-rbind.fill(gtex_meta,pre_meta_temp)
}

write.csv(gtex_meta,'gtex_meta.csv',row.names = FALSE)
tcga_gtex <- read.table(gzfile('TcgaTargetGTEX_phenotype.txt.gz'),sep = '\t',header = TRUE)
tcga_gtex <- tcga_gtex[which(tcga_gtex$X_study=='TCGA'| tcga_gtex$X_study =='GTEX'),]
write.table(unique(tcga_gtex[,c(2,7)]),'TCGA_GTEX_TISSUE.txt',sep = '\t',quote = FALSE,row.names = FALSE,col.names = FALSE)
