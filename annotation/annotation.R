setwd('/localdisk/home/s2172876/project/annotation')
#--------------------------------------------------------------------------------------#
# Step 1. exclude single cell RNA-Seq
#--------------------------------------------------------------------------------------#
# load project description
pro_description <- read.csv("recount3_selection_2022-05-11 21_54_01.csv")
# find out single cell project
single_cell_row <- (unique(c(
    which(grepl("single cell", pro_description$study_abstract, ignore.case = TRUE)),
    which(grepl("single-cell", pro_description$study_abstract, ignore.case = TRUE)),
    which(grepl("scrna-seq", pro_description$study_abstract, ignore.case = TRUE)),
    which(grepl("scrna", pro_description$study_abstract, ignore.case = TRUE)),
    which(grepl("single cell", pro_description$study_title, ignore.case = TRUE)),
    which(grepl("single-cell", pro_description$study_title, ignore.case = TRUE)),
    which(grepl("scrna-seq", pro_description$study_title, ignore.case = TRUE)),
    which(grepl("scrna", pro_description$study_title, ignore.case = TRUE))
)))
# how many projects are single cell RNA-Seq?
length(single_cell_row) ## 412
# exclude single cell project
pro_description <- pro_description[-single_cell_row, ]
# write.csv(pro_description,'recount3_nonSC.csv',row.names=FALSE)
table(pro_description$project_home)
# geo experiment
geo_project <- pro_description[pro_description$project_home == "data_sources/sra", ]

# load file path
path <- read.csv("recount3_metadata_files.csv")
# geo sample file path
geo_path <- path[grepl("/sra/", path$project_meta, fixed = TRUE), ]
# extract project ID
library(stringr)
# geo_path$ProjectID <- str_match(geo_path$project_meta, "sra\\.sra\\.(.*?)\\.MD\\.gz")[, 2]
# exclude scRNA
geo_path <- geo_path[which(geo_path$ProjectID %in% geo_project$project), ]

#--------------------------------------------------------------------------------------#
# Step 2. generate sample description dataframe
#--------------------------------------------------------------------------------------#
# define pre_annotation dataframe
pre_anno <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(pre_anno) <- c(
    "Project ID", "External ID",
    "Sample attributes"
)
for (i in 1:nrow(geo_path)) {
    print(paste0("Downloading metadata for ", samp_id))
    samp_url <- geo_path[i, 1]
    samp_id <- geo_path[i, 6]
    download_command <- paste("wget -O", samp_id, samp_url)
    try(system(download_command), silent = TRUE)
    # Did it download successfully?
    if (samp_id %in% list.files()) {
        pre_meta_temp <- read.csv(gzfile(samp_id), sep = "\t")
        # Delete file
        system(paste0("rm ", samp_id))
        try(pre_anno <- rbind(pre_anno, pre_meta_temp[, c("study", "external_id", "sample_attributes")]), silent = TRUE)
    }
}
pre_anno[,3] <- str_replace_all(pre_anno[,3],'Barrett\\?\\?\\?','Barrett\'s')
write.table(pre_anno, "pre_re_annotation.csv", sep = "\t", row.names = FALSE, quote = FALSE)
dim(pre_anno)

#--------------------------------------------------------------------------------------#
# Step 3. Keyword Extraction
#--------------------------------------------------------------------------------------#

# Tissue
Sys.time()
tissue <- str_match((pre_anno$sample_attributes), "(?i).*(Tissue-type|tissue type|tissue|tissue_type|organism part);;(.*?)(\\||$)")[, 3]
Sys.time()
save.image(file = "/localdisk/home/s2172876/project/annotation/Annotation.RData")

# Cell line
Sys.time()
cellline <- str_match(pre_anno$sample_attributes, "(?i).*(cell line|cell_line|cell line background);;(.*?)(\\||$)")[, 3]
Sys.time()
save.image(file = "/localdisk/home/s2172876/project/annotation/Annotation.RData")

# Cell type
Sys.time()
celltype <- str_match(pre_anno$sample_attributes, "(?i).*(cell type|cell_type|CellType);;(.*?)(\\||$)")[, 3]
Sys.time()
save.image(file = "/localdisk/home/s2172876/project/annotation/Annotation.RData")

# Condition
Sys.time()
condition <- str_match(pre_anno$sample_attributes, "(?i).*(status|disease staging|disease status|disease|disease state|infection agent|infection|health state|mutational status|health state|histological type|sample comment|Description|developmental stage|condition);;(.*?)(\\||$)")[, 3]
Sys.time()
save.image(file = "/localdisk/home/s2172876/project/annotation/Annotation.RData")

# treatment
Sys.time()
treatment <- str_match(pre_anno$sample_attributes, "(?i).*(agent|treatment|treated with|intervention);;(.*?)(\\||$)")[, 3]
Sys.time()
save.image(file = "/localdisk/home/s2172876/project/annotation/Annotation.RData")

# time/dose
Sys.time()
treatment_duration <- str_match(pre_anno$sample_attributes, "(?i).*(dose|time|intervention duration);;(.*?)(\\||$)")[, 3]
Sys.time()
save.image(file = "/localdisk/home/s2172876/project/annotation/Annotation.RData")

# genotype
Sys.time()
genotype <- str_match(pre_anno$sample_attributes, "(?i).*(genotype|genotype/variation);;(.*?)(\\||$)")[, 3]
Sys.time()
save.image(file = "/localdisk/home/s2172876/project/annotation/Annotation.RData")

# subtype
Sys.time()
subtype <- str_match(pre_anno$sample_attributes, "(?i).*(Subtype|disease_stage);;(.*?)(\\||$)")[, 3]
Sys.time()
save.image(file = "/localdisk/home/s2172876/project/annotation/Annotation.RData")

# Combine all columns
sample_annotation <- cbind(pre_anno$study, pre_anno$external_id, tissue, cellline, celltype, condition, treatment, treatment_duration,genotype,subtype)
colnames(sample_annotation)[1:2] <- c('study','external_id')
# remove donor information
sample_annotation[,3] <- str_replace_all(sample_annotation[,3],"(?i)donor(#)*( )*(_)*\\d+(;)*","")
sample_annotation[,3] <- str_replace_all(sample_annotation[,3],"(?i)(\\.|\\_|\\#|, repeat |rep|replicate \\#|Control-rep|control replicate|control |control)\\d{1,5}(\\.){0,1}$","")
sample_annotation[,3] <- str_replace_all(sample_annotation[,3] ,"(?i)[Patient|Donor|rep]( ){0,2}\\d{1,5}( ){0,2}","")
sample_annotation[,3] <- str_replace_all(sample_annotation[,3],"(?i)Repeat [A-Z]\\+[A-Z]","Repeat")
sample_annotation[,3] <- str_replace_all(sample_annotation[,3],"(?i) [A-Ea-e]$","")
sample_annotation[,3] <- str_replace_all(sample_annotation[,3],"(?i)(\\(){0,1}Replicate(\\_|-#| #|#)[0-9]{1,3}(\\)){0,1}$","")
sample_annotation[,3] <- str_replace_all(sample_annotation[,3],"(?i) \\({0,1}Replicate [0-9]{1,3}\\){0,1}$","")
sample_annotation[,3] <- str_replace_all(sample_annotation[,3],"(?i)KiaZcKO-\\d$","KiaZcKO")
sample_annotation[,3] <- str_replace_all(sample_annotation[,3],"(?i)control\\_\\d","control")
sample_annotation[,3] <- str_replace_all(sample_annotation[,3],"(?i)control\\-\\d$","control")
sample_annotation[,3] <- str_replace_all(sample_annotation[,3],"(?i)KiaZcKO\\-\\d$","KiaZcKO")
sample_annotation[,4] <- str_replace_all(sample_annotation[,4],"(?i)Biological replicate \\d{1,5}","")
sample_annotation[,5] <- str_replace_all(sample_annotation[,5],"donor(#)*( )*(_)*\\d+(;)*","")
sample_annotation[,5] <- str_replace_all(sample_annotation[,5],"(?i)replicate \\d$","")
sample_annotation[,6] <- str_replace_all(sample_annotation[,6],"donor(#)*( )*(_)*\\d+(;)*","")
sample_annotation[,6] <- str_replace_all(sample_annotation[,6],"donation(#)*( )*(_)*\\d+(;)*","")
sample_annotation[,6] <- str_replace_all(sample_annotation[,6],"(?i).Biological Replicate \\d$","")
sample_annotation[,6] <- str_replace_all(sample_annotation[,6],"(?i).Technical Replicate \\d$","")
sample_annotation[,6] <- str_replace_all(sample_annotation[,6],"(?i).Replicate( )*\\d$","")
sample_annotation[,6] <- str_replace_all(sample_annotation[,6],"(?i)sequence replicate( )*\\d.","")
sample_annotation[,6] <- str_replace_all(sample_annotation[,6],"(?i)of participant( )*\\d{1,5}","")
sample_annotation[,6] <- str_replace_all(sample_annotation[,6],"(?i), ([a-zA-Z])+ mouse","")
sample_annotation[,6] <- str_replace_all(sample_annotation[,6],"(?i)[patient |person ]\\d{1,3}","")
sample_annotation[,7] <- str_replace_all(sample_annotation[,7],"(?i)\\({0,1}(Biological Replicate |replicate |replicate)(\\d|[A-Z])\\){0,1}$","")
sample_annotation[,7] <- str_replace_all(sample_annotation[,7],"repeat\\-\\d","")
sample_annotation[,7] <- str_replace_all(sample_annotation[,7],"replicate\\_\\d","")
sample_annotation[,7] <- str_replace_all(sample_annotation[,7],"replicate\\_\\d","")
sample_annotation[,7] <- str_replace_all(sample_annotation[,7],"(?i)(control)\\-\\d$","control")
sample_annotation[,7] <- str_replace_all(sample_annotation[,7],"(?i)(infection)\\-\\d$","infection")
sample_annotation[,7] <- str_replace_all(sample_annotation[,7],"(?i)control \\d$","control")
sample_annotation[,7] <- str_replace_all(sample_annotation[,7],"(?i)Silence \\d$","Silence")
sample_annotation[,7] <- str_replace_all(sample_annotation[,7],"(?i)knock-down( ){0,1}\\d$","knock-down")
sample_annotation[,7] <- str_replace_all(sample_annotation[,7],"(?i)knock-out( ){0,1}\\d$","knock-out")
sample_annotation[,9] <- str_replace_all(sample_annotation[,9],"replicate \\d","replicate")
sample_annotation[,9] <- str_replace_all(sample_annotation[,9],"\\(replicate\\d\\)","")
sample_annotation[,9] <- str_replace_all(sample_annotation[,9],"control \\d$","")
sample_annotation[,9] <- str_replace_all(sample_annotation[,9],"(?i)control#\\d$","Control")
sample_annotation[,10] <- str_replace_all(sample_annotation[,10],"(?i)replicate(_| )[0-9]{1,3}","")
for(i in 3:ncol(sample_annotation)){
    sample_annotation[,i] <- str_replace_all(sample_annotation[,i],"#\\d","")
}
sample_annotation[is.na(sample_annotation)] <- ''
# write.table(sample_annotation,'sample_annotation.csv',sep = '\t',quote = FALSE,row.names = FALSE)
# Annotation assessment
# samples fail to be annotated
length(which(rowSums(sample_annotation=='') == 8)) / length(pre_anno$sample_attributes)
# studies fail to annotate
length(unique(sample_annotation[which(rowSums(sample_annotation=='') == 8),1])) # 724 out of 18357 studies fails to annotate ---> 874
length(which(rowSums(sample_annotation=='') == 8)) # 22840

# Only one annotation attribute sample
length(which(rowSums(sample_annotation=='') == 7)) # 87668 samples ---> 141615
length(unique(sample_annotation[which(rowSums(sample_annotation=='') == 7),1])) #2642 studies ---> 4931

# Two annotation attributes
length(which(rowSums(sample_annotation=='') == 6)) # 172972 samples---> 172022
length(unique(sample_annotation[which(rowSums(sample_annotation=='') == 6),1])) #6604 studies --->7329

# Three annotation attributes
length(which(rowSums(sample_annotation=='') == 5)) # 140111 samples---> 100930
length(unique(sample_annotation[which(rowSums(sample_annotation=='') == 5),1])) #6830 studies ---> 4571

# Four annotation attributes
length(which(rowSums(sample_annotation=='') == 4)) # 46392 samples---> 33806
length(unique(sample_annotation[which(rowSums(sample_annotation=='') == 4),1])) #1999 studies ---> 1256

# Five annotation attributes
length(which(rowSums(sample_annotation=='') == 3)) # 4932 samples---> 3616
length(unique(sample_annotation[which(rowSums(sample_annotation=='') == 3),1])) #312 studies ---> 206

# More than two annotation attributes
length(which(rowSums(sample_annotation=='') <= 6)) # 311170

# Assessment of keyword extraction
annotate_assess <- data.frame(study = c(length(unique(sample_annotation[which(rowSums(sample_annotation=='') == 8),1])),
                                        length(unique(sample_annotation[which(rowSums(sample_annotation=='') == 7),1])),
                                        length(unique(sample_annotation[which(rowSums(sample_annotation=='') == 6),1])),
                                        length(unique(sample_annotation[which(rowSums(sample_annotation=='') == 5),1])),
                                        length(unique(sample_annotation[which(rowSums(sample_annotation=='') == 4),1])),
                                        length(unique(sample_annotation[which(rowSums(sample_annotation=='') == 3),1])),
                                        length(unique(sample_annotation[which(rowSums(sample_annotation=='') == 2),1])),
                                        length(unique(sample_annotation[which(rowSums(sample_annotation=='') == 1),1]))),
                            sample = c(length(which(rowSums(sample_annotation=='') == 8)),
                                        length(which(rowSums(sample_annotation=='') == 7)),
                                        length(which(rowSums(sample_annotation=='') == 6)),
                                        length(which(rowSums(sample_annotation=='') == 5)),
                                        length(which(rowSums(sample_annotation=='') == 4)),
                                        length(which(rowSums(sample_annotation=='') == 3)),
                                        length(which(rowSums(sample_annotation=='') == 2)),
                                        length(which(rowSums(sample_annotation=='') == 1))))
annotate_assess$`Number of annotated attribute(s)` <- as.integer(rownames(annotate_assess))-1
library(ggplot2)
library(ggbreak)
library(patchwork)

annotate_study <- ggplot(annotate_assess, aes(`Number of annotated attribute(s)` ,study)) +
        geom_col()+ scale_x_continuous(labels = seq(0,7), breaks = seq(0,7))+
        scale_y_cut(breaks=c(100, 2100), which=c(1, 3), scales=c(3, 0.5))+
        xlab("Number of annotated attribute(s)") + ylab("Number of studies") +
        ggtitle("Number of attributes annotated \n\t\tfor studies")
annotate_sample <- ggplot(annotate_assess, aes(`Number of annotated attribute(s)` ,sample)) +
        geom_col()+ scale_x_continuous(labels = seq(0,7), breaks = seq(0,7))+
        xlab("Number of annotated attribute(s)") + ylab("Number of samples") +
        scale_y_cut(breaks=c(200, 5000), which=c(1, 3), scales=c(3, 0.5)) +
        ggtitle("Number of attributes annotated \n\t\tfor samples")
pdf("annotation_assessment.pdf")
annotate_sample
dev.off()
write.csv(sample_annotation,'sample_annotation.csv',row.names=FALSE)
#--------------------------------------------------------------------------------------#
# Step 4. Grouping
#--------------------------------------------------------------------------------------#
# Comebind same group 
# ignore samples with too long condition
sample_annotation <- sample_annotation[which(str_count(sample_annotation[,6] ,"\\w+") < 30),]
Group <- aggregate(external_id ~ study+ tissue + cellline + celltype + condition + treatment + treatment_duration + genotype + subtype, sample_annotation, FUN = paste, collapse=",")
dim(Group) # 60505

Group <- Group[order(Group$study),]
# filter group contains only one sample
Group <- Group[which(grepl( ',', Group[,10], fixed = TRUE)),]
# projects that only contain one sample in each group
length(unique(Group[-which(grepl( ',', Group[,10], fixed = TRUE)),1]))
dim(Group) # 45587
Group <- Group[-which(rowSums(Group=='') >= 7),]
dim(Group) # 34455

