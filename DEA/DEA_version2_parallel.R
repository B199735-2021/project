setwd("/localdisk/home/s2172876/project/DEA")
library(recount3)
library(stringr)
library(DESeq2)
library(BiocParallel)
library(tibble)
library(edgeR)
library(doParallel)  
library(BiocParallel) 
register(MulticoreParam(10))
no_cores <- 10
study_id <- read.csv('../annotation/pre_re_annotation.csv',sep = '\t')[,c(1,2)]
experiment <- read.csv('../annotation/Contrasts_group_description',sep = '\t')
abstract <- read.csv('/localdisk/home/s2172876/project/DEA/recount3_selection_2022-05-26 17_28_49.csv')[,c(1,3)]
study_id <-  merge(study_id,abstract,by.x = 'study',by.y = 'project')

# Numbering the experiments
experiment$ExpID <- paste0('Exp',rownames(experiment))
projects <- rbind(
    recount3::available_projects("human"),
    recount3::available_projects("mouse")
)

registerDoParallel(cores=no_cores)  
cl <- makeCluster(no_cores, type="FORK")  
# If two adjacent comparison groups originate from the same experiment, there is no need to re-download the expression profile.
# define null for the first experiment
pre_study_id <- ''
getPrimeNumbers <- function(i){
     try({
        # study id
        control_study <- unique(study_id[which(study_id[,2] %in% str_split(experiment$Control.Sample.ID[i],',')[[1]]  ),1])
        treatment_study <- unique(study_id[which(study_id[,2] %in% str_split(experiment$Treatment.Sample.ID[i],',')[[1]]),1])
        if(length(control_study)==1&length(treatment_study)==1){
            if(control_study==treatment_study){
                # sample id
                control_sample <- experiment$Control.Sample.ID[i]
                treatment_sample <- experiment$Treatment.Sample.ID[i]
                # ignore comparisons which only contain one sample
                if(length(str_split(control_sample,",")[[1]])>1 & length(str_split(treatment_sample,",")[[1]])>1){
                    if(length(unique(study_id[which(study_id$study==control_study),3]))==1){
                        ## Create a RSE object at the gene level
                        if(control_study!=pre_study_id){
                            rse_gene <- create_rse_manual(control_study,organism =  unique(study_id[which(study_id$study==control_study),3]))
                            ## Scale the counts using the AUC
                            assays(rse_gene)$counts <- transform_counts(rse_gene)
                        }
                        
                    ## DEA 
                    # select sample columns and make it into DESeq Object
                    Expr <- as.data.frame(assays(rse_gene)$counts[,c(which(colnames(assays(rse_gene)$counts) %in% str_split(control_sample,",")[[1]]),
                                        which(colnames(assays(rse_gene)$counts) %in% str_split(treatment_sample,",")[[1]]))])

                    coldata <- data.frame(row.names=colnames(Expr), 
                                        factor = c(rep('control',length(str_split(control_sample,",")[[1]])),rep('treatment',length(str_split(treatment_sample,",")[[1]]))))
                    dds <- DESeqDataSetFromMatrix(countData = Expr,
                                                colData = coldata,
                                                design = ~factor)

                    # Apply filters to remove genes that do not pass quality filters
                    # first filter based upon counts
                    keep <- rowSums(counts(dds)) > 1
                    dds <- dds[keep,]
                    nrow(dds)
                    # second filter-require >10 in >= 2 samples
                    keep <- rowSums(counts(dds) >= 10) >= 2
                    dds <- dds[keep,]
                    nrow(dds)

                    # size factor estimation for normalisation
                    dds <- estimateSizeFactors(dds)
                    sizeFactors(dds)

                    # Defferential gene expression
                    deg <- DESeq(dds,parallel=TRUE)
                    res <- results(deg, contrast = c('factor', 'treatment', 'control'))
                    deg_res <- data.frame(res)
                    # deg filter
                    deg_res <- deg_res[which((abs(deg_res$log2FoldChange)>1) & (deg_res$padj<0.05)),]
                    # Convert row names into first column 
                    # It is possible to have duplicate differentially expressed genes (the same gene is differently expressed in different experiment), 
                    # and dataframe does not allow duplicate row names. 
                    # To avoid errors being reported, the row names are converted to the first column.
                    deg_res <- tibble::rownames_to_column(deg_res, "Gene")
                    deg_res$experimentID <- rep(experiment[i,9],nrow(deg_res))
                    pre_study_id <- control_study
                    file_name <- paste0("./DEA_Result/DEA_result_",experiment[i,10],".csv")
                    write.table(deg_res,file_name,sep = '\t',quote = FALSE,row.names = FALSE)
                    return(dim(deg_res[1])[1])
                    rm(deg_res)
                }
            }
        }
        }
        }
    )
}


t1 <- Sys.time()
#parLapply(cl, 58756:58858, getPrimeNumbers)
result <- foreach(i=940:dim(experiment)[1]) %dopar% getPrimeNumbers(i)  
t2 <- Sys.time()
print(t2-t1)
result
stopCluster(cl)  
# save.image(file = "./DEA.RData")
