setwd("/localdisk/home/s2172876/project/DEA")
# library('devtools')
# install_github('rsait/ORdensity')
library(recount3)
library(stringr)
library(DESeq2)
library("BiocParallel")
library(tibble)
register(MulticoreParam(10))
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

# define output dataframe
outDEG_df <- data.frame(matrix(ncol = 8, nrow = 0))
# define summary dataframe
summary_DEG <- data.frame(matrix(ncol = 2, nrow = 0))
# create file folder DEA_Result
mainDir <- getwd()
subDir <- "DEA_Result"
dir.create(file.path(mainDir, subDir), showWarnings = FALSE)


# counter to mark result file
k = 1
# If two adjacent comparison groups originate from the same experiment, there is no need to re-download the expression profile.
# define null for the first experiment
pre_study_id <- ''
x = 1
t1 <- Sys.time()
# start to loop and analyze
for(i in 58758){
    try({
        # study id
        control_study <- unique(study_id[which(study_id[,2] %in% str_split(experiment$Control.Sample.ID[i],',')[[1]]  ),1])
        treatment_study <- unique(study_id[which(study_id[,2] %in% str_split(experiment$Treatment.Sample.ID[i],',')[[1]]),1])
        if (control_study!=treatment_study){
            next
        }
        # sample id
        control_sample <- experiment$Control.Sample.ID[i]
        treatment_sample <- experiment$Treatment.Sample.ID[i]
        # ignore comparisons which only contain one sample
        if(length(str_split(control_sample,",")[[1]])>1 & length(str_split(treatment_sample,",")[[1]])>1){
            if(control_study!=pre_study_id){
                ## Create a RSE object at the gene level
                if(length(unique(study_id[which(study_id$study==control_study),3]))!=1){
                    next
                }
                rse_gene <- create_rse_manual(control_study,organism =  unique(study_id[which(study_id$study==control_study),3]))
                ## Scale the counts using the AUC
                assays(rse_gene)$counts <- transform_counts(rse_gene)
                x = x+1
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
            outDEG_df <- rbind(outDEG_df,deg_res)
            summary_DEG <- rbind(summary_DEG,c(experiment[i,9],nrow(deg_res)))
            pre_study_id <- control_study
            file_name <- paste0("./DEA_Result/DEA_result_",k,".csv")
            if(nrow(outDEG_df)>1000000){
                write.table(outDEG_df,file_name,sep = '\t',quote = FALSE,row.names = FALSE)
                k = k+1
                # empty outDEG_df
                outDEG_df <- data.frame(matrix(ncol = 8, nrow = 0))
            }
        }
        }
    )
    print('Finish')
}

t2 <- Sys.time()
t2-t1 # Time difference of 2.314128 mins
x
dim(summary_DEG);dim(outDEG_df)
colnames(summary_DEG) <- c('ExpID','DEGnum')
file_name <- paste0("DEA_Result/DEA_result_",k,".csv")
write.table(outDEG_df,file_name,sep = '\t',quote = FALSE,row.names = FALSE)
write.table(summary_DEG,"./summary_DEG.csv",sep = '\t',quote = FALSE,row.names = FALSE)
save.image(file = "./DEA.RData")










control_study <- unique(study_id[which(study_id[,2] %in% str_split(experiment$Control.Sample.ID[i],',')[[1]]  ),1])
        treatment_study <- unique(study_id[which(study_id[,2] %in% str_split(experiment$Treatment.Sample.ID[i],',')[[1]]),1])
      
        # sample id
        control_sample <- experiment$Control.Sample.ID[i]
        treatment_sample <- experiment$Treatment.Sample.ID[i]
        # ignore comparisons which only contain one sample
        if(length(str_split(control_sample,",")[[1]])>1 & length(str_split(treatment_sample,",")[[1]])>1){
            if(control_study!=pre_study_id){
                ## Create a RSE object at the gene level
                if(length(unique(study_id[which(study_id$study==control_study),3]))!=1){
                    next
                }
                rse_gene <- create_rse_manual(control_study,organism =  unique(study_id[which(study_id$study==control_study),3]))
                ## Scale the counts using the AUC
                assays(rse_gene)$counts <- transform_counts(rse_gene)
                x = x+1
            }
            ## DEA 
            # select sample columns and make it into DESeq Object
            Expr <- as.data.frame(assays(rse_gene)$counts[,c(which(colnames(assays(rse_gene)$counts) %in% str_split(control_sample,",")[[1]]),
                                which(colnames(assays(rse_gene)$counts) %in% str_split(treatment_sample,",")[[1]]))])

            coldata <- data.frame(row.names=colnames(Expr), 
                                factor = c(rep('0',length(str_split(control_sample,",")[[1]])),rep('1',length(str_split(treatment_sample,",")[[1]]))))
            coldata
            dds <- DESeqDataSetFromMatrix(countData = Expr,
                                        colData = coldata,
                                        design = ~factor)
            deg <- DESeq(dds,parallel=TRUE)
            res <- results(deg, contrast = c('factor', '1', '0'))
            deg_res <- data.frame(res)

dds <- makeExampleDESeqDataSet(n=100,m=12)
dds$genotype <- factor(rep(rep(c("I","II"),each=3),2))
dds$condition
design(dds) <- ~ genotype + condition + genotype:condition
dds <- DESeq(dds)
resultsNames(dds)
# the condition effect for genotype I (the main effect)
results(dds, contrast=c("condition","B","A"))
# the condition effect for genotype II
# this is, by definition, the main effect *plus* the interaction term
# (the extra condition effect in genotype II compared to genotype I).
results(dds, list( c("condition_B_vs_A","genotypeII.conditionB") ))
# the interaction term, answering: is the condition effect *different* across genotypes?
results(dds, name="genotypeII.conditionB")