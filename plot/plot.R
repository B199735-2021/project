setwd('/localdisk/home/s2172876/project/result_plot')
# Result Figure1
nonSC <- read.csv('../annotation/recount3_nonSC.csv')
study_source <- as.data.frame(table(nonSC$project_home))
study_source$Var1 <- c("GTEx","SRA Bulk RNA-Seq","TCGA")
study_source <- rbind(study_source,c('SRA single cell RNA-Seq',412))
colnames(study_source) <- c('Type/Source', 'Num')
study_source$Num <- as.integer(study_source$Num)
study_source <- study_source[c(2,4,3,1),]
library(ggplot2)
library(ggbreak)
library(wesanderson)
mycolors <- c("#0F6444", "#3182BD", "#8576AB", "#D68484")
png('Result_figure1_project_type.png')
ggplot(study_source, aes(`Type/Source` ,Num)) +
        geom_col(aes(fill=`Type/Source`))+ scale_y_cut(breaks=c(50,500),which=c(1, 3), scales=c(3, 1))+
        xlab("Project type/source") + ylab("Number of studies") + ylim(0, 19000)+
        ggtitle("Original source and type of projects in recount3") +
        geom_text(aes(label=Num),vjust=1,color = "white",size=4)+
        theme_bw()+ 
        theme(plot.title = element_text(hjust = 0.5),legend.position="bottom")+
        scale_fill_manual(values = mycolors)
dev.off()

# Result Figure2
sample_annotation <- read.csv('../annotation/sample_annotation.csv',sep = '\t')
sample_annotation <- sample_annotation[which(sample_annotation$study!=''),]
sample_annotation[is.na(sample_annotation)] <- ''
dim(sample_annotation)
# Assessment of keyword extraction
annotate_assess <- data.frame(sample = c(length(which(rowSums(sample_annotation=='') == 8)),
                                        length(which(rowSums(sample_annotation=='') == 7)),
                                        length(which(rowSums(sample_annotation=='') == 6)),
                                        length(which(rowSums(sample_annotation=='') == 5)),
                                        length(which(rowSums(sample_annotation=='') == 4)),
                                        length(which(rowSums(sample_annotation=='') == 3)),
                                        length(which(rowSums(sample_annotation=='') == 2)),
                                        length(which(rowSums(sample_annotation=='') == 1)),
                                        length(which(rowSums(sample_annotation=='') == 0))))
annotate_assess$`Number of annotated attribute(s)` <- as.integer(rownames(annotate_assess))-1
annotate_assess

annotate_assess[,2] <- as.factor(annotate_assess[,2])
mycolors <- c("#0F6444", "#3182BD", "#8576AB", "#481010","#af6729",
              "#857c7c","#aa2a2a","#28e454f0","#3803e5")
png('Result_figure2_annotation_assessment.png')
ggplot(annotate_assess, aes(`Number of annotated attribute(s)` ,sample)) +
        geom_col(aes(fill=`Number of annotated attribute(s)`))+ scale_y_cut(breaks=c(5,700,5000),which=c(1, 3), scales=c(3, 1))+
        xlab("Annotated attribute(s)") +ylab('')+  ylim(0, 200000)+
        ggtitle("Number of attributes annotated") +
        geom_text(aes(label=sample),vjust=1,color = "white",size=4,fontface = "bold")+
        theme_bw()+ 
        theme(plot.title = element_text(hjust = 0.5),legend.position="none")+
        scale_fill_manual(values = mycolors)
dev.off()

study = c(length(unique(sample_annotation[which(rowSums(sample_annotation=='') == 8),1])),
                                        length(unique(sample_annotation[which(rowSums(sample_annotation=='') == 7),1])),
                                        length(unique(sample_annotation[which(rowSums(sample_annotation=='') == 6),1])),
                                        length(unique(sample_annotation[which(rowSums(sample_annotation=='') == 5),1])),
                                        length(unique(sample_annotation[which(rowSums(sample_annotation=='') == 4),1])),
                                        length(unique(sample_annotation[which(rowSums(sample_annotation=='') == 3),1])),
                                        length(unique(sample_annotation[which(rowSums(sample_annotation=='') == 2),1])),
                                        length(unique(sample_annotation[which(rowSums(sample_annotation=='') == 1),1])))

# Result Figure : summary result (number of comparison)
comp <- c(system('ls ~/project_backup/DEA_Result/singlefactor_twolevel/ | wc -l',intern = TRUE),
 system('cat ~/project_backup/DEA_Result/summary_1.txt | wc -l',intern = TRUE),
 system('ls ~/project_backup/DEA_Result/multifactor | wc -l',intern = TRUE),
 system('ls ~/project_backup/DEA_Result/TCGA | wc -l',intern = TRUE),
 system('ls ~/project_backup/DEA_Result/GTEx | wc -l',intern = TRUE))
comp <- as.integer(comp)
comp <- as.data.frame(cbind(number = comp, type = c('single-factor \n two levels (SRA)',
                                      'single-factor \n multilevel (SRA)',
                                      'multifactorial \n (SRA)',
                                      'TCGA','GTEx')))
comp$type <- factor(comp$type, levels = c('single-factor \n two levels (SRA)',
                                      'single-factor \n multilevel (SRA)',
                                      'multifactorial \n (SRA)',
                                      'TCGA','GTEx'))
comp$number <- as.integer(comp$number)
comp$number[2] <- comp$number[2]-1
mycolors <- c("#0F6444", "#3182BD", "#8576AB", "#481010","#af6729")

ggplot(comp, aes(type ,number)) +
        geom_col(aes(fill=type)) + 
        scale_y_cut(breaks=c(50,250,5000),which=c(1, 3), scales=c(3, 1))+
        xlab("Comparison type") + ylab("Number of comparisons") + 
        ggtitle("Comparison type/source") +
        geom_text(aes(label=number),vjust=1,color = "white",size=5,,fontface = "bold")+
        ylim(0, 11000)+
        theme_bw()+ 
        theme(plot.title = element_text(hjust = 0.65),text = element_text(size = 20),
        legend.position="non")+
        scale_fill_manual(values = mycolors)
ggsave('Result_figure3_comparison_type.png', height = 9 , width = 11)
