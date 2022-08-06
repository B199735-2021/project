# assign project ID to different type of comparison
setwd("/localdisk/home/s2172876/project/DEA")
study_id <- read.csv('../annotation/pre_re_annotation.csv',sep = '\t')[,c(1,2)]
experiment <- read.csv('../annotation/group.txt',sep = '\t')
# identify factors
single_factor_two_levels <- vector()
single_factor_multi_levels <- vector()
two_factor <- vector()
three_factor <- vector()
more_than_three_factor <- vector()
tissue_specific <- vector()
for(i in unique(experiment$study)){
    summary_level<-vector()
    for(j in 1:9){
    summary_level <-c(summary_level,length(table(experiment[which(experiment$study == i),j])))
    }
    if(length(which(summary_level !=1))==1){#single factor
        if(summary_level[(which(summary_level !=1))]==2){
            single_factor_two_levels <- c(single_factor_two_levels,i)
        }else{
            single_factor_multi_levels <- c(single_factor_multi_levels,i)
        }
    } 
    if(length(which(summary_level !=1))==2){#2 factors
        two_factor <- c(two_factor,i)
        if(summary_level[2]!=1){
            tissue_specific <- c(tissue_specific,i)
        }
    }
    if(length(which(summary_level !=1))==3){#3 factors
        three_factor <- c(three_factor,i)
        if(summary_level[2]!=1){
            tissue_specific <- c(tissue_specific,i)
        }
    }
    if(length(which(summary_level !=1))>3){#more than3 factors
        more_than_three_factor <- c(more_than_three_factor,i)
        if(summary_level[2]!=1){
            tissue_specific <- c(tissue_specific,i)
        }
    }
}
length(tissue_specific)

factors_summary <- data.frame(Number=c(length(single_factor_two_levels),length(single_factor_multi_levels),length(two_factor),length(three_factor),length(more_than_three_factor)))
factors_summary$Type <- factor(c("single factor (two levels)","single factor (multiple levels)",'two factors','three factors','more than three factors'),
                                levels = c('single factor (two levels)','single factor (multiple levels)','two factors','three factors','more than three factors'))
factors_summary

# write.table(single_factor_two_levels,'single_factor_two_levels.txt',sep = '\t',row.names = FALSE,col.names = FALSE,quote = FALSE)
# write.table(single_factor_multi_levels,'single_factor_multi_levels.txt',sep = '\t',row.names = FALSE,col.names = FALSE,quote = FALSE)
# write.table(two_factor,'two_factor.txt',sep = '\t',row.names = FALSE,col.names = FALSE,quote = FALSE)
# write.table(three_factor,'three_factor.txt',sep = '\t',row.names = FALSE,col.names = FALSE,quote = FALSE)
# write.table(more_than_three_factor,'more_than_three_factor.txt',sep = '\t',row.names = FALSE,col.names = FALSE,quote = FALSE)
# write.table(tissue_specific,'tissue_specific.txt',sep = '\t',row.names = FALSE,col.names = FALSE,quote = FALSE)
