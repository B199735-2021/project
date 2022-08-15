setwd('/home/s2172876/mirna')
library('ampir')
mirlist <- read_faa(file = 'mature.fa.gz')
mirlist$motif <- reverseComplement(substring(mirlist[,2],2,7))
library(spgs)
library(stringr)
library(data.table)
library(readr)
# reverse complement

mirlist$motif <- toupper(mirlist$motif)


mirlist$id <- sapply(str_split(mirlist[,1],' '), "[[", 2)
mirlist$mir <- sapply(str_split(mirlist[,1],' '), "[[", 1)
mirlist$organism <- paste(sapply(str_split(mirlist[,1],' '), "[[", 3),sapply(str_split(mirlist[,1],' '), "[[", 4))

motif_5UTR <- data.frame(matrix(ncol = 6))
topN = 5
for(i in 1:as.integer(system('ls /home/s2172876/mirna/motif_result/size_6/sylamer_result_5UTR | wc -l',TRUE))){
    file = paste0('/home/s2172876/mirna/motif_result/size_6/sylamer_result_5UTR/Exp',i,'_5UTR.out')
    # Read table of sylamer
    tab <- cbind(0,data.frame(fread(file,header = TRUE), row.names=1))
    # ordering sylamer tab
    mycols <- 1:ncol(tab) # columns to consider in the ordering
    tab <- tab[order(apply(abs(tab[,mycols]), 1, max), decreasing=TRUE),] # ordering table based in max enrichment
    #   top5 motif
    motif_5UTR[i,] <- c(paste('Exp',i,sep = ''),rownames(tab)[1:topN])
    print(i)
}

motif_3UTR <- data.frame(matrix(ncol = 6))
for(i in 1:as.integer(system('ls /home/s2172876/mirna/motif_result/size_6/sylamer_result_3UTR | wc -l',TRUE))){
    file = paste0('/home/s2172876/mirna/motif_result/size_6/sylamer_result_3UTR/Exp',i,'_3UTR.out')
    # Read table of sylamer
    tab <- cbind(0,data.frame(fread(file,header = TRUE), row.names=1))
    #   tab <- abs(tab)
    # ordering sylamer tab
    mycols <- 1:ncol(tab) # columns to consider in the ordering
    tab <- tab[order(apply(abs(tab[,mycols]), 1, max), decreasing=TRUE),] # ordering table based in max enrichment
    #   top5 motif
    motif_3UTR[i,] <- c(paste('Exp',i,sep = ''),rownames(tab)[1:topN])
    print(i)
}

motif_5UTR_summary <- data.frame( ExpID = rep(motif_5UTR[,1],5),
                                  motif = c(motif_5UTR[,2],motif_5UTR[,3],motif_5UTR[,4],motif_5UTR[,5],motif_5UTR[,6]),
                                  type = rep('Five Prime Untranslated Region',length(motif_5UTR[,1])*5))
motif_5UTR_summary_agg <- aggregate(ExpID ~ motif + type, motif_5UTR_summary, FUN = paste, collapse=",")
library(stringr)
motif_5UTR_summary_agg$count <- str_count(motif_5UTR_summary_agg[,3],',')
motif_5UTR_summary_agg$count <- motif_5UTR_summary_agg$count+1

motif_3UTR_summary <- data.frame( ExpID = rep(motif_3UTR[,1],5),
                                  motif = c(motif_3UTR[,2],motif_3UTR[,3],motif_3UTR[,4],motif_3UTR[,5],motif_3UTR[,6]),
                                  type = rep('Three Prime Untranslated Region',length(motif_3UTR[,1])*5))
motif_3UTR_summary_agg <- aggregate(ExpID ~ motif + type, motif_3UTR_summary, FUN = paste, collapse=",")
motif_3UTR_summary_agg$count <- str_count(motif_3UTR_summary_agg[,3],',')
motif_3UTR_summary_agg$count <- motif_3UTR_summary_agg$count+1
motif_summary_agg <- rbind(motif_3UTR_summary_agg,motif_5UTR_summary_agg)

# write.table(motif_summary_agg,'motif_summary_agg.txt',sep = '\t',row.names = FALSE,quote = FALSE)

exp_summary <- read.csv('../DEG_summary_new_fix',sep = '\t')[,c(1,2)]
project <- read.csv('../annotation/recount3_selection_2022-05-11 21_54_01.csv')
exp_summary <- merge(exp_summary,project,by.x = 'ProjectID',by.y = 'project')
motif_5UTR_species <- merge(motif_5UTR,exp_summary[,c(2,3)],by.x = 'X1',by.y = 'ExpID')
motif_3UTR_species <- merge(motif_3UTR,exp_summary[,c(2,3)],by.x = 'X1',by.y = 'ExpID')

# 5'UTR motif frequency
# human
human_motif_5UTR <- as.vector(as.matrix(motif_5UTR_species[which(motif_5UTR_species[,7]=='human'),2:6]))
human_motif_5UTR <- human_motif_5UTR[which(human_motif_5UTR!='')]
human_mir_5UTR <- as.data.frame(table(human_motif_5UTR[which(human_motif_5UTR %in% mirlist[which(mirlist$organism=='Homo sapiens'),'motif'])]))
human_nonmir_5UTR <- as.data.frame(table(human_motif_5UTR[-which(human_motif_5UTR %in% mirlist[which(mirlist$organism=='Homo sapiens'),'motif'])]))

# mouse
mouse_motif_5UTR <- as.vector(as.matrix(motif_5UTR_species[which(motif_5UTR_species[,7]=='mouse'),2:6]))
mouse_motif_5UTR <- mouse_motif_5UTR[which(mouse_motif_5UTR!='')]
mouse_mir_5UTR <- as.data.frame(table(mouse_motif_5UTR[which(mouse_motif_5UTR %in% mirlist[which(mirlist$organism=='Mus musculus'),'motif'])]))
mouse_nonmir_5UTR <- as.data.frame(table(mouse_motif_5UTR[-which(mouse_motif_5UTR %in% mirlist[which(mirlist$organism=='Mus musculus'),'motif'])]))

nrow(human_mir_5UTR)+nrow(mouse_mir_5UTR)
nrow(human_nonmir_5UTR)+nrow(mouse_nonmir_5UTR)

# 3'UTR motif frequency
# human
human_motif_3UTR <- as.vector(as.matrix(motif_3UTR_species[which(motif_3UTR_species[,7]=='human'),2:6]))
human_motif_3UTR <- human_motif_3UTR[which(human_motif_3UTR!='')]
human_mir_3UTR <- as.data.frame(table(human_motif_3UTR[which(human_motif_3UTR %in% mirlist[which(mirlist$organism=='Homo sapiens'),'motif'])]))
human_nonmir_3UTR <- as.data.frame(table(human_motif_3UTR[-which(human_motif_3UTR %in% mirlist[which(mirlist$organism=='Homo sapiens'),'motif'])]))

# mouse
mouse_motif_3UTR <- as.vector(as.matrix(motif_3UTR_species[which(motif_3UTR_species[,7]=='mouse'),2:6]))
mouse_motif_3UTR <- mouse_motif_3UTR[which(mouse_motif_3UTR!='')]
mouse_mir_3UTR <- as.data.frame(table(mouse_motif_3UTR[which(mouse_motif_3UTR %in% mirlist[which(mirlist$organism=='Mus musculus'),'motif'])]))
mouse_nonmir_3UTR <- as.data.frame(table(mouse_motif_3UTR[-which(mouse_motif_3UTR %in% mirlist[which(mirlist$organism=='Mus musculus'),'motif'])]))

nrow(human_mir_3UTR)+nrow(mouse_mir_3UTR)
nrow(human_nonmir_3UTR)+nrow(mouse_nonmir_3UTR)

# summary table
motif_summary_compare <- data.frame(matrix(ncol=4)) 
motif_summary_compare[1,] <- c(nrow(human_mir_5UTR)+nrow(mouse_mir_5UTR),
                               sum(human_mir_5UTR$Freq) + sum(mouse_mir_5UTR$Freq),
                               nrow(human_mir_3UTR)+nrow(mouse_mir_3UTR),
                               sum(human_mir_3UTR$Freq) + sum(mouse_mir_3UTR$Freq))

colnames(motif_summary_compare) <- c('5\'UTR miRNA-related motif','5\'UTR miRNA-related motif frequency',
                                     '3\'UTR miRNA-related motif','3\'UTR miRNA-related motif frequency')


motif_5UTR_cutoff8 <- data.frame(matrix(ncol = 6))
for(i in 1:as.integer(system('ls /home/s2172876/mirna/motif_result/size_6/sylamer_result_5UTR | wc -l',TRUE))){
    file = paste0('/home/s2172876/mirna/motif_result/size_6/sylamer_result_5UTR/Exp',i,'_5UTR.out')
    # Read table of sylamer
    tab <- cbind(0,data.frame(fread(file,header = TRUE), row.names=1))
    # ordering sylamer tab
    mycols <- 1:ncol(tab) # columns to consider in the ordering
    tab <- tab[order(apply(abs(tab[,mycols]), 1, max), decreasing=TRUE),] # ordering table based in max enrichment
    # filter
    if(length(which(apply(abs(tab[1:5,2:(ncol(tab)-1)]), 1, max)<8))){
        m <- c(names(which(apply(abs(tab[1:5,2:(ncol(tab)-1)]), 1, max)>8)),
               rep('',5-length(which(apply(abs(tab[1:5,2:(ncol(tab)-1)]), 1, max)>8))))
    }else{
        m <- names(which(apply(abs(tab[1:5,2:(ncol(tab)-1)]), 1, max)>8))
    }
    motif_5UTR_cutoff8[i,] <- c(paste('Exp',i,sep = ''),m)
    print(i)
    print(m)
}
motif_3UTR_cutoff8 <- data.frame(matrix(ncol = 6))
for(i in 1:as.integer(system('ls /home/s2172876/mirna/motif_result/size_6/sylamer_result_3UTR | wc -l',TRUE))){
    file = paste0('/home/s2172876/mirna/motif_result/size_6/sylamer_result_3UTR/Exp',i,'_3UTR.out')
    # Read table of sylamer
    tab <- cbind(0,data.frame(fread(file,header = TRUE), row.names=1))
    # ordering sylamer tab
    mycols <- 1:ncol(tab) # columns to consider in the ordering
    tab <- tab[order(apply(abs(tab[,mycols]), 1, max), decreasing=TRUE),] # ordering table based in max enrichment
    # filter
    if(length(which(apply(abs(tab[1:5,2:(ncol(tab)-1)]), 1, max)<8))){
        m <- c(names(which(apply(abs(tab[1:5,2:(ncol(tab)-1)]), 1, max)>8)),
               rep('',5-length(which(apply(abs(tab[1:5,2:(ncol(tab)-1)]), 1, max)>8))))
    }else{
        m <- names(which(apply(abs(tab[1:5,2:(ncol(tab)-1)]), 1, max)>8))
    }
    motif_3UTR_cutoff8[i,] <- c(paste('Exp',i,sep = ''),m)
    print(i)
    print(m)
}

motif_5UTR_species <- merge(motif_5UTR_cutoff8,exp_summary[,c(2,3)],by.x = 'X1',by.y = 'ExpID')
motif_3UTR_species <- merge(motif_3UTR_cutoff8,exp_summary[,c(2,3)],by.x = 'X1',by.y = 'ExpID')

# 5'UTR motif frequency
# human
human_motif_5UTR <- as.vector(as.matrix(motif_5UTR_species[which(motif_5UTR_species[,7]=='human'),2:6]))
human_motif_5UTR <- human_motif_5UTR[which(human_motif_5UTR!='')]
human_mir_5UTR <- as.data.frame(table(human_motif_5UTR[which(human_motif_5UTR %in% mirlist[which(mirlist$organism=='Homo sapiens'),'motif'])]))
human_nonmir_5UTR <- as.data.frame(table(human_motif_5UTR[-which(human_motif_5UTR %in% mirlist[which(mirlist$organism=='Homo sapiens'),'motif'])]))

# mouse
mouse_motif_5UTR <- as.vector(as.matrix(motif_5UTR_species[which(motif_5UTR_species[,7]=='mouse'),2:6]))
mouse_motif_5UTR <- mouse_motif_5UTR[which(mouse_motif_5UTR!='')]
mouse_mir_5UTR <- as.data.frame(table(mouse_motif_5UTR[which(mouse_motif_5UTR %in% mirlist[which(mirlist$organism=='Mus musculus'),'motif'])]))
mouse_nonmir_5UTR <- as.data.frame(table(mouse_motif_5UTR[-which(mouse_motif_5UTR %in% mirlist[which(mirlist$organism=='Mus musculus'),'motif'])]))

nrow(human_mir_5UTR)+nrow(mouse_mir_5UTR)
nrow(human_nonmir_5UTR)+nrow(mouse_nonmir_5UTR)

# 3'UTR motif frequency
# human
human_motif_3UTR <- as.vector(as.matrix(motif_3UTR_species[which(motif_3UTR_species[,7]=='human'),2:6]))
human_motif_3UTR <- human_motif_3UTR[which(human_motif_3UTR!='')]
human_mir_3UTR <- as.data.frame(table(human_motif_3UTR[which(human_motif_3UTR %in% mirlist[which(mirlist$organism=='Homo sapiens'),'motif'])]))
human_nonmir_3UTR <- as.data.frame(table(human_motif_3UTR[-which(human_motif_3UTR %in% mirlist[which(mirlist$organism=='Homo sapiens'),'motif'])]))

# mouse
mouse_motif_3UTR <- as.vector(as.matrix(motif_3UTR_species[which(motif_3UTR_species[,7]=='mouse'),2:6]))
mouse_motif_3UTR <- mouse_motif_3UTR[which(mouse_motif_3UTR!='')]
mouse_mir_3UTR <- as.data.frame(table(mouse_motif_3UTR[which(mouse_motif_3UTR %in% mirlist[which(mirlist$organism=='Mus musculus'),'motif'])]))
mouse_nonmir_3UTR <- as.data.frame(table(mouse_motif_3UTR[-which(mouse_motif_3UTR %in% mirlist[which(mirlist$organism=='Mus musculus'),'motif'])]))
nrow(human_mir_3UTR)+nrow(mouse_mir_3UTR)
nrow(human_nonmir_3UTR)+nrow(mouse_nonmir_3UTR)
motif_summary_compare[2,] <- c(nrow(human_mir_5UTR)+nrow(mouse_mir_5UTR),
                               sum(human_mir_5UTR$Freq) + sum(mouse_mir_5UTR$Freq),
                               nrow(human_mir_3UTR)+nrow(mouse_mir_3UTR),
                               sum(human_mir_3UTR$Freq) + sum(mouse_mir_3UTR$Freq))





library(spgs)
library(stringr)
library(data.table)
library(readr)
# reverse complement
mirlist$motif <- reverseComplement(substring(mirlist[,2],2,8))
mirlist$motif <- toupper(mirlist$motif)


motif_5UTR_7mer <- data.frame(matrix(ncol = 6))
topN = 5
for(i in 1:as.integer(system('ls /home/s2172876/mirna/motif_result/size_7/sylamer_result_5UTR | wc -l',TRUE))){
    file = paste0('/home/s2172876/mirna/motif_result/size_7/sylamer_result_5UTR/Exp',i,'_5UTR.out')
    # Read table of sylamer
    tab <- cbind(0,data.frame(fread(file,header = TRUE), row.names=1))
    # ordering sylamer tab
    mycols <- 1:ncol(tab) # columns to consider in the ordering
    tab <- tab[order(apply(abs(tab[,mycols]), 1, max), decreasing=TRUE),] # ordering table based in max enrichment
    #   top5 motif
    motif_5UTR_7mer[i,] <- c(paste('Exp',i,sep = ''),rownames(tab)[1:topN])
    print(i)
}

motif_3UTR_7mer <- data.frame(matrix(ncol = 6))
for(i in 1:as.integer(system('ls /home/s2172876/mirna/motif_result/size_7/sylamer_result_3UTR | wc -l',TRUE))){
    file = paste0('/home/s2172876/mirna/motif_result/size_7/sylamer_result_3UTR/Exp',i,'_3UTR.out')
    # Read table of sylamer
    tab <- cbind(0,data.frame(fread(file,header = TRUE), row.names=1))
    #   tab <- abs(tab)
    # ordering sylamer tab
    mycols <- 1:ncol(tab) # columns to consider in the ordering
    tab <- tab[order(apply(abs(tab[,mycols]), 1, max), decreasing=TRUE),] # ordering table based in max enrichment
    #   top5 motif
    motif_3UTR_7mer[i,] <- c(paste('Exp',i,sep = ''),rownames(tab)[1:topN])
    print(i)
}
motif_5UTR_species <- merge(motif_5UTR_7mer,exp_summary[,c(2,3)],by.x = 'X1',by.y = 'ExpID')
motif_3UTR_species <- merge(motif_3UTR_7mer,exp_summary[,c(2,3)],by.x = 'X1',by.y = 'ExpID')

# 5'UTR motif frequency
# human
human_motif_5UTR <- as.vector(as.matrix(motif_5UTR_species[which(motif_5UTR_species[,7]=='human'),2:6]))
human_motif_5UTR <- human_motif_5UTR[which(human_motif_5UTR!='')]
human_mir_5UTR <- as.data.frame(table(human_motif_5UTR[which(human_motif_5UTR %in% mirlist[which(mirlist$organism=='Homo sapiens'),'motif'])]))
human_nonmir_5UTR <- as.data.frame(table(human_motif_5UTR[-which(human_motif_5UTR %in% mirlist[which(mirlist$organism=='Homo sapiens'),'motif'])]))

# mouse
mouse_motif_5UTR <- as.vector(as.matrix(motif_5UTR_species[which(motif_5UTR_species[,7]=='mouse'),2:6]))
mouse_motif_5UTR <- mouse_motif_5UTR[which(mouse_motif_5UTR!='')]
mouse_mir_5UTR <- as.data.frame(table(mouse_motif_5UTR[which(mouse_motif_5UTR %in% mirlist[which(mirlist$organism=='Mus musculus'),'motif'])]))
mouse_nonmir_5UTR <- as.data.frame(table(mouse_motif_5UTR[-which(mouse_motif_5UTR %in% mirlist[which(mirlist$organism=='Mus musculus'),'motif'])]))

nrow(human_mir_5UTR)+nrow(mouse_mir_5UTR)
nrow(human_nonmir_5UTR)+nrow(mouse_nonmir_5UTR)

# 3'UTR motif frequency
# human
human_motif_3UTR <- as.vector(as.matrix(motif_3UTR_species[which(motif_3UTR_species[,7]=='human'),2:6]))
human_motif_3UTR <- human_motif_3UTR[which(human_motif_3UTR!='')]
human_mir_3UTR <- as.data.frame(table(human_motif_3UTR[which(human_motif_3UTR %in% mirlist[which(mirlist$organism=='Homo sapiens'),'motif'])]))
# mouse
mouse_motif_3UTR <- as.vector(as.matrix(motif_3UTR_species[which(motif_3UTR_species[,7]=='mouse'),2:6]))
mouse_motif_3UTR <- mouse_motif_3UTR[which(mouse_motif_3UTR!='')]
mouse_mir_3UTR <- as.data.frame(table(mouse_motif_3UTR[which(mouse_motif_3UTR %in% mirlist[which(mirlist$organism=='Mus musculus'),'motif'])]))
(sum(human_mir_3UTR$Freq) + sum(mouse_mir_3UTR$Freq))/(length(human_motif_3UTR)+length(mouse_motif_3UTR))
nrow(human_mir_3UTR)+nrow(mouse_mir_3UTR)
motif_summary_compare[3,] <- c(nrow(human_mir_5UTR)+nrow(mouse_mir_5UTR),
                               sum(human_mir_5UTR$Freq) + sum(mouse_mir_5UTR$Freq),
                               nrow(human_mir_3UTR)+nrow(mouse_mir_3UTR),
                               sum(human_mir_3UTR$Freq) + sum(mouse_mir_3UTR$Freq),
                               nrow(human_nonmir_5UTR)+nrow(mouse_nonmir_5UTR),
                               sum(human_nonmir_5UTR$Freq) + sum(mouse_nonmir_5UTR$Freq),
                               nrow(human_nonmir_3UTR)+nrow(mouse_nonmir_3UTR),
                               sum(human_nonmir_3UTR$Freq) + sum(mouse_nonmir_3UTR$Freq))


motif_5UTR_7mer_cutoff8 <- data.frame(matrix(ncol = 6))
for(i in 1:as.integer(system('ls /home/s2172876/mirna/motif_result/size_7/sylamer_result_5UTR | wc -l',TRUE))){
    file = paste0('/home/s2172876/mirna/motif_result/size_7/sylamer_result_5UTR/Exp',i,'_5UTR.out')
    # Read table of sylamer
    tab <- cbind(0,data.frame(fread(file,header = TRUE), row.names=1))
    # ordering sylamer tab
    mycols <- 1:ncol(tab) # columns to consider in the ordering
    tab <- tab[order(apply(abs(tab[,mycols]), 1, max), decreasing=TRUE),] # ordering table based in max enrichment
    # filter
    if(length(which(apply(abs(tab[1:5,2:(ncol(tab)-1)]), 1, max)<8))){
        m <- c(names(which(apply(abs(tab[1:5,2:(ncol(tab)-1)]), 1, max)>8)),
               rep('',5-length(which(apply(abs(tab[1:5,2:(ncol(tab)-1)]), 1, max)>8))))
    }else{
        m <- names(which(apply(abs(tab[1:5,2:(ncol(tab)-1)]), 1, max)>8))
    }
    motif_5UTR_7mer_cutoff8[i,] <- c(paste('Exp',i,sep = ''),m)
    print(i)
    print(m)
}


motif_3UTR_7mer_cutoff8 <- data.frame(matrix(ncol = 6))
for(i in 1:as.integer(system('ls /home/s2172876/mirna/motif_result/size_7/sylamer_result_3UTR | wc -l',TRUE))){
    file = paste0('/home/s2172876/mirna/motif_result/size_7/sylamer_result_3UTR/Exp',i,'_3UTR.out')
    # Read table of sylamer
    tab <- cbind(0,data.frame(fread(file,header = TRUE), row.names=1))
    # ordering sylamer tab
    mycols <- 1:ncol(tab) # columns to consider in the ordering
    tab <- tab[order(apply(abs(tab[,mycols]), 1, max), decreasing=TRUE),] # ordering table based in max enrichment
    # filter
    if(length(which(apply(abs(tab[1:5,2:(ncol(tab)-1)]), 1, max)<8))){
        m <- c(names(which(apply(abs(tab[1:5,2:(ncol(tab)-1)]), 1, max)>8)),
               rep('',5-length(which(apply(abs(tab[1:5,2:(ncol(tab)-1)]), 1, max)>8))))
    }else{
        m <- names(which(apply(abs(tab[1:5,2:(ncol(tab)-1)]), 1, max)>8))
    }
    motif_3UTR_7mer_cutoff8[i,] <- c(paste('Exp',i,sep = ''),m)
    print(i)
    print(m)
}

motif_5UTR_species <- merge(motif_5UTR_7mer_cutoff8,exp_summary[,c(2,3)],by.x = 'X1',by.y = 'ExpID')
motif_3UTR_species <- merge(motif_3UTR_7mer_cutoff8,exp_summary[,c(2,3)],by.x = 'X1',by.y = 'ExpID')

# 5'UTR motif frequency
# human
human_motif_5UTR <- as.vector(as.matrix(motif_5UTR_species[which(motif_5UTR_species[,7]=='human'),2:6]))
human_motif_5UTR <- human_motif_5UTR[which(human_motif_5UTR!='')]
human_mir_5UTR <- as.data.frame(table(human_motif_5UTR[which(human_motif_5UTR %in% mirlist[which(mirlist$organism=='Homo sapiens'),'motif'])]))
human_nonmir_5UTR <- as.data.frame(table(human_motif_5UTR[-which(human_motif_5UTR %in% mirlist[which(mirlist$organism=='Homo sapiens'),'motif'])]))

# mouse
mouse_motif_5UTR <- as.vector(as.matrix(motif_5UTR_species[which(motif_5UTR_species[,7]=='mouse'),2:6]))
mouse_motif_5UTR <- mouse_motif_5UTR[which(mouse_motif_5UTR!='')]
mouse_mir_5UTR <- as.data.frame(table(mouse_motif_5UTR[which(mouse_motif_5UTR %in% mirlist[which(mirlist$organism=='Mus musculus'),'motif'])]))
mouse_nonmir_5UTR <- as.data.frame(table(mouse_motif_5UTR[-which(mouse_motif_5UTR %in% mirlist[which(mirlist$organism=='Mus musculus'),'motif'])]))

nrow(human_mir_5UTR)+nrow(mouse_mir_5UTR)
nrow(human_nonmir_5UTR)+nrow(mouse_nonmir_5UTR)

# 3'UTR motif frequency
# human
human_motif_3UTR <- as.vector(as.matrix(motif_3UTR_species[which(motif_3UTR_species[,7]=='human'),2:6]))
human_motif_3UTR <- human_motif_3UTR[which(human_motif_3UTR!='')]
human_mir_3UTR <- as.data.frame(table(human_motif_3UTR[which(human_motif_3UTR %in% mirlist[which(mirlist$organism=='Homo sapiens'),'motif'])]))
human_nonmir_3UTR <- as.data.frame(table(human_motif_3UTR[-which(human_motif_3UTR %in% mirlist[which(mirlist$organism=='Homo sapiens'),'motif'])]))

# mouse
mouse_motif_3UTR <- as.vector(as.matrix(motif_3UTR_species[which(motif_3UTR_species[,7]=='mouse'),2:6]))
mouse_motif_3UTR <- mouse_motif_3UTR[which(mouse_motif_3UTR!='')]
mouse_mir_3UTR <- as.data.frame(table(mouse_motif_3UTR[which(mouse_motif_3UTR %in% mirlist[which(mirlist$organism=='Mus musculus'),'motif'])]))
mouse_nonmir_3UTR <- as.data.frame(table(mouse_motif_3UTR[-which(mouse_motif_3UTR %in% mirlist[which(mirlist$organism=='Mus musculus'),'motif'])]))

nrow(human_mir_3UTR)+nrow(mouse_mir_3UTR)
nrow(human_nonmir_3UTR)+nrow(mouse_nonmir_3UTR)

motif_summary_compare[4,] <- c(nrow(human_mir_5UTR)+nrow(mouse_mir_5UTR),
                               sum(human_mir_5UTR$Freq) + sum(mouse_mir_5UTR$Freq),
                               nrow(human_mir_3UTR)+nrow(mouse_mir_3UTR),
                               sum(human_mir_3UTR$Freq) + sum(mouse_mir_3UTR$Freq),
                               nrow(human_nonmir_5UTR)+nrow(mouse_nonmir_5UTR),
                               sum(human_nonmir_5UTR$Freq) + sum(mouse_nonmir_5UTR$Freq),
                               nrow(human_nonmir_3UTR)+nrow(mouse_nonmir_3UTR),
                               sum(human_nonmir_3UTR$Freq) + sum(mouse_nonmir_3UTR$Freq))

row.names(motif_summary_compare) <- c('6mer','6mer cutoff 8','7mer','7mer cutoff 8')
# save.image(file = "mirna.RData")
library(ggplot2)
library(ggbreak)
library(wesanderson)
# positions <- mir_motif[c(1,6,11,16,2,7,12,17,
#                          3,8,13,18,4,9,14,19,
#                          5,10,15,20),1]
motif_summary_compare_stack <- stack(motif_summary_compare)
motif_summary_compare_stack$Type <- rep(row.names(motif_summary_compare),4)
colnames(motif_summary_compare_stack)[2] <- 'Position'

motif_summary_compare_stack[,2] <- str_replace_all(motif_summary_compare_stack[,2],'miRNA-related motif','miRNA response element')
p1 <- ggplot(motif_summary_compare_stack[c(1:4,9:12),], aes(fill=Position, y=values, x= Type)) +
    geom_bar(position="dodge", stat="identity") +
    geom_col(width = 0.6, position = "dodge")+ ylim(0,2500)+
    # scale_y_cut(breaks=c(105,1000)) +
    geom_text(aes(label=values),vjust=-0.6,color = "black",size=4,position = position_dodge(width = .95),fontface = "bold")+
    theme_bw()+  theme(plot.title = element_text(hjust = 0.5),legend.position="bottom")+
    # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    xlab('') + ylab('The number of unique miRNA response element')+
    scale_fill_manual(values = mycolors[c(2,1)])#+labs(tag = "A")

p2 <- ggplot(motif_summary_compare_stack[c(5:8,13:16),], aes(fill=Position, y=values, x= Type)) +
    geom_bar(position="dodge", stat="identity") +
    geom_col(width = 0.6, position = "dodge")+ ylim(0,35000)+
    # scale_y_cut(breaks=c(6000,20000)) +
    geom_text(aes(label=values),vjust=-0.6,color = "black",size=4,position = position_dodge(width = .95),fontface = "bold")+
    theme_bw()+  theme(plot.title = element_text(hjust = 0.5),legend.position="bottom")+
    # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    xlab('') + ylab('miRNA response element occurance')+
    scale_fill_manual(values = mycolors[c(2,1)])#+labs(tag = "B")

mir_motif <- rbind(human_mir_3UTR[order(human_mir_3UTR$Freq,decreasing=TRUE),][1:15,],
                   mouse_mir_3UTR[order(mouse_mir_3UTR$Freq,decreasing=TRUE),][1:15,],
                   human_mir_5UTR[order(human_mir_5UTR$Freq,decreasing=TRUE),][1:15,],
                   mouse_mir_5UTR[order(mouse_mir_5UTR$Freq,decreasing=TRUE),][1:15,])
mir_motif$Type <- c(rep('human 3\'UTR',15),
                    rep('mouse 3\'UTR',15),
                    rep('human 5\'UTR',15),
                    rep('mouse 5\'UTR',15))
mir_motif[,1] <- as.character(mir_motif[,1])
mir_motif[c(1:15,31:45),1] <- paste(mir_motif[c(1:15,31:45),1] ,'(hsa)',sep='')
mir_motif[c(16:30,46:60),1] <- paste(mir_motif[c(16:30,46:60),1],'(mmu)',sep='')
mir_motif <- mir_motif[order(mir_motif$Type),]
mir_motif$order <-rep(1:15,4)
mir_motif[c(1:15,31:45),1]
p3 <- ggplot(mir_motif[c(1:10,31:40),], aes(x=order , y=Freq, group=Type)) +
    geom_line(aes(linetype=Type),color=mycolors[2])+
    geom_point(color=mycolors[2])+ylab('miRNA-related motif occurance')+xlab('Top 10 motifs')+
    theme_bw()+  theme(plot.title = element_text(hjust = 0.5),legend.position="bottom")+
    scale_x_continuous(breaks=seq(1,15,1)) +
    geom_hline(yintercept=41, linetype="dashed", color = "red")+
    annotate(geom="text", x=9.5, y=180, label="y=41",color="red")+
    annotate(geom="text", x=1.8, y=1650, label=mir_motif[1,1],color='black')+ 
    annotate(geom="text", x=1.8, y=1265, label=mir_motif[31,1],color='black')
    # scale_color_manual(values=mycolors[c(1,3)])#+labs(tag = "C")


p4 <- ggplot(mir_motif[c(16:25,46:55),], aes(x=order , y=Freq, group=Type)) +
    geom_line(aes(linetype=Type),color=mycolors[1])+
    geom_point(color=mycolors[1])+ylab('miRNA-related motif occurance')+xlab('Top 15 motifs')+
    theme_bw()+  theme(plot.title = element_text(hjust = 0.5),legend.position="bottom")+
    scale_x_continuous(breaks=seq(1,15,1))+
    geom_hline(yintercept=4, linetype="dashed", color = "red")+  
    annotate(geom="text", x=9.5, y=500, label="y=4",color="red")+
    annotate(geom="text", x=1.9, y=3700, label=mir_motif[46,1],color='black')+ 
    annotate(geom="text", x=1.9, y=1260, label=mir_motif[16,1],color='black')
    # scale_color_manual(values=mycolors[c(1,3)])#+labs(tag = "D")
library(ggpubr)
ggarrange(p1, p2, p3, p4, 
          labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2)


# exclude CG rich
motif_5UTR_species <- merge(motif_5UTR_cutoff8,exp_summary[,c(2,3)],by.x = 'X1',by.y = 'ExpID')
motif_3UTR_species <- merge(motif_3UTR_cutoff8,exp_summary[,c(2,3)],by.x = 'X1',by.y = 'ExpID')

# 5'UTR motif frequency
# human
human_motif_5UTR <- as.vector(as.matrix(motif_5UTR_species[which(motif_5UTR_species[,7]=='human'),2:6]))
human_motif_5UTR <- human_motif_5UTR[which(human_motif_5UTR!='')]
human_mir_5UTR <- as.data.frame(table(human_motif_5UTR[which(human_motif_5UTR %in% mirlist[which(mirlist$organism=='Homo sapiens'),'motif'])]))
human_nonmir_5UTR <- as.data.frame(table(human_motif_5UTR[-which(human_motif_5UTR %in% mirlist[which(mirlist$organism=='Homo sapiens'),'motif'])]))

# mouse
mouse_motif_5UTR <- as.vector(as.matrix(motif_5UTR_species[which(motif_5UTR_species[,7]=='mouse'),2:6]))
mouse_motif_5UTR <- mouse_motif_5UTR[which(mouse_motif_5UTR!='')]
mouse_mir_5UTR <- as.data.frame(table(mouse_motif_5UTR[which(mouse_motif_5UTR %in% mirlist[which(mirlist$organism=='Mus musculus'),'motif'])]))
mouse_nonmir_5UTR <- as.data.frame(table(mouse_motif_5UTR[-which(mouse_motif_5UTR %in% mirlist[which(mirlist$organism=='Mus musculus'),'motif'])]))

nrow(human_mir_5UTR)+nrow(mouse_mir_5UTR)
nrow(human_nonmir_5UTR)+nrow(mouse_nonmir_5UTR)

# 3'UTR motif frequency
# human
human_motif_3UTR <- as.vector(as.matrix(motif_3UTR_species[which(motif_3UTR_species[,7]=='human'),2:6]))
human_motif_3UTR <- human_motif_3UTR[which(human_motif_3UTR!='')]
human_mir_3UTR <- as.data.frame(table(human_motif_3UTR[which(human_motif_3UTR %in% mirlist[which(mirlist$organism=='Homo sapiens'),'motif'])]))
human_nonmir_3UTR <- as.data.frame(table(human_motif_3UTR[-which(human_motif_3UTR %in% mirlist[which(mirlist$organism=='Homo sapiens'),'motif'])]))

# mouse
mouse_motif_3UTR <- as.vector(as.matrix(motif_3UTR_species[which(motif_3UTR_species[,7]=='mouse'),2:6]))
mouse_motif_3UTR <- mouse_motif_3UTR[which(mouse_motif_3UTR!='')]
mouse_mir_3UTR <- as.data.frame(table(mouse_motif_3UTR[which(mouse_motif_3UTR %in% mirlist[which(mirlist$organism=='Mus musculus'),'motif'])]))
mouse_nonmir_3UTR <- as.data.frame(table(mouse_motif_3UTR[-which(mouse_motif_3UTR %in% mirlist[which(mirlist$organism=='Mus musculus'),'motif'])]))
nrow(human_mir_3UTR)+nrow(mouse_mir_3UTR)
nrow(human_nonmir_3UTR)+nrow(mouse_nonmir_3UTR)


human_mir_5UTR$status <- rep('non',nrow(human_mir_5UTR))
human_mir_5UTR$Var1 <- as.vector(human_mir_5UTR$Var1)
human_mir_5UTR <- human_mir_5UTR[which(human_mir_5UTR$Var1!=''),]
for (i in 1:nrow(human_mir_5UTR)) {
  C = lengths(regmatches(human_mir_5UTR[i,1], gregexpr("C", human_mir_5UTR[i,1])))
  G = lengths(regmatches(human_mir_5UTR[i,1], gregexpr("G", human_mir_5UTR[i,1])))
  GC <- C+G
  CpG <- lengths(regmatches(human_mir_5UTR[i,1], gregexpr("CG", human_mir_5UTR[i,1])))
  if(CpG!=0){
    exp_1 = C*G/6
    exp_2 = ((C+G)/2)^2/6
    if(GC>=3 & (CpG/exp_1>0.6|(CpG/exp_2>0.6))){
      human_mir_5UTR[i,3] <- 'CpG motif'
    }
  }
}
human_mir_5UTR[grepl("TAAAT", human_mir_5UTR$Var1,TRUE),3] <- 'ARE' 
sum(human_mir_5UTR[which(human_mir_5UTR[,3]=='non'),2])

mouse_mir_5UTR$status <- rep('non',nrow(mouse_mir_5UTR))
mouse_mir_5UTR$Var1 <- as.vector(mouse_mir_5UTR$Var1)
mouse_mir_5UTR <- mouse_mir_5UTR[which(mouse_mir_5UTR$Var1!=''),]
for (i in 1:nrow(mouse_mir_5UTR)) {
  C = lengths(regmatches(mouse_mir_5UTR[i,1], gregexpr("C", mouse_mir_5UTR[i,1])))
  G = lengths(regmatches(mouse_mir_5UTR[i,1], gregexpr("G", mouse_mir_5UTR[i,1])))
  GC <- C+G
  CpG <- lengths(regmatches(mouse_mir_5UTR[i,1], gregexpr("CG", mouse_mir_5UTR[i,1])))
  if(CpG!=0){
    exp_1 = C*G/6
    exp_2 = ((C+G)/2)^2/6
    if(GC>=3 & (CpG/exp_1>0.6|(CpG/exp_2>0.6))){
      mouse_mir_5UTR[i,3] <- 'CpG motif'
    }
  }
}
mouse_mir_5UTR[grepl("ATTTA", mouse_mir_5UTR$Var1,TRUE),3] <- 'ARE' 
sum(mouse_mir_5UTR[which(mouse_mir_5UTR[,3]=='non'),2])
sum(human_mir_5UTR[which(human_mir_5UTR[,3]=='non'),2]) + sum(mouse_mir_5UTR[which(mouse_mir_5UTR[,3]=='non'),2])
nrow(human_mir_5UTR[which(human_mir_5UTR[,3]=='non'),]) + nrow(mouse_mir_5UTR[which(mouse_mir_5UTR[,3]=='non'),])



human_mir_3UTR$status <- rep('non',nrow(human_mir_3UTR))
human_mir_3UTR$Var1 <- as.vector(human_mir_3UTR$Var1)
human_mir_3UTR <- human_mir_3UTR[which(human_mir_3UTR$Var1!=''),]
for (i in 1:nrow(human_mir_3UTR)) {
  C = lengths(regmatches(human_mir_3UTR[i,1], gregexpr("C", human_mir_3UTR[i,1])))
  G = lengths(regmatches(human_mir_3UTR[i,1], gregexpr("G", human_mir_3UTR[i,1])))
  GC <- C+G
  CpG <- lengths(regmatches(human_mir_3UTR[i,1], gregexpr("CG", human_mir_3UTR[i,1])))
  if(CpG!=0){
    exp_1 = C*G/6
    exp_2 = ((C+G)/2)^2/6
    if(GC>=3 & (CpG/exp_1>0.6|(CpG/exp_2>0.6))){
      human_mir_3UTR[i,3] <- 'CpG motif'
    }
  }
}
human_mir_3UTR[grepl("TAAAT", human_mir_3UTR$Var1,TRUE),3] <- 'ARE' 
sum(human_mir_3UTR[which(human_mir_3UTR[,3]=='non'),2])


mouse_mir_3UTR
mouse_mir_3UTR$status <- rep('non',nrow(mouse_mir_3UTR))
mouse_mir_3UTR$Var1 <- as.vector(mouse_mir_3UTR$Var1)
mouse_mir_3UTR <- mouse_mir_3UTR[which(mouse_mir_3UTR$Var1!=''),]
for (i in 1:nrow(mouse_mir_3UTR)) {
  C = lengths(regmatches(mouse_mir_3UTR[i,1], gregexpr("C", mouse_mir_3UTR[i,1])))
  G = lengths(regmatches(mouse_mir_3UTR[i,1], gregexpr("G", mouse_mir_3UTR[i,1])))
  GC <- C+G
  CpG <- lengths(regmatches(mouse_mir_3UTR[i,1], gregexpr("CG", mouse_mir_3UTR[i,1])))
  if(CpG!=0){
    exp_1 = C*G/6
    exp_2 = ((C+G)/2)^2/6
    if(GC>=3 & (CpG/exp_1>0.6|(CpG/exp_2>0.6))){
      mouse_mir_3UTR[i,3] <- 'CpG motif'
    }
  }
}
mouse_mir_3UTR[grepl("TAAAT", mouse_mir_3UTR$Var1,TRUE),3] <- 'ARE'
sum(mouse_mir_3UTR[which(mouse_mir_3UTR[,3]=='non'),2]) + sum(human_mir_3UTR[which(human_mir_3UTR[,3]=='non'),2])

# mouse_mir_3UTR,human_mir_3UTR,human_mir_5UTR,mouse_mir_5UTR
mirmotif_type <- as.data.frame(matrix(ncol = 3,nrow = 9))
mirmotif_type[,1] <- c(sum(human_mir_5UTR[which(human_mir_5UTR$status=='CpG motif'),2]),
  sum(human_mir_5UTR[which(human_mir_5UTR$status=='non'),2]),
  sum(mouse_mir_5UTR[which(mouse_mir_5UTR$status=='CpG motif'),2]),
  sum(mouse_mir_5UTR[which(mouse_mir_5UTR$status=='non'),2]),
  sum(human_mir_3UTR[which(human_mir_3UTR$status=='CpG motif'),2]),
  sum(human_mir_3UTR[which(human_mir_3UTR$status=='non'),2]),
  sum(mouse_mir_3UTR[which(mouse_mir_3UTR$status=='CpG motif'),2]),
  sum(mouse_mir_3UTR[which(mouse_mir_3UTR$status=='non'),2]),
  sum(mouse_mir_3UTR[which(mouse_mir_3UTR$status=='ARE'),2]))
mirmotif_type[,2] <- c(rep(c('CpG','Other'),4),'ARE')
mirmotif_type[,3] <- c(rep('human 5\'UTR',2),
                       rep('mouse 5\'UTR',2),
                       rep('human 3\'UTR',2),
                       rep('mouse 3\'UTR',3))
colnames(mirmotif_type) <- c('Num','Type','Region')
library(ggplot2)
library(dplyr)

# https://r-charts.com/part-whole/pie-chart-labels-outside-ggplot2/
# Get the positions
df2 <- mirmotif_type[1:2,1:2] %>% 
  mutate(csum = rev(cumsum(rev(Num))), 
         pos = Num/2 + lead(csum, 1),
         pos = if_else(is.na(pos), Num/2, pos))
df2$per <- round(df2$Num/sum(df2$Num),3)
pie1 <- ggplot(df2, aes(x = "" , y = Num, fill = fct_inorder(Type))) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  geom_label_repel(data = df2,
                   aes(y = pos, label = paste0(per*100, "%")),
                   size = 4.5, nudge_x = 1, show.legend = FALSE) +
  guides(fill = guide_legend(title = "Category")) +
  theme(plot.title = element_text(size=6)) +
  ggtitle("Predicted human 5'UTR \n miRNA response element") +
  theme_void()+
  theme(plot.title = element_text(hjust = 0.5))

df2 <- mirmotif_type[3:4,1:2] %>% 
  mutate(csum = rev(cumsum(rev(Num))), 
         pos = Num/2 + lead(csum, 1),
         pos = if_else(is.na(pos), Num/2, pos))
df2$per <- round(df2$Num/sum(df2$Num),3)
pie2 <- ggplot(df2, aes(x = "" , y = Num, fill = fct_inorder(Type))) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  geom_label_repel(data = df2,
                   aes(y = pos, label = paste0(per*100, "%")),
                   size = 4.5, nudge_x = 1, show.legend = FALSE) +
  guides(fill = guide_legend(title = "Category")) +
  theme(plot.title = element_text(size=6)) +
  ggtitle("Predicted mouse 5'UTR \n miRNA response element") +
  theme_void()+
  theme(plot.title = element_text(hjust = 0.5))

df2 <- mirmotif_type[5:6,1:2] %>% 
  mutate(csum = rev(cumsum(rev(Num))), 
         pos = Num/2 + lead(csum, 1),
         pos = if_else(is.na(pos), Num/2, pos))
df2$per <- round(df2$Num/sum(df2$Num),3)
pie3 <- ggplot(df2, aes(x = "" , y = Num, fill = fct_inorder(Type))) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  geom_label_repel(data = df2,
                   aes(y = pos, label = paste0(per*100, "%")),
                   size = 4.5, nudge_x = 1, show.legend = FALSE) +
  guides(fill = guide_legend(title = "Category")) +
  theme(plot.title = element_text(size=6)) +
  ggtitle("Predicted human 3'UTR \n miRNA response element") +
  theme_void()+
  theme(plot.title = element_text(hjust = 0.5))


# Get the positions
mirmotif_type[8,1] <-  mirmotif_type[8,1] +mirmotif_type[9,1]
df2 <- mirmotif_type[7:8,1:2] %>% 
  mutate(csum = rev(cumsum(rev(Num))), 
         pos = Num/2 + lead(csum, 1),
         pos = if_else(is.na(pos), Num/2, pos))

df2$per <- round(df2$Num/sum(df2$Num),3)

pie4 <- ggplot(df2, aes(x = "" , y = Num, fill = fct_inorder(Type))) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  geom_label_repel(data = df2,
                   aes(y = pos, label = paste0(per*100, "%")),
                   size = 4.5, nudge_x = 1, show.legend = FALSE) +
  guides(fill = guide_legend(title = "Category")) +
  theme(plot.title = element_text(size=6)) +
  ggtitle("Predicted mouse 3'UTR \n miRNA response element") +
  theme_void()+
  theme(plot.title = element_text(hjust = 0.5))

ggarrange(pie1, pie2, pie3, pie4, 
          labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2)



# human_nonmir_5UTR$status <- rep('non',nrow(human_nonmir_5UTR))
# human_nonmir_5UTR$Var1 <- as.vector(human_nonmir_5UTR$Var1)
# human_nonmir_5UTR <- human_nonmir_5UTR[which(human_nonmir_5UTR$Var1!=''),]
# for (i in 1:nrow(human_nonmir_5UTR)) {
#   C = lengths(regmatches(human_nonmir_5UTR[i,1], gregexpr("C", human_nonmir_5UTR[i,1])))
#   G = lengths(regmatches(human_nonmir_5UTR[i,1], gregexpr("G", human_nonmir_5UTR[i,1])))
#   GC <- C+G
#   CpG <- lengths(regmatches(human_nonmir_5UTR[i,1], gregexpr("CG", human_nonmir_5UTR[i,1])))
#   if(CpG!=0){
#     exp_1 = C*G/6
#     exp_2 = ((C+G)/2)^2/6
#     if(GC>=3 & (CpG/exp_1>0.6|(CpG/exp_2>0.6))){
#       human_nonmir_5UTR[i,3] <- 'CpG motif'
#     }
#   }
# }
# human_nonmir_5UTR[grepl("TAAAT", human_nonmir_5UTR$Var1,TRUE),3] <- 'ARE'
# sum(human_nonmir_5UTR[which(human_nonmir_5UTR[,3]=='non'),2])

p5_df <- as.data.frame(matrix(ncol = 2,nrow = 2))
p5_df[,1] <- c(nrow(human_mir_5UTR[human_mir_5UTR$status=='non',])+nrow(mouse_mir_5UTR[mouse_mir_5UTR$status=='non',]),
               nrow(human_mir_3UTR[human_mir_3UTR$status=='non',])+nrow(mouse_mir_3UTR[mouse_mir_3UTR$status=='non',]))
p5_df[,2] <- c('5\'UTR','3\'UTR')
p5_df[,3] <- rep('miRNA-related motif',2)
colnames(p5_df) <- c('Num','Region','Typ')
p5 <-ggplot(p5_df, aes(fill=Region, y=Num, x= Region)) +
  geom_bar(position="dodge", stat="identity") +
  geom_col(width = 0.6, position = "dodge")+
  # scale_y_cut(breaks=c(15,100)) +
  geom_text(aes(label=Num),vjust=-0.5,color = "black",size=4,fontface = "bold")+
  theme_bw()+  theme(plot.title = element_text(hjust = 0.5),legend.position="bottom")+
  # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab('') + ylab('The number of unique miRNA-related motif')+
  scale_fill_manual(values = mycolors[c(2,1)]) +labs(tag = "A")
p6_df <- data.frame(Num = c(sum(mirmotif_type[which(mirmotif_type$Type=='Other'),c(1,3)][1:2,1]),
                            sum(mirmotif_type[which(mirmotif_type$Type=='Other'),c(1,3)][3:4,1])),
                    Region = c('5\'UTR','3\'UTR'),Type=rep('miRNA-related motif',2))

p6 <- ggplot(p6_df, aes(fill=Region, y=Num, x= Region)) +
  geom_bar(position="dodge", stat="identity") +
  geom_col(width = 0.6, position = "dodge")+ ylim(0,8800)+
  geom_text(aes(label=Num),vjust=-0.5,color = "black",size=4,fontface = "bold")+
  theme_bw()+  theme(plot.title = element_text(hjust = 0.5),legend.position="bottom")+
  # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab('') + ylab('miRNA-related motif occurance')+
  scale_fill_manual(values = mycolors[c(2,1)])+labs(tag = "B")
library(patchwork)
p5+p6


# exclude CG rich (7mer)
motif_5UTR_species <- merge(motif_5UTR_7mer_cutoff8,exp_summary[,c(2,3)],by.x = 'X1',by.y = 'ExpID')
motif_3UTR_species <- merge(motif_3UTR_7mer_cutoff8,exp_summary[,c(2,3)],by.x = 'X1',by.y = 'ExpID')
mirlist$motif <- reverseComplement(substring(mirlist[,2],2,8))
mirlist$motif <- toupper(mirlist$motif)

# 5'UTR motif frequency
# human
human_motif_5UTR <- as.vector(as.matrix(motif_5UTR_species[which(motif_5UTR_species[,7]=='human'),2:6]))
human_motif_5UTR <- human_motif_5UTR[which(human_motif_5UTR!='')]
human_mir_5UTR <- as.data.frame(table(human_motif_5UTR[which(human_motif_5UTR %in% mirlist[which(mirlist$organism=='Homo sapiens'),'motif'])]))
human_nonmir_5UTR <- as.data.frame(table(human_motif_5UTR[-which(human_motif_5UTR %in% mirlist[which(mirlist$organism=='Homo sapiens'),'motif'])]))

# mouse
mouse_motif_5UTR <- as.vector(as.matrix(motif_5UTR_species[which(motif_5UTR_species[,7]=='mouse'),2:6]))
mouse_motif_5UTR <- mouse_motif_5UTR[which(mouse_motif_5UTR!='')]
mouse_mir_5UTR <- as.data.frame(table(mouse_motif_5UTR[which(mouse_motif_5UTR %in% mirlist[which(mirlist$organism=='Mus musculus'),'motif'])]))
mouse_nonmir_5UTR <- as.data.frame(table(mouse_motif_5UTR[-which(mouse_motif_5UTR %in% mirlist[which(mirlist$organism=='Mus musculus'),'motif'])]))

nrow(human_mir_5UTR)+nrow(mouse_mir_5UTR)
nrow(human_nonmir_5UTR)+nrow(mouse_nonmir_5UTR)

# 3'UTR motif frequency
# human
human_motif_3UTR <- as.vector(as.matrix(motif_3UTR_species[which(motif_3UTR_species[,7]=='human'),2:6]))
human_motif_3UTR <- human_motif_3UTR[which(human_motif_3UTR!='')]
human_mir_3UTR <- as.data.frame(table(human_motif_3UTR[which(human_motif_3UTR %in% mirlist[which(mirlist$organism=='Homo sapiens'),'motif'])]))
human_nonmir_3UTR <- as.data.frame(table(human_motif_3UTR[-which(human_motif_3UTR %in% mirlist[which(mirlist$organism=='Homo sapiens'),'motif'])]))

# mouse
mouse_motif_3UTR <- as.vector(as.matrix(motif_3UTR_species[which(motif_3UTR_species[,7]=='mouse'),2:6]))
mouse_motif_3UTR <- mouse_motif_3UTR[which(mouse_motif_3UTR!='')]
mouse_mir_3UTR <- as.data.frame(table(mouse_motif_3UTR[which(mouse_motif_3UTR %in% mirlist[which(mirlist$organism=='Mus musculus'),'motif'])]))
mouse_nonmir_3UTR <- as.data.frame(table(mouse_motif_3UTR[-which(mouse_motif_3UTR %in% mirlist[which(mirlist$organism=='Mus musculus'),'motif'])]))
nrow(human_mir_3UTR)+nrow(mouse_mir_3UTR)
nrow(human_nonmir_3UTR)+nrow(mouse_nonmir_3UTR)


human_mir_5UTR$status <- rep('non',nrow(human_mir_5UTR))
human_mir_5UTR$Var1 <- as.vector(human_mir_5UTR$Var1)
human_mir_5UTR <- human_mir_5UTR[which(human_mir_5UTR$Var1!=''),]
for (i in 1:nrow(human_mir_5UTR)) {
  C = lengths(regmatches(human_mir_5UTR[i,1], gregexpr("C", human_mir_5UTR[i,1])))
  G = lengths(regmatches(human_mir_5UTR[i,1], gregexpr("G", human_mir_5UTR[i,1])))
  GC <- C+G
  CpG <- lengths(regmatches(human_mir_5UTR[i,1], gregexpr("CG", human_mir_5UTR[i,1])))
  if(CpG!=0){
    exp_1 = C*G/6
    exp_2 = ((C+G)/2)^2/6
    if(GC>=3 & (CpG/exp_1>0.6|(CpG/exp_2>0.6))){
      human_mir_5UTR[i,3] <- 'CpG motif'
    }
  }
}
human_mir_5UTR[grepl("TAAAT", human_mir_5UTR$Var1,TRUE),3] <- 'ARE' 
sum(human_mir_5UTR[which(human_mir_5UTR[,3]=='non'),2])

mouse_mir_5UTR$status <- rep('non',nrow(mouse_mir_5UTR))
mouse_mir_5UTR$Var1 <- as.vector(mouse_mir_5UTR$Var1)
mouse_mir_5UTR <- mouse_mir_5UTR[which(mouse_mir_5UTR$Var1!=''),]
for (i in 1:nrow(mouse_mir_5UTR)) {
  C = lengths(regmatches(mouse_mir_5UTR[i,1], gregexpr("C", mouse_mir_5UTR[i,1])))
  G = lengths(regmatches(mouse_mir_5UTR[i,1], gregexpr("G", mouse_mir_5UTR[i,1])))
  GC <- C+G
  CpG <- lengths(regmatches(mouse_mir_5UTR[i,1], gregexpr("CG", mouse_mir_5UTR[i,1])))
  if(CpG!=0){
    exp_1 = C*G/6
    exp_2 = ((C+G)/2)^2/6
    if(GC>=3 & (CpG/exp_1>0.6|(CpG/exp_2>0.6))){
      mouse_mir_5UTR[i,3] <- 'CpG motif'
    }
  }
}
mouse_mir_5UTR[grepl("ATTTA", mouse_mir_5UTR$Var1,TRUE),3] <- 'ARE' 
sum(mouse_mir_5UTR[which(mouse_mir_5UTR[,3]=='non'),2])
sum(human_mir_5UTR[which(human_mir_5UTR[,3]=='non'),2]) + sum(mouse_mir_5UTR[which(mouse_mir_5UTR[,3]=='non'),2])
nrow(human_mir_5UTR[which(human_mir_5UTR[,3]=='non'),]) + nrow(mouse_mir_5UTR[which(mouse_mir_5UTR[,3]=='non'),])



human_mir_3UTR$status <- rep('non',nrow(human_mir_3UTR))
human_mir_3UTR$Var1 <- as.vector(human_mir_3UTR$Var1)
human_mir_3UTR <- human_mir_3UTR[which(human_mir_3UTR$Var1!=''),]
for (i in 1:nrow(human_mir_3UTR)) {
  C = lengths(regmatches(human_mir_3UTR[i,1], gregexpr("C", human_mir_3UTR[i,1])))
  G = lengths(regmatches(human_mir_3UTR[i,1], gregexpr("G", human_mir_3UTR[i,1])))
  GC <- C+G
  CpG <- lengths(regmatches(human_mir_3UTR[i,1], gregexpr("CG", human_mir_3UTR[i,1])))
  if(CpG!=0){
    exp_1 = C*G/6
    exp_2 = ((C+G)/2)^2/6
    if(GC>=3 & (CpG/exp_1>0.6|(CpG/exp_2>0.6))){
      human_mir_3UTR[i,3] <- 'CpG motif'
    }
  }
}
human_mir_3UTR[grepl("ATTTA", human_mir_3UTR$Var1,TRUE),3] <- 'ARE' 
sum(human_mir_3UTR[which(human_mir_3UTR[,3]=='non'),2])


mouse_mir_3UTR
mouse_mir_3UTR$status <- rep('non',nrow(mouse_mir_3UTR))
mouse_mir_3UTR$Var1 <- as.vector(mouse_mir_3UTR$Var1)
mouse_mir_3UTR <- mouse_mir_3UTR[which(mouse_mir_3UTR$Var1!=''),]
for (i in 1:nrow(mouse_mir_3UTR)) {
  C = lengths(regmatches(mouse_mir_3UTR[i,1], gregexpr("C", mouse_mir_3UTR[i,1])))
  G = lengths(regmatches(mouse_mir_3UTR[i,1], gregexpr("G", mouse_mir_3UTR[i,1])))
  GC <- C+G
  CpG <- lengths(regmatches(mouse_mir_3UTR[i,1], gregexpr("CG", mouse_mir_3UTR[i,1])))
  if(CpG!=0){
    exp_1 = C*G/6
    exp_2 = ((C+G)/2)^2/6
    if(GC>=3 & (CpG/exp_1>0.6|(CpG/exp_2>0.6))){
      mouse_mir_3UTR[i,3] <- 'CpG motif'
    }
  }
}
mouse_mir_3UTR[grepl("TAAAT", mouse_mir_3UTR$Var1,TRUE),3] <- 'ARE'
sum(mouse_mir_3UTR[which(mouse_mir_3UTR[,3]=='non'),2]) + sum(human_mir_3UTR[which(human_mir_3UTR[,3]=='non'),2])

# mouse_mir_3UTR,human_mir_3UTR,human_mir_5UTR,mouse_mir_5UTR
mirmotif_type <- as.data.frame(matrix(ncol = 3,nrow = 9))
mirmotif_type[,1] <- c(sum(human_mir_5UTR[which(human_mir_5UTR$status=='CpG motif'),2]),
                       sum(human_mir_5UTR[which(human_mir_5UTR$status=='non'),2]),
                       sum(mouse_mir_5UTR[which(mouse_mir_5UTR$status=='CpG motif'),2]),
                       sum(mouse_mir_5UTR[which(mouse_mir_5UTR$status=='non'),2]),
                       sum(human_mir_3UTR[which(human_mir_3UTR$status=='CpG motif'),2]),
                       sum(human_mir_3UTR[which(human_mir_3UTR$status=='non'),2]),
                       sum(mouse_mir_3UTR[which(mouse_mir_3UTR$status=='CpG motif'),2]),
                       sum(mouse_mir_3UTR[which(mouse_mir_3UTR$status=='non'),2]),
                       sum(mouse_mir_3UTR[which(mouse_mir_3UTR$status=='ARE'),2]))
mirmotif_type[,2] <- c(rep(c('CpG','Other'),4),'ARE')
mirmotif_type[,3] <- c(rep('human 5\'UTR',2),
                       rep('mouse 5\'UTR',2),
                       rep('human 3\'UTR',2),
                       rep('mouse 3\'UTR',3))
colnames(mirmotif_type) <- c('Num','Type','Region')

p7_df <- as.data.frame(matrix(ncol = 2,nrow = 2))
p7_df[,1] <- c(nrow(human_mir_5UTR[human_mir_5UTR$status=='non',])+nrow(mouse_mir_5UTR[mouse_mir_5UTR$status=='non',]),
               nrow(human_mir_3UTR[human_mir_3UTR$status=='non',])+nrow(mouse_mir_3UTR[mouse_mir_3UTR$status=='non',]))
p7_df[,2] <- c('5\'UTR','3\'UTR')
p7_df[,3] <- rep('miRNA-related motif',2)
colnames(p7_df) <- c('Num','Region','Typ')
p7 <-ggplot(p7_df, aes(fill=Region, y=Num, x= Region)) +
  geom_bar(position="dodge", stat="identity") +
  geom_col(width = 0.6, position = "dodge")+
  # scale_y_cut(breaks=c(15,100)) +
  geom_text(aes(label=Num),vjust=-0.5,color = "black",size=4,fontface = "bold")+
  theme_bw()+  theme(plot.title = element_text(hjust = 0.5),legend.position="bottom")+
  # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab('') + ylab('The number of unique miRNA-related motif')+
  scale_fill_manual(values = mycolors[c(2,1)]) +labs(tag = "A")
p8_df <- data.frame(Num = c(sum(mirmotif_type[which(mirmotif_type$Type=='Other'),c(1,3)][1:2,1]),
                            sum(mirmotif_type[which(mirmotif_type$Type=='Other'),c(1,3)][3:4,1])),
                    Region = c('5\'UTR','3\'UTR'),Type=rep('miRNA-related motif',2))

p8 <- ggplot(p8_df, aes(fill=Region, y=Num, x= Region)) +
  geom_bar(position="dodge", stat="identity") +
  geom_col(width = 0.6, position = "dodge")+ ylim(0,4000)+
  geom_text(aes(label=Num),vjust=-0.5,color = "black",size=4,fontface = "bold")+
  theme_bw()+  theme(plot.title = element_text(hjust = 0.5),legend.position="bottom")+
  # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab('') + ylab('miRNA-related motif occurance')+
  scale_fill_manual(values = mycolors[c(2,1)])+labs(tag = "B")
library(patchwork)
p7+p8
remove_cpg <- data.frame(rbindlist(list(p5_df,p6_df,p7_df,p8_df)))
remove_cpg$Source <- c(rep('6mer cutoff8',2),
                       rep('6mer cutoff8 \n frequency',2),
                       rep('7mer cutoff8',2),
                       rep('7mer cutoff8 \n frequency',2))
ggplot(remove_cpg[c(1,2,5,6),], aes(fill=Region, y=Num, x= Source)) +
  geom_bar(position="dodge", stat="identity") +
  geom_col(width = 0.6, position = "dodge")+ 
  geom_text(aes(label=Num),vjust=-0.6,color = "black",size=4,position = position_dodge(width = .95),fontface = "bold")+
  theme_bw()+  theme(plot.title = element_text(hjust = 0.5),legend.position="bottom")+
  # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab('') + ylab('Number of unique miRNA response element')+
  scale_fill_manual(values = mycolors[c(2,1)]) +
ggplot(remove_cpg[c(3,4,7,8),], aes(fill=Region, y=Num, x= Source)) +
  geom_bar(position="dodge", stat="identity") +
  geom_col(width = 0.6, position = "dodge")+ 
  geom_text(aes(label=Num),vjust=-0.6,color = "black",size=4,position = position_dodge(width = .95),fontface = "bold")+
  theme_bw()+  theme(plot.title = element_text(hjust = 0.5),legend.position="bottom")+
  # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab('') + ylab('miRNA response element frequency')+
  scale_fill_manual(values = mycolors[c(2,1)])


library(ggplot2)
library(dplyr)

# https://r-charts.com/part-whole/pie-chart-labels-outside-ggplot2/
# Get the positions
df2 <- mirmotif_type[1:2,1:2] %>% 
  mutate(csum = rev(cumsum(rev(Num))), 
         pos = Num/2 + lead(csum, 1),
         pos = if_else(is.na(pos), Num/2, pos))
df2$per <- round(df2$Num/sum(df2$Num),3)
pie5 <- ggplot(df2, aes(x = "" , y = Num, fill = fct_inorder(Type))) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  geom_label_repel(data = df2,
                   aes(y = pos, label = paste0(per*100, "%")),
                   size = 4.5, nudge_x = 1, show.legend = FALSE) +
  guides(fill = guide_legend(title = "Category")) +
  theme(plot.title = element_text(size=6)) +
  ggtitle("Predicted human 5'UTR \n miRNA response element") +
  theme_void()+
  theme(plot.title = element_text(hjust = 0.5))

df2 <- mirmotif_type[3:4,1:2] %>% 
  mutate(csum = rev(cumsum(rev(Num))), 
         pos = Num/2 + lead(csum, 1),
         pos = if_else(is.na(pos), Num/2, pos))
df2$per <- round(df2$Num/sum(df2$Num),3)
pie6 <- ggplot(df2, aes(x = "" , y = Num, fill = fct_inorder(Type))) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  geom_label_repel(data = df2,
                   aes(y = pos, label = paste0(per*100, "%")),
                   size = 4.5, nudge_x = 1, show.legend = FALSE) +
  guides(fill = guide_legend(title = "Category")) +
  theme(plot.title = element_text(size=6)) +
  ggtitle("Predicted mouse 5'UTR \n miRNA response element") +
  theme_void()+
  theme(plot.title = element_text(hjust = 0.5))

df2 <- mirmotif_type[5:6,1:2] %>% 
  mutate(csum = rev(cumsum(rev(Num))), 
         pos = Num/2 + lead(csum, 1),
         pos = if_else(is.na(pos), Num/2, pos))
df2$per <- round(df2$Num/sum(df2$Num),3)
pie7 <- ggplot(df2, aes(x = "" , y = Num, fill = fct_inorder(Type))) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  geom_label_repel(data = df2,
                   aes(y = pos, label = paste0(per*100, "%")),
                   size = 4.5, nudge_x = 1, show.legend = FALSE) +
  guides(fill = guide_legend(title = "Category")) +
  theme(plot.title = element_text(size=6)) +
  ggtitle("Predicted human 3'UTR \n miRNA response element") +
  theme_void()+
  theme(plot.title = element_text(hjust = 0.5))


# Get the positions
mirmotif_type[8,1] <-  mirmotif_type[8,1] +mirmotif_type[9,1]
df2 <- mirmotif_type[7:8,1:2] %>% 
  mutate(csum = rev(cumsum(rev(Num))), 
         pos = Num/2 + lead(csum, 1),
         pos = if_else(is.na(pos), Num/2, pos))

df2$per <- round(df2$Num/sum(df2$Num),3)

pie8 <- ggplot(df2, aes(x = "" , y = Num, fill = fct_inorder(Type))) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  geom_label_repel(data = df2,
                   aes(y = pos, label = paste0(per*100, "%")),
                   size = 4.5, nudge_x = 1, show.legend = FALSE) +
  guides(fill = guide_legend(title = "Category")) +
  theme(plot.title = element_text(size=6)) +
  ggtitle("Predicted mouse 3'UTR \n miRNA response element") +
  theme_void()+
  theme(plot.title = element_text(hjust = 0.5))

ggarrange(pie5, pie6, pie7, pie8, 
          labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2)



# cancer vs normal comparison
exp_human <- exp[which(exp$organism=='human'),]
cancer_motif_human <- motif_3UTR_cutoff8[which(motif_3UTR_cutoff8$X1 %in% as.vector(exp_human[
  unique(c(which(grepl('cancer',tolower(as.vector(exp_human$study_abstract)), fixed=TRUE)),
           which(grepl('tumor',tolower(as.vector(exp_human$study_abstract)), fixed=TRUE)),
           which(grepl('carcinoma',tolower(as.vector(exp_human$study_abstract)), fixed=TRUE)),
           which(grepl('neoplasm',tolower(as.vector(exp_human$study_abstract)), fixed=TRUE)),
           which(grepl('malignan',tolower(as.vector(exp_human$study_abstract)), fixed=TRUE))
           )),7])),]
# cancer_motif_human <- motif_3UTR_7mer_cutoff8[which(motif_3UTR_7mer_cutoff8$X1 %in% as.vector(exp_human[
#   unique(c(which(grepl('cancer',tolower(as.vector(exp_human$study_abstract)), fixed=TRUE)),
#            which(grepl('tumor',tolower(as.vector(exp_human$study_abstract)), fixed=TRUE)),
#            which(grepl('carcinoma',tolower(as.vector(exp_human$study_abstract)), fixed=TRUE)),
#            which(grepl('neoplasm',tolower(as.vector(exp_human$study_abstract)), fixed=TRUE)),
#            which(grepl('malignan',tolower(as.vector(exp_human$study_abstract)), fixed=TRUE))
#   )),7])),]
exp_mouse <- exp[which(exp$organism=='mouse'),]
cancer_motif_mouse <- motif_3UTR_cutoff8[which(motif_3UTR_cutoff8$X1 %in% as.vector(exp_mouse[
  unique(c(which(grepl('cancer',tolower(as.vector(exp_mouse$study_abstract)), fixed=TRUE)),
           which(grepl('tumor',tolower(as.vector(exp_mouse$study_abstract)), fixed=TRUE)),
           which(grepl('carcinoma',tolower(as.vector(exp_mouse$study_abstract)), fixed=TRUE)),
           which(grepl('neoplasm',tolower(as.vector(exp_mouse$study_abstract)), fixed=TRUE)),
           which(grepl('malignan',tolower(as.vector(exp_mouse$study_abstract)), fixed=TRUE))
  )),7])),]
# cancer_motif_mouse <- motif_3UTR_7mer_cutoff8[which(motif_3UTR_7mer_cutoff8$X1 %in% as.vector(exp_mouse[
#   unique(c(which(grepl('cancer',tolower(as.vector(exp_mouse$study_abstract)), fixed=TRUE)),
#            which(grepl('tumor',tolower(as.vector(exp_mouse$study_abstract)), fixed=TRUE)),
#            which(grepl('carcinoma',tolower(as.vector(exp_mouse$study_abstract)), fixed=TRUE)),
#            which(grepl('neoplasm',tolower(as.vector(exp_mouse$study_abstract)), fixed=TRUE)),
#            which(grepl('malignan',tolower(as.vector(exp_mouse$study_abstract)), fixed=TRUE))
#   )),7])),]
human_mir_cancer <- as.data.frame(table(c(cancer_motif_human$X2,cancer_motif_human$X3,cancer_motif_human$X4,
                                          cancer_motif_human$X5,cancer_motif_human$X6))[
                                            which(names(table(c(cancer_motif_human$X2,cancer_motif_human$X3,
                                                                cancer_motif_human$X4,cancer_motif_human$X5,
                                                                cancer_motif_human$X6))) %in% 
                                                    c(mirlist[mirlist$organism=='Homo sapiens','motif']))])
mouse_mir_cancer <- as.data.frame(table(c(cancer_motif_mouse$X2,cancer_motif_mouse$X3,cancer_motif_mouse$X4,cancer_motif_mouse$X5,cancer_motif_mouse$X6))[
  which(names(table(c(cancer_motif_mouse$X2,cancer_motif_mouse$X3,cancer_motif_mouse$X4,cancer_motif_mouse$X5,cancer_motif_mouse$X6)))
        %in% c(mirlist[mirlist$organism=='Mus musculus','motif']))])

mir_cancer <- rbind(human_mir_cancer,mouse_mir_cancer)
mir_cancer$organism <- c(rep('human',nrow(human_mir_cancer)),
                         rep('mouse',nrow(mouse_mir_cancer)))
mir_cancer <- merge(mir_cancer,
      mirlist[which((mirlist$organism=='Homo sapiens')|(mirlist$organism=='Mus musculus')),c(4,6)],
      by.x = 'Var1',by.y='motif')
mir_cancer$mirorganism <- ifelse(substring(mir_cancer$mir,1,3)=='hsa', 'human', 'mouse')
mir_cancer <- mir_cancer[mir_cancer$organism==mir_cancer$mirorganism,1:4]
mir_cancer <- aggregate(mir ~ Var1+ Freq + organism , mir_cancer, FUN = paste, collapse=",")
mir_cancer <- mir_cancer[order(mir_cancer$Freq,decreasing = TRUE),]
mir_cancer$order <- as.integer(seq(1,nrow(mir_cancer)))
# remove ARE

ggplot(mir_cancer[1:15,], aes(x=order , y=Freq)) +
  geom_line()+
  geom_point()+ ylab('miRNA-related motif occurance')+xlab('Top 15 cancer-related motifs')+
  theme_bw()+  theme(plot.title = element_text(hjust = 0.5),legend.position="bottom")+
  scale_x_continuous(breaks=seq(1,15,1))+
  annotate(geom="text", x=2, y=400, label=mir_cancer[1,1],color=mycolors[1])
