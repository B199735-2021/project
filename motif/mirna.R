setwd('/home/s2172876/mirna')
library('ampir')
mirlist <- read_faa(file = 'mature.fa.gz')
mirlist$seed <- substring(mirlist[,2],2,7)
library(spgs)
library(stringr)
library(data.table)
library(readr)
# reverse complement
mirlist$motif <- gsub('U','T',mirlist$seed)
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


# # plot
# library(ggplot2)
# mir_motif <- rbind(human_mir_3UTR[order(human_mir_3UTR$Freq,decreasing=TRUE),][1:5,],
#                    mouse_mir_3UTR[order(mouse_mir_3UTR$Freq,decreasing=TRUE),][1:5,],
#                    human_mir_5UTR[order(human_mir_5UTR$Freq,decreasing=TRUE),][1:5,],
#                    mouse_mir_5UTR[order(mouse_mir_5UTR$Freq,decreasing=TRUE),][1:5,])
# mir_motif$Type <- c(rep('human 3\'UTR',5),
#                     rep('mouse 3\'UTR',5),
#                     rep('human 5\'UTR',5),
#                     rep('mouse 5\'UTR',5))
# mir_motif[,1] <- as.character(mir_motif[,1])
# mir_motif[c(1:5,11:15),1] <- paste(mir_motif[c(1:5,11:15),1],'(hsa)',sep='')
# mir_motif[c(6:10,16:20),1] <- paste(mir_motif[c(6:10,16:20),1],'(mmu)',sep='')
# library(ggplot2)
# library(ggbreak)
# library(wesanderson)
# positions <- mir_motif[c(1,6,11,16,2,7,12,17,  
#                          3,8,13,18,4,9,14,19,
#                          5,10,15,20),1]
# mycolors <- c("#0F6444", "#3182BD", "#8576AB", "#2c966d")
# ggplot(mir_motif, aes(fill=Type, y=Freq, x= Var1)) + 
#     # geom_bar(position="dodge", stat="identity") +
#     geom_col(width = 0.6, position = "dodge")+
#     scale_y_cut(breaks=c(80,300,4000)) +
#     theme_bw()+ ylim(0, 7000)+
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
#     scale_x_discrete(limits = positions) + xlab('Motif') +
#     scale_fill_manual(values = mycolors)
# 
# # piechart
# (sum(human_mir_3UTR$Freq) + sum(mouse_mir_3UTR$Freq)) / sum(table(motif_3UTR_summary$motif[which(motif_3UTR_summary$motif!='')]))
# (sum(human_mir_5UTR$Freq) + sum(mouse_mir_5UTR$Freq)) / sum(table(motif_5UTR_summary$motif[which(motif_5UTR_summary$motif!='')]))
# 
# ggplot(data, aes(x="", y=prop, fill=group)) +
#     geom_bar(stat="identity", width=1, color="white") +
#     coord_polar("y", start=0) +
#     theme_void() + 
#     theme(legend.position="none") +
#     geom_text(aes(y = ypos, label = group), color = "white", size=6) +
#     scale_fill_brewer(palette="Set1")



# 7mer
mirlist$seed <- substring(mirlist[,2],2,8)
library(spgs)
library(stringr)
library(data.table)
library(readr)
# reverse complement
mirlist$motif <- gsub('U','T',mirlist$seed)
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
p1 <- ggplot(motif_summary_compare_stack[c(1:4,9:12),], aes(fill=Position, y=values, x= Type)) +
    geom_bar(position="dodge", stat="identity") +
    geom_col(width = 0.6, position = "dodge")+
    scale_y_cut(breaks=c(105,1000)) +
    geom_text(aes(label=values),vjust=1,color = "white",size=4,position = position_dodge(width = .95),fontface = "bold")+
    theme_bw()+  theme(plot.title = element_text(hjust = 0.5),legend.position="bottom")+
    # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    xlab('') + ylab('The number of unique miRNA-related motif')+
    scale_fill_manual(values = mycolors)#+labs(tag = "A")

p2 <- ggplot(motif_summary_compare_stack[c(5:8,13:16),], aes(fill=Position, y=values, x= Type)) +
    geom_bar(position="dodge", stat="identity") +
    geom_col(width = 0.6, position = "dodge")+
    scale_y_cut(breaks=c(6000,20000)) +
    geom_text(aes(label=values),vjust=1,color = "white",size=4,position = position_dodge(width = .95),fontface = "bold")+
    theme_bw()+  theme(plot.title = element_text(hjust = 0.5),legend.position="bottom")+
    # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    xlab('') + ylab('miRNA-related motif occurance')+
    scale_fill_manual(values = mycolors)#+labs(tag = "B")

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
p3 <- ggplot(mir_motif[c(1:15,31:45),], aes(x=order , y=Freq, group=Type)) +
    geom_line(aes(color=Type))+
    geom_point(aes(color=Type))+ylab('miRNA-related motif occurance')+xlab('Top 15 motifs')+
    theme_bw()+  theme(plot.title = element_text(hjust = 0.5),legend.position="bottom")+
    scale_x_continuous(breaks=seq(1,15,1))+
    geom_hline(yintercept=11, linetype="dashed", color = "red")+
    annotate(geom="text", x=15, y=30, label="y=11",color="red")+
    annotate(geom="text", x=2.5, y=370, label=mir_motif[1,1],color=mycolors[1])+ # CCCAGG
    annotate(geom="text", x=2.5, y=330, label=mir_motif[31,1],color=mycolors[3])+# CTTCTT
    scale_color_manual(values=mycolors[c(1,3)])#+labs(tag = "C")


p4 <- ggplot(mir_motif[c(16:30,46:60),], aes(x=order , y=Freq, group=Type)) +
    geom_line(aes(color=Type))+
    geom_point(aes(color=Type))+ylab('miRNA-related motif occurance')+xlab('Top 15 motifs')+
    theme_bw()+  theme(plot.title = element_text(hjust = 0.5),legend.position="bottom")+
    scale_x_continuous(breaks=seq(1,15,1))+
    geom_hline(yintercept=3, linetype="dashed", color = "red")+  
    annotate(geom="text", x=15, y=300, label="y=3",color="red")+
    annotate(geom="text", x=2.6, y=6800, label=mir_motif[46,1],color=mycolors[3])+ # CCCAGG
    annotate(geom="text", x=2.7, y=2400, label=mir_motif[16,1],color=mycolors[1])+# CTTCTT
    scale_color_manual(values=mycolors[c(1,3)])#+labs(tag = "D")
library(ggpubr)
ggarrange(p1, p2, p3, p4, 
          labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2)

library(patchwork)
ps1 <- p1+p2
ps2 <- p3+p4
ps1/ps2
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
mouse_mir_5UTR[grepl("TAAAT", mouse_mir_5UTR$Var1,TRUE),3] <- 'ARE' 
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
  ggtitle("Predicted human 5'UTR \npotential miRNA-related motifs") +
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
  ggtitle("Predicted mouse 5'UTR \npotential miRNA-related motifs") +
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
  ggtitle("Predicted human 3'UTR \npotential miRNA-related motifs") +
  theme_void()+
  theme(plot.title = element_text(hjust = 0.5))


# Get the positions
df2 <- mirmotif_type[7:9,1:2] %>% 
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
  ggtitle("Predicted mouse 3'UTR \npotential miRNA-related motifs") +
  theme_void()+
  theme(plot.title = element_text(hjust = 0.5))

ggarrange(pie1, pie2, pie3, pie4, 
          labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2)

# dev.off()

# 
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
p5 <- ggplot(p5_df, aes(fill=Region, y=Num, x= Region)) +
  geom_bar(position="dodge", stat="identity") +
  geom_col(width = 0.6, position = "dodge")+
  scale_y_cut(breaks=c(15,100)) +
  geom_text(aes(label=Num),vjust=1,color = "white",size=4,position = position_dodge(width = .95),fontface = "bold")+
  theme_bw()+  theme(plot.title = element_text(hjust = 0.5),legend.position="bottom")+
  # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab('') + ylab('The number of unique miRNA-related motif')+
  scale_fill_manual(values = mycolors[c(2,1)]) +labs(tag = "A")
p6_df <- data.frame(Num = c(sum(mirmotif_type[which(mirmotif_type$Type=='Other'),c(1,3)][1:2,1]),
                            sum(mirmotif_type[which(mirmotif_type$Type=='Other'),c(1,3)][3:4,1])),
                    Region = c('5\'UTR','3\'UTR'),Type=rep('miRNA-related motif',2))

p6 <- ggplot(p6_df, aes(fill=Region, y=Num, x= Region)) +
  geom_bar(position="dodge", stat="identity") +
  geom_col(width = 0.6, position = "dodge")+
  scale_y_cut(breaks=c(30,100)) +  
  geom_text(aes(label=Num),vjust=1,color = "white",size=4,position = position_dodge(width = .95),fontface = "bold")+
  theme_bw()+  theme(plot.title = element_text(hjust = 0.5),legend.position="bottom")+
  # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab('') + ylab('miRNA-related motif occurance')+
  scale_fill_manual(values = mycolors[c(2,1)])+labs(tag = "B")
library(patchwork)
p5+p6

human_nonmir_3UTR
mouse_nonmir_3UTR
human_nonmir_5UTR
mouse_nonmir_5UTR



exp_summary <- read.csv('../DEG_summary_new_fix',sep = '\t')[,c(1,2,3)]
project <- read.csv('../annotation/recount3_selection_2022-05-11 21_54_01.csv')
sample <- read.csv('../annotation/sample_annotation.csv',sep='\t',stringsAsFactors=FALSE)
exp <- merge(project,exp_summary,by.x = 'project',by.y='ProjectID')

exp_CCCAGG <- exp[which(exp$ExpID %in% c(motif_3UTR_cutoff8[which(motif_3UTR_cutoff8[,2]=='CCCAGG'),1],
                                         motif_3UTR_cutoff8[which(motif_3UTR_cutoff8[,3]=='CCCAGG'),1],
                                         motif_3UTR_cutoff8[which(motif_3UTR_cutoff8[,4]=='CCCAGG'),1],
                                         motif_3UTR_cutoff8[which(motif_3UTR_cutoff8[,5]=='CCCAGG'),1],
                                         motif_3UTR_cutoff8[which(motif_3UTR_cutoff8[,6]=='CCCAGG'),1])),]
exp_non_CCCAGG <- exp[-which(exp$ExpID %in% c(motif_3UTR_cutoff8[which(motif_3UTR_cutoff8[,2]=='CCCAGG'),1],
                                             motif_3UTR_cutoff8[which(motif_3UTR_cutoff8[,3]=='CCCAGG'),1],
                                             motif_3UTR_cutoff8[which(motif_3UTR_cutoff8[,4]=='CCCAGG'),1],
                                             motif_3UTR_cutoff8[which(motif_3UTR_cutoff8[,5]=='CCCAGG'),1],
                                             motif_3UTR_cutoff8[which(motif_3UTR_cutoff8[,6]=='CCCAGG'),1])),]
# sample <- sample[which(sample$study %in% exp_CCGGAA$project),]
library("tm")
library("SnowballC")
library("wordcloud")
library("RColorBrewer")
kw <- function(text_input){
  text_input <- tolower(text_input)
  text_input <- str_replace_all(text_input,'(cells|wild|type|/ml|-/-|untreated|wildetype|derived|line|dmso|cell|hours|none|knockout|vehicle)','')
  docs <- Corpus(VectorSource(text_input))
  # inspect(docs)
  
  toSpace <- content_transformer(function (x , pattern ) gsub(pattern, " ", x))
  docs <- tm_map(docs, toSpace, "@")
  docs <- tm_map(docs, toSpace, "\\|")
  docs <- tm_map(docs, toSpace, "\\(")
  docs <- tm_map(docs, toSpace, "\\)")
  docs <- tm_map(docs, toSpace, ";")
  docs <- tm_map(docs, toSpace, ",")
  
  
  # Convert the text to lower case
  docs <- tm_map(docs, content_transformer(tolower))
  docs <- tm_map(docs, toSpace, "hrs")
  docs <- tm_map(docs, toSpace, "days")
  docs <- tm_map(docs, toSpace, "um")
  docs <- tm_map(docs, toSpace, "ng")
  docs <- tm_map(docs, toSpace, "wild")
  docs <- tm_map(docs, toSpace, "tissue")
  # Remove numbers
  # docs <- tm_map(docs, removeNumbers)
  # Remove english common stopwords
  docs <- tm_map(docs, removeWords, stopwords("english"))
  # Remove your own stop word
  # specify your stopwords as a character vector
  docs <- tm_map(docs, removeWords, c("blabla1", "blabla2"))
  # Remove punctuations
  # docs <- tm_map(docs, removePunctuation)
  # Eliminate extra white spaces
  # docs <- tm_map(docs, stripWhitespace)
  # Text stemming
  # docs <- tm_map(docs, stemDocument)
  
  dtm <- TermDocumentMatrix(docs)
  m <- as.matrix(dtm)
  v <- sort(rowSums(m),decreasing=TRUE)
  d <- data.frame(word = names(v),freq=v)
  return(d)
}


# set.seed(1234)
# wordcloud(words = d1$word, freq = d1$freq, min.freq = 1,
#           max.words=200, random.order=FALSE, rot.per=0.35, 
#           colors=brewer.pal(8, "Dark2"))

# CCGGAA
exp_CCGGAA <- exp[which(exp$ExpID %in% c(motif_5UTR_cutoff8[which(motif_5UTR_cutoff8[,2]=='CCGGAA'),1],
                                         motif_5UTR_cutoff8[which(motif_5UTR_cutoff8[,3]=='CCGGAA'),1],
                                         motif_5UTR_cutoff8[which(motif_5UTR_cutoff8[,4]=='CCGGAA'),1],
                                         motif_5UTR_cutoff8[which(motif_5UTR_cutoff8[,5]=='CCGGAA'),1],
                                         motif_5UTR_cutoff8[which(motif_5UTR_cutoff8[,6]=='CCGGAA'),1])),]
exp_non_CCGGAA <- exp[-which(exp$ExpID %in% c(motif_5UTR_cutoff8[which(motif_5UTR_cutoff8[,2]=='CCGGAA'),1],
                                              motif_5UTR_cutoff8[which(motif_5UTR_cutoff8[,3]=='CCGGAA'),1],
                                              motif_5UTR_cutoff8[which(motif_5UTR_cutoff8[,4]=='CCGGAA'),1],
                                              motif_5UTR_cutoff8[which(motif_5UTR_cutoff8[,5]=='CCGGAA'),1],
                                              motif_5UTR_cutoff8[which(motif_5UTR_cutoff8[,6]=='CCGGAA'),1])),]
text_CCGGAA <- as.vector(exp_CCGGAA[,8])
text_CCGGAA <- gsub('_',' ',text_CCGGAA)
text_non_CCGGAA <- as.vector(exp_non_CCGGAA[,8])
text_non_CCGGAA <- gsub('_',' ',text_non_CCGGAA)
d5 <- kw(text_CCGGAA)
d6 <- kw(text_non_CCGGAA)

pvalue<-vector()
for (i in d5$word) {
  mat <- matrix(c(d5[which(d5$word==i),2],d6[which(d6$word==i),2],
                  sum(d5[which(d5$word!=i),2]),sum(d6[which(d6$word!=i),2])),nrow = 2)
  pvalue <- c(pvalue ,fisher.test(mat, alternative = "greater")$p.value)
}
d5$pvalue <- pvalue

d5_plot <- d5[which(d5$pvalue<0.01),]
d5_plot <- d5_plot[1:10,]
d5_plot$word <- as.vector(d5_plot$word)
d5_plot[order(d5_plot$freq,decreasing = TRUE),1]
# plot: dot plot
CCGGAA_dot <- ggplot(data = d5_plot, aes(x = freq, y = word, 
                           color =pvalue, size = freq)) + 
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +
  theme_bw() + 
  ylab("Word") + scale_y_discrete(limits=d5_plot[order(d5_plot$freq,decreasing = FALSE),1])+
  xlab("Frequency") +ggtitle('miRNA-associated motif \n CCGGAA key word')+
  theme(plot.title = element_text(hjust = 0.5,size = 13))


# CGGAAG
exp_CGGAAG <- exp[which(exp$ExpID %in% c(motif_5UTR_cutoff8[which(motif_5UTR_cutoff8[,2]=='CGGAAG'),1],
                                         motif_5UTR_cutoff8[which(motif_5UTR_cutoff8[,3]=='CGGAAG'),1],
                                         motif_5UTR_cutoff8[which(motif_5UTR_cutoff8[,4]=='CGGAAG'),1],
                                         motif_5UTR_cutoff8[which(motif_5UTR_cutoff8[,5]=='CGGAAG'),1],
                                         motif_5UTR_cutoff8[which(motif_5UTR_cutoff8[,6]=='CGGAAG'),1])),]
exp_non_CGGAAG <- exp[-which(exp$ExpID %in% c(motif_5UTR_cutoff8[which(motif_5UTR_cutoff8[,2]=='CGGAAG'),1],
                                              motif_5UTR_cutoff8[which(motif_5UTR_cutoff8[,3]=='CGGAAG'),1],
                                              motif_5UTR_cutoff8[which(motif_5UTR_cutoff8[,4]=='CGGAAG'),1],
                                              motif_5UTR_cutoff8[which(motif_5UTR_cutoff8[,5]=='CGGAAG'),1],
                                              motif_5UTR_cutoff8[which(motif_5UTR_cutoff8[,6]=='CGGAAG'),1])),]
text_CGGAAG <- as.vector(exp_CGGAAG[,8])
text_CGGAAG <- gsub('_',' ',text_CGGAAG)
text_non_CGGAAG <- as.vector(exp_non_CGGAAG[,8])
text_non_CGGAAG <- gsub('_',' ',text_non_CGGAAG)
d7 <- kw(text_CGGAAG)
d8 <- kw(text_non_CGGAAG)

pvalue<-vector()
for (i in d7$word) {
  mat <- matrix(c(d7[which(d7$word==i),2],d8[which(d8$word==i),2],
                  sum(d7[which(d7$word!=i),2]),sum(d8[which(d8$word!=i),2])),nrow = 2)
  pvalue <- c(pvalue, fisher.test(mat, alternative = "greater")$p.value)

}
d7$pvalue <- pvalue

d7_plot <- d7[which(d7$pvalue<0.05),]
d7_plot <- d7_plot[1:10,]
d7_plot$word <- as.vector(d7_plot$word)
d7_plot[order(d7_plot$freq,decreasing = TRUE),1]
# plot: dot plot
CGGAAG_dot <- ggplot(data = d7_plot, aes(x = freq, y = word, 
                                         color =pvalue, size = freq)) + 
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +
  theme_bw() + 
  ylab("Word") + scale_y_discrete(limits=d7_plot[order(d7_plot$freq,decreasing = FALSE),1])+
  xlab("Frequency") +ggtitle('miRNA-associated motif \n CGGAAG key word')+
  theme(plot.title = element_text(hjust = 0.5,size = 13))
ggarrange(CCGGAA_dot, CGGAAG_dot, 
          labels = c("A", "B"),
          ncol = 2, nrow = 1)
