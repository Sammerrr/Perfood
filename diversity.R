rm(list=ls())
setwd("/home/samer/projects/perfood/diet_modelling/Rscripts/")
library(vegan)
library(stringr)
mic.table <- read.table("/home/samer/projects/perfood/diet_modelling/preliminary_data/AGORA_mic.table",
                        sep = "\t", header = T, row.names = 1)
mic.table <- mic.table[rowSums(mic.table)>0,]
diversity <- as.data.frame(diversity(mic.table, MARGIN = 2, index = "shannon"))
names(diversity) <- "shanon"
chao1 <- as.data.frame(t(estimateR(t(mic.table))))
diversity$chao1 <- chao1$S.chao1
diversity$richness <- specnumber(t(mic.table))
diversity$simpson <- diversity(t(mic.table), index = "simpson")


diversity$sample.id <- str_remove(rownames(diversity), "X")
write.table(diversity, "/home/samer/projects/perfood/diet_modelling/preliminary_data/diversity.tsv",
            sep = "\t", row.names = T, col.names = T)
