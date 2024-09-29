rm(list=ls())  
groups <- read.csv("/home/samer/projects/perfood/diet_modelling/results/KW/groups.tsv", sep = "\t")
rownames(groups) <- groups$id
flux_diab <- read.csv("/home/samer/projects/perfood/diet_modelling/results/KW/KW_flux_vs_diab.tsv", sep = "\t")
rownames(flux_diab) <- flux_diab$id
flux_diab$Group <- groups[rownames(flux_diab),"Group"]

host_to_mic_ind <- grepl("host", flux_diab$metabolites)&flux_diab$consumed>flux_diab$produced
mic_to_host_ind <- grepl("host", flux_diab$metabolites)&flux_diab$consumed<flux_diab$produced
mic_ind <- !grepl("host", flux_diab$metabolites)

inds <- as.data.frame(cbind (host_to_mic_ind&flux_diab$Diabetic.Healthy == 1,
                             host_to_mic_ind&flux_diab$Diabetic.Healthy == -1,
                             mic_to_host_ind&flux_diab$Diabetic.Healthy == 1,
                             mic_to_host_ind&flux_diab$Diabetic.Healthy == -1,
                             mic_ind&flux_diab$Diabetic.Healthy == 1,
                             mic_ind&flux_diab$Diabetic.Healthy == -1
                             ))
  names(inds) <- c("host_to_mic_d", "host_to_mic_h", "mic_to_host_d", "mic_to_host_h",
                   "mic_d", "nic_h")
  
  i=1
  for (i in 1: ncol(inds)){
    groups <- as.data.frame(table(flux_diab$Group[inds[,i]]))
    groups <- groups[order(groups$Freq, decreasing = T),]
    if (length(groups)!=0){
      names(groups) <- c("Category", "Number_of_metabolites")
      all_sum <- sum(groups$Number_of_metabolites)
      groups <- rbind(as.matrix(groups), c("Sum", all_sum))
      write.csv(groups, paste0("/home/samer/projects/perfood/diet_modelling/results/KW/met_groups/",
                               names(inds)[i], ".csv"), row.names = F)
    }
  }
  
  
