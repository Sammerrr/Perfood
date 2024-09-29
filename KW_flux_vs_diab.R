rm(list=ls())
library(pgirmess)
library(stringr)
pheno_flux <- read.table("/home/samer/projects/perfood/diet_modelling/sim_data/flux_pheno.csv", header = T)
pheno_flux <- pheno_flux[!pheno_flux$remove,]
pheno_flux <- pheno_flux[!duplicated(pheno_flux$patient_id),]

pheno_flux$ailment_diabetes_mellitus <- str_replace(pheno_flux$ailment_diabetes_mellitus, "True", "Diabetic")
pheno_flux$ailment_diabetes_mellitus <- str_replace(pheno_flux$ailment_diabetes_mellitus, "PreDiab", "Healthy")
pheno_flux$ailment_diabetes_mellitus <- str_replace(pheno_flux$ailment_diabetes_mellitus, "False", "Healthy")
sum(pheno_flux$fast_glucose>125&pheno_flux$ailment_diabetes_mellitus!="Diabetic")
pheno_flux$ailment_diabetes_mellitus[pheno_flux$fast_glucose>125]
ind <- grep("EX_cpd", names(pheno_flux))
flux_to_corr <- pheno_flux[,ind]
flux_to_corr <- flux_to_corr/pheno_flux$growth
names(flux_to_corr) <- str_replace(names(flux_to_corr), "-c0","")
names(flux_to_corr) <- str_replace(names(flux_to_corr), "_e0","")
names(flux_to_corr) <- str_replace(names(flux_to_corr), "EX_","")
flux_to_corr  <- Filter(function(x) sd(x) != 0, flux_to_corr)
num_in_diab <- sapply(flux_to_corr, function(x)sum(x!=0&pheno_flux$ailment_diabetes_mellitus=="Diabetic"))
ratio_in_diab <- num_in_diab/sum(pheno_flux$ailment_diabetes_mellitus=="Diabetic")
flux_to_corr <- flux_to_corr[,num_in_diab>0.5]
prod <- sapply(flux_to_corr, function(x)sum(x>0))
cons <- sapply(flux_to_corr, function(x)sum(x<0))
pvalue <- as.data.frame(t(matrix(NA,6,ncol(flux_to_corr))))
colnames(pvalue) = c("id", "P.value", "chi_squared", "num_samps", "Diabetic-Healthy", "log2 Folds")
pvalue$id <- colnames(flux_to_corr)
healthy_ind <- pheno_flux$ailment_diabetes_mellitus=="Healthy"
diab_ind <- pheno_flux$ailment_diabetes_mellitus=="Diabetic"

i=1
for (i in 1:ncol(flux_to_corr)){
  current_compar <- kruskal.test(abs(flux_to_corr[,i]) ~ pheno_flux$ailment_diabetes_mellitus)
  pvalue$P.value[i] <- current_compar[["p.value"]]
  pvalue$chi_squared[i] <- current_compar[["statistic"]][["Kruskal-Wallis chi-squared"]]
  pvalue$num_samps[i] <- nrow(flux_to_corr)-sum(is.na(flux_to_corr[,i])|flux_to_corr[,i]==0)
  median_diab <- abs(median(flux_to_corr[diab_ind,i]))
  median_healthy <- abs(median(flux_to_corr[healthy_ind,i]))
  pvalue[i, "Diabetic-Healthy"] = sign(median_diab - median_healthy)
  if(median_diab*median_healthy!=0){
    pvalue[i, "log2 Folds"] = log2(median_diab/median_healthy)
  }
  else{
    pvalue[i, "log2 Folds"] = max(c(median_diab, median_healthy))
  }
}
pvalue$fdr_p.value <- p.adjust(pvalue$P.value, method = "fdr", n = nrow(pvalue))
segn_pvalue <- pvalue[pvalue$fdr_p.value<=0.05,]
segn_pvalue <- segn_pvalue[segn_pvalue$num_samps>nrow(pheno_flux)/2,]
pvalue <-pvalue[order(pvalue$fdr_p.value),]

metabolite_names <- read.csv("/home/samer/projects/gapseq_models/all_metabolites_14.11.2022.tsv", sep = "\t",
                             row.names = 1, header = T)
segn_pvalue <- merge(segn_pvalue, metabolite_names, by.x = "id", by.y = "ids"  )
segn_pvalue <-segn_pvalue[order(segn_pvalue$fdr_p.value),]

segn_pvalue[,"consumed"] <- cons[segn_pvalue$id]
segn_pvalue[,"produced"] <- prod[segn_pvalue$id]

write.table(segn_pvalue, "/home/samer/projects/perfood/diet_modelling/results/KW/KW_flux_vs_diab.tsv",
            sep = "\t", row.names = F, col.names = T)
