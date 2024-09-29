rm(list=ls())
library(stringr)
pheno_flux <- read.table("/home/samer/projects/perfood/diet_modelling/sim_data/flux_pheno.csv", header = T)
pheno_flux <- pheno_flux[!pheno_flux$remove,]
pheno_flux <- pheno_flux[pheno_flux$ailment_diabetes_mellitus=="True",]
ind <- grep("EX_cpd", names(pheno_flux))
flux_to_corr <- pheno_flux[,ind]
names(flux_to_corr) <- str_replace(names(flux_to_corr), "-c0","")
names(flux_to_corr) <- str_replace(names(flux_to_corr), "_e0","")
names(flux_to_corr) <- str_replace(names(flux_to_corr), "EX_","")
flux_to_corr  <- Filter(function(x) sd(x) != 0, flux_to_corr)
pvalue <- as.data.frame(t(matrix(NA,5,ncol(flux_to_corr))))
colnames(pvalue) = c("id", "P.value", "chi_squared", "num_samps", "direction")
pvalue$id <- colnames(flux_to_corr)
ind <- pheno_flux$med_anti_diabetics=="True"
for (i in 1:ncol(flux_to_corr)){
  current_compar <- kruskal.test(flux_to_corr[,i] ~ pheno_flux$med_anti_diabetics)
  pvalue$P.value[i] <- current_compar[["p.value"]]
  pvalue$chi_squared[i] <- current_compar[["statistic"]][["Kruskal-Wallis chi-squared"]]
  pvalue$num_samps[i] <- nrow(flux_to_corr)-sum(is.na(flux_to_corr[,i])|flux_to_corr[,i]==0)
  median_diab <- median(flux_to_corr[ind,i])
  median_normal <- median(flux_to_corr[-ind,i])
  pvalue$direction[i] <- sign(median_diab - median_normal)
}
pvalue$fdr_p.value <- p.adjust(pvalue$P.value, method = "fdr", n = nrow(pvalue))
pvalue <-pvalue[order(pvalue$fdr_p.value),]
segn_pvalue <- pvalue[pvalue$fdr_p.value<=0.05,]
metabolite_names <- read.csv("/home/samer/projects/gapseq_models/all_metabolites_14.11.2022.tsv", sep = "\t",
                             row.names = 1, header = T)
segn_pvalue <-segn_pvalue[order(segn_pvalue$fdr_p.value),]
write.table(segn_pvalue, "/home/samer/projects/perfood/diet_modelling/results/KW_flux_vs_medication.tsv",
            sep = "\t", row.names = F, col.names = T)