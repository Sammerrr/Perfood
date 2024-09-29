rm(list=ls())
library(pgirmess)
library(stringr)
pheno_flux <- read.table("/home/samer/projects/perfood/diet_modelling/sim_data/flux_pheno.csv", header = T)
pheno_flux <- pheno_flux[!pheno_flux$remove,]
pheno_flux <- pheno_flux[!duplicated(pheno_flux$patient_id),]
rownames(pheno_flux) <- paste0("X", pheno_flux$sample.id)
mic <- read.csv("/home/samer/projects/perfood/diet_modelling/preliminary_data/AGORA_mic.table", sep = "\t")
mic <- as.data.frame(t(mic))

common <- intersect(rownames(pheno_flux), rownames(mic))
pheno_flux <- pheno_flux[common,]
mic <- mic[common,]

pheno_flux$ailment_diabetes_mellitus <- str_replace(pheno_flux$ailment_diabetes_mellitus, "True", "Diabetic")
pheno_flux$ailment_diabetes_mellitus <- str_replace(pheno_flux$ailment_diabetes_mellitus, "PreDiab", "Healthy")
pheno_flux$ailment_diabetes_mellitus <- str_replace(pheno_flux$ailment_diabetes_mellitus, "False", "Healthy")
sum(pheno_flux$fast_glucose>125&pheno_flux$ailment_diabetes_mellitus!="Diabetic")
pheno_flux$ailment_diabetes_mellitus[pheno_flux$fast_glucose>125]
mic_to_corr  <- Filter(function(x) sd(x) != 0, mic)

num_in_diab <- sapply(mic_to_corr, function(x)sum(x!=0&pheno_flux$ailment_diabetes_mellitus=="Diabetic"))
ratio_in_diab <- num_in_diab/sum(pheno_flux$ailment_diabetes_mellitus=="Diabetic")
mic_to_corr <- mic_to_corr[,num_in_diab>0.5]
prod <- sapply(mic_to_corr, function(x)sum(x>0))
cons <- sapply(mic_to_corr, function(x)sum(x<0))
pvalue <- as.data.frame(t(matrix(NA,6,ncol(mic_to_corr))))
colnames(pvalue) = c("id", "P.value", "chi_squared", "num_samps", "Diabetic-Healthy", "log2 Folds")
pvalue$id <- colnames(mic_to_corr)
healthy_ind <- pheno_flux$ailment_diabetes_mellitus=="Healthy"
diab_ind <- pheno_flux$ailment_diabetes_mellitus=="Diabetic"

i=10
for (i in 1:ncol(mic_to_corr)){
  current_compar <- kruskal.test(abs(mic_to_corr[,i]) ~ pheno_flux$ailment_diabetes_mellitus)
  pvalue$P.value[i] <- current_compar[["p.value"]]
  pvalue$chi_squared[i] <- current_compar[["statistic"]][["Kruskal-Wallis chi-squared"]]
  pvalue$num_samps[i] <- nrow(mic_to_corr)-sum(is.na(mic_to_corr[,i])|mic_to_corr[,i]==0)
  median_diab <- abs(median(mic_to_corr[diab_ind,i]))
  median_healthy <- abs(median(mic_to_corr[healthy_ind,i]))
  pvalue[i, "Diabetic-Healthy"] = sign(median_diab - median_healthy)
  if(median_diab*median_healthy!=0){
    pvalue[i, "log2 Folds"] = log2(median_diab/median_healthy)
  }
  else{
    pvalue[i, "log2 Folds"] = max(c(median_diab, median_healthy))
  }
}
pvalue$fdr_p.value <- p.adjust(pvalue$P.value, method = "fdr", n = nrow(pvalue))
segn_pvalue <- pvalue[pvalue$fdr_p.value<=0.05&pvalue$`Diabetic-Healthy`!=0,]
segn_pvalue <- segn_pvalue[segn_pvalue$num_samps>nrow(pheno_flux)/4,]
pvalue <-pvalue[order(pvalue$fdr_p.value),]
segn_pvalue <-segn_pvalue[order(segn_pvalue$fdr_p.value),]


write.table(segn_pvalue, "/home/samer/projects/perfood/diet_modelling/results/KW/KW_mic_vs_diab.tsv",
            sep = "\t", row.names = F, col.names = T)