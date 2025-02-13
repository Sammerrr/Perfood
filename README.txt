Scripts for the statistical analysis:

- KW_flux_vs_diab.R:
Kruskal–Wallis test to compare the predicted metabolic fluxes between two groups of participants,
(Diabetic and healthy).

- KW_flux_vs_medication.R:
Kruskal–Wallis test to compare the predicted metabolic fluxes between two groups of diabetic patients,
taking diabetic medication and not taking diabetic medication).

- KW_mic_vs_diab.R:
Kruskal–Wallis test to compare the gut microbiome species abundance between two groups of participants,
(Diabetic and healthy).

- diversity.R:
estimate alpha diversity measures of the gut microbial commuinities

- make_KW_metabolic_groups.R:
To count the number of metabolies per metabolic group that are segnificantly increased or decreased in diabetic patients. 

###################
Result datafiles:

- Pvalues.zip:
A collection of statistical analysis results containing:

-- Abundances_vs_diabetes.tsv:
The statistical analysis results of the comparasion in speices abundance between diabetic patients and healthy individuals.

-- Ecological_relation_ratio_vs_diabetes.tsv:
The statistical analysis results of the comparasion in ecological relation ratio (log10 transformed) between diabetic patients and healthy individuals.

-- Ecological_relation_ratio_vs_glucose.tsv
The statistical analysis results of the association between ecological relation ratio (log10 transformed) and blood glucose levels.

-- enrichment.tsv
The enrichment analysis results results of pathway enrichment in reactions that are associated with diabetes.

-- Fluxes_vs_diabetes.tsv
The statistical analysis results of the comparasion in predicted metabolic flux rates (log10 transformed) between diabetic patients and healthy individuals.

-- Fluxes_vs_glucouse.tsv
The statistical analysis results of the association between the predicted metabolic flux rates and blood glucose levels.

-- reacts_vs_gluc.tsv
The statistical analysis results of the association between the predicted internal reaction rates and blood glucose levels.

- community_growth.csv:
The predicted growth of the gut bacterial communities using flux balance analysis.

- eco_mat.csv:
A matrix of the predicted 6 types of ecological interactions between all possible pairs of bacterial species found in the gut microbial communities.

- flux_among_bacteria.csv:
The predicted metabolic fluxe exchange within the microbial communities (among the member species).

- flux_with_host.csv:
The predicted metabolic fluxe exchange between gut microbial communities and the human host.

- internal_reacts1.csv, internal_reacts2.csv, internal_reacts3.csv, internal_reacts4.csv:
The predicted internal reactions within the different bacterial species summed up for each of the microbial communities.

- relation_ratios_log.csv:
the ratio between ecological interaction frequencies for each of the microbioal communities (log10 transformed)