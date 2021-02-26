#Pool quantiles together and plot resulting histograms for OOA model
#Input the snakemake parameters
population_pairs = snakemake@params$paired_values
print(population_pairs)
#Grab the results across the replicates and the different GWAS pairings
LD_Results_File      = list.files(path = "results/UKBB", pattern = sprintf("%s_LD_quantile*", population_pairs), full.names = T)
AF_Results_File      = list.files(path = "results/UKBB", pattern = sprintf("%s_AF_quantile*", population_pairs), full.names = T)


LD_Pooled = lapply(LD_Results_File,readRDS,.GlobalEnv)
AF_Pooled = lapply(AF_Results_File,readRDS,.GlobalEnv)

#Write out the pooled quantiles
saveRDS(LD_Pooled,sprintf("results/UKBB/%s_LD_pooled_quantile_counts.RData", population_pairs), version = 2)
saveRDS(AF_Pooled,sprintf("results/UKBB/%s_AF_pooled_quantile_counts.RData", population_pairs), version = 2)
