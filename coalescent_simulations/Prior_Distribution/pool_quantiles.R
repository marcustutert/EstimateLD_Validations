#Pool quantiles together and plot resulting histograms

#Input the snakemake parameters
divergence_time = snakemake@params$divergence_time

#Grab the results across the replicates
LD_Results_File      = list.files(path = "results/pop_split", pattern = sprintf("split_%s_LD_quantile", divergence_time), full.names = T)
print(LD_Results_File)
AF_Results_File      = list.files(path = "results/pop_split", pattern = sprintf("split_%s_AF_quantile", divergence_time), full.names = T)

LD_Pooled = lapply(LD_Results_File,readRDS,.GlobalEnv)
AF_Pooled = lapply(AF_Results_File,readRDS,.GlobalEnv)

#Write out the pooled quantiles
saveRDS(LD_Pooled,sprintf("results/pop_split/split_%s_LD_pooled_quantile_counts.RData", divergence_time), version = 2)
saveRDS(AF_Pooled,sprintf("results/pop_split/split_%s_AF_pooled_quantile_counts.RData", divergence_time), version = 2)



