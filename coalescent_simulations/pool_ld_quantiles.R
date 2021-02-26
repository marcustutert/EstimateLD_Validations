#Pool quantiles together and plot resulting histograms

#Input the snakemake parameters
divergence_time = snakemake@params$divergence_time

#Grab the results across the replicates
LD_Results_File      = list.files(path = "results", pattern = sprintf("split_%s_LD", divergence_time), full.names = T)
AF_Results_File      = list.files(path = "results", pattern = sprintf("split_%s_AF", divergence_time), full.names = T)

LD_Pooled = lapply(LD_Results_File,readRDS,.GlobalEnv)
AF_Pooled = lapply(AF_Results_File,readRDS,.GlobalEnv)

#Write out the pooled quantiles
saveRDS(LD_Pooled,sprintf("results/split_%s_LD_pooled_quantile_counts.RData", divergence_time), version = 2)
saveRDS(AF_Pooled,sprintf("results/split_%s_AF_pooled_quantile_counts.RData", divergence_time), version = 2)



