#Generating weights from prior
#Split into jobs in order to improve speed of the inference
library(data.table)
source("/well/mcvean/mtutert/thesis/coalescent_coverage/helper_functions.R")
file                               = snakemake@params$replicates
chunk                              = snakemake@params$chunk
divergence_time_fst_parameter      = snakemake@params$paired_values #Do this across all divergence/Fst values (grid of graphs in the end)

print(divergence_time_fst_parameter)
divergence_time = strsplit(divergence_time_fst_parameter, "_")[[1]][1]
fst_parameter   = strsplit(divergence_time_fst_parameter, "_")[[1]][2]
print(divergence_time)
print(fst_parameter)

#Import the reference & gwas panels based on the divergence time and the replicate number (done in parallel through snakemake)
ref      = as.matrix(fread(sprintf("/well/mcvean/mtutert/thesis/coalescent_coverage/prior_draws/msprime_data/population_split/Ref_panel_replicate_%s_split_%s.csv", file, divergence_time), header = T))
gwas     = as.matrix(fread(sprintf("/well/mcvean/mtutert/thesis/coalescent_coverage/prior_draws/msprime_data/population_split/GWAS_panel_replicate_%s_split_%s.csv", file, divergence_time), header = T))

#Perform filtering (removing non-segregating and low freq SNPs)
res      = msprime_gwas_sumstats(gwas_haplotypes = gwas, reference_haplotypes = ref)
gwas     = res[[1]]
ref      = res[[2]]
nSamples = 500

#Write out GWAS & REF (matched) tables
write.table(gwas,sprintf("msprime_data/population_split/matched_panels/GWAS_panel_replicate_%s_split_%s.csv",file,divergence_time),quote = F,col.names = T, row.names = F)
write.table(ref,sprintf("msprime_data/population_split/matched_panels/Ref_panel_replicate_%s_split_%s.csv",file,divergence_time),quote = F,col.names = T, row.names = F)

#Back out the correct Fst given divergence time
nhaps_ref     = nrow(ref)
nhaps_gwas    = nrow(gwas)
nsnps         = ncol(ref)
effective_fst = as.numeric(fst_parameter) #Note that these values have already been backed out of the graph

#Create matrix & array to store AF and LD results
AF_Inferred_Results     =  matrix(data = NA, nrow = nsnps, ncol = nSamples)
LD_Results              =  array(data = NA, dim = c(nsnps, nsnps, nSamples))

#Loop across to get the draws I need
for (i in 1:nSamples) {
  print(i)
  #Draw from gamma_quantiled_weights nhaps times
  gamma_draw               = rgamma(n = nhaps_ref, shape =  1/( nhaps_ref * ( effective_fst / (1-effective_fst))), scale = ( nhaps_ref * (effective_fst/(1-effective_fst))))
  #Extend into matrix
  weight_matrix            = matrix(rep(gamma_draw,nsnps), ncol = nsnps)
  #Normalize matrix
  norm_weight_matrix       = weight_matrix/colSums(weight_matrix)[col(weight_matrix)]
  ####### Remove Ascertainment Bias
  #Ask if the weight matrix we are generating will
  AF_Inferred_Results[,i]  = colSums(ref*norm_weight_matrix)
  cov                      = cov.wt(ref,norm_weight_matrix[,1],cor = TRUE, method = "ML")
  LD_Results[,,i]          = cov$cor
}

#Save object in format Fst_#_Replicate_#_Chunk_#
saveRDS(object = LD_Results, file = sprintf("./results/pop_split/%s_split_LD_Replicate_%s_Chunk_%s.RData", divergence_time_fst_parameter, file, chunk), version = 2)
saveRDS(object = AF_Inferred_Results, file = sprintf("./results/pop_split/%s_split_AF_Replicate_%s_Chunk_%s.RData",divergence_time_fst_parameter, file, chunk), version = 2)

print("Done Drawing Weights")


