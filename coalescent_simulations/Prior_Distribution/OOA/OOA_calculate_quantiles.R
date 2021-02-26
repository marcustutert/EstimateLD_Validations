#####OOA_calculate_quantiles######
#Script to take the Rdata objec¢ts generated across all replicates and chunks and summarize the
library(data.table)
#install.packages("abind")
library(abind)
#Input the snakemake parameters
replicate            = snakemake@params$replicates
population_pairs     = snakemake@params$paired_values

#Get populations (GWAS & Reference)
split_pops = strsplit(population_pairs, "_")
gwas_pop   = split_pops[[1]][1]
ref_pop    = split_pops[[1]][3]
fst        = split_pops[[1]][5]

#Read in all the chunks together for each replicate (and each divergence time) for the Allele Frequencies
AF_list              = do.call('rbind', lapply(list.files(path = "results/OOA", pattern = sprintf("%s_GWAS_%s_Ref_%s_AF_Replicate_%s_Chunk",gwas_pop,ref_pop,fst,replicate),full.names = T), readRDS))

#Read in all chunks, per replicate (for each divergence time) for the LD
LD_file_list         = list.files(path = "results/OOA", pattern = sprintf("%s_GWAS_%s_Ref_%s_LD_Replicate_%s_Chunk",gwas_pop,ref_pop,fst,replicate), full.names = T)
LD_Array             = do.call('abind', lapply(LD_file_list, readRDS))
print(dim(LD_Array))
#Get the quantiles
nquantiles           = 100
LD_quantile_matrx    = apply(LD_Array, 1:2, quantile, prob = c(seq(0,1,length = nquantiles)), na.rm = T, type = 6)

#Read in the actual GWAS matrix and calculate the LD
GWAS    = as.matrix(fread(sprintf("msprime_data/OOA/matched_panels/%s_GWAS_matched_to_%s_replicate_%s.csv", gwas_pop, ref_pop, replicate, header = T)))
GWAS_LD = cov.wt(GWAS,cor = TRUE, method = "ML")
GWAS_LD = GWAS_LD$cor
GWAS_AF = colMeans(GWAS)
nsnps   = ncol(GWAS_LD)

#Loop through rows and columns in M to find which index of the array matches up closest with element M[i,j]
LD_results = matrix(data = NA, nrow = nsnps, ncol = nsnps)
for (i in 1:nsnps) {
  print(i)
  for (j in i:nsnps) {
    true_value         = GWAS_LD[i,j]
    #Subset Q to the ith and jth element (vector of nquantiles)
    quantiles          = LD_quantile_matrx[,i,j]
    LD_results[i,j]    = (which.min(abs(quantiles-true_value)))
  }
}

# #Get the AF quantiles across all nsnps
AF_results            = c()
#Store indexes of the SNPs to remove (precision errors)
#Do this per sample as a list, then uniquify
# AF_remove_index       = list()
# af_boundary_tolerance = 1e-5
# ##### Filter SNPs on the Boundary #######
# for (i in 1:(dim(AF_list)[2])) {
#   #Loop through alll columns (Samples)
#   AF_remove_index[[i]] = which(AF_list[,i] > (1 -  af_boundary_tolerance))
# }
#Get unique subset of SNPs
#unique_remove_AF_Index = unique(unlist(AF_remove_index))
all_snps      = seq(1:nsnps)
#if (length(unique_remove_AF_Index)>0) {
#  filtered_snps = all_snps[-unique_remove_AF_Index]}
#if (length(unique_remove_AF_Index) == 0) {
#  filtered_snps = all_snps
#}

#Generate the quantiles for AF across columns
AF_quantiles = apply(AF_list,1, quantile, prob = c(seq(0,1,length = nquantiles)), na.rm = T, type = 6)
print(dim(AF_quantiles))
for (i in all_snps) { #Only sort through SNPs segreating (note filter above)
  #Get the quantiled data
  #data      = AF_list[i,]
  #quantiles = split(data, cut(data, quantile(data, prob = 0:nquantiles / nquantiles, names = FALSE), include.lowest = TRUE))
  #intervals = names(quantiles)
  #Grep out the second element of each name to get the overlapping sets
  #boundaries = c()
  #SNP is too close to boundary
  #for (j in 1:length(intervals)) {
  #  boundaries[j]= as.numeric(substr(sub(".*,", "", intervals[j]) ,1,nchar(sub(".*,", "", intervals[j]) )-1))
  #}
  #Now look for where in the data does our true value lie (to fit the quantiles)
  pi_pop                 = GWAS_AF[i]
  AF_results[i]          = which.min(abs(AF_quantiles[,i] - pi_pop))

}
saveRDS(LD_results,sprintf("results/OOA/%s_GWAS_%s_Ref_%s_LD_quantile_counts_replicate_%s.RData", gwas_pop, ref_pop, fst, replicate), version = 2)
saveRDS(AF_results,sprintf("results/OOA/%s_GWAS_%s_Ref_%s_AF_quantile_counts_replicate_%s.RData", gwas_pop, ref_pop, fst, replicate), version = 2)


