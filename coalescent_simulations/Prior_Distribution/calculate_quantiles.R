#Script to take the Rdata objects generated across all replicates and chunks and summarize the
library(data.table)
#install.packages("abind")
library(abind)
#Input the snakemake parameters
replicate                     = snakemake@params$replicates
divergence_time_fst_parameter = snakemake@params$paired_values #Do this across all divergence/Fst values (grid of graphs in the end)

divergence_time = strsplit(divergence_time_fst_parameter, "_")[[1]][1]
fst_parameter   = strsplit(divergence_time_fst_parameter, "_")[[1]][2]

#Read in all the chunks together for each replicate (and each paired fst-divergence time) for the Allele Frequencies
AF_list              = do.call('rbind', lapply(list.files(path = "results/pop_split", pattern = sprintf("^%s_%s_split_AF_Replicate_%s_Chunk",divergence_time,fst_parameter,replicate),full.names = T), readRDS))

#Read in all chunks, per replicate (for each divergence time) for the LD
LD_file_list         = list.files(path = "results/pop_split", pattern = sprintf("^%s_%s_split_LD_Replicate_%s_Chunk",divergence_time,fst_parameter,replicate),full.names = T)
LD_Array             = do.call('abind', lapply(LD_file_list, readRDS))
#Convert to a matrix

#Get the quantiles
nquantiles           = 100
LD_quantile_matrx    = apply(LD_Array, 1:2, quantile, prob = c(seq(0,1,length = nquantiles)), na.rm	= TRUE)

#Read in the actual GWAS matrix and calculate the LD
GWAS    = as.matrix(fread(sprintf("msprime_data/population_split/matched_panels/GWAS_panel_replicate_%s_split_%s.csv", replicate, divergence_time), header = T))
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
AF_results           = c()

for (i in 1:nsnps) {
  #Get the quantiled data
  data      = AF_list[i,]
  quantiles = split(data, cut(data, quantile(data, prob = 0:nquantiles / nquantiles, names = FALSE), include = TRUE))
  intervals = names(quantiles)
  #Grep out the second element of each name to get the overlapping sets
  boundaries = c()
  for (j in 1:length(intervals)) {
    boundaries[j]= as.numeric(substr(sub(".*,", "", intervals[j]) ,1,nchar(sub(".*,", "", intervals[j]) )-1))
  }
  #Now look for where in the data does our true value lie
  pi_pop                 = GWAS_AF[i]
  AF_results[i]          = which.min(abs(boundaries - pi_pop))
}
saveRDS(LD_results,sprintf("results/pop_split/%s_split_LD_quantile_counts_replicate_%s.RData", divergence_time_fst_parameter, replicate), version = 2)
saveRDS(AF_results,sprintf("results/pop_split/%s_split_AF_quantile_counts_replicate_%s.RData", divergence_time_fst_parameter, replicate), version = 2)

print("Done Calculating Quantiles")


