#Now do the sumstat imputation across one set of SNPs
set.seed(5223)
#Source R functions that we will use
source("/well/mcvean/mtutert/thesis/distr_ld/helper_functions.R")
#Read in snakemake params
replicates = snakemake@params$replicates
#Read in the reference table and the gwas sumstats
ref      = read.table(sprintf("msprime_data/filtered_panels/filtered_ref_%s",replicates),header = T)
sumstats = read.table(sprintf("msprime_data/filtered_panels/filtered_gwas_%s",replicates),header = T)
#Create list of imputed SNPs to use, 50% of the data
nsnps              = nrow(sumstats)
imputed_snps_index = sample(1:nsnps, 0.5*nsnps, replace = F) #Note that this refers to the INDEX and NOT the RSID

#Extract the LD from reference panel
LD                 = LD_Matrix(ref)

#Perform the sumstat imputation
results = sumstat_impute(typed_snps         = setdiff(1:nsnps,imputed_snps_index),
                         untyped_snps_index = imputed_snps_index,
                         genotyped_sumstats = sumstats[-imputed_snps_index,],
                         imputed_sumstats   = sumstats[imputed_snps_index,],
                         LD                 = LD)

#Write out the results to /results directory as an Rdata object
saveRDS(results,sprintf("results/sumstat_imputed_ref_%s",replicates), version =  2)
