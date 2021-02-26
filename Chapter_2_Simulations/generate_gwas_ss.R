#Generate GWAS sumstats under the null model
source("/well/mcvean/mtutert/thesis/distr_ld/helper_functions.R")

#Note that this only needs to be done a single time
#Read in the GWAS panel
gwas_panel = read.csv("msprime_data/gwas_test.csv", header = T)
results = msprime_gwas_sumstats(gwas_haplotypes = gwas_panel)

#Write out: GWAS sumstats, the GWAS panel (matched filtered), Reference Panel (matched filtered)
write.table(results[[2]], "msprime_data/gwas_sumstats", quote = F, col.names = T, row.names = F)
