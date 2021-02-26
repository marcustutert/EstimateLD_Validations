#Script to find which SNPs need to be filtered from the data
replicates = snakemake@params$replicates
#Read in the reference panel data (from msprime--in the form of a CSV)
ref = read.csv(sprintf("msprime_data/ref_%s.csv",replicates), header = T)
#Add the column names to the reference panel data
colnames(ref) = seq(1:ncol(ref))
#Find the SNPs that have <1% MAF
filtered_snps_ref = which(colMeans(ref) < 0.01 | colMeans(ref) > 0.99)
#Write out these the reference panel as a TABLE and the filtered SNPs
write.table(ref, sprintf("msprime_data/filtered_panels/ref_%s",replicates),quote = F, row.names = F, col.names = T) #With column names
write.table(filtered_snps_ref, sprintf("msprime_data/filtered_panels/filtered_snps_ref_%s",replicates),quote = F, row.names = F, col.names = F)
