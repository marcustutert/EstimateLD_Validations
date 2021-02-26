#Use the filtered data
library(data.table)
replicates = snakemake@params$replicates
#Read in the reference panel
#Read in all the filtered SNPs
joint_snps_filtered_list = c()
snp_filtered_list = list.files(path = "msprime_data/filtered_panels", pattern = "filtered", full.names = T)
for (i in 1:length(snp_filtered_list)) {
  snps_filtered = fread(snp_filtered_list[i], header = F)
  joint_snps_filtered_list = append(joint_snps_filtered_list,snps_filtered)
}
joint_snps_filtered_list = unique(sort(unlist(joint_snps_filtered_list)))
#Write this out to check for mistakes (will be over written each time per parallel scripts but thats ok)
write.table(joint_snps_filtered_list, "msprime_data/filtered_panels/joint_snps_filtered_list", quote = F, row.names = F, col.names = F)
#Read in the reference panel
ref           = as.matrix(fread(sprintf("msprime_data/filtered_panels/ref_%s", replicates), header = T))
#Add colnames
colnames(ref) = as.character(seq(1:ncol(ref)))
#Read in the GWAS data
sumstats       = fread("msprime_data/gwas_sumstats", header = T)
#Find all the SNPs in the dataset
snps_in_gwas = unlist(sumstats[,1])
#Read in the original reference panel (just to see how many SNPs there are)
nsnps = ncol(ref)
#Find SNPs not in GWAS sumstats that we need to add to the remove list for the ref panel
snps_not_in_gwas = seq(1:nsnps)[-unlist(snps_in_gwas)]
#Remove SNPs not in GWAS sumstats (filtered out) and SNPs jointly not in any of the other reference panels (filtered out)
ref = ref[,-c(snps_not_in_gwas)]
#Now we need to remove from BOTH the ref and the GWAS sumstats, the joint_snps_filtered list SNPs (note that they correspond to the RSID now)
ref = ref[,-which(joint_snps_filtered_list %in% colnames(ref))]
#Remove the SNPs from the GWAS data now
sumstats = sumstats[which(unlist(sumstats[,1]) %in% colnames(ref)),]
write.table(ref,sprintf("msprime_data/filtered_panels/filtered_ref_%s", replicates), quote = F, col.names = T, row.names = F)
write.table(sumstats,sprintf("msprime_data/filtered_panels/filtered_gwas_%s",replicates), quote = F, col.names = T, row.names = F)
#Now we can remove any of the SNPs that are still in the GWAS data that are NOT in the ref data
print(all(unlist(sumstats[,1]) == colnames(ref))) #Check if I didn't fuck anything up this time

