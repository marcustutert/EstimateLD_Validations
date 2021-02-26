################################
#Parameters to change (Possibly)
library(data.table)
################################
nsnps            = 50
total_regions    = 10 #This I will need to change manually!
#####################
regions              = as.numeric(snakemake@params$regions)
super_pop            = snakemake@params$populations

#Split up chr22 into the chunks that we need
regions_total    = seq(1,nsnps*total_regions+1,nsnps)
#Pick the regions we need by subsetting the bim file rsids to the 1:500, 501:100.....
ukbb_bim         = read.table("prior_draw_analysis_panels/ukbb_genotype_qc_wba_chr22.bim", header = F)
#Extract the rsids that we need for each section
ukbb_rsid_region = ukbb_bim[regions_total[regions]:regions_total[regions + 1],2]
#Write out these rsids
write.table(ukbb_rsid_region, sprintf("prior_draw_analysis_panels/regions/rsid_region_%s", regions), quote = F, row.names = F, col.names = F)
#Now we want to extract from both the 1000G population & from the UKBB data this set of SNPs
#Note that we will need to do this from the PHASED HAPLOTYPE data in BGEN format
#And we will need to filter out the non WBA ancestry individuals at the same time
system(sprintf("/well/mcvean/mtutert/software/plink2 --bgen prior_draw_analysis_panels/ukb_hap_chr22_v2.bgen ref-first --keep prior_draw_analysis_panels/wba_sample_qc_ids --sample prior_draw_analysis_panels/ukbb_imp_v3.sample --extract prior_draw_analysis_panels/regions/rsid_region_%s --export haps --out prior_draw_analysis_panels/regions/ukbb_region_%s_matched_1000G_%s --oxford-single-chr 22 -max-alleles 2", regions, regions, super_pop))
Sys.sleep(20)
#Now we want to extract these same SNPs from the 1000G population that we are looking at
system(sprintf("/well/mcvean/mtutert/software/plink2 --vcf prior_draw_analysis_panels/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --keep prior_draw_analysis_panels/sample_pops/%s_samples --extract prior_draw_analysis_panels/regions/rsid_region_%s --export haps --out prior_draw_analysis_panels/regions/1000G_%s_region_%s --max-alleles 2", super_pop, regions, super_pop, regions))
Sys.sleep(20)
#Now read in the UKBB panel
ref_panel  = fread(sprintf("prior_draw_analysis_panels/regions/1000G_%s_region_%s.haps", super_pop, regions), header = F)
#Extract the rsids
ref_rsid   = ref_panel[,2]
#Add the rsid's as column names (trust me, this makes things WAY easier down the road)
rownames(ref_panel) = unlist(ref_rsid)
#Split off the naming notation for Oxford haps format and convert to normal panel
ref_panel  = t(ref_panel[,-c(1:6)])
#First read in the 1000G haplotype panel and the UKBB panel into R
ukbb_panel = fread(as.matrix(sprintf("prior_draw_analysis_panels/regions/ukbb_region_%s_matched_1000G_%s.haps", regions, super_pop), header = F, stringsAsFactors = F))
#Extract the rsids
ukbb_rsid   = (ukbb_panel[,2])
#Split off the naming notation for Oxford haps format and convert to normal panel
ukbb_panel  = t(ukbb_panel[,-c(1:6)])

####CHANGE THIS SO WE LOOK AT SAME SNPs across ALL SIMULATIONS!
remove_ukbb_panel = which(!unlist(ukbb_rsid) %in% unlist(ref_rsid))
remove_kg_panel   = which(!unlist(ref_rsid) %in% unlist(ukbb_rsid))

#Remove variants from the panel
ukbb_panel = ukbb_panel[,-c(remove_ukbb_panel)]
ref_panel  = ref_panel[,-c(remove_kg_panel)]

if (ncol(ref_panel) != ncol(ukbb_panel)) {
  #Filtering has gone wrong!
  break
}

#Now we also wish to remove any of the non-segregating variants in the reference panel:
high_fixed_freq_ref_snps = which(colMeans(ref_panel) == 1)
low_fixed_freq_ref_snps  = which(colMeans(ref_panel) == 0)
remove_ref_snps_index    = unique(c(high_fixed_freq_ref_snps,low_fixed_freq_ref_snps))
#Remove these SNPs from both the ref & gwas panel
if (length(remove_ref_snps_index)>0) {
  ukbb_panel = ukbb_panel[,-c(remove_ref_snps_index)]
  ref_panel  = ref_panel[,-c(remove_ref_snps_index)]
}

if (ncol(ref_panel) != ncol(ukbb_panel)) {
  #Filtering has gone wrong!
  break
}

#Write out the tables to rescomp
write.table(ref_panel, sprintf("prior_draw_analysis_panels/regions/matched_panels/1000G_%s_region_%s_matched_ukbb_ref_panel.haps", super_pop, regions), quote = F, row.names = F, col.names = T)
write.table(ukbb_panel,sprintf("prior_draw_analysis_panels/regions/matched_panels/ukbb_region_%s_matched_1000G_%s.haps", regions, super_pop), quote = F, row.names = F, col.names = F)
