#Read in SNP QC file
for (i in 1:22) {
  snp_qcs             = fread(sprintf("/well/ukbb-wtchg/v3/imputation/ukb_mfi_chr%s_v3.txt",i), header = FALSE)
  #Filter based on QC criteria
  high_qual_info_rsid = snp_qcs[which(snp_qcs$V6 > 0.01 & snp_qcs$V7 > 0.8),]$V2
  #Write out table of variants
  write.table(high_qual_info_rsid, sprintf("/well/mcvean/ukbb12788/mtutert/impute_snp_qc/neale_variant_rsid_qc_chr%s",i),quote = F, row.names = F, col.names = F)
}
