#Function that will read in some UKBB summary statistics (genotyped) and then build a reference panel using 1000G (EUR) then output the MAF & the LD
ukbb_analysis_inference = function(name_of_sumstats_file,
                                   BOLT                          = TRUE,
                                   quanitative_scale_convert     = 0.22,
                                   KG_Reference_Panel_Population = "EUR"){

  #### CURRENTLY ONLY WORKS WITH A SINGLE CHROMOSOME #####
  #Load in the library first
  install.packages("/well/mcvean/mtutert/myPackages/InferLD_0.1.0.tar.gz")
  library(InferLD)
  #Read in some list of summary statistics (from BOLT?)
  ukbb_sumstats = read.table(sprintf("/well/mcvean/ukbb12788/mtutert/%s",name_of_sumstats_file), header = T)
  #Cut down to the first 1000 summary statistics (NEED A BETTER WAY TO DO THIS)
  ukbb_sumstats = ukbb_sumstats[1:200,]
  #Assume BOLT summary statistics
  if (BOLT == TRUE) {
  #Extract out the SE's and convert to the quantitative scale
  #Run our inference method using the 1000G results, as a result, we need to download the haplotypes using PLINK
  #Extract the populations of interest
  KG_pop_identifiers         = read.table("/well/mcvean/mtutert/1000_Genomes/VCF/integrated_call_samples_v3.20130502.ALL.panel", header = T)
  KG_subset_pop_samples      = KG_pop_identifiers[which(KG_pop_identifiers$super_pop == KG_Reference_Panel_Population),1]
  write.table(KG_subset_pop_samples,"KG_subset_pop_samples", quote = F, row.names = F, col.names = F)
  #Use PLINK2 to extract out the specific haplotype segments we are using
  gwas_rsid                  = ukbb_sumstats$SNP
  write.table(gwas_rsid,"gwas_rsid", quote = F, row.names = F, col.names = F)
  system("/well/mcvean/mtutert/software/plink2 --max-alleles 2 --vcf /well/mcvean/mtutert/1000_Genomes/VCF/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz --keep KG_subset_pop_samples --extract gwas_rsid --export haps --out 1000G_test1")
  ref_hap_panel = t(read.table("1000G_test1.haps", header = F))
  #Remove SNPs in GWAS SumStats that are NOT in 1000G and those non-segreating
  #Extract rsids from 1000G file
  KG_rsid                 = ref_hap_panel[2,]
  KG_rsid_missing_indexes = which(!(gwas_rsid%in% KG_rsid))
  ref_hap_panel           = ref_hap_panel[6:length(ref_hap_panel[,1]),]
  nrow                    = dim(ref_hap_panel)[1]
  ncol                    = dim(ref_hap_panel)[2]
  ref_hap_panel           = matrix(as.numeric(as.matrix(ref_hap_panel)),nrow = nrow, ncol = ncol)
  KG_rsid_fixed           = which(colMeans(ref_hap_panel) == 1)
  ukbb_sumstats           = ukbb_sumstats[-c(KG_rsid_missing_indexes,KG_rsid_fixed),]
  ref_hap_panel           = ref_hap_panel[,-c(KG_rsid_fixed)]
  gwas_observed_SE        = ukbb_sumstats$SE/0.17

  library(matrixStats)
  inference_results_1000G = LD_from_GSHMM(ref_panel_haplotypes = (ref_hap_panel),
                                        Fst                  = 0.09,
                                        betas                = FALSE,
                                        alpha                = 1e1,
                                        nSamples             = 5,
                                        recomb_rate          = 1e-100,
                                        weights_resolution   = 20,
                                        likelihood_toggle    = TRUE,
                                        se_observed          = gwas_observed_SE,
                                        LD_Infer             = FALSE,
                                        genetic_map          = FALSE,
                                        chain_likelihood     = TRUE,
                                        nChains              = 1,
                                        recombination        = FALSE)
  #Allele Freq Results#
  allele_freq_results  = inference_results_1000G[[1]][[4]]
  total_updates        = ncol(allele_freq_results)
  allele_freq_mean     = rowMeans(allele_freq_results[,(.9*total_updates):total_updates])
  }
}
