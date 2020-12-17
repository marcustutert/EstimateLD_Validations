
#Function that will read in some UKBB summary statistics (genotyped) and then build a reference panel using 1000G (EUR) then output the MAF & the LD
ukbb_analysis_inference = function(name_of_sumstats_file,
                                   BOLT                          = TRUE,
                                   quanitative_scale_convert     = 0.22,
                                   KG_Reference_Panel_Population = c("AFR"),
                                   imputation                    = TRUE){

  #### CURRENTLY ONLY WORKS WITH A SINGLE CHROMOSOME #####
  #Load in the library first
  install.packages("/well/mcvean/mtutert/myPackages/InferLD_0.1.0.tar.gz")
  library(data.table)
  library(InferLD)
  library(matrixStats)

  #Read in some list of summary statistics (from BOLT?)
  ukbb_sumstats = fread(sprintf("/well/mcvean/ukbb12788/mtutert/%s",name_of_sumstats_file), header = T)
  #Cut down to the first nsnps_genotyped summary statistics
  nsnps_genotyped = 1000
  nSamples        = 5

  if (BOLT == TRUE) {
    #Extract out the SE's and convert to the quantitative scale
    #Run our inference method using the 1000G results, as a result, we need to download the haplotypes using PLINK
    #Extract the populations of interest
    KG_pop_identifiers         = read.table("/well/mcvean/mtutert/1000_Genomes/VCF/integrated_call_samples_v3.20130502.ALL.panel", header = T)
    samples_ll                 = KG_pop_identifiers[which(KG_pop_identifiers$super_pop == KG_Reference_Panel_Population),1:3]
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
    ukbb_sumstats           = ukbb_sumstats[-c(KG_rsid_missing_indexes),]
    ref_hap_panel           = ref_hap_panel[6:length(ref_hap_panel[,1]),]
    colnames(ref_hap_panel) = KG_rsid
    nrow                    = dim(ref_hap_panel)[1]
    ncol                    = dim(ref_hap_panel)[2]
    haps_names              = rep("haps",nrow)
    ref_hap_panel           = matrix(as.numeric(as.matrix(ref_hap_panel)),nrow = nrow, ncol = ncol,dimnames = list(haps_names,colnames(ref_hap_panel)))
    KG_rsid_fixed           = which(colMeans(ref_hap_panel) > 0.98 | colMeans(ref_hap_panel) < 0.02)
    ref_hap_panel           = ref_hap_panel[,-c(KG_rsid_fixed)]
    #Check to see which rsids are left to remove in the SumStats
    rsid_difference         = which(!ukbb_sumstats$SNP %in%  colnames(ref_hap_panel))
    ukbb_sumstats           = ukbb_sumstats[-c(rsid_difference),]
    gwas_observed_SE        = ukbb_sumstats$SE/0.17
    if(length(gwas_observed_SE)!=ncol(ref_hap_panel)){print("WARNING!")}
    if(length(setdiff(ukbb_sumstats$SNP,colnames(ref_hap_panel))) != 0){print("rsids are mismatched: ABORT!")}
    write.table(ukbb_sumstats$SNP[1:nsnps_genotyped],"ld_rsid",quote = F, row.names = F, col.names = F)
        #Take KG_pop_identifiers to create a table of ID/subpopulation/population
    #Use this for the heatmap
    heatmap_id_input = KG_pop_identifiers[which(KG_pop_identifiers$super_pop == KG_Reference_Panel_Population),1:3]
    write.table(heatmap_id_input,"heatmap_id_input")
    nhaps = ncol(ref_hap_panel)
    inference_results_1000G = LD_from_GSHMM(ref_panel_haplotypes = (ref_hap_panel[,1:nsnps_genotyped]),
                                          fst                  = 0.1,
                                          betas                = FALSE,
                                          alpha                = 1e3,
                                          nSamples             = nSamples,
                                          recomb_rate          = 1e-300,
                                          weights_resolution   = 5,
                                          likelihood_toggle    = TRUE,
                                          se_observed          = (gwas_observed_SE^2)[1:nsnps_genotyped],
                                          LD_Infer             = TRUE,
                                          genetic_map          = FALSE,
                                          chain_likelihood     = TRUE,
                                          nChains              = 1,
                                          recombination        = FALSE,
                                          case_control_constant = 122848)
    HW_Array              = inference_results_1000G$Gibbs_Array
    saveRDS(HW_Array,"HW_test")
    allele_freq_results  = inference_results_1000G$inferred_af_given_weights
    ll                   = inference_results_1000G$log_likelihood
    total_updates        = ncol(allele_freq_results)
    allele_freq_mean     = rowMeans(allele_freq_results[,(.9*total_updates):total_updates])
    true_inferred_r2 = (summary(lm(allele_freq_mean~ukbb_sumstats$A1FREQ[1:nsnps_genotyped]))$r.squared)
    true_reference_r2 = (summary(lm(colMeans(ref_hap_panel[,1:nsnps_genotyped])~ukbb_sumstats$A1FREQ[1:nsnps_genotyped]))$r.squared)
    print(sprintf("True to inferred correlation (r2) is %s",true_inferred_r2 ))
    print(sprintf("True to reference correlation (r2) is %s",true_reference_r2 ))
    jpeg(file = "af_accuracy.jpeg")
    plot(ukbb_sumstats$A1FREQ[1:nsnps_genotyped],allele_freq_mean)
    dev.off()
    jpeg(file = "af_reference_accuracy.jpeg")
    plot(ukbb_sumstats$A1FREQ[1:nsnps_genotyped],colMeans(ref_hap_panel[,1:nsnps_genotyped]))
    dev.off()
    jpeg(file = "ll_plot.jpeg")
    plot(ll[2:length(ll)])
    dev.off()

    #Write a function which has as input a reference panel, and outputs likelihoods
    likelihood_output = function(ref_panel = ref_hap_panel[,1:nsnps_genotyped],
                                 gwas_observed = (gwas_observed_SE^2)[1:nsnps_genotyped],
                                 weight_matrix,
                                 noise = 1e3,
                                 case_control_constant = 122848){
      #Normalize the weight matrix
      weight_matrix_sum           = (colSums(weight_matrix))
      weight_matrix_normalized    = weight_matrix/weight_matrix_sum
      allele_freqs                = colSums(ref_panel*weight_matrix_normalized)
      SE                          = (1/(case_control_constant*(allele_freqs*(1-allele_freqs))))
      ratio                       = gwas_observed/SE
      ll                          = sum(dgamma(x = ratio, shape = noise, rate = noise,log = T))
      return(ll)
    }

    #Get likelihoods, by creating a weight matrix for the different cases
    #Note that we don't have haps (yet for the ll)
    true_ll = sum(log(dgamma(x = rep(1,length(allele_freq_mean)), shape = 1e3, rate = 1e3)))
    print(sprintf("True likelihood (GWAS Panel) is %s",true_ll))
    print(sprintf("Final Inferred Likelihood is %s",ll[length(ll)]))
    #Equal weights across the panel
    equal_weight_matrix = matrix(data = 1, nrow = nrow(ref_hap_panel),ncol = nsnps_genotyped)
    equal_weight_ll     = likelihood_output(weight_matrix = equal_weight_matrix)
    print(sprintf("Ref likelihood is %s",equal_weight_ll))

    #Only EUR weights, first find index of all EUR indivs
    EUR_indivs = which(samples_ll$super_pop == "EUR")
    #Create weight matrix
    EUR_weight_matrix = matrix(data = 0, nrow = nrow(ref_hap_panel),ncol = nsnps_genotyped)
    EUR_weight_matrix[EUR_indivs,] = 1
    #print(EUR_weight_matrix[,1])
    EUR_weight_ll     = likelihood_output(weight_matrix = EUR_weight_matrix)
    print(sprintf("EUR likelihood is %s",EUR_weight_ll))

    #Only EUR weights, first find index of all EUR indivs
    AFR_indivs = which(samples_ll$super_pop == "AFR")
    #Create weight matrix
    AFR_weight_matrix = matrix(data = 0, nrow = nrow(ref_hap_panel),ncol = nsnps_genotyped)
    AFR_weight_matrix[AFR_indivs,] = 1
    AFR_weight_ll     = likelihood_output(weight_matrix = AFR_weight_matrix)
    print(sprintf("AFR likelihood is %s",AFR_weight_ll))

    #######IMPUTATION OF IMPUTED SNPS ########

    #Find range from the first SNP in our reference panel to our last (this will be what we span)
    first_snp_genotyped = colnames(ref_hap_panel)[1]
    last_snp_genotyped  = colnames(ref_hap_panel)[nsnps_genotyped]

    #Use PLINK2 to extract the AF from this range
    #system(sprintf("/well/mcvean/mtutert/software/plink2 --max-alleles 2 --freq --from %s --to %s --bgen /well/mcvean/ukbb12788/mtutert/impute_snp_qc/ukbb_imputed_qc_wba_chr1.bgen ref-first  --out impute_freqs",first_snp_genotyped,last_snp_genotyped))

    #Extract the rsids from the AF
    AF = read.table("impute_freqs.afreq",header = T)
    AF_rsid = AF[,2] #rsids from imputed data
    write.table(AF_rsid,"AF_rsid", quote = F, row.names = F, col.names = F)

    #Match up the SNPs
    system("/well/mcvean/mtutert/software/plink2 --max-alleles 2 --vcf /well/mcvean/mtutert/1000_Genomes/VCF/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz --keep KG_subset_pop_samples --extract AF_rsid --export haps --out 1000G_imputed")
    ref_hap_panel = t(read.table("1000G_imputed.haps", header = F))

    #We will have some SNPs in the UKBB imputed dataset that we will not have in the 1000G subset reference panel: remove those in imputed dataset

    KG_rsid                 = ref_hap_panel[2,]
    KG_rsid_missing_indexes = which(!(AF_rsid%in% KG_rsid))
    AF                      = AF[-c(KG_rsid_missing_indexes),] #Remove SNPs
    ref_hap_panel           = ref_hap_panel[6:length(ref_hap_panel[,1]),]
    colnames(ref_hap_panel) = KG_rsid
    nrow                    = dim(ref_hap_panel)[1]
    ncol                    = dim(ref_hap_panel)[2]
    haps_names              = rep("haps",nrow)
    ref_hap_panel           = matrix(as.numeric(as.matrix(ref_hap_panel)),nrow = nrow, ncol = ncol,dimnames = list(haps_names,colnames(ref_hap_panel)))
    KG_rsid_fixed           = which(colMeans(ref_hap_panel) > 0.98 | colMeans(ref_hap_panel) < 0.02) #Remove frequency of some variants from reference panel
    ref_hap_panel           = ref_hap_panel[,-c(KG_rsid_fixed)]
    rsid_difference         = which(!AF[,2] %in% colnames(ref_hap_panel))
    AF                      = AF[-c(rsid_difference),]
    #Remove duplicated variants
    AF = AF[!duplicated(AF[,2]), ]
    if(length(AF[,2])!=ncol(ref_hap_panel)){print("WARNING!")}
    if(length(setdiff(AF[,2],colnames(ref_hap_panel))) != 0){print("rsids are mismatched: ABORT!")}
    nsnps_imputed = ncol(ref_hap_panel)

    #Take the weights that we have calculated (For now it'll just be the last sample-->Deal with burn in later?)
    HW_Array              = inference_results_1000G$Gibbs_Array
    nhaps                 = nrow(ref_hap_panel)
    HW_Per_Sample_Array   = HW_Array[,,(nSamples-1)*nhaps]
    #print(HW_Per_Sample_Array)
    single_column_weights = HW_Per_Sample_Array[,1] #First column of weights
    imputed_weights       = matrix(rep(single_column_weights,nsnps_imputed),nrow = length(single_column_weights),ncol =nsnps_imputed,byrow = FALSE)
    imputed_af            = colSums(imputed_weights * ref_hap_panel)
    print(summary(lm(colMeans(ref_hap_panel)~AF[,5]))$r.squared)
    print(summary(lm(imputed_af~AF[,5]))$r.squared)

    jpeg(file = "af_imputed_accuracy_ref.jpeg")
    plot(colMeans(ref_hap_panel),AF[,5])
    dev.off()

    jpeg(file = "af_imputed_accuracy_inferred.jpeg")
    plot(imputed_af,AF[,5])
    dev.off()

    #Save Data
    LD = inference_results_1000G$LD_Array
    saveRDS(LD, "LD")
    print("DONE")
  }
}
