
#Function that will read in some UKBB summary statistics (genotyped) and then build a reference panel using 1000G (EUR) then output the MAF & the LD
ukbb_analysis_inference = function(name_of_sumstats_file,
                                   BOLT                          = TRUE,
                                   quanitative_scale_convert     = 0.22,
                                   KG_Reference_Panel_Population = c("EUR","AFR"),
                                   imputation                    = TRUE){

  #### CURRENTLY ONLY WORKS WITH A SINGLE CHROMOSOME #####
  #Load in the library first
  install.packages("/well/mcvean/mtutert/myPackages/InferLD_0.1.0.tar.gz")
  library(data.table)
  library(InferLD)
  #Read in some list of summary statistics (from BOLT?)
  ukbb_sumstats = fread(sprintf("/well/mcvean/ukbb12788/mtutert/%s",name_of_sumstats_file), header = T)
  #Cut down to the first 1000 summary statistics (NEED A BETTER WAY TO DO THIS)
  #Assume BOLT summary statistics
  nsnps = 50
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
    library(matrixStats)
    #Take KG_pop_identifiers to create a table of ID/subpopulation/population
    #Use this for the heatmap
    heatmap_id_input = KG_pop_identifiers[which(KG_pop_identifiers$super_pop == KG_Reference_Panel_Population),1:3]
    write.table(heatmap_id_input,"heatmap_id_input")
    inference_results_1000G = LD_from_GSHMM(ref_panel_haplotypes = (ref_hap_panel[,1:nsnps]),
                                          fst                  = 0.1,
                                          betas                = FALSE,
                                          alpha                = 1e3,
                                          nSamples             = 10,
                                          recomb_rate          = 1e-300,
                                          weights_resolution   = 5,
                                          likelihood_toggle    = TRUE,
                                          se_observed          = (gwas_observed_SE^2)[1:nsnps],
                                          LD_Infer             = FALSE,
                                          genetic_map          = FALSE,
                                          chain_likelihood     = TRUE,
                                          nChains              = 1,
                                          recombination        = FALSE)
    #Allele Freq Results#
    allele_freq_results  = inference_results_1000G[[3]]
    ll                   = inference_results_1000G[[3]]
    total_updates        = ncol(allele_freq_results)
    allele_freq_mean     = rowMeans(allele_freq_results[,(.9*total_updates):total_updates])
    true_inferred_r2 = (summary(lm(allele_freq_mean~ukbb_sumstats$A1FREQ[1:nsnps]))$r.squared)
    true_reference_r2 = (summary(lm(colMeans(ref_hap_panel[,1:nsnps])~ukbb_sumstats$A1FREQ[1:nsnps]))$r.squared)
    print(sprintf("True to inferred correlation (r2) is %s",true_inferred_r2 ))
    print(sprintf("True to reference correlation (r2) is %s",true_reference_r2 ))
    jpeg(file = "af_accuracy.jpeg")
    plot(ukbb_sumstats$A1FREQ[1:nsnps],allele_freq_mean)
    dev.off()
    jpeg(file = "af_reference_accuracy.jpeg")
    plot(ukbb_sumstats$A1FREQ[1:nsnps],colMeans(ref_hap_panel[,1:nsnps]))
    dev.off()
    jpeg(file = "ll_plot.jpeg")
    plot(ll[2:length(ll)])
    dev.off()

    #Write a function which has as input a reference panel, and outputs likelihoods
    likelihood_output = function(ref_panel = ref_hap_panel[,1:nsnps],
                                 gwas_observed = (gwas_observed_SE^2)[1:nsnps],
                                 weight_matrix,
                                 noise = 1e3){
      #Normalize the weight matrix
      weight_matrix_sum           = (colSums(weight_matrix))
      weight_matrix_normalized    = weight_matrix/weight_matrix_sum
      allele_freqs                = colSums(ref_panel*weight_matrix_normalized)
      SE                          = (1/(122848*(allele_freqs*(1-allele_freqs))))
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
    equal_weight_matrix = matrix(data = 1, nrow = nrow(ref_hap_panel),ncol = nsnps)
    equal_weight_ll     = likelihood_output(weight_matrix = equal_weight_matrix)
    print(sprintf("Ref likelihood is %s",equal_weight_ll))

    #Only EUR weights, first find index of all EUR indivs
    EUR_indivs = which(samples_ll$super_pop == "EUR")
    #Create weight matrix
    EUR_weight_matrix = matrix(data = 0, nrow = nrow(ref_hap_panel),ncol = nsnps)
    EUR_weight_matrix[EUR_indivs,] = 1
    #print(EUR_weight_matrix[,1])
    EUR_weight_ll     = likelihood_output(weight_matrix = EUR_weight_matrix)
    print(sprintf("EUR likelihood is %s",EUR_weight_ll))

    #Only EUR weights, first find index of all EUR indivs
    AFR_indivs = which(samples_ll$super_pop == "AFR")
    #Create weight matrix
    AFR_weight_matrix = matrix(data = 0, nrow = nrow(ref_hap_panel),ncol = nsnps)
    AFR_weight_matrix[AFR_indivs,] = 1
    AFR_weight_ll     = likelihood_output(weight_matrix = AFR_weight_matrix)
    print(sprintf("AFR likelihood is %s",AFR_weight_ll))
    Gibbs_Array = inference_results_1000G[[1]]
    #print(seq(1,dim(Gibbs_Array)[2],dim(Gibbs_Array)[1]))
    #print("here")
    #HW_Per_Sample_Array = Gibbs_Array[,seq(1,dim(Gibbs_Array)[2],dim(Gibbs_Array)[1])]
    #saveRDS(ref_hap_panel[,1:nsnps], "ref_hap_panel")
    #saveRDS(object = HW_Per_Sample_Array,file = "HW_Per_Sample_Array")
    LD = inference_results_1000G[[2]]
    saveRDS(LD, "LD")
    saveRDS(object = Gibbs_Array,file = "HW_Per_Sample_Array")
    saveRDS(object = (gwas_observed_SE^2)[1:nsnps], "gwas_se")
    #saveRDS(object = ref_hap_panel[,1:nsnps],file = "ref_hap_panel")
    print("DONE")
  }
}
