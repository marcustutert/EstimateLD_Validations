# Helper Functions ----
match_rsids = function(ref_hap_panel,
                       sumstats,
                       AF_filter = c(0.01,0.99),
                       impute = FALSE){

#Read in reference panel
print(ref_hap_panel)
ref_hap_panel = t(read.table(ref_hap_panel,header = F))
if (impute == TRUE) {
  #Replace colnames in imputed file
  colnames(sumstats)[2] = "SNP"
}
KG_rsid                 = ref_hap_panel[2,]
KG_rsid_missing_indexes = which(!(sumstats$SNP %in% KG_rsid))
sumstats                = sumstats[-c(KG_rsid_missing_indexes),]
ref_hap_panel           = ref_hap_panel[6:length(ref_hap_panel[,1]),]
colnames(ref_hap_panel) = KG_rsid
nrow                    = dim(ref_hap_panel)[1]
ncol                    = dim(ref_hap_panel)[2]
haps_names              = rep("haps",nrow)
ref_hap_panel           = matrix(as.numeric(as.matrix(ref_hap_panel)),nrow = nrow, ncol = ncol,dimnames = list(haps_names,colnames(ref_hap_panel)))
KG_rsid_fixed           = which(colMeans(ref_hap_panel) > AF_filter[2] | colMeans(ref_hap_panel) < AF_filter[1])
ref_hap_panel           = ref_hap_panel[,-c(KG_rsid_fixed)]
rsid_difference         = which(!sumstats$SNP %in%  colnames(ref_hap_panel))
sumstats                = sumstats[-c(rsid_difference),]
if (impute == TRUE) {
  return(list(sumstats,ref_hap_panel))
}
else{
  sumstats$SE                = (sumstats$SE/0.17)^2 #DONT HARD CODE THIS IN
return(list(sumstats,ref_hap_panel))
  }
}

#Write a function which has as input a reference panel, and outputs likelihoods
likelihood_output = function(ref_hap_panel = ref_hap_panel,
                             gwas_observed,
                             weight_matrix,
                             noise                 = 1e3,
                             case_control_constant = 122848){
  #Normalize the weight matrix
  weight_matrix_sum           = (colSums(weight_matrix))
  weight_matrix_normalized    = weight_matrix/weight_matrix_sum
  allele_freqs                = colSums(ref_hap_panel*weight_matrix_normalized)
  SE                          = (1/(case_control_constant*(allele_freqs*(1-allele_freqs))))
  ratio                       = gwas_observed/SE
  ll                          = sum(dgamma(x = ratio, shape = noise, rate = noise,log = T))
  return(ll)
}

#Perform AF imputation ••Note only works with constants weights across genome atm
AF_Imputation = function(ref_hap_panel_string)
{
  ref_hap_panel     = read.table(sprintf("./analysis/reference_panel/%s_panel",ref_hap_panel_string))
  #Find range from the first SNP in our reference panel to our last (this will be what we span)
  first_snp_genotyped = colnames(ref_hap_panel)[1]
  last_snp_genotyped  = colnames(ref_hap_panel)[ncol(ref_hap_panel)]

  #Use PLINK2 to extract the AF from this range
  #system(sprintf("/well/mcvean/mtutert/software/plink2 --max-alleles 2 --freq --from %s --to %s --bgen /well/mcvean/ukbb12788/mtutert/impute_snp_qc/ukbb_imputed_qc_wba_chr1.bgen ref-first  --out ./temp/impute_freqs",first_snp_genotyped,last_snp_genotyped))

  #Extract the rsids from the AF
  AF = read.table("./temp/impute_freqs.afreq", header = T, stringsAsFactors = F)
  AF_rsid =AF[,2] #rsids from imputed data
  write.table(AF_rsid,"./temp/AF_rsid", quote = F, row.names = F, col.names = F)

  #Match up the SNPs
  system(sprintf("/well/mcvean/mtutert/software/plink2 --max-alleles 2 --pfile /well/mcvean/mtutert/1000_Genomes/1000G_chr1 --keep ./analysis/pop_sample_tables/%s_pop_samples --extract ./temp/AF_rsid --export haps --out ./analysis/imputed_reference_panel/imputed_%s_panel",ref_hap_panel_string,ref_hap_panel_string))

  #Perform rsid matching using above helper function
  results               = match_rsids(ref_hap_panel = sprintf("./analysis/imputed_reference_panel/imputed_%s_panel.haps",ref_hap_panel_string),
                                      sumstats = AF,
                                      impute = TRUE)

  ref_hap_panel_imputed = results[[2]]
  AF                    = results[[1]]
  #Remove duplicates (only happens with imputed data(?))
  AF = AF[!duplicated(AF[,2]), ]
  return(list(ref_hap_panel_imputed,AF))
}

#Create list of KG sample IDs
KG_samples_split_by_pop = function(pop_list,
                                   reference_panel = "1000G"){
  if (reference_panel == "1000G") {
    #Read in 1000G reference sample file
    KG_pop_identifiers         = read.table("/well/mcvean/mtutert/1000_Genomes/VCF/integrated_call_samples_v3.20130502.ALL.panel", header = T, stringsAsFactors = F)
  }
  #Loop through all populations in the pop_list (likely super populations)
  for (i in 1:length(pop_list)) {
    KG_pop_table                = KG_pop_identifiers[which(KG_pop_identifiers$super_pop == pop_list[[i]]),1]
    write.table(KG_pop_table, sprintf("./analysis/pop_sample_tables/%s_pop_samples", pop_list[[i]]), quote = F, row.names = F, col.names = F)
  }
}

#Function that gets the ll per reference panel (assuming equal weights on each population)
#Can add in fancier admixture later on!
ll_per_pop_panel = function(ref_hap_panel,
                            gwas_observed){

    #Create weight matrix that uses only the base reference panel
    weight_matrix                  = matrix(data = 1, nrow = nrow(ref_hap_panel), ncol = ncol(ref_hap_panel))
    pop_weight_ll                  = likelihood_output(weight_matrix = weight_matrix,
                                                       ref_hap_panel = ref_hap_panel,
                                                       gwas_observed = gwas_observed)
    return(pop_weight_ll)
}

#Function that takes in inference results and plots the results against the reference panel chosen
plot_genotyped_markers_vs_reference = function(){

}

#Function that takes in the inference results and plots: the LL curves, the reference panel line and the zoomed in inset graph portion
plot_ll_curves = function(ll,
                          population,
                          ref_hap_panel,
                          gwas_ss){

  # Font Options
  t         = list(size = 10)
  ll_graph  = plot_ly()

  #Add in first graph
  ll_graph  = add_trace(ll_graph,
                        x =~(1:length(ll)),
                        y = ~ll,
                        name = sprintf("%s Inferred Weights",population))

  #Add max-likelihood line (with base reference panel--no adjustment to weights)
  ref_hap_panel   = as.matrix(read.table(sprintf("./analysis/reference_panel/%s_panel",population),header = T))
  sumstats        = read.table(sprintf("./analysis/reference_panel/sumstats_matched_to_%s_panel",population), header = T)
  likelihood_line = ll_per_pop_panel(ref_hap_panel = ref_hap_panel,
                                     gwas_observed = gwas_ss)
  #Add in reference panel line
  ll_graph        = ll_graph  %>% add_segments(ll_graph,
                                               x    = 0,
                                               xend = length(ll),
                                               y    = likelihood_line,
                                               yend = likelihood_line,
                                               line = list(dash = "dash"),
                                               name = sprintf("%s Reference Panel",population))

  #Add in post burn-in (zoomed in graph)
  ll_graph  <- add_trace(ll_graph,
                         x          = ~seq(from = .8*length(ll), to = length(ll)),
                         y          = ~ll[(.8*length(ll)):length(ll)],
                         xaxis      = 'x2',
                         yaxis      = 'y2',
                         mode       = 'markers',
                         showlegend = F,
                         marker     = list(color = "#1f77b4"),
                         title      = "Post Burn In Convergence")

  #Add in design of showing zoomed in box plot
  ll_graph  = layout(ll_graph,
                     shapes = list(list(type = "rect",
                                        line = list(color = "black"),
                                        x0   = .9*length(ll),
                                        x1   = length(ll),
                                        xref = "x",
                                        y0   = ll[length(ll)]+ll[length(ll)]*1.3,
                                        y1   = ll[length(ll)]-ll[length(ll)]*1.3,
                                        yref = "y")))
  #Layout of Plot
  ll_graph  <- layout(ll_graph,
                      title  = sprintf("Log Likelihood Convergence (%s)",population),
                      xaxis  = list(title = "Update Number (Per Haplotype Sweep)"),
                      yaxis  = list(title = "Log Likelihood"),
                      xaxis2 = list(title = "Update Number", domain = c(0.7, 1.0), anchor='y2'),
                      yaxis2 = list(title = "Log Likelihood", domain = c(0.2, 0.4), anchor='x2'),
                      font = t,
                      margin = list(b = 100, l = 100))
  return(ll_graph)
}

# Analysis ----

#Function that will read in some UKBB summary statistics (genotyped) and then build a reference panel using 1000G (EUR)
#Then output the MAF & the LD and use the LD to do the imputation
ukbb_analysis_inference = function(name_of_sumstats_file         = "test_genotyped_chr1",
                                   name_of_imputed_sumstats      = "imputed_test_chr1",
                                   BOLT                          = TRUE,
                                   quanitative_scale_convert     = 0.22,
                                   KG_Reference_Panel_Population = c("EUR"),
                                   imputation                    = TRUE,
                                   nsnps                         = 100)
{

  # Load Libraries ----
  # Load in package (update if it need be)
  install.packages("/well/mcvean/mtutert/myPackages/InferLD_0.1.0.tar.gz")
  library(InferLD)

  #Load in other associated libraries
  library(data.table)
  library(matrixStats)
  library(plotly)
  library(webshot)

  #Read in table of summary statistics
  ukbb_sumstats   = fread(sprintf("/well/mcvean/ukbb12788/mtutert/%s",name_of_sumstats_file), header = T)
  nsnps           = nsnps
  ukbb_sumstats   = ukbb_sumstats[1:nsnps,]

  if (BOLT == TRUE) {

    #Prepare data for inference ----

    #Create list of all 1000G Populations in a list
    KG_super_pop_names = list("EAS",
                              "SAS",
                              "AFR",
                              "EUR",
                              "AMR")

    #Write out list of samples split by super population in 1000G
    KG_samples_split_by_pop(pop_list = KG_super_pop_names)
    gwas_rsid = ukbb_sumstats$SNP #SNPs in the GWAS
    write.table(gwas_rsid,"./temp/gwas_rsid", quote = F, row.names = F, col.names = F)

    #Loop through ref panels ----

    #Go through all 1000G populations to use as a reference panel
    #Will add in options to mix & match later on!!!!
    for (i in 1:length(KG_super_pop_names)) {
      system(sprintf("/well/mcvean/mtutert/software/plink2 --max-alleles 2 --pfile /well/mcvean/mtutert/1000_Genomes/1000G_chr1 --keep ./analysis/pop_sample_tables/%s_pop_samples --extract ./temp/gwas_rsid --export haps --out  ./analysis/reference_panel/1000G_%s --silent",KG_super_pop_names[[i]],KG_super_pop_names[[i]]))
      results       = match_rsids(ref_hap_panel = sprintf("./analysis/reference_panel/1000G_%s.haps", KG_super_pop_names[[i]]),sumstats = ukbb_sumstats)
      sumstats      = results[[1]]
      ref_hap_panel = results[[2]]
      #Write the matched sumstats/ref_hap_panel
      write.table(x = sumstats,      file = sprintf("./analysis/reference_panel/sumstats_matched_to_%s_panel",KG_super_pop_names[[i]]),quote = F,row.names = F, col.names = T)
      write.table(x = ref_hap_panel, file = sprintf("./analysis/reference_panel/%s_panel",KG_super_pop_names[[i]]),quote = F, row.names = F, col.names = T)

      #Inference Parameters ----
      nSamples        = 5
      alpha           = 1e3
      fst             = 0.1
    }
    #   #Run Genotype Inference (InferLD) ----
      inference_results_1000G = LD_from_GSHMM(ref_panel_haplotypes   = ref_hap_panel,
                                              fst                    = fst,
                                              betas                  = FALSE,
                                              alpha                  = alpha,
                                              nSamples               = nSamples,
                                              recomb_rate            = 1e-300,
                                              weights_resolution     = 5,
                                              likelihood_toggle      = TRUE,
                                              se_observed            = sumstats$SE,
                                              LD_Infer               = TRUE,
                                              genetic_map            = FALSE,
                                              chain_likelihood       = TRUE,
                                              nChains                = 1,
                                              recombination          = FALSE,
                                              case_control_constant  = 122848)
      saveRDS(inference_results_1000G,file = sprintf("./analysis/results/%s_inference_results",KG_super_pop_names[[i]]))
    # }

    fig          = list()
    fig_genotype = list()
    for (i in 1:length(KG_super_pop_names)) { #Do this across reference panels

      #Plot Likelihoods ----
      #For each reference panel plot:
      #The likelihood curve, a zoomed in version of the curve & a likelihood line of 'base' reference panel
      results  = readRDS(sprintf("./analysis/results/%s_inference_results", KG_super_pop_names[[i]]))
      ll       = results$log_likelihood

      ll_graph = plot_ll_curves(ll            = ll ,
                                population    = KG_super_pop_names[[i]],
                                ref_hap_panel = ref_hap_panel,
                                gwas_ss       = sumstats$SE)

      #Export plot
      export(ll_graph, file = sprintf("./analysis/figures/%s_ll.jpeg",KG_super_pop_names[[i]]))
      #Genotyped SNPs r2 & Graphs ----

      allele_freq_results   = results$inferred_af_given_weights
      HW_Array              = results$Gibbs_Array
      total_updates         = ncol(allele_freq_results)
      allele_freq_mean      = rowMeans(allele_freq_results[,(.9*total_updates):total_updates])
      true_inferred_r2      = (summary(lm(allele_freq_mean~sumstats$A1FREQ))$r.squared)
      true_reference_r2     = (summary(lm(colMeans(ref_hap_panel)~sumstats$A1FREQ))$r.squared)

      #Plot these results

      fig_genotype[[i]] = plot_ly()
      fig_genotype[[i]] <- plot_ly(data = iris)
      fig_genotype[[i]] <- plot_ly(x = ~sumstats$A1FREQ, y = ~colMeans(ref_hap_panel), name = "reference", mode = 'markers',type='scatter')
      fig_genotype[[i]] <- fig_genotype[[i]] %>% add_trace(fig_genotype[[i]], y = ~allele_freq_mean, name = "inferred",mode = 'markers',type='scatter',evaluate = TRUE)
      fig_genotype[[i]] <- plotly_build(fig_genotype[[i]])
    }
    genotype_accuracy_plot <- subplot(fig_genotype[[1]], fig_genotype[[2]],fig_genotype[[3]],fig_genotype[[4]],fig_genotype[[5]])
    export(genotype_accuracy_plot, file = "./analysis/figures/genotype_accuracy_plot.jpeg")

      #AF Imputation -----

      AF                    = AF_Imputation(KG_super_pop_names[[i]])[[2]]
      ref_hap_panel_imputed = AF_Imputation(KG_super_pop_names[[i]])[[1]]

      #Take the weights that we have calculated (For now it'll just be the last sample-->Deal with burn in later?)
      nhaps                 = nrow(ref_hap_panel_imputed)
      nsnps_imputed         = ncol(ref_hap_panel_imputed)
      #Filter down to per sample
      HW_Per_Sample_Array   = HW_Array[,,(nSamples-1)*nhaps]
      single_column_weights = HW_Per_Sample_Array[,1] #First column of weights
      imputed_weights       = matrix(rep(single_column_weights,nsnps_imputed),nrow = length(single_column_weights),ncol = nsnps_imputed,byrow = FALSE)
      imputed_af            = colSums(imputed_weights * ref_hap_panel_imputed)

      inferred_r2_imputed = (summary(lm(imputed_af~AF[,5]))$r.squared)
      print(inferred_r2_imputed)


    #First get max possible likelihood (assuming we match allele frequencies exactly)
    #true_ll = sum(log(dgamma(x = rep(1,length(allele_freq_mean)), shape = alpha, rate = alpha)))

    #Bind all likelihoods together


    #Get r2 measures (for genotyped SNPs across all possible reference panels): EUR/AFR/EUR+AFR
    #This requires the true genotyped AF's (which I can grab from bolt) & the reference panel AF's which I can output with PLINK?
    write.table(colnames(ref_hap_panel),"./temp/rsid_r2_genotyped",quote = F, row.names = F, col.names = F)
    #Use PLINK to subset the reference panel to these values
    system("/well/mcvean/mtutert/software/plink2 --max-alleles 2 --pfile /well/mcvean/mtutert/1000_Genomes/1000G_chr1 --keep  ./temp/EUR_pop_samples  --freq --out ./temp/EUR_AF --keep-allele-order")
    EUR_AF = read.table("./temp/EUR_AF.afreq", header = T)
    #Get correlation
    genotyped_AF_EUR_match     = sumstats$A1FREQ[sumstats$SNP %in% EUR_AF[,2]]
    EUR_AF_genotyped           = EUR_AF[EUR_AF[,2] %in% sumstats$SNP,5]
    EUR_genotyped_r2           = (summary(lm(EUR_AF_genotyped~genotyped_AF_EUR_match))$r.squared)
    system("/well/mcvean/mtutert/software/plink2 --max-alleles 2 --pfile /well/mcvean/mtutert/1000_Genomes/1000G_chr1 --keep  ./temp/AFR_pop_samples --freq --out ./temp/AFR_AF --keep-allele-order")
    AFR_AF = read.table("./temp/AFR_AF.afreq", header = T)
    genotyped_AF_AFR_match     = sumstats$A1FREQ[sumstats$SNP %in% AFR_AF[,2]]
    AFR_AF_genotyped           = AFR_AF[AFR_AF[,2] %in% sumstats$SNP,5]
    AFR_genotyped_r2 = (summary(lm(AFR_AF_genotyped~genotyped_AF_AFR_match))$r.squared)
    system("/well/mcvean/mtutert/software/plink2 --max-alleles 2 --pfile /well/mcvean/mtutert/1000_Genomes/1000G_chr1 --keep  ./temp/EUR_AFR_pop_samples --freq --out ./temp/EUR_AFR_AF --keep-allele-order")
    EUR_AFR_AF = read.table("./temp/EUR_AFR_AF.afreq", header = T)
    genotyped_AF_EUR_AFR_match = sumstats$A1FREQ[sumstats$SNP %in% EUR_AFR_AF[,2]]
    EUR_AFR_AF_genotyped           = EUR_AFR_AF[EUR_AFR_AF[,2] %in% sumstats$SNP,5]
    EUR_AFR_genotyped_r2 = (summary(lm(EUR_AFR_AF_genotyped~genotyped_AF_EUR_AFR_match))$r.squared)

    ##Do the same thing as above but this time for the imputed data
    #First it makes sense to pre-download all the allele frequencies with PLINK
    write.table(colnames(ref_hap_panel_imputed),"./temp/imputed_rsid",quote = F,col.names = F,row.names = F)
    system("/well/mcvean/mtutert/software/plink2 --max-alleles 2 --pfile  /well/mcvean/ukbb12788/mtutert/impute_snp_qc/ukbb_chr1_pgen --extract ./temp/imputed_rsid --freq --out ./temp/bgen_AF --keep-allele-order")
    imputed_af_EUR_AFR              = read.table("./temp/bgen_AF.afreq", header = T)
    imputed_af_EUR_AFR              = imputed_af_EUR_AFR[!(duplicated(imputed_af_EUR_AFR[,2])),]
    reference_imputed_af_EUR_AFR    = EUR_AFR_AF[EUR_AFR_AF[,2] %in% (imputed_af_EUR_AFR[,2]),]
    EUR_AFR_imputed_r2              = (summary(lm(reference_imputed_af_EUR_AFR[,5]~imputed_af_EUR_AFR[,5]))$r.squared)

    imputed_af_AFR              = read.table("./temp/bgen_AF.afreq", header = T)
    imputed_af_AFR              = imputed_af_AFR[!(duplicated(imputed_af_AFR[,2])),]
    reference_imputed_af_AFR    = AFR_AF[AFR_AF[,2] %in% (imputed_af_AFR[,2]),]
    AFR_imputed_r2              = (summary(lm(reference_imputed_af_AFR[,5]~imputed_af_AFR[,5]))$r.squared)

    imputed_af_EUR              = read.table("./temp/bgen_AF.afreq", header = T)
    imputed_af_EUR              = imputed_af_EUR[!(duplicated(imputed_af_EUR[,2])),]
    reference_imputed_af_EUR    = EUR_AF[EUR_AF[,2] %in% (imputed_af_EUR[,2]),]
    EUR_imputed_r2          = (summary(lm(reference_imputed_af_EUR[,5]~imputed_af_EUR[,5]))$r.squared)

    #Graphs ----

    #Plot Genotyped Data
    x <- list(
      title = "True AF (from GWAS)",
      titlefont = F
    )
    y <- list(
      title = "Comparison AF",
      titlefont = F)

    AF_Comparison <- plot_ly(data = iris)
    AF_Comparison <- AF_Comparison %>% add_trace(y=~allele_freq_mean, x=~sumstats$A1FREQ, name = 'Inferred AF',mode = 'markers')
    AF_Comparison <- AF_Comparison %>% add_trace(y=~1-EUR_AF_genotyped, x=~genotyped_AF_EUR_match, name = 'EUR Reference Panel Genotyped AF',mode = 'markers')
    AF_Comparison <- AF_Comparison %>% add_trace(y=~1-AFR_AF_genotyped, x=~genotyped_AF_AFR_match, name = 'AFR Reference Panel Genotyped AF',mode = 'markers')
    AF_Comparison <- AF_Comparison %>% add_trace(y=~1-EUR_AFR_AF_genotyped, x=~genotyped_AF_EUR_AFR_match, name = 'EUR & AFR Reference Panel Genotyped AF',mode = 'markers')
    export(AF_Comparison, file = "./figures/AF_Comparison.jpeg")

    #Plot Log_likelihood
    x <- list(
      title = "Update Number",
      titlefont = F
    )
    y <- list(
      title = "Log Likelihood",
      titlefont = F)
    ll_plot <- plot_ly(data = iris, x = ~seq(1:length(ll[2:length(ll)])), y = ~ll[2:length(ll)])
    ll_plot <- ll_plot %>% layout(xaxis = x, yaxis = y)
    export(ll_plot, file = "./figures/ll.jpeg")

    #Plot imputed data
    x <- list(
      title = "True Imputed AF",
      titlefont = F
    )
    y <- list(
      title = "Comparison Imputed AF",
      titlefont = F)

    AF_Imputed_Comparison <- plot_ly(data = iris)
    AF_Imputed_Comparison <- AF_Imputed_Comparison %>% add_trace(x = ~imputed_af, y= ~AF[,5], name = 'Inferred Imputed AF',mode = 'markers')
    AF_Imputed_Comparison <- AF_Imputed_Comparison %>% add_trace(x=~1-reference_imputed_af_EUR[,5], y=~imputed_af_EUR[,5], name = 'EUR Reference Panel Imputed AF',mode = 'markers')
    AF_Imputed_Comparison <- AF_Imputed_Comparison %>% add_trace(x=~1-reference_imputed_af_AFR[,5], y=~imputed_af_AFR[,5], name = 'AFR Reference Panel Imputed AF', mode = 'markers')
    AF_Imputed_Comparison <- AF_Imputed_Comparison %>% add_trace(x=~1-reference_imputed_af_EUR_AFR[,5],y=~imputed_af_EUR_AFR[,5], name = 'EUR & AFR Reference Panel Imputed AF', mode = 'markers')
    AF_Imputed_Comparison <- AF_Imputed_Comparison %>% layout(xaxis = x, yaxis = y)
    export(AF_Imputed_Comparison, file = "./figures/AF_Imputed_Comparison.jpeg")
    #Output to Table ----
    #Combine likelihood results
    results            = cbind(true_ll,
                               ll[length(ll)],
                               EUR_weight_ll,
                               AFR_weight_ll,
                               equal_weight_ll,
                               inferred_r2_imputed,
                               true_inferred_r2,
                               EUR_genotyped_r2,
                               AFR_genotyped_r2,
                               EUR_AFR_genotyped_r2,
                               EUR_AFR_imputed_r2,
                               AFR_imputed_r2,
                               EUR_imputed_r2)
    results_name       = c("True Likelihood",
                          "Inferred Max Likelkhood",
                          "EUR Likelihood",
                          "AFR Likelihood",
                          "EUR & AFR Likelihood",
                          "Inferred Imputed r2 (Gibbs)",
                          "Inferred Genotype r2 (Gibbs)",
                          "EUR Genotyped r2",
                          "AFR Genotyped r2",
                          "EUR & AFR Genotyped r2",
                          "EUR & AFR Imputed r2",
                          "AFR Imputed r2",
                          "EUR Imputed r2")
    colnames(results) = results_name
    write.table(results, "likelihood_table",quote = F, col.names = T, row.names = F)

    #SumStats Imputation ----

    #Extract from BOLT "Imputed" File
    #Then we search for these positions in the imputed BOLT file and choose inclusive set of them
    #Find where first and last SNP is in the genotyped SNPs
    ukbb_sumstats_imputed = as.data.frame(fread(name_of_imputed_sumstats,header = T))
    first_snp_genotyped = ukbb_sumstats$BP[1]
    last_snp_genotyped  = ukbb_sumstats$BP[length(ukbb_sumstats$BP)]

    #Cut down the BOLT Imputed file to this range
    ukbb_sumstats_imputed = subset(ukbb_sumstats_imputed, BP>first_snp_genotyped & BP<last_snp_genotyped)

    write.table(ukbb_sumstats_imputed$SNP, "./temp/ss_imp_af",quote = F,col.names = F, row.names = F)

    #Now we need to match this up with our reference panel in the region!
    system("/well/mcvean/mtutert/software/plink2 --max-alleles 2 --pfile /well/mcvean/mtutert/1000_Genomes/1000G_chr1 --keep ./temp/KG_subset_pop_samples --extract ./temp/ss_imp_af --export haps --out ./temp/1000G_BOLT_Imputed_Test")

    #Do the matching across the SNPs to the reference panel
    results = match_rsids(ref_hap_panel = "./temp/1000G_BOLT_Imputed_Test.haps", sumstats = ukbb_sumstats_imputed)
    bolt_imputed_ref_hap_panel = results[[2]]
    ukbb_sumstats_imputed      = results[[1]]

    #Now what we need to do is extend the weights from our results across the length of this imputed
    HW_Per_Sample_Array   = HW_Array[,,((nSamples-1)*nrow(bolt_imputed_ref_hap_panel))]
    single_column_weights = HW_Per_Sample_Array[,1] #First column of weights
    imputed_weights       = matrix(rep(single_column_weights,ncol(bolt_imputed_ref_hap_panel)),nrow = length(single_column_weights),ncol = ncol(bolt_imputed_ref_hap_panel),byrow = FALSE)

    #Perform the LD inference on this set of weights
    Qx_Array              = markovian_flow(HW_Matrix    = imputed_weights,
                                          start         = 1,
                                          end           = ncol(imputed_weights),
                                          recombination = F)

    LD_imputed            = LD_flow_expectation(HW_Matrix         = imputed_weights,
                                                ref_allele_matrix = bolt_imputed_ref_hap_panel,
                                                Qx_array          = Qx_Array,
                                                recombination     = F)
    saveRDS(LD_imputed,"LD_imputed")

    #LD Results ----

    #Get the LD across the individual level data
    #First thing we we will need to do is to convert the BGEN files to PLINK (keeping PHASE)
    write.table(colnames(bolt_imputed_ref_hap_panel),"./temp/indiv_rsid",quote = F,row.names = F,col.names = F)
    write.table(cbind(EUR_pop_samples,EUR_pop_samples),"./temp/EUR_pop_LD",quote = F,row.names = F,col.names = F)
    write.table(cbind(AFR_pop_samples,AFR_pop_samples),"./temp/AFR_pop_LD",quote = F,row.names = F,col.names = F)
    write.table(cbind(EUR_AFR_pop_samples,EUR_AFR_pop_samples),"./temp/EUR_AFR_pop_LD",quote = F,row.names = F,col.names = F)

    system("/apps/well/plink/1.90b3/plink --vcf /well/mcvean/mtutert/1000_Genomes/VCF/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz   --r square --out ./temp/EUR_REF_LD --extract ./temp/indiv_rsid --keep-allele-order --keep ./temp/EUR_pop_LD")
    #system("/apps/well/plink/1.90b3/plink --vcf /well/mcvean/mtutert/1000_Genomes/VCF/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz   --r square --out ./temp/AFR_REF_LD --extract ./temp/indiv_rsid --keep-allele-order --keep ./temp/AFR_pop_LD")
    #system("/apps/well/plink/1.90b3/plink --vcf /well/mcvean/mtutert/1000_Genomes/VCF/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz   --r square --out ./temp/EUR_AFR_REF_LD --extract ./temp/indiv_rsid --keep-allele-order --keep ./temp/EUR_AFR_pop_LD")
    system("/apps/well/plink/1.90b3/plink --bfile /well/mcvean/ukbb12788/mtutert/impute_snp_qc/plink/ukbb_imputed_plink_chr1 --r square --out ./temp/true_LD --extract ./temp/indiv_rsid --keep-allele-order")

    true_LD     = (as.matrix(read.table("./temp/true_LD.ld")))
    true_LD     = true_LD[upper.tri(true_LD)]
    inferred_LD = (LD_imputed)
    inferred_LD = inferred_LD[upper.tri(inferred_LD)]
    EUR_LD      = (as.matrix(read.table("./temp/EUR_REF_LD.ld")))
    EUR_LD      = EUR_LD[upper.tri(EUR_LD)]
    AFR_LD      = (as.matrix(read.table("./temp/AFR_REF_LD.ld")))
    AFR_LD      = AFR_LD[upper.tri(AFR_LD)]
    EUR_AFR_LD  = (as.matrix(read.table("./temp/EUR_AFR_REF_LD.ld")))
    EUR_AFR_LD  = EUR_AFR_LD[upper.tri(EUR_AFR_LD)]

    x <- list(
      title = "True LD (In-Sample)",
      titlefont = F
    )
    y <- list(
      title = "Comparison LD",
      titlefont = F)

    LD_Plot <- plot_ly(data = iris)
    LD_Plot <- LD_Plot %>% add_trace(x=~ c(true_LD),y = ~c(EUR_LD), name = 'EUR Reference LD', mode = 'markers')
    LD_Plot <- LD_Plot %>% add_trace(x=~ c(true_LD),y = ~c(inferred_LD), name = 'Inferred LD',mode = 'markers')
    # LD_Plot <- LD_Plot %>% add_trace(y = ~c(AFR_LD), name = 'AFR Reference LD', mode = 'markers')
    # LD_Plot <- LD_Plot %>% add_trace(y = ~c(EUR_AFR_LD), name = 'EUR & AFR Reference LD', mode = 'markers')
    # LD_Plot <- LD_Plot %>% layout(xaxis = x, yaxis = y)
    export(LD_Plot, file = "./figures/LD_Plot.jpeg")

  }
}
