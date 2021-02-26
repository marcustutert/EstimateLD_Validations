####Functions that will be used to simulate Data used in the notebooks #####
#This function uses the 1000G data located in /Data and exports the
simulate_genetic_drift = function(Fst,
                                  nhaps,
                                  nsnps,
                                  noise,
                                  weights_resolution = 5,
                                  KG_Genome_ref = FALSE){

  #Set number of haplotypes and variants to a small number so we can do testing and debugging

  path = dirname(getwd()) #Set path relative to testing folder
  #path <- dirname(getwd())

  #Create reference panel based on binomial distribution
  ref_hap_panel_filtered = matrix(data = rbinom(nhaps*nsnps, size = c(0,1), prob = 0.5), nrow = nhaps, ncol = nsnps)

  if (KG_Genome_ref == TRUE) {
    #Load in 1000G data from ./Test/Data
    ref_hap_panel     = as.matrix(read.table(sprintf("/Users/marcustutert/Desktop/Oxford_Dphil/InferLD_Validations/Data/ref_panel_filtered", path),
                                             header = TRUE))
    #Subset to nhaps and nsnps size to make testing easier
    ref_hap_panel_filtered     = ref_hap_panel[ 1:nhaps, 1:nsnps,drop = F ]
  }

  #Get AF of the variants in the reference panel
  if (nsnps == 1) {
    ref_variants_af = mean(ref_hap_panel_filtered)
  }
  else{
    ref_variants_af   = colMeans(ref_hap_panel_filtered)
  }

  # #Generate weight state space from gamma quantiled weight distribution (discretized)
  Gamma_Quantiled_Weights = qgamma(p = (1:weights_resolution)/(1+weights_resolution), shape = 1 / ( nhaps * (Fst/(1-Fst))), rate  = 1/( nhaps * (Fst/(1-Fst))) )
  Gamma_Quantiled_Weights = Gamma_Quantiled_Weights[is.finite(Gamma_Quantiled_Weights) & Gamma_Quantiled_Weights > 0 ]
  # #Create a GWAS matrix of weights (no recombination)
  # #Each row will be the same vector of weights
  gwas_haplotype_matrix = matrix(data =  sample(x       = Gamma_Quantiled_Weights,
                                                size    = nrow(ref_hap_panel_filtered),
                                                replace = TRUE),
                                 nrow = nrow(ref_hap_panel_filtered),
                                 ncol = ncol(ref_hap_panel_filtered),
                                 byrow = FALSE)

  gwas_sum                                  = colSums(gwas_haplotype_matrix)
  gwas_haplotype_matrix_norm                = sweep(gwas_haplotype_matrix, 2, gwas_sum, FUN = '/')
  #Get GWAS AF
  gwas_variants_af         = colSums(ref_hap_panel_filtered * gwas_haplotype_matrix_norm)

  #Filter out SNPs that are above/below a threshhold (0-100% exclusive)
  if (nsnps>1) {
    filter_ref_snp_index_1       = which(colMeans(ref_hap_panel_filtered) > 0.99)
    filter_ref_snp_index_2       = which(colMeans(ref_hap_panel_filtered) < 0.01)
    filter_gwas_snp_index_1      = which(gwas_variants_af > 0.99)
    filter_gwas_snp_index_2      = which(gwas_variants_af < 0.01)
    filter_snp_index_merged = unique(c(filter_ref_snp_index_1,filter_ref_snp_index_2,filter_gwas_snp_index_1,filter_gwas_snp_index_2))

    if (length(filter_snp_index_merged)>0) {
      ref_hap_panel_filtered     = ref_hap_panel_filtered[,-filter_snp_index_merged]
      gwas_variants_af           = gwas_variants_af[-filter_snp_index_merged]
      ref_variants_af            = ref_variants_af[-filter_snp_index_merged]
      gwas_haplotype_matrix      = gwas_haplotype_matrix[,filter_snp_index_merged]
    }
  }

  observed_gwas_se  = 1/(gwas_variants_af*(1-gwas_variants_af)) * rgamma(length(gwas_variants_af),shape = noise,rate = noise)
  ll                = sum(dgamma(observed_gwas_se/observed_gwas_se,shape = noise, rate = noise, log = T))
  data = list(observed_gwas_se,
              gwas_variants_af,
              ref_hap_panel_filtered,
              gwas_haplotype_matrix,
              ll)
  names(data) <- c("Observed_GWAS_SE", "GWAS_AF", "ref_panel", "gwas_matix", "ll")
  return(data)
}

#This function takes in the allele frequencies at each nSample point in the chain and sees how close they are
#To the true allele frequencies

allele_frequency_accuracy_across_chains = function(allele_frequencies, #This should be a matrix (or an array?)
                                                   true_gwas_af){        #When using simulated data)
  nsnps    = length(true_gwas_af)
  nSamples = ncol(allele_frequencies)
  p <- plot_ly()
  ma <- function(x, n = 2){filter(x, rep(1 / n, n), sides = 2)}
  for(k in 1:nsnps) {
    allele_frequency_error = abs(true_gwas_af[k]-allele_frequencies[k,])
    allele_frequency_error_average = stats::filter(allele_frequency_error, rep(1/2,2))
    p <- add_trace(p, x=seq(1:nSamples), y=allele_frequency_error_average, type="scatter", mode="lines",evaluate = TRUE)
  }

  p <- p %>% layout(title = 'Allele Frequency Error across Gibbs (Full) Sweeps',
                        xaxis = list(title = 'nSamples (Full Sweeps)'),
                        yaxis = list(title = 'Abs Difference to True AF (Rolling Avg)'))
  return(p)
}

r2_true_inferred_af_across_chains = function(inferred_allele_freq, #Matrix across chains
                                             true_allele_freq,
                                             graphic = FALSE,
                                             nhaps,
                                             nSamples){    #Vector of AFs
  nSamples = nSamples
  r2       = c()
  for (i in (1:((nSamples-1)*nhaps+1))) {
    r2[i] = summary(lm(inferred_allele_freq[,i]~true_allele_freq))$r.squared
  }
  if (graphic == FALSE) {
    return(r2)
  }
  else{
  #Create plot
  fig <- plot_ly(data = iris, x = ~seq(1:nSamples), y = ~r2)
  fig <- fig %>% layout(title = 'Correlation of True to Inferred Allele Frequency Across Sweeps',
                    xaxis = list(title = 'nSamples (Full Sweeps)'),
                    yaxis = list(title = 'r2 correlation (Inferred~True)'))
  return(fig)
  }
}
