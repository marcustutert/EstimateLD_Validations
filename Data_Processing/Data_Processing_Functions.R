####Functions that will be used to simulate Data used in the notebooks #####

#This function uses the 1000G data located in /Data and exports the
simulate_genetic_drift = function(Fst,
                                  nhaps,
                                  nsnps){

  #Set number of haplotypes and variants to a small number so we can do testing and debugging
  nhaps = nhaps #Haplotypes in the reference panel
  nsnps = nsnps #Variants in the reference panel

  ## path = dirname(getwd()) #Set path relative to testing folder
  #path <- dirname(getwd())

  #Load in 1000G data from ./Test/Data
  ref_hap_panel     = as.matrix(read.table(sprintf('%s/Data/ref_panel_filtered',path),
                                           header = TRUE))

  #Subset to nhaps and nsnps size to make testing easier
  ref_hap_panel     = ref_hap_panel[ 1:nhaps, 1:nsnps ]

  #Get AF of the variants in the reference pnale
  ref_variants_af   = colMeans(ref_hap_panel)
  Fst               = Fst

  #Assume that allele frequencies are the same in the reference panel and the GWAS panel
  sample_weights           = rgamma(n = nhaps, shape =  1 / ( nhaps * (Fst/(1-Fst))), scale =  1 / ( nhaps * (Fst/(1-Fst))))
  normalize_sample_weights = sample_weights/(sum(sample_weights))
  gwas_variants_af         = colSums(ref_hap_panel * normalize_sample_weights)

  #Filter out SNPs that are above/below a threshhold (0-100% exclusive)
  filter_ref_snp_index_1       = which(colMeans(ref_hap_panel) > 0.97)
  filter_ref_snp_index_2       = which(colMeans(ref_hap_panel) < 0.03)
  filter_gwas_snp_index_1      = which(gwas_variants_af > 0.97)
  filter_gwas_snp_index_2      = which(gwas_variants_af < 0.03)
  filter_snp_index_merged = c(filter_ref_snp_index_1,filter_ref_snp_index_2,filter_gwas_snp_index_1,filter_gwas_snp_index_2)
  if (length(filter_snp_index_merged)>0) {
    ref_hap_panel     = ref_hap_panel[,-filter_snp_index_merged]
    gwas_variants_af  = gwas_variants_af[-filter_snp_index_merged]
    ref_variants_af   = ref_variants_af[-filter_snp_index_merged]
  }
  observed_gwas_se  = 1/(gwas_variants_af*(1-gwas_variants_af))
  return(list(observed_gwas_se,gwas_variants_af,ref_hap_panel))
}

#This function takes in the allele frequencies at each nSample point in the chain and sees how close they are
#To the true allele frequencies

allele_frequency_accuracy_across_chains = function(allele_frequencies, #This should be a matrix (or an array?)
                                                   true_gwas_af){        #When using simulated data)
  nsnps    = length(true_gwas_af)
  nSamples = ncol(allele_frequencies)
  p <- plot_ly()
  ma <- function(x, n = 5){filter(x, rep(1 / n, n), sides = 2)}
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
                                             graphic = FALSE){    #Vector of AFs
  nSamples = ncol(inferred_allele_freq)
  r2       = c()
  for (i in 1:nSamples) {
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
