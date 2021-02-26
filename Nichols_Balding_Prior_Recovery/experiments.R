####Recovery of Nichols Balding Prior under different parameter values####
install.packages("~/Desktop/Oxford_Dphil/InferLD_0.1.0.tar.gz")
library(InferLD)
#Simulate reference panel & GWAS under generative model
results = simulate_genetic_drift(Fst                = 1e-1,
                                 nhaps              = 500,
                                 nsnps              = 5,
                                 weights_resolution = 5,
                                 noise              = 100)
ref_hap_panel = results$ref_panel
se_observed   = results$Observed_GWAS_SE
#Inference Parameters ----
nSamples        = 50
alpha           = 100
fst             = 0.01

#Run Genotype Inference (InferLD) ----
inference_results_1000G = LD_from_GSHMM(ref_panel_haplotypes   = ref_hap_panel,
                                        fst                    = 0.01,
                                        betas                  = FALSE,
                                        alpha                  = 100,
                                        nSamples               = nSamples,
                                        recomb_rate            = 1e-300,
                                        weights_resolution     = 5,
                                        likelihood_toggle      = FALSE,
                                        se_observed            = se_observed,
                                        LD_Infer               = FALSE,
                                        genetic_map            = FALSE,
                                        chain_likelihood       = TRUE,
                                        nChains                = 1,
                                        recombination          = FALSE,
                                        case_control_constant  = 1,
                                        BurnIn                 = TRUE)

#Want these densities to compare with Nichols and Balding Model
density1 = density(inference_results_1000G$inferred_af_given_weights[8,]) #Empirical Density

#Nichols Balding Model parametrization
c = 1/Fst-1
alpha = c*0.25
beta  = c*(1-0.25)
density2 = density(rbeta(n = 5e4, shape1 = alpha, shape2 = beta)) #N-B Density

#Plot the Densities of Empirical vs Observed
p <- plot_ly(x = ~density1$x, y = ~density1$y, type = 'scatter', mode = 'lines', name = 'Empirical Density', fill = 'tozeroy') %>%
  add_trace(x = ~density2$x, y = ~density2$y, name = 'Theoretical Density', fill = 'tozeroy') %>%
  layout(xaxis = list(title = 'Pi_Pop'),
         yaxis = list(title = 'Density'))

show(p)
