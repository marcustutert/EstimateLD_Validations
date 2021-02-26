file = snakemake@params[[1]][1]
nsnps = 10
source("/well/mcvean/mtutert/thesis/coalescent_coverage/helper_functions.R")
library(data.table)
#install.packages("/well/mcvean/mtutert/myPackages/InferLD_0.1.0.tar.gz")
library(InferLD)
#Read in coalescent simulations
gwas = as.matrix(fread(sprintf("./msprime_data/population_split/GWAS_panel_replicate_%s_split_0.csv",file), header = T))
ref  = as.matrix(fread(sprintf("./msprime_data/population_split/GWAS_panel_replicate_%s_split_0.csv",file), header = T))

#Generate the sumstats with glm() in R and perform necessary filtering
res           = msprime_gwas_sumstats(gwas_haplotypes = gwas, reference_haplotypes = ref)
se            = res[[1]]
gwas          = res[[2]]
ref           = res[[3]]
write.table(x =    gwas,
            file      = sprintf("./msprime_data/population_split/matched_panels/GWAS_panel_replicate_%s_split_0.csv", file),
            quote     = F,
            row.names = F,
            col.names = T)
case_control_constant = 2*(nrow(gwas)*0.5*(1-0.5))
#Run the inference
print("Starting Inference")
inference_results_1000G = LD_from_GSHMM(ref_panel_haplotypes   = ref[,1:nsnps],
                                        fst                    = 0.01,
                                        betas                  = FALSE,
                                        alpha                  = 100,
                                        nSamples               = 3,
                                        recomb_rate            = 1e-300,
                                        weights_resolution     = 5,
                                        likelihood_toggle      = FALSE,
                                        se_observed            = se[1:nsnps,3]^2,
                                        LD_Infer               = FALSE,
                                        genetic_map            = FALSE,
                                        chain_likelihood       = TRUE,
                                        nChains                = 1,
                                        recombination          = FALSE,
                                        case_control_constant  = case_control_constant,
                                        BurnIn                 = TRUE)
saveRDS(inference_results_1000G, sprintf("./inference_results/inference_results_rep_%s",file))
