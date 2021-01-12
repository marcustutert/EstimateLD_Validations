  # Helper Functions ----
  match_rsids = function(ref_hap_panel,
                         sumstats,
                         AF_filter = c(0.01,0.99),
                         impute = FALSE){

  #Read in reference panel and tranpose it
  ref_hap_panel = t(read.table(ref_hap_panel,header = F))
  if (impute == TRUE) {
     #Replace colnames in BOLT imputed file name
     #This is done because in the UKBB imputed AF, the colname is 'ID' NOT 'SNP'
     colnames(sumstats)[2] = "SNP"
  }

  #Extract the rsids from the .haps file
  KG_rsid                 = ref_hap_panel[2,]

  #Find SNPs in UKBB sumstats that are not in 1000G
  KG_rsid_missing_indexes = which(!(sumstats$SNP %in% KG_rsid))
  sumstats                = sumstats[-c(KG_rsid_missing_indexes),]

  #Add column  &row names
  ref_hap_panel           = ref_hap_panel[6:length(ref_hap_panel[,1]),]
  colnames(ref_hap_panel) = KG_rsid
  nrow                    = dim(ref_hap_panel)[1]
  ncol                    = dim(ref_hap_panel)[2]
  haps_names              = rep("haps",nrow)
  ref_hap_panel           = matrix(as.numeric(as.matrix(ref_hap_panel)),nrow = nrow, ncol = ncol,dimnames = list(haps_names,colnames(ref_hap_panel)))

  #Filter based on AF in UKBB/1000G
  KG_rsid_fixed           = which(colMeans(ref_hap_panel) > AF_filter[2] | colMeans(ref_hap_panel) < AF_filter[1])
  ref_hap_panel           = ref_hap_panel[,-c(KG_rsid_fixed)]
  rsid_difference         = which(!sumstats$SNP %in%  colnames(ref_hap_panel))
  sumstats                = sumstats[-c(rsid_difference),]

  if (impute == TRUE) {
    #sumstats$SE                = (sumstats$SE/0.17)^2 #DONT HARD CODE THIS IN
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
                               noise,
                               case_control_constant = 122848){

    #Normalize the weight matrix
    weight_matrix_sum           = (colSums(weight_matrix))
    weight_matrix_normalized    = weight_matrix/weight_matrix_sum
    allele_freqs                = colSums(ref_hap_panel*weight_matrix_normalized)
    SE                          = 1/(2*(case_control_constant*(allele_freqs*(1-allele_freqs))))
    ratio                       = gwas_observed/SE
    ll                          = sum(dgamma(x = ratio, shape = noise, rate = noise,log = T))
    return(ll)
  }

  #Perform AF imputation: Note only works with constants weights across genome atm
  AF_Imputation = function(ref_hap_panel_string)
  {

    #Read in reference panel based on the population
    ref_hap_panel       = read.table(sprintf("./analysis/genotyped_reference_panel/%s_panel",ref_hap_panel_string),header = T)
    #Find range from the first SNP in our reference panel to our last (this will be what weeights we span)
    first_snp_genotyped = colnames(ref_hap_panel)[1]
    last_snp_genotyped  = colnames(ref_hap_panel)[ncol(ref_hap_panel)]

    #Use PLINK2 to compute the AF from across entire CHR in the UKBB BGENs
    #system(sprintf("/well/mcvean/mtutert/software/plink2 --max-alleles 2 --freq --from %s --to %s --bgen /well/mcvean/ukbb12788/mtutert/impute_snp_qc/ukbb_imputed_qc_wba_chr1.bgen ref-first  --out ./temp/impute_freqs",first_snp_genotyped,last_snp_genotyped))

    #Extract the rsids from the AF
    AF      = fread("./analysis/temp_files/UKBB_chr1_allele_freqs.afreq", header = T, stringsAsFactors = F)
    AF      = AF[which(AF[,2] == first_snp_genotyped):which(AF[,2] == last_snp_genotyped),] #rsids from imputed data (first to last SNP)
    AF_rsid = AF[,2]
    write.table(AF_rsid,"./analysis/temp_files/AF_rsid_UKBB_chr1_imputed_rsid", quote = F, row.names = F, col.names = F)

    #Match up the SNPs
    system(sprintf("/well/mcvean/mtutert/software/plink2 --max-alleles 2 --pfile /well/mcvean/mtutert/1000_Genomes/1000G_chr1 --keep ./analysis/pop_sample_tables/%s_pop_samples --extract ./analysis/temp_files/AF_rsid_UKBB_chr1_imputed_rsid --export haps --out ./analysis/imputed_reference_panel/imputed_%s_panel --silent",ref_hap_panel_string,ref_hap_panel_string))

    #Match the rsids in the imputed reference panel, with the values in the sumstats and according to our AF filter
    results               = match_rsids(ref_hap_panel = sprintf("./analysis/imputed_reference_panel/imputed_%s_panel.haps",ref_hap_panel_string),
                                        sumstats = AF,
                                        impute = TRUE)

    ref_hap_panel_imputed = results[[2]]
    AF                    = results[[1]]
    #Remove duplicates (only happens with imputed data(?))
    AF = AF[!duplicated(AF[,2]), ]

    #Write results to imputed reference panel directory
    write.table(ref_hap_panel_imputed, file = sprintf("./analysis/imputed_reference_panel/imputed_reference_panel_matched_%s", ref_hap_panel_string), quote = F, col.names = F, row.names = F)
    write.table(AF[,2], file = sprintf("./analysis/imputed_reference_panel/rsid_imputed_matched_%s", ref_hap_panel_string), quote = F, col.names = F, row.names = F)
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
                              gwas_observed,
                              noise){

      #Create weight matrix that uses only the base reference panel
      weight_matrix                  = matrix(data = 1, nrow = nrow(ref_hap_panel), ncol = ncol(ref_hap_panel))
      pop_weight_ll                  = likelihood_output(weight_matrix = weight_matrix,
                                                         ref_hap_panel = ref_hap_panel,
                                                         gwas_observed = gwas_observed,
                                                         noise)
      return(pop_weight_ll)
  }

  #Function that takes in inference results and plots the results against the reference panel chosen
  plot_genotyped_markers_vs_reference = function(allele_freq_results,
                                                 Gibbs_Array,
                                                 ref_hap_panel,
                                                 sumstats,
                                                 population){

    #Read in data from results file
    total_updates         = ncol(allele_freq_results)
    allele_freq_mean      = rowMeans(allele_freq_results[,(.9*total_updates):total_updates])

    #Calculate correlation to true AF (from BOLT)
    true_inferred_r2      = (summary(lm(allele_freq_mean~sumstats$A1FREQ))$r.squared)
    true_reference_r2     = (summary(lm(colMeans(ref_hap_panel)~sumstats$A1FREQ))$r.squared)

    #Plot these results
    fig_genotype = plot_ly()
    fig_genotype <- plot_ly(data = iris)
    fig_genotype <- plot_ly(x = ~sumstats$A1FREQ, y = ~colMeans(ref_hap_panel),     name = sprintf("%s Reference Panel (r2 %.3f)", population,true_reference_r2), mode = 'markers',type='scatter')
    fig_genotype <- fig_genotype %>% add_trace(fig_genotype, y = ~allele_freq_mean, name = sprintf("Inferred %s Panel (r2 %.3f)", population,true_inferred_r2) ,mode = 'markers',type='scatter',evaluate = TRUE)
    fig_genotype <- plotly_build(fig_genotype)

    #Layout of Plot
    fig_genotype  <- layout(fig_genotype,
                            title  = sprintf("Genotyped Markers Comparison to %s Panel",population),
                            xaxis  = list(title = "True Allele Frequency (BOLT)"),
                            yaxis  = list(title = "Comparison Allele Frequency"))

    return(fig_genotype)
  }

  plot_imputed_markers_vs_reference = function(population,
                                               HW_Array){

    #Perform the AF imputation and matching
    AF_matching           = AF_Imputation(population)
    AF                    = unlist(AF_matching[[2]][,5])
    ref_hap_panel_imputed = AF_matching[[1]]

    #Take the weights that we have calculated (For now it'll just be the last sample-->Deal with burn in later?)
    nhaps                 = nrow(ref_hap_panel_imputed)
    nsnps_imputed         = ncol(ref_hap_panel_imputed)

    #Filter down to per sample basis (prev: per haplotype basis) to do the imputed weights
    HW_Per_Sample_Array   = HW_Array[,,dim(HW_Array)[3]]
    single_column_weights = HW_Per_Sample_Array[,1] #First column of weights
    imputed_weights       = matrix(rep(single_column_weights,nsnps_imputed),nrow = length(single_column_weights),ncol = nsnps_imputed,byrow = FALSE)
    inferred_imputed_af   = colSums(imputed_weights * ref_hap_panel_imputed)
    #Get correlation of true/inferred AF across imputed data
    inferred_r2_imputed  = (summary(lm(inferred_imputed_af~AF))$r.squared)
    reference_r2_imputed = (summary(lm(colMeans(ref_hap_panel_imputed)~AF))$r.squared)

    fig_impute = plot_ly()
    fig_impute <- plot_ly(data = iris)
    fig_impute <- plot_ly(x = ~AF, y = ~colMeans(ref_hap_panel_imputed), name = sprintf("%s Reference Panel (r2 %.3f)", population, reference_r2_imputed), mode = 'markers',type='scatter')
    fig_impute <- fig_impute %>% add_trace(fig_impute, y = ~inferred_imputed_af, name = sprintf("Inferred %s Panel (r2 %.3f)", population,inferred_r2_imputed) ,mode = 'markers',type='scatter',evaluate = TRUE)
    fig_impute <- plotly_build(fig_impute)

    #Layout of Plot
    fig_impute  <- layout(fig_impute,
                          title  = sprintf("Imputed Markers Comparison to %s Panel",population),
                          xaxis  = list(title = "True Allele Frequency (BOLT)"),
                          yaxis  = list(title = "Comparison Allele Frequency"))

    return(fig_impute)
  }

  #Function that takes in the inference results and plots: the LL curves, the reference panel line and the zoomed in inset graph portion
  plot_ll_curves = function(ll,
                            population,
                            ref_hap_panel,
                            gwas_ss,
                            noise){

    # Font Options
    t         = list(size = 10)
    ll_graph  = plot_ly()

    #Add in first graph
    ll_graph  = add_trace(ll_graph,
                          x =~(1:length(ll)),
                          y = ~ll,
                          name = sprintf("%s Inferred Weights",population))

    #Add max-likelihood line (with base reference panel--no adjustment to weights)
    ref_hap_panel   = as.matrix(read.table(sprintf("./analysis/genotyped_reference_panel/%s_panel",population),header = T))
    sumstats        = read.table(sprintf("./analysis/genotyped_reference_panel/sumstats_matched_to_%s_panel",population), header = T)
    likelihood_line = ll_per_pop_panel(ref_hap_panel = ref_hap_panel,
                                       gwas_observed = gwas_ss,
                                       noise)
    #Add in reference panel line
    ll_graph        = ll_graph  %>% add_segments(ll_graph,
                                                 x    = 0,
                                                 xend = length(ll),
                                                 y    = likelihood_line,
                                                 yend = likelihood_line,
                                                 line = list(dash = "dash"),
                                                 name = sprintf("%s Reference Panel",population))

    #Add in max likelihood line
    max_ll          = sum(dgamma(rep(1,ncol(ref_hap_panel)),shape = noise, rate = noise))
    ll_graph        = ll_graph  %>% add_segments(ll_graph,
                                                 x    = 0,
                                                 xend = length(ll),
                                                 y    = max_ll,
                                                 yend = max_ll,
                                                 line = list(dash = "dash"),
                                                 name = "Max Likelihood (AFs matched)")


    #Add in post burn-in (zoomed in graph)
    # ll_graph  <- add_trace(ll_graph,
    #                        x          = ~seq(from = .8*length(ll), to = length(ll)),
    #                        y          = ~ll[(.8*length(ll)):length(ll)],
    #                        xaxis      = 'x2',
    #                        yaxis      = 'y2',
    #                        mode       = 'markers',
    #                        showlegend = F,
    #                        marker     = list(color = "#1f77b4"),
    #                        title      = "Post Burn In Convergence")

    #Layout of Plot
    ll_graph  <- layout(ll_graph,
                        title  = sprintf("Log Likelihood Convergence (%s)",population),
                        xaxis  = list(title = "Update Number (Per Haplotype Sweep)"),
                        yaxis  = list(title = "Log Likelihood"),
                        #xaxis2 = list(title = "Update Number", domain = c(0.7, 1.0), anchor='y2'),
                        #yaxis2 = list(title = "Log Likelihood", domain = c(0.2, 0.4), anchor='x2'),
                        font = t,
                        margin = list(b = 100, l = 100))

    #Add in design of showing zoomed in box plot
    # ll_graph  = layout(ll_graph,
    #                    shapes = list(list(type = "rect",
    #                                       line = list(color = "black"),
    #                                       x0   = .9*length(ll),
    #                                       x1   = length(ll),
    #                                       xref = "x",
    #                                       y0   = ll[length(ll)]+ll[length(ll)]*1.3,
    #                                       y1   = ll[length(ll)]-ll[length(ll)]*1.3,
    #                                       yref = "y")))
    return(ll_graph)
  }

  #Function that will calculate LD across 1000G populations and compare to our inferred method and true underlying LD
  plot_ld_genotyped_markers = function(population,
                                       HW_Array,
                                       Inferred_LD){

    #SNPs should be already matched up from previous steps (when performing inference) so we can just read in the rsids that we need
    sumstats       = read.table(sprintf("./analysis/genotyped_reference_panel/sumstats_matched_to_%s_panel", population), header = T)
    genotyped_rsid = sumstats$SNP
    #Write out this file
    write.table(x = genotyped_rsid, file = sprintf("./analysis/ld_results/genotyped_%s_rsid", population), quote = F, row.names = F, col.names = F)

    #Calculate the LD at the genotyped sites using both the reference panel (1000G) and with the individual level data (plink bfiles from imputed UKBB)
    #Generate individual level LD (from the individual level GWAS data)
    system(sprintf("/apps/well/plink/1.90b3/plink --bfile /well/mcvean/ukbb12788/mtutert/genotyped_qc_wba/ukbb_genotype_qc_wba_chr1 --r square --out ./analysis/ld_results/genotyped_ld_indiv_ukbb_data_matched_to_%s --extract ./analysis/ld_results/genotyped_%s_rsid --keep-allele-order --silent", population, population))

    #Get the LD across 1000G population specified in the function and across specific rsid's
    system(sprintf("/apps/well/plink/1.90b3/plink  --bfile /well/mcvean/mtutert/1000_Genomes/plink/1000G_chr1 --r square --out ./analysis/ld_results/genotyped_%s_ld_reference --extract ./analysis/ld_results/genotyped_%s_rsid --keep-allele-order --keep-fam ./analysis/pop_sample_tables/%s_pop_samples", population, population, population))

    #Use the LD we inferred at the LAST sample from the
    Inferred_LD = Inferred_LD[,,dim(Inferred_LD)[3]]

    #Read in the LD (reference, true & inferred) and downsample to plot points
    Reference_LD  = as.matrix(read.table(sprintf("./analysis/ld_results/genotyped_%s_ld_reference.ld",population),header = F))
    True_LD       = as.matrix(read.table(sprintf("./analysis/ld_results/genotyped_ld_indiv_ukbb_data_matched_to_%s.ld",population),header = F))

    #Plot the resulting graphs
    #Down-sample the LD (only if points are >>>>)
    if (length(c(True_LD))>1e4) {
      index             = seq(1:length(c(True_LD)))
      downsampled_index = sample(x = index, size = 1e4, replace = F)
    }
    else{
      downsampled_index = seq(1:length(c(True_LD)))
    }

    #Get r2 correlations
    reference_true_r2 = summary(lm(c(Reference_LD)~c(True_LD)))$r.squared
    imputed_true_r2   = summary(lm(c(Inferred_LD)~c(True_LD)))$r.squared

    #Plot results
    fig_genotyped_LD <- plot_ly(data = iris, x = ~ c(True_LD)[downsampled_index])
    fig_genotyped_LD <- fig_genotyped_LD %>% add_trace(y = ~ c(Reference_LD)[downsampled_index], name = sprintf('Reference LD (r2 %.3f)', reference_true_r2), mode = 'markers')
    fig_genotyped_LD <- fig_genotyped_LD %>% add_trace(y = ~c(Inferred_LD)[downsampled_index], name = sprintf('Inferred LD (r2 %.3f)', imputed_true_r2), mode = 'markers')

    fig_genotyped_LD  <- layout(fig_genotyped_LD,
                                title  = sprintf("Genotyped (Training Model) SNPs Pairwise LD Comparison to %s Panel",population),
                                xaxis  = list(title = "True LD (LD r)"),
                                yaxis  = list(title = "Comparison LD (LD r)"))
    return(fig_genotyped_LD)
  }

  #Function that will calculate LD across the 1000G populations and compare that to our inferred method and true underlying LD
  #This will be done at IMPUTED SNPs (validation set)
  plot_ld_imputed_markers = function(population,
                                     HW_Array){

    #First match up the rsids between 1000G panel and the Imputed SNPs used in BOLT
    AF_matching           = AF_Imputation(population)
    rsid                  = unlist(AF_matching[[2]][,2])
    write.table(rsid,sprintf("./analysis/ld_results/rsid_imputed_%s",population),quote = F,row.names = F,col.names = F)
    ref_hap_panel_imputed = AF_matching[[1]] #This is the imputed reference panel (1000G panel with all the imputed (matched) SNPs)

    #Take the weights that we have calculated (For now it'll just be the last sample-->Deal with burn in later?)
    nhaps                 = nrow(ref_hap_panel_imputed)
    nsnps_imputed         = ncol(ref_hap_panel_imputed)

    #Get the post burnin LD
    HW_Per_Sample_Array   = HW_Array[,,dim(HW_Array)[3]]
    single_column_weights = HW_Per_Sample_Array[,1] #Extract first column of weights

    imputed_weights       = matrix(rep(single_column_weights,nsnps_imputed),nrow = length(single_column_weights),ncol = nsnps_imputed,byrow = FALSE)
    Qx_Array              = markovian_flow(HW_Matrix    = imputed_weights,
                                           start         = 1,
                                           end           = ncol(imputed_weights),
                                           recombination = F)

    LD_imputed            = LD_flow_expectation(HW_Matrix         = imputed_weights,
                                                ref_allele_matrix = ref_hap_panel_imputed,
                                                Qx_array          = Qx_Array,
                                                recombination     = F)
    #Write out imputed inferred LD results
    write.table(LD_imputed,sprintf("./analysis/ld_results/LD_imputed_from_inference_matched_%s", population),quote = F, row.names = F, col.names = F)

    #Generate individual level LD (from the individual level GWAS data)
    system(sprintf("/apps/well/plink/1.90b3/plink --bfile /well/mcvean/ukbb12788/mtutert/impute_snp_qc/plink/ukbb_imputed_plink_chr1 --r square --out ./analysis/ld_results/ld_indiv_ukbb_data_matched_to_%s --extract ./analysis/ld_results/rsid_imputed_%s --keep-allele-order --silent", population, population))

    #Get the LD across 1000G population specified in the function and across specific rsid's
    system(sprintf("/apps/well/plink/1.90b3/plink --bfile /well/mcvean/mtutert/1000_Genomes/plink/1000G_chr1 --r square --out ./analysis/ld_results/%s_ld_reference --extract ./analysis/ld_results/rsid_imputed_%s --keep-allele-order --keep-fam ./analysis/pop_sample_tables/%s_pop_samples --silent", population, population, population))

    #Read in the LD (reference, true & inferred) and downsample to plot points
    Reference_LD  = as.matrix(read.table(sprintf("./analysis/ld_results/%s_ld_reference.ld",population),header = F))
    True_LD       = as.matrix(read.table(sprintf("./analysis/ld_results/ld_indiv_ukbb_data_matched_to_%s.ld",population),header = F))
    Inferred_LD   = LD_imputed

    #Downsample the LD
    index             = seq(1:length(c(True_LD)))
    downsampled_index = sample(x = index, size = 1e4, replace = F)

    #Get r2 correlations
    reference_true_r2 = summary(lm(c(Reference_LD)~c(True_LD)))$r.squared
    imputed_true_r2   = summary(lm(c(Inferred_LD)~c(True_LD)))$r.squared

    #PLot results
    fig_imputed_LD <- plot_ly(data = iris, x = ~ c(True_LD)[downsampled_index])
    fig_imputed_LD <- fig_imputed_LD %>% add_trace(y = ~ c(Reference_LD)[downsampled_index], name = sprintf('Reference LD (r2 %.3f)',reference_true_r2) ,mode = 'markers')
    fig_imputed_LD <- fig_imputed_LD %>% add_trace(y = ~c(Inferred_LD)[downsampled_index], name = sprintf('Inferred LD (r2 %.3f)',imputed_true_r2) , mode = 'markers')

    #Add layout
    fig_imputed_LD  <- layout(fig_imputed_LD,
                              title  = sprintf("Imputed (Training Model) SNPs Pairwise LD Comparison",population),
                              xaxis  = list(title = "True LD (LD r)"),
                              yaxis  = list(title = "Comparison LD (LD r)"))
    #Build figure
    fig_imputed_LD = plotly_build(fig_imputed_LD)
    return(fig_imputed_LD)
  }

  #Function that will actually perform the sumstat imputation
  sumstat_impute = function(genotyped_sumstats,
                            imputed_sumstats,
                            LD)
  {
    #Remove from the genotyped_sumstats file, SNPs that don't match at all
    remove_snps_index  = which(!genotyped_sumstats$SNP %in% imputed_sumstats$SNP )
    genotyped_sumstats = genotyped_sumstats[-c(remove_snps_index),]
    #First we need to find which SNPs, based on their indexes, were typed or imputed
    typed_snps_index   = which(imputed_sumstats$SNP %in% genotyped_sumstats$SNP)

    untyped_snps_index = seq(1:length(imputed_sumstats$SNP))
    LD_typed_untyped   = LD[typed_snps_index,untyped_snps_index] #LD of typed - untyped pairs
    inv_LD_typed       = solve( LD[typed_snps_index,typed_snps_index] + 0.001*diag(length(typed_snps_index))) #inverse of LD of typed SNPs
    W                  = LD[untyped_snps_index, typed_snps_index] %*% inv_LD_typed #these are the weights that turn typed to imputed
    infos              = as.numeric(rowSums(W * LD[untyped_snps_index,typed_snps_index])) #info measures per each untyped
    z.imp              = (W %*% (genotyped_sumstats$BETA/genotyped_sumstats$SE^0.5))/sqrt(infos) #use scaling 1/sqrt(infos) to get var=1 for z-scores
    true_z             = -imputed_sumstats$BETA/(imputed_sumstats$SE/0.17) #Flip beta alleles
    return(list(z.imp,true_z))
  }

  #Plotting the resulting sumstat imputation comparison
  plot_sumstat_imputation_comparison = function(population,
                                                Inferred_LD)
  {
    #Read in the imputed LD (which should be of the same length as imputed zscores)
    imputed_reference_panel                    = read.table(sprintf("./analysis/imputed_reference_panel/imputed_reference_panel_matched_%s", population), header = F)
    #Read in the matched reference panel GENOTYPED UKBB sumstats
    genotyped_matched_population_ukbb_sumstats = fread(sprintf("./analysis/genotyped_reference_panel/sumstats_matched_to_%s_panel", population), header = T)
    genotyped_zscores                          = genotyped_matched_population_ukbb_sumstats$BETA/genotyped_matched_population_ukbb_sumstats$SE^2
    #Read in the imputed z scores
    imputed_matched_population_ukbb_sumstats   = fread("imputed_test_chr1", stringsAsFactors = F, header = T)

    #Calculate the LD in the reference panel & BGEN files (ie the imputed data)
    system(sprintf("/apps/well/plink/1.90b3/plink --bfile /well/mcvean/ukbb12788/mtutert/impute_snp_qc/plink/ukbb_imputed_plink_chr1 --r square --out ./analysis/results/ld_results/ld_indiv_ukbb_data_matched_to_%s --extract ./analysis/imputed_reference_panel/rsid_imputed_matched_%s --keep-allele-order --silent", population, population))
    system(sprintf("/apps/well/plink/1.90b3/plink --bfile /well/mcvean/mtutert/1000_Genomes/plink/1000G_chr1 --r square --out ./analysis/results/ld_results/%s_ld_reference --extract ./analysis/imputed_reference_panel/rsid_imputed_matched_%s --keep-allele-order --keep-fam ./analysis/pop_sample_tables/%s_pop_samples --silent", population, population, population))

    #Read in the LD (including the imputed LD)
    Inferred_LD  = as.matrix(read.table(sprintf("./analysis/ld_results/LD_imputed_from_inference_matched_%s", population), header = F))
    True_LD      = as.matrix(read.table(sprintf("./analysis/ld_results/ld_indiv_ukbb_data_matched_to_%s.ld", population), header = F))
    Reference_LD = as.matrix(read.table(sprintf("./analysis/ld_results/%s_ld_reference.ld", population), header = F))

    #Read in the imputed matched zscores
    rsid         = read.table(sprintf("./analysis/ld_results/rsid_imputed_%s", population), header = F, stringsAsFactors = F)
    rsid         = rsid$V1
    imputed_matched_population_ukbb_sumstats = imputed_matched_population_ukbb_sumstats[which(imputed_matched_population_ukbb_sumstats$SNP %in% rsid),]
    imputed_matched_population_ukbb_sumstats = imputed_matched_population_ukbb_sumstats %>% distinct(SNP, .keep_all = TRUE)  #First we will get the sumstats from the reference LD

    #Perform the SumStatimp
    sumstat_results = sumstat_impute(genotyped_sumstats = genotyped_matched_population_ukbb_sumstats,
                                     imputed_sumstats   = imputed_matched_population_ukbb_sumstats,
                                     LD                 = True_LD)
    imputed_z_indiv_ld       = sumstat_results[[1]]
    true_z                   = sumstat_results[[2]]
    indiv_ld_r2              = summary(lm(imputed_z_indiv_ld~true_z))$r.squared

    sumstat_results = sumstat_impute(genotyped_sumstats = genotyped_matched_population_ukbb_sumstats,
                                     imputed_sumstats   = imputed_matched_population_ukbb_sumstats,
                                     LD                 = Inferred_LD)
    imputed_z_inferred       = sumstat_results[[1]]
    true_z                   = sumstat_results[[2]]
    inferred_ld_r2           = summary(lm(imputed_z_inferred~true_z))$r.squared

    sumstat_results = sumstat_impute(genotyped_sumstats = genotyped_matched_population_ukbb_sumstats,
                                     imputed_sumstats   = imputed_matched_population_ukbb_sumstats,
                                     LD                 = Reference_LD)
    imputed_z_reference      = sumstat_results[[1]]
    true_z                   = sumstat_results[[2]]
    reference_ld_r2          = summary(lm(imputed_z_reference~true_z))$r.squared
    print(length(true_z))

    #Do the plotting!
    #PLot results
    fig_ss_imputed_LD <- plot_ly(data = iris, x = ~c(true_z))
    fig_ss_imputed_LD <- fig_ss_imputed_LD %>% add_trace(y = ~c(imputed_z_indiv_ld), name = sprintf('Indiv LD (r2 %.3f)',indiv_ld_r2) ,mode = 'markers')
    fig_ss_imputed_LD <- fig_ss_imputed_LD %>% add_trace(y = ~c(imputed_z_inferred), name = sprintf('Inferred LD (r2 %.3f)',inferred_ld_r2) , mode = 'markers')
    fig_ss_imputed_LD <- fig_ss_imputed_LD %>% add_trace(y = ~c(imputed_z_reference), name = sprintf('Reference LD (r2 %.3f)',reference_ld_r2) , mode = 'markers')

    #Add layout
    fig_ss_imputed_LD  <- layout(fig_ss_imputed_LD,
                              title  = sprintf("GWAS Summary Statistic Imputation Comparison",population),
                              xaxis  = list(title = "True Zscores (BOLT)"),
                              yaxis  = list(title = "Comparison Zscores (Imputed)"))

    return((fig_ss_imputed_LD))
  }


  # Analysis ----
  ukbb_analysis_inference = function(name_of_sumstats_file         = "test_genotyped_chr1",
                                     name_of_imputed_sumstats      = "imputed_test_chr1",
                                     BOLT                          = TRUE,
                                     quanitative_scale_convert     = 0.22,
                                     imputation                    = TRUE,
                                     model_snps                    = 100)
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
    library(corpcor)
    library(dplyr)

    #Filter model SNPs ----
    #Read in table of summary statistics
    ukbb_sumstats   = fread(sprintf("/well/mcvean/ukbb12788/mtutert/%s",name_of_sumstats_file), header = T)

    #Trivial Case with just X number of model_snps chosen
    if (length(model_snps) == 1) {
      nsnps           = model_snps
      ukbb_sumstats   = ukbb_sumstats[1:nsnps,]
    }

    #This case is if we want to extract a specific region in this chromosome
    if (grepl(":", model_snps, fixed = TRUE)) {
      start_position = model_snps
      end_postion    = model_snps
    }

    #This case is if we have a specific list of rsids we are interested in
    #Rsids will be a series of comma seperated values
    if (length(model_snps) > 1) {
      print("TO DO")
    }

    if (BOLT == TRUE) {

      #Prepare data for inference ----

      #Create list of all 1000G Populations in a list (this is what we can expand)
      KG_super_pop_names = list(c("EUR","AFR"))

      #Write out list of samples split by super population in 1000G
      #KG_samples_split_by_pop(pop_list = KG_super_pop_names)
      gwas_rsid = ukbb_sumstats$SNP #SNPs in the GWAS
      write.table(gwas_rsid,"./analysis/temp_files/gwas_rsid", quote = F, row.names = F, col.names = F)

      fig_imputed_markers_accuracy            = list()
      fig_genotyped_markers_accuracy          = list()
      fig_ld_imputed_accuracy                 = list()
      fig_ld_genotyped_accuracy               = list()
      fig_sumstat_imputation_accuracy         = list()

      #Loop through ref panels ----

      for (i in 1:length(KG_super_pop_names)) {

        #Diverse Reference Panels ----
        if (length(KG_super_pop_names[[i]])>1) {
          #Split up panel into component subpopulations
          pop_A            = KG_super_pop_names[[i]][1]
          pop_B            = KG_super_pop_names[[i]][2]
          pop_A_samples    = read.table(sprintf("./analysis/pop_sample_tables/%s_pop_samples", pop_A), header = F)
          pop_B_samples    = read.table(sprintf("./analysis/pop_sample_tables/%s_pop_samples", pop_B), header = F)
          diverse_samples  = rbind(pop_A_samples,pop_B_samples)
          diverse_name     = paste(sprintf("%s",pop_A), sep = "_and_",sprintf("%s",pop_B))
          write.table(file = sprintf("./analysis/pop_sample_tables/%s_pop_samples",diverse_name), x = diverse_samples, quote = F, row.names = F, col.names = F)
          #Replace name in KG_super_pop_names with the name pasted name (pop_A & pop_B)
          KG_super_pop_names[[i]] = diverse_name
          }
        system(sprintf("/well/mcvean/mtutert/software/plink2 --max-alleles 2 --pfile /well/mcvean/mtutert/1000_Genomes/1000G_chr1 --keep ./analysis/pop_sample_tables/%s_pop_samples --extract ./analysis/temp_files/gwas_rsid --export haps --out  ./analysis/genotyped_reference_panel/1000G_%s --silent", KG_super_pop_names[[i]],KG_super_pop_names[[i]]))
        results       = match_rsids(ref_hap_panel = sprintf("./analysis/genotyped_reference_panel/1000G_%s.haps", KG_super_pop_names[[i]]),sumstats = ukbb_sumstats)
        sumstats      = results[[1]]
        ref_hap_panel = results[[2]]
        #Write the matched sumstats/ref_hap_panel
        write.table(x = sumstats,      file = sprintf("./analysis/genotyped_reference_panel/sumstats_matched_to_%s_panel",KG_super_pop_names[[i]]),quote = F,row.names = F, col.names = T)
        write.table(x = ref_hap_panel, file = sprintf("./analysis/genotyped_reference_panel/%s_panel",KG_super_pop_names[[i]]),quote = F, row.names = F, col.names = T)

        #Inference Parameters ----
        nSamples        = 5
        alpha           = 100
        fst             = 0.001

        #Run Genotype Inference (InferLD) ----
        # inference_results_1000G = LD_from_GSHMM(ref_panel_haplotypes   = ref_hap_panel,
        #                                         fst                    = fst,
        #                                         betas                  = FALSE,
        #                                         alpha                  = alpha,
        #                                         nSamples               = nSamples,
        #                                         recomb_rate            = 1e-300,
        #                                         weights_resolution     = 10,
        #                                         likelihood_toggle      = TRUE,
        #                                         se_observed            = sumstats$SE,
        #                                         LD_Infer               = TRUE,
        #                                         genetic_map            = FALSE,
        #                                         chain_likelihood       = TRUE,
        #                                         nChains                = 1,
        #                                         recombination          = FALSE,
        #                                         case_control_constant  = 122848,
        #                                         BurnIn                 = TRUE)
        #saveRDS(inference_results_1000G,file = sprintf("./analysis/inference_results/%s_inference_results",KG_super_pop_names[[i]]))

        #Plot Likelihoods ----
        #For each reference panel plot:
        #The likelihood curve, the max likelihood (matched) and the reference panel likelihood
        results  = readRDS(sprintf("./analysis/inference_results/%s_inference_results", KG_super_pop_names[[i]]))

        ll       = results$log_likelihood

        ll_graph = plot_ll_curves(ll            = ll ,
                                population    = KG_super_pop_names[[i]],
                                ref_hap_panel = ref_hap_panel,
                                gwas_ss       = sumstats$SE,
                                noise         = alpha)

        #Export plot (per population)
        export(ll_graph, file = sprintf("./analysis/ll/%s_ll.jpeg",KG_super_pop_names[[i]]))

        #Genotyped Markers Accuracy----
        fig_genotyped_markers_accuracy[[i]] = plot_genotyped_markers_vs_reference(allele_freq_results = results$inferred_af_given_weights ,
                                                                                  Gibbs_Array         = results$Gibbs_Array,
                                                                                  ref_hap_panel       = ref_hap_panel,
                                                                                  sumstats            = sumstats,
                                                                                  population          = KG_super_pop_names[[i]])

        export(fig_genotyped_markers_accuracy[[i]], file = sprintf("./analysis/genotyped_reference_panel/genotyped_markers_accuracy_%s_panel.jpeg", KG_super_pop_names[[i]]))
        #
        #Imputed Markers Accuracy -----
        fig_imputed_markers_accuracy[[i]]  = plot_imputed_markers_vs_reference(population = KG_super_pop_names[[i]], HW_Array = results$Gibbs_Array)
        export(fig_imputed_markers_accuracy[[i]], file = sprintf("./analysis/imputed_reference_panel/imputed_markers_accuracy_%s_panel.jpeg", KG_super_pop_names[[i]]))

        #Genotyped LD Accuracy ----

        fig_ld_genotyped_accuracy[[i]] = plot_ld_genotyped_markers(population  = KG_super_pop_names[[i]],
                                                                   HW_Array    = results$Gibbs_Array,
                                                                   Inferred_LD = results$LD_Array)
        export(fig_ld_genotyped_accuracy[[i]], file = sprintf("./analysis/ld_results/genotyped_ld_accuracy_%s_panel.jpeg", KG_super_pop_names[[i]]))

        #Imputed LD Accuracy ----
        fig_ld_imputed_accuracy[[i]] = plot_ld_imputed_markers(KG_super_pop_names[[i]],
                                                               results$Gibbs_Array)
        export(fig_ld_imputed_accuracy[[i]], file = sprintf("./analysis/ld_results/imputed_ld_accuracy_%s_panel.jpeg", KG_super_pop_names[[i]]))

        #SumStats Imputation ----
        fig_sumstat_imputation_accuracy[[i]] = plot_sumstat_imputation_comparison(population  = KG_super_pop_names[[i]],
                                                                                  Inferred_LD = results$Gibbs_Array)
        export(fig_sumstat_imputation_accuracy[[i]], file = sprintf("./analysis/sumstat_imputation_results/sumstat_imputation_accuracy_to_%s_Panel.jpeg", KG_super_pop_names[[i]]))

        # Generating Subplots ----

        ######MAKE THIS SHIT SQUARE SOMEHOW

        #genotype_accuracy_plot <- subplot(fig_genotype[[1]], fig_genotype[[2]],fig_genotype[[3]],fig_genotype[[4]],fig_genotype[[5]]) %>% layout(scene = list(aspectration=list(x=1,y=1)))
        #export(genotype_accuracy_plot, file = "./analysis/figures/genotype_accuracy_plot.jpeg")

        #fig_impute_plot <- subplot(fig_impute[[1]], fig_impute[[2]],fig_impute[[3]],fig_impute[[4]],fig_impute[[5]])
        #export(fig_impute_plot, file = "./analysis/figures/imputation_accuracy_plot.jpeg")

        #fig_ld_imputed_subplot <- subplot(fig_impute[[1]], fig_impute[[2]],fig_impute[[3]],fig_impute[[4]],fig_impute[[5]])
        #export(fig_ld_imputed[[1]], file = sprintf("./analysis/results/ld_results/ld_results_%s.jpeg",KG_super_pop_names[[i]]))

        #fig_ld_imputed_subplot <- subplot(fig_impute[[1]], fig_impute[[2]],fig_impute[[3]],fig_impute[[4]],fig_impute[[5]])
        #export(fig_ld_imputed[[1]], file = sprintf("./analysis/results/ld_results/ld_results_%s.jpeg",KG_super_pop_names[[i]])
      }
    }
  }
