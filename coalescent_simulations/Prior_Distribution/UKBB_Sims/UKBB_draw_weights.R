#Generating weights from prior
#Split into jobs in order to improve speed of the inference
library(data.table)
library(abind)
#source("match_ukbb.R")
regions              = snakemake@params$regions
chunk                = snakemake@params$chunk
population_pairs     = snakemake@params$paired_values

#Get population (reference)
split_pops = strsplit(population_pairs, "_")
ref_pop    = split_pops[[1]][1]
fst        = split_pops[[1]][2]

#MATCHING PREVIOUSLY DONE IN THE LAST RULE
#Read in the tables we just generated (matched by ref super pop)
ref  = as.matrix(fread(sprintf("prior_draw_analysis_panels/regions/matched_panels/1000G_%s_region_%s_matched_ukbb_ref_panel.haps",ref_pop,regions), header = T))
gwas = as.matrix(fread(sprintf("prior_draw_analysis_panels/regions/matched_panels/ukbb_region_%s_matched_1000G_%s.haps",regions, ref_pop), header = T))

nhaps_ref     = nrow(ref)
nhaps_gwas    = nrow(gwas)
nsnps         = ncol(ref)
effective_fst = as.numeric(fst)


#This filter determines what the skew of that data is
maf_filter = 1/nhaps_gwas
nSamples   = 100

#Results matrix & array initially filled with NAs
AF_Inferred_Results = matrix(data = NA, nrow = nsnps, ncol = nSamples)
LD_Inferred_Results = array(data = NA, dim = c(nsnps,nsnps,nSamples))

#Check if we have any more NAs in the AF result matrix AND the LD Array
while(any(is.na(LD_Inferred_Results))) {
  #Just for sanity purposes check how many more NAs are left to fill
  print(sum(is.na(AF_Inferred_Results)))
  print(sum(is.na(LD_Inferred_Results)))
  #Draw from gamma_quantiled_weights nhaps times
  gamma_draw               = rgamma(n = nhaps_ref, shape =  1/( nhaps_ref * ( effective_fst / (1-effective_fst))), scale = ( nhaps_ref * (effective_fst/(1-effective_fst))))
  #Extend into matrix
  weight_matrix            = matrix(rep(gamma_draw,nsnps), ncol = nsnps)
  #Normalize matrix
  norm_weight_matrix       = weight_matrix/colSums(weight_matrix)[col(weight_matrix)]

  ####### Remove Ascertainment Bias
  #Check which indexes pass the filtering step for AF (ie: rare variants)

  #Check which variants pass this filter
  af_passed_filter_index    = which(colSums(ref*norm_weight_matrix) > maf_filter & colSums(ref*norm_weight_matrix) < (1-maf_filter))

  #Also check to see if any of the variants will need to be removed if we have already populated the full matrix
  filled_variants_af        = which(apply(AF_Inferred_Results, 1, function(x) all(!is.na(x))))

  #If no variants passed check, move to next weight draw
  if (length(af_passed_filter_index) == 0) {
    print("No SNPs pass the AF Filter")
  }

  #Some variants pass MAF filter
  else{

    #Check if the remaining filtered SNPs have snp-snp LD r that is Nan and therefore we need to skip this matrix
    if (any(is.nan((cov.wt(ref,norm_weight_matrix[,1],cor = TRUE, method = "ML")$cor)[af_passed_filter_index,af_passed_filter_index]))) {
      print("LD Matrix with given weights generates NaN")
    }

    else{

      #Check if matrix is full of NAs (this would be when we start the weight draws)
      if (all(is.na(AF_Inferred_Results))) {
        firstNonNA = rep(1,length(af_passed_filter_index))
      }

      else{

        ##################AF RESEULTS#######################

        #Matrices pass both the LD and AF filters, so we can store the values in the AF matrix and weight array accordingly
        #Ask where in the AF results matrix is the first minimum data point, aka the first non NA value
        firstNonNA = c()
        #Loop through all the variants which have passed filter for this weight matrix draw
        for (i in c(1:nsnps)[af_passed_filter_index]) {
          #For these variants where is the first  NA value
          NonNAindex    = which(is.na(AF_Inferred_Results[i,]))
          #Append the minimum value
          firstNonNA    = append(firstNonNA, min(NonNAindex))
        }
      }

      #Replace any infinities with 1
      firstNonNA[!is.finite(firstNonNA)] <- 1

      #Populate AFs in the results matrix
      for (i in 1:length(af_passed_filter_index)) {
        #If we have filled up that row in the matrix then we skip it again
        AF_Inferred_Results[af_passed_filter_index[i],firstNonNA[i]] = (colSums(ref*norm_weight_matrix)[af_passed_filter_index])[i]
      }

      ##################LD RESEULTS#######################

      #Pre-calculate the LD to store later depending on index location
      cov                    = cov.wt(ref,norm_weight_matrix[,1], cor = TRUE, method = "ML")
      LD_Results             = as.array(cov$cor)

      #First we ask which matrix indexes are filled with NA in array
      #Have to loop through the array across all SNPs that have passed the filter
      #This will now be a matrix of NonNA indexes across all SNPs
      firstNonNA = matrix(data = 1, nrow = nsnps, ncol = nsnps )

      for (i in c(1:nsnps)[af_passed_filter_index]) {
        for (j in c(1:nsnps)[af_passed_filter_index]) {
          NonNAindex         = which(is.na(LD_Inferred_Results[i,j,]))
          firstNonNA[i,j]    = min(NonNAindex)
        }
      }

      #Replace any infinities with 1
      firstNonNA[!is.finite(firstNonNA)] <- 1

      #Populate LD in the results matrix
      for (i in 1:length(af_passed_filter_index)) {
        for (j in 1:length(af_passed_filter_index)) {
          LD_Inferred_Results[af_passed_filter_index[i],af_passed_filter_index[j],firstNonNA[af_passed_filter_index[i],af_passed_filter_index[j]]] = LD_Results[af_passed_filter_index[i],af_passed_filter_index[j]]
        }
      }
    }
  }
}

print("here")
#Save object in format Fst_#_Replicate_#_Chunk_#
saveRDS(object = LD_Inferred_Results, file = sprintf("results/UKBB/%s_%s_LD_Region_%s_Chunk_%s.RData",ref_pop,fst,regions,chunk), version = 2)
saveRDS(object = AF_Inferred_Results, file = sprintf("results/UKBB/%s_%s_AF_Region_%s_Chunk_%s.RData",ref_pop,fst,regions,chunk), version = 2)
