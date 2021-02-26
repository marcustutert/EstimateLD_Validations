#Plot OOA results
#Plot results as a grid across the results
library(plotly)

#Read in results file (for both the LD and the Allele Frequencies)
AF      = list()
LD      = list()

#We will split up the analysis into GWAS/REF and REF/GWAS --> Talk to Luke about this!
#Fst_values = c(fst1,fst2,fst3) #Sprintf this in later
#Read in all the pooled quantile data (both the AFs and the LD)
paired_values = c("CHB_GWAS_CEU_Ref_0.15",  "YRI_GWAS_CEU_Ref_0.15", "YRI_GWAS_CHB_Ref_0.15",
                  "CHB_GWAS_CEU_Ref_0.25",  "YRI_GWAS_CEU_Ref_0.25", "YRI_GWAS_CHB_Ref_0.25",
                  "CHB_GWAS_CEU_Ref_0.35",  "YRI_GWAS_CEU_Ref_0.35", "YRI_GWAS_CHB_Ref_0.35")

par(mfrow = c(3, 3))
#Plot AFs
for (i in 1:length(paired_values)) {
  file = sprintf("~/Desktop/%s_AF_pooled_quantile_counts.RData",paired_values[i])
  print(file)
  AF[[i]] = readRDS(file)
  #Grab the data
  af_data  = c((AF[[i]]))
  #Plot QQ plot compared to uniform
  qqplot(unlist(af_data),randu$x*100,
         xlab = "",
         ylab = "")
  abline(0,1)
}

par(mfrow = c(3, 3))
#Plot AFs
for (i in 1:length(paired_values)) {
  file = sprintf("~/Desktop/%s_LD_pooled_quantile_counts.RData",paired_values[i])
  print(file)
  #Grab data and clean
  LD[[i]] = readRDS(file)
  ld_data  = unlist(LD[[i]])
  ld_data  = ld_data[!is.na(ld_data)]
  #Plot QQ plot compared to uniform
  qqplot(unlist(ld_data),randu$x*100,
         xlab = "",
         ylab = "")
  abline(0,1)
}
