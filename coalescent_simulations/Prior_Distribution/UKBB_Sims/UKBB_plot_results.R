#Plot OOA results
#Plot results as a grid across the results
library(plotly)

#Read in results file (for both the LD and the Allele Frequencies)
AF      = list()
LD      = list()

#We will split up the analysis into GWAS/REF and REF/GWAS --> Talk to Luke about this!
#Fst_values = c(fst1,fst2,fst3) #Sprintf this in later
#Read in all the pooled quantile data (both the AFs and the LD)
paired_values = c("EUR_0.0001", "EUR_0.05",  "EUR_0.1", "EUR_0.25",
                  "AFR_0.0001", "AFR_0.05",  "AFR_0.1", "AFR_0.25",
                  "AMR_0.0001", "AMR_0.05",  "AMR_0.1", "AMR_0.25",
                  "SAS_0.0001", "SAS_0.05",  "SAS_0.1", "SAS_0.25",
                  "EAS_0.0001", "EAS_0.05",  "EAS_0.1", "EAS_0.25",)

par(mfrow = c(4,5))
#Plot AFs
for (i in 1:length(paired_values)) {
  file = sprintf("~/Desktop/%s_AF_pooled_quantile_counts.RData",paired_values[i])
  print(file)
  AF[[i]] = readRDS(file)
  #Grab the data
  af_data  = c((AF[[i]]))
  #Plot QQ plot compared to uniform
  qqplot(randu$x*100,unlist(af_data),
         xlab = sprintf("Fst Prior %s", strsplit(paired_values, "_")[[i]][2]),
         ylab = sprintf("1000G %s Panel", strsplit(paired_values, "_")[[i]][1]))
  abline(0,1)
}

par(mfrow = c(4,5))
#Plot AFs
for (i in 1:length(paired_values)) {
  file = sprintf("~/Desktop/%s_LD_pooled_quantile_counts.RData",paired_values[i])
  print(file)
  #Grab data and clean
  LD[[i]] = readRDS(file)
  ld_data  = unlist(LD[[i]])
  ld_data  = ld_data[!is.na(ld_data)]
  #Plot QQ plot compared to uniform
  qqplot(randu$x*100,unlist(ld_data),
         xlab = sprintf("Fst Prior %s", strsplit(paired_values, "_")[[i]][2]),
         ylab = sprintf("1000G %s Panel", strsplit(paired_values, "_")[[i]][1]))
  abline(0,1)
}
