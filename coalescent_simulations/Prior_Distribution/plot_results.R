#Plot results as a grid across the results
library(plotly)

#Read in results file (across all divergence times)
AF      = list()
LD      = list()

#Read in all the pooled quantile data (both the AFs and the LD)
paired_values = c("0_0.0012",    "0_0.0052",   "0_0.0211",   "0_0.0506",
                  "10_0.0012",   "10_0.0052",  "10_0.0211",  "10_0.0506",
                  "50_0.0012",   "50_0.0052",  "50_0.0211",  "50_0.0506",
                  "125_0.0012",  "125_0.0052", "125_0.0211", "125_0.0506")


par(mfrow = c(4, 4))
#Plot AFs
for (i in 1:length(paired_values)) {
  file = sprintf("~/Desktop/results_split/%s_split_AF_pooled_quantile_counts.RData",paired_values[i])
  print(file)
  AF[[i]] = readRDS(file)
  #Grab the data
  af_data  = c((AF[[i]]))
  #Plot QQ plot compared to uniform
  qqplot(randu$x*100,unlist(af_data),
         xlab = sprintf("", strsplit(paired_values, "_")[[i]][1]),
         ylab = sprintf("", strsplit(paired_values, "_")[[i]][2]))
  abline(0,1)
}

par(mfrow = c(4, 4))
#Plot AFs
for (i in 1:length(paired_values)) {
  #Read in file
  file = sprintf("~/Desktop/results_split/%s_split_LD_pooled_quantile_counts.RData",paired_values[i])
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
