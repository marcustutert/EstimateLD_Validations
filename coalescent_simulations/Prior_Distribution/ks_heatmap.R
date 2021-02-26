#KS-Test
library(stats)
split_times = c("0","10","50","125")
#Calculate the KS metrics across the pooled LD and pooled AFs
#Read results into this list
AF         = list()
LD         = list()
KS_results = c()
for (i in 1:length(split_times)) {
  #Read in RData (LD & AFs)
  LD[[i]]    = readRDS(sprintf("results/split_%s_LD_pooled_quantile_counts.RData",split_times[i]))
  AF[[i]]    = readRDS(sprintf("results/split_%s_AF_pooled_quantile_counts.RData",split_times[i]))
  #Run Ks test
  KS_results = ks.test(c(unlist(AF[[i]])),"punif",min = 1, max = 100)
  print(KS_results)
}

#Now we calculate the KS-test

