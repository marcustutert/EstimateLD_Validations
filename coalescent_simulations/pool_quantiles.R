#Data to summarize the results
#Need to wait until all replicates inference has finished before running this script (figure out how to do this in SNAKEMAKE)
pooled_quantiles = c()
file_names = list.files(path = "/well/mcvean/mtutert/thesis/coalescent_coverage/inference_results/", full.names = T)
for (i in 1:length(file_names)) {
  #Load in the resulting RData
  print(file_names[i])
  results = readRDS(file_names[i])
  pooled_quantiles <- append(pooled_quantiles, results)
}
saveRDS(pooled_quantiles, "/well/mcvean/mtutert/thesis/coalescent_coverage/inference_results/pooled_quantiles")
