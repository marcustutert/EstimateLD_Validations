#Pool results
r2 = c()

sumstat_files = list.files(path = "~/Desktop/", pattern = "sumstat_imputed_ref", full.names = T)

#Loop through files and get r2s
for (i in 1:38) {
  data  = readRDS(sumstat_files[i])
  print(sumstat_files[i])
  print(data[[1]][1])
  r2[i] = summary(lm(data[[1]]~data[[2]]))$r.squared
}


