#Scratch
library(data.table)
marcus_sumstats = fread("~/Desktop/imputed_hypertension_sumstats", header = T)
#plot(-log10(marcus_sumstats$P_BOLT_LMM_INF))
#neale_sumstats = read.table("~/Downloads/disease_HYPERTENSION_DIAGNOSED.txt", header = T)
merged = merge(x = neale_sumstats,y= marcus_sumstats,by = "SNP")
plot(-log10(merged$P_BOLT_LMM.x),-log10(merged$P_BOLT_LMM.y))


