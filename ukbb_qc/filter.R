library(data.table)
library(ggplot2)
#Load sample QC data
#Set directory
dir = ("/well/mcvean/ukbb12788/mtutert")
sqc = fread(sprintf("%s/ukb_sqc_v2.txt",dir))

#Check PCs
sqc[, in.PC.outlier.check:=
      het.missing.outliers==0 &
      excluded.from.kinship.inference == 0 &
      excess.relatives == 0 &
      used.in.pca.calculation == 1 &
      in.white.British.ancestry.subset == 1]

#Filter on PCs (7 SDs away from every dimension)
PCs = paste0("PC",1:6)
sds_threshold = 7

distance_sds_squared = function(x, sub){
  ( abs(mean(x[sub])-x) / sd(x[sub]) )^2
}

dists_uni_sdsq = sqc[, lapply(.SD, distance_sds_squared, sub = in.PC.outlier.check==TRUE), .SD=PCs]
dists_multi_sdsq = rowSums(dists_uni_sdsq)
sqc[, PCoutlier :=  dists_multi_sdsq > sds_threshold^2]

#Prepare link between sample QC and phenotype data, add eid's
link_path = sprintf("%s/ukbb_v2.fam",dir)
link = fread(link_path)

nrow(link)==nrow(sqc)
sqc[,eid:=link[,1]]
#Add phenotype data (includes ethnicity)
pheno = fread(sprintf("%s/ukb7749.csv",dir))
included_ethnicities = c(1,1001,1002,1003) #See https://biobank.ctsu.ox.ac.uk/crystal/coding.cgi?id=1001
missing_ethnicities = c(-3,-1,NA)
pheno$eid = as.numeric(as.character(pheno$eid))
sqcpheno = merge(sqc, pheno[,.(eid,`21000-0.0`)], by="eid", all.x=TRUE, all.y=FALSE)

#Filter samples in SQC file
filtered_sqcpheno = sqcpheno[which(sqcpheno$het.missing.outliers==0 &
                                     sqcpheno$excluded.from.kinship.inference == 0 &
                                     sqcpheno$excess.relatives == 0 &
                                     sqcpheno$used.in.pca.calculation == 1 &
                                     sqcpheno$PCoutlier == FALSE &
                                     sqcpheno$`21000-0.0` %in% c(included_ethnicities, missing_ethnicities)),]
write.table(filtered_sqcpheno, sprintf("%s/filtered_sqcpheno",dir),quote = F, row.names = F, col.names = T)
#Read in full sample table, and then subset
ukbb_sample              = read.table(sprintf("%s/ukbb_imp_v3.sample",dir), header = TRUE)
ukbb_sample_neale_filter = ukbb_sample[which(filtered_sqcpheno$eid %in% ukbb_sample$ID_1), ]
write.table(ukbb_sample_neale_filter, sprintf("%s/ukbb_sample_neale_filter",dir),quote = F, row.names = F, col.names = T)
#Convert to a FAM file to do the filtering on UKBB plink files
ukbb_plink_id_neale_filter = ukbb_sample_neale_filter[2:length(ukbb_sample_neale_filter[,1]),]
write.table(ukbb_plink_id_neale_filter,sprintf("%s/ukbb_plink_id_neale_filter",dir),quote = F, row.names = F, col.names = F)
