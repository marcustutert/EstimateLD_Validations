####RUN ALL SCRIPTS FROM #######
##/well/mcvean/ukbb12788/mtutert###

#Function which performs Ben Neale Sample QCing to white british ancestry
sample_qc_generate_clean_wba_id = function(){
  library(data.table)
  #Load sample QC data from UKBB
  sqc = fread("ukb_sqc_v2.txt")

  ### Check PCs ###
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
  ############

  #Prepare link between sample QC and phenotype data, add eid's
  link_path = c("ukbb_v2.fam")
  link = fread(link_path)
  #Check we have the correct number of samples
  nrow(link)==nrow(sqc)
  #Add in the eids via the fam file to the sqc table
  sqc[,eid:=link[,1]]

  #Add phenotype data (includes ethnicity)
  pheno = fread("ukb7749.csv")
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
  write.table(filtered_sqcpheno, "filtered_sqcpheno",quote = F, row.names = F, col.names = T)
  #These are the filtered EIDS^

  #Read in full sample table, and then subset
  ukbb_sample              = read.table("ukbb_imp_v3.sample", header = TRUE)
  ukbb_sample_neale_filter = ukbb_sample[which(filtered_sqcpheno$eid %in% ukbb_sample$ID_1), ]
  write.table(ukbb_sample_neale_filter, "ukbb_sample_neale_filter",quote = F, row.names = F, col.names = T)
  #Convert to a FAM file to do the filtering on UKBB plink files
  ukbb_plink_id_neale_filter = ukbb_sample_neale_filter[2:length(ukbb_sample_neale_filter[,1]),]
  write.table(ukbb_plink_id_neale_filter,"ukbb_plink_id_neale_filter",quote = F, row.names = F, col.names = F)
}


#Function that creates (at least) one FAM file for the phenotypes used downstream in BOLT analysis
generate_fam_for_pheno = function(ICD,ukbb_population){}
#######Linking up phenotypes of interest and the SQC#####
pheno = fread("ukb27864.csv") #Read in CSV
#Find all instances of the I10
pheno_ICD = pheno[,594:806]
cases_hypertension    = as.numeric(pheno$eid[which(apply(pheno_ICD, 1, function(pheno_ICD) any(pheno_ICD %in% c("I10"))))])
controls_hypertension = as.numeric(pheno$eid[!which(apply(pheno_ICD, 1, function(pheno_ICD) any(pheno_ICD %in% c("I10"))))])
#Get eIDS of the patients with hypertension
cases_hypertension_eid              = cbind(cases_hypertension,rep(1,length(cases_hypertension)))
controls_hypertension_eid           = cbind(controls_hypertension,rep(0,length(controls_hypertension)))
ukbb_eids                           = read.table("ukbb_imputed_neale_filter.sample",skip = 1) #Load in clean eIDS, skip mock-up row
fam                                 = read.table("/well/mcvean/ukbb12788/mtutert/genotyped_qc_wba/ukbb_genotype_qc_wba_chr1.fam")
cases_hypertension_eid              = as.data.table(cases_hypertension_eid)
controls_hypertension_eid           = as.data.table(controls_hypertension_eid)
colnames(cases_hypertension_eid)    = c("ID","PHENO")
colnames(controls_hypertension_eid) = c("ID","PHENO")
fam$V6                              = cases_hypertension_eid$PHENO[match(fam$V1,cases_hypertension_eid$ID)]
fam[is.na(fam)] <- 0
write.table(fam, "hypertension.fam",quote = F, row.names = F, col.names = F)
####Filter the V2 Data (Genotyped), such that it only takes BN'd sample QCs in the PLINK files #####
###BOLT-LMM#####
/apps/well/bolt-lmm/2.3.2/bolt --bed ukbb_qc_chr21.bed --bim ukbb_qc_chr21.bim --fam hypertension.fam --phenoUseFam --numThreads 8 --statsFile=test --geneticMapFile /apps/well/bolt-lmm/2.3.2/tables/genetic_map_hg19_withX.txt.gz  --LDscoresFile /apps/well/bolt-lmm/2.3.2/tables/LDSCORE.1000G_EUR.tab.gz
