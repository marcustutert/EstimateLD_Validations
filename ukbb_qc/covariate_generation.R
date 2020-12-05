#Generating BNeale Covariates (sex + age + age^2 + sex*age + sex*age2)
#Create a covariate file (BNeale_Covariate) which has as column:
#FID/IID (will be the same in UKBB)/sex/age/age^2/sex*age/sex*age2)

#Read in filtered_sqcpheno from Chris script
library(data.table)
fam = fread("hypertension.fam")

pheno           = fread("ukb7749.csv") #Add in age and sex column, match with eIDs
pheno$eid       = as.numeric(as.character(pheno$eid))
#Add age and sex column to covariates
#Create a file which has
column_header   = c("FID","IID","sex","age","age^2","sex*age","sex*age^2")
pheno_filtered  = pheno[which( pheno$eid %in% unlist(fam[,1])),]
#Get each of the columns necessary
age             = as.numeric(pheno_filtered$'21003-0.0')
sex             = as.numeric(pheno_filtered$'31-0.0')
eid             = pheno_filtered$eid
covariate       = cbind(eid,eid,sex,age,age^2,sex*age,sex*age^2)
covariate       = rbind(column_header,covariate)
#Write out covariate file
write.table(covariate,"ukbb_neale_covariate", quote = F,row.names = F, col.names = F)
