#Funcion that reads in a *haplotype* panel (generated from msprime) and reads out GWAS sumstats (SEs) form a GLM
msprime_gwas_sumstats = function(gwas_haplotypes){
  #Add rsids to make the filtering way easy
  rsids            = paste("rsid",seq(1:ncol(gwas_haplotypes)),sep = "")
  colnames(gwas_haplotypes)      = rsids

  #Split GWAS population into cases & controls
  #First convert GWAS haplotypes -> genotypes
  #gwas_genotypes          = t(sapply(seq(1,nrow(gwas_haplotypes),by=2),function(i) colSums(gwas_haplotypes[i:(i+1),])))
  #Remove duplicate haplotypes!
  #gwas_genotypes          = gwas_genotypes[!duplicated(gwas_genotypes), ]
  set.seed(143)
  status_index            = rbinom(seq(1:nrow(gwas_haplotypes)),1,0.5)

  #Remove everything in the GWAS panel that is outside of AF range
  #qc_variants_ref_min              = which(colMeans(reference_haplotypes)<0.01)
  #qc_variants_ref_max              = which(colMeans(reference_haplotypes)>0.99)
  qc_variants_gwas_max             = which(colMeans(gwas_haplotypes)<0.01)
  qc_variants_gwas_min             = which(colMeans(gwas_haplotypes)>0.99)

  filtered_snps            = unique(c(#qc_variants_ref_min,
                                      #qc_variants_ref_max,
                                      qc_variants_gwas_min,
                                      qc_variants_gwas_max))

  if (length(filtered_snps) > 0) {
    # #Remove these variants
    gwas_haplotypes           = gwas_haplotypes[,-filtered_snps]
    #reference_haplotypes      = reference_haplotypes[,-filtered_snps]
  }
   nsnps = ncol(gwas_haplotypes)
   op<-cbind(1:nsnps, rep(0, nsnps), rep(0, nsnps));
   colnames(op)<-c("SNP", "BETA", "SE");
   for (i in 1:nsnps) {
   flush.console()
   glm<-summary(glm(status_index~gwas_haplotypes[,i], family="binomial"));
   op[i,2:3]<-glm$coef[2,1:2];
   }
  op[,1] = seq(1:length(rsids))[-filtered_snps]
  return(list(gwas_haplotypes,op))
}

#Function that will actually perform the sumstat imputation
sumstat_impute = function(typed_snps,
                          untyped_snps_index,
                          genotyped_sumstats,
                          imputed_sumstats,
                          LD)
{
  LD_typed_untyped   = LD[typed_snps,untyped_snps_index] #LD of typed - untyped pairs
  print(LD_typed_untyped)
  inv_LD_typed       = solve( LD[typed_snps,typed_snps] + 0.001*diag(length(typed_snps))) #inverse of LD of typed SNPs
  W                  = LD[untyped_snps_index, typed_snps] %*% inv_LD_typed #these are the weights that turn typed to imputed
  infos              = as.numeric(rowSums(W * LD[untyped_snps_index,typed_snps])) #info measures per each untyped
  z.imp              = (W %*% (genotyped_sumstats[,2]/genotyped_sumstats[,3]^0.5))/sqrt(infos) #use scaling 1/sqrt(infos) to get var=1 for z-scores
  true_z             = -imputed_sumstats[,2]/(imputed_sumstats[,3]) #Flip beta alleles
  return(list(z.imp,true_z))
}


#This function calculates LD between all SNP's
LD_Matrix = function(haplotypes){

  #Takes in haplotype matrix (dim of haps x snps) to calculate various LD metrics (r here)
  haplotypes = t(haplotypes)
  pAB = haplotypes %*% t( haplotypes ) / ncol( haplotypes)
  pA  = rowMeans( haplotypes )
  pB  = rowMeans( haplotypes )
  D  = pAB - pA %*% t( pB )
  denA = sqrt( 1 / pA / ( 1 - pA ) )
  denB = sqrt( 1 / pB / ( 1 - pB ) )
  r  = t( D * denA ) * denB
  return(r)
}



