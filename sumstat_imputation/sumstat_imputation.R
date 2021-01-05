########################################################
#Function to do SS imputation from reference hap panel
########################################################

impute_ss = function(ref_hap_panel,
                     sites.typed,
                     betas,
                     sd.betas,
                     case_control_constant,
                     min.AF     = 0.01,
                     shrink.fac = 0.001,
                     LD_Matrix  = FALSE,
                     LD         = NULL){

  #Parameters
  nsnps = ncol(ref_hap_panel)
  nhaps = nrow(ref_hap_panel)

  if (LD_Matrix == FALSE) { #Perform summary statistic imputation using the reference panel
    LD<-cor(ref_hap_panel)
  }

  if (LD_Matrix == TRUE) {  #Use LD matrix from inferred method
    LD = LD
  }

  #Use a small amount of shrinkage to make panel invertible / correct for LD overestimation
  ii = seq(1:nsnps)
  im<-array(ii,c(nsnps,nsnps))
  del<-exp(-abs(im - t(im))*shrink.fac) #Use linear alg to shrink non-diagonals

  #Get normalisation for Sigma
  ref_freq = colMeans(ref_hap_panel)
  gt.var   = sqrt(2*ref_freq*(1-ref_freq));
  vv       = gt.var %*% t(gt.var);
  R        = LD*del/vv;

  # #Figure out which sites are typed and untyped
  set.typed   = sites.typed
  set.untyped = setdiff(1:ncol(ref_hap_panel), sites.typed);

  #Solve correlation matrix
  inv.sig.typed = pseudoinverse(R[set.typed,set.typed], tol = 1e-10)
  W             = R[set.untyped,set.typed]
  expected_z    = W %*% inv.sig.typed %*% (betas/sd.betas)
  return(expected_z) #Return expected z scores
}
