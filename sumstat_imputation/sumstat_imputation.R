########################################################
#Function to do SS imputation from reference hap panel
########################################################

impute.ss<-function(ref.hap, sites.typed, betas, sd.betas, N, phi, min.AF=0.01,
                    shrink.fac=0.001, show.plot=TRUE){

  #Parameters
  l<-ncol(ref.hap);

  #First construct LD matrix
  rij<-cor(ref.hap);

  #Use a small amount of shrinkage to make it invertible / correct for LD overestimation
  ii<-1:l;
  im<-array(ii,c(l,l));
  del<-exp(-abs(im - t(im))*shrink.fac);

  #Get normalisation for Sigma
  fi<-apply(ref.hap, 2, mean);
  gt.var<-sqrt(2*fi*(1-fi));
  vv<- gt.var %*% t(gt.var);
  sig<-rij*del/vv;

  set.typed<-sort(sites.typed);
  set.untyped<-setdiff(1:l, sites.typed);

  #Perform mean imputation
  inv.sig.typed<-solve(sig[set.typed,set.typed]);
  sig12<-sig[set.untyped,set.typed];
  e.beta1<- (sig12 %*% inv.sig.typed) %*% (betas);
  sig11<-sig[set.untyped, set.untyped];
  cond.sig11<- sig11 - (sig12 %*% inv.sig.typed) %*% t(sig12);
  cond.sd1 <- sqrt(diag(cond.sig11)/(N*phi*(1-phi)));

  zs.typed<-abs(betas)/sd.betas;
  sd.betas1<-sqrt(diag(sig)[set.untyped]/(N*phi*(1-phi)));
  zs.untyped<-abs(e.beta1)/sd.betas1;

  if (show.plot) {
    plot(sites.typed, zs.typed, col="darkblue", pch=19, xlim=c(1,l), xlab="Position", ylab="|Z|-score");
    w.show<-which(fi[set.untyped]>min.AF & fi[set.untyped]<(1-min.AF));
    points(set.untyped[w.show], zs.untyped[w.show], col="lightblue", pch=19);
    segments(set.untyped[w.show], zs.untyped[w.show]-2*cond.sd1[w.show]/sd.betas1[w.show],
             set.untyped[w.show], zs.untyped[w.show]+2*cond.sd1[w.show]/sd.betas1[w.show], col=grey(0.75));
  }

  #Return object
  op<-array(0,c(l,7));
  colnames(op)<-c("Site", "Ref.freq", "Typed", "Beta",
                  "SD.Beta", "Cond.SD.Beta", "Z-score");
  op[,1]<-1:l;
  op[,2]<-fi;
  op[set.typed,3]<-1;
  op[set.typed,4]<-betas;
  op[set.typed,5]<-sd.betas;
  op[set.typed,7]<-abs(betas)/sd.betas;
  op[set.untyped,4]<-e.beta1;
  op[set.untyped,5]<-sd.betas1;
  op[set.untyped,6]<-cond.sd1;
  op[set.untyped,7]<-(e.beta1)/sd.betas1;

  return(op);
}
