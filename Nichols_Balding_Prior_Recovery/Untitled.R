#Simulate some allele freqs
nruns = 1e7
fq<-runif(nruns, 0.001, 0.999);
low_fq = fq[which(fq < 0.01)]
#Shape parameter for simulation of Obs/Exp VAR(BETA) ratio
alpha.ratio.base<-100;
alpha.var<-alpha.ratio.base*fq*(1-fq);
ratio.var<-rgamma(n = nruns , shape=alpha.var, rate=alpha.var);

#Calculate density of each simulated ratio
llk.obs<-dgamma(ratio.var, shape=alpha.var, rate=alpha.var, log=T);
hist(llk.obs)
print(exp(mean(llk.obs)))

low_fq  = fq[which(fq < 0.01)]
high_fq = fq[which(fq > 0.99)]
rare_fq = c(low_fq,high_fq)
alpha.ratio.base<-100;
alpha.var<-alpha.ratio.base*2*rare_fq*(1-rare_fq);
ratio.var<-rgamma(length(rare_fq), shape=alpha.var, rate=alpha.var);

#Calculate density of each simulated ratio
llk.obs<-dgamma(ratio.var, shape=alpha.var, rate=alpha.var, log=T);
hist(llk.obs)
print(exp(mean(llk.obs)))
