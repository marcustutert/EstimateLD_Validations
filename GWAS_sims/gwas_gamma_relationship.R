library("statmod");

#Look at distribution of SE values from logistic model

n.gwas<-1e4;
n.sim<-1000;
phi<-0.3;
y<-rbinom(n.gwas, 1, phi);

op<-array(0, c(n.sim, 3));
colnames(op)<-c("FQ", "beta", "se.beta");
for (i in 1:n.sim) {
  fq<-runif(1,0.02, 0.98);
  x<-rbinom(n.gwas, 2, fq);
  y<-sample(y);
  lg<-summary(glm(y~x, family="binomial"));
  op[i,]<-c(mean(x)/2, lg$coef[2,1:2]);
}
e.se<-sqrt(1/(n.gwas*phi*(1-phi))*1/(2*op[,1]*(1-op[,1])));

plot(e.se, op[,3]);

#Plot graphics in plotly

abline(0,1);
rr<-op[,3]^2/e.se^2;
hist(rr);
mean(rr);
sd(rr);
cat("\nImplied shape parameter for Gamma distribution for VAR = ", mean(rr)/var(rr));

plot(op[,1], rr, xlab="Frequency", ylab="Obs Var / Exp Var");

#Suggests that lower-frequency variants have lower shape parameter.

#Obtain implied gamma shape parameter as a function of allele frequency
nq<-10;
qq<-seq(0,1,l=11);
fi<-findInterval(op[,1], qq);
gm.op<-array(0, c(nq, 4));
colnames(gm.op)<-c("Mean.FQ", "Mean.R", "Med.R", "Shape");
for (i in 1:nq) {
  wi<-which(fi==i);
  gm.op[i,]<-c(mean(op[wi,1]), mean(rr[wi]), median(rr[i]), mean(rr[wi])/var(rr[wi]));
}

print(gm.op);

#Look at Heterozygosity versus shape parameter
plot(2*gm.op[,1]*(1-gm.op[,1]), gm.op[,4], xlab="Het", ylab="Shape parameter");
abline(0,1);

#Looks linear - i.e. variance of ratio (inverse of shape parameter) goes with x(1-x)


