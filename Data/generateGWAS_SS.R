
#Simulate GWAS SS

library("gtools");

ref_hap.fn<-"./Data/ref1.txt";
fst<-0.2;
N.gwas<-1e4;
phi<-0.5;
shp.noise<-10000;


ref.hap<-as.matrix(read.table(ref_hap.fn, as.is=T, head=F, colClasses="numeric"));
fq.ref<-colMeans(ref.hap);
N.ref<-nrow(ref.hap);
l.ref<-ncol(ref.hap);

#Simulate WTS
alpha<-1/N.ref * (1-fst)/fst;

#Highly skewed
set.seed(234);
wts<-rdirichlet(1, rep(alpha, N.ref));


#Simulate GWAS popn
gwas<-array(0, c(N.gwas, l.ref));
ii<-sample(N.ref, 2*N.gwas, p=wts, rep=T);
for (i in 1:N.gwas) gwas[i,]<-ref.hap[ii[2*i-1],]+ref.hap[ii[2*i],];
y<-rbinom(N.gwas, 1, phi);
op<-cbind(1:l.ref, rep(0, l.ref), rep(0, l.ref));
colnames(op)<-c("Pos", "Beta", "se.Beta");
for (i in 1:l.ref) {
	slm<-summary(glm(y~gwas[,i], family="binomial"));
	op[i,2:3]<-slm$coef[2,1:2];
}
fq.gwas<-colMeans(gwas)/2;

write.table(op, file="./Data/sim1.txt", col=T, row=F, quote=F, sep="\t");





