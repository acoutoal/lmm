gwas_lmm_demo <- function(){
  library(lme4)
  source("lmm_eigen.r")
  source("gen_sim.r")
  
  #' Simulate genotype data with 
  #' 1000 SNPs, 
  #' 300 individuals in each one of the 4 sub-populations, 
  #' 100 ancestry informative SNPs, 
  #' 100 causal SNPs, 
  #' Ancestry informative SNPs have d allele frequency differences among sub-populations, 
  #' Vg is the variance of the genetic effects, 
  #' Vs is the variance of the ancestry-related effects, and  
  #' Ve is the variance of the error term . 
  data = gen_sim(m = 1000, k = 4, n = 300, Ns = 100, Nc = 100, d = c(1.4,1.2,1,0.7), Vg = 2, Ve = 4 , Vs = 40 )
  #detach(data)
  attach(data)
  
  #' Show histograms of the phenotype Y, additive genetic component Yc, ancestry component Yz
  layout(matrix(1:3))
  hist(y, main="observed phenotype")
  hist(yz, main="ancestry-associated effect"); 
  hist(yc, main="additive genetic effect"); 
  
  #' Standard deviation of the genetic component
  print(sd(yc));
  #' Standard deviation of the ancestry component
  print(sd(yz))
  
  readline(prompt="Press [enter] to continue")
  
  #' Calculate genotype PCs using SVD 
  m.g=svd(geno)
  PC.svd =  m.g$u %*% diag(m.g$d)  
  
  
  #' Demonstrate population structure
  #' Plot scores (PC1 and PC2)
  layout(matrix(1))
  plot(PC.svd[,1:2], main="Population structure",cex.main=.8)
  points(PC.svd[group==1,1:2],col="blue")
  points(PC.svd[group==2,1:2],col="red")
  points(PC.svd[group==3,1:2],col="green")
  points(PC.svd[group==4,1:2],col="purple")
  
  readline(prompt="Press [enter] to continue")
  
  #' Run a GWAS 
  GWAS=matrix(nrow = m,ncol = 3)
  for(i in 1:m){
    GWAS[i,]=coef(summary(lm(y ~ geno[,i])))[2,-3]
  }
  
  #' Genomic inflation
  #' Calculate lambda
  chisq  = qchisq(GWAS[,3],1,lower.tail=FALSE)
  (lambda = median(chisq)/qchisq(.5,1,lower.tail=FALSE))
  
  #' QQ plot of GWAS P values
  qqplot(y= -log10(GWAS[,3]),x=-log10(runif(m)) )
  abline(a=0,b=1)
  usr <- par( "usr" )
  text(x = usr[2]-1,usr[3]+2,labels = sprintf("Inflation lambda=%1.2f",lambda))
  
  readline(prompt="Press [enter] to continue")
  
  #' Manhattan plot of GWAS P values 
  #' Structural and causal SNPs highlighted in red and blue  
  plot(-log10(GWAS[,3]),main = "Manhattan plot of GWAS without correction for pop structure",cex.main=.8)
  points(SNPc,-log10(GWAS[SNPc,3]),col="blue")
  points(SNPs,-log10(GWAS[SNPs,3]),col="red")
  abline(h=-log10(0.05/m),lty=2)
  legend(x="topright",legend = c("Ancestry SNP","Causal SNP"),col=c("red","blue"),pch=c(1,1))
  
  
  readline(prompt="Press [enter] to continue")
  
  #' Genetic relatedness matrix 
  #' Covariance and correlation among individuals 
  sgeno = scale(geno)
  GRM   = sgeno %*% t(sgeno) / (n*k-1)               #Covariance
  GRMs  = sgeno[,SNPs] %*% t(sgeno[,SNPs]) / (n*k-1) #Correlation

  image(GRMs,main="Relatedness matrix",col=rev(rainbow(5)))
  
  #' Ancestry informative PCs
  #' Demonstrate that PCA capture ancestry information 
  #' Calculate correlation of the first 10 PC with phenotype (Y) and with 
  #' ancestry component of the phenotype (Yc)   
  Ry.pc=cor(y,PC.svd)
  Ryz.pc=cor(yz,PC.svd)
  #' Ancestry correlations
  print(abs(Ryz.pc[1:10]))
  #' Phenotypic correlations
  print(abs(Ry.pc[1:10]))
  
  readline(prompt="Press [enter] to continue")
  
  #' Run a GWAS conditionining on the first 3 PCs
  GWAS.pca=matrix(nrow = m,ncol = 3)
  for(i in 1:m){
    GWAS.pca[i,]=coef(summary(lm(y ~ geno[,i] + PC.svd[,1:3])))[2,-3]
  }
  
  #' Genomic inflation
  #' Calculate lambda
  chisq  = qchisq(GWAS.pca[,3],1,lower.tail=FALSE)
  (lambda = median(chisq)/qchisq(.5,1,lower.tail=FALSE))
  
  #' QQ plot GWAS P values
  qqplot(y= -log10(GWAS.pca[,3]),x=-log10(runif(m)),main = "GWAS corrected for ancestry informative PCs",cex.main=.8)
  abline(a=0,b=1)
  usr <- par( "usr" )
  text(x = usr[2]-1,usr[3]+2,labels = sprintf("Inflation lambda=%1.2f",lambda))
  
  readline(prompt="Press [enter] to continue")
  
  #' Manhattan plot of GWAS P values 
  #' Structural and causal SNPs highlighted in red and blue  
  plot(-log10(GWAS.pca[,3]),main = "Manhattan plot for GWAS corrected for ancestry informative PCs",cex.main=.8)
  points(SNPc,-log10(GWAS.pca[SNPc,3]),col="blue")
  points(SNPs,-log10(GWAS.pca[SNPs,3]),col="red")
  abline(h=-log10(0.05/m),lty=2)
  legend(x="topright",legend = c("Ancestry SNP","Causal SNP"),col=c("red","blue"),pch=c(1,1))
  
  readline(prompt="Press [enter] to continue")

  #' Impact of PCA adjustment on ancestry informative SNPs
  #' GWAS P-values for the ancestry informative SNPs before and after correction for PC
  print(quantile(GWAS[SNPs,3]))
  print(quantile(GWAS.pca[SNPs,3]))

  
  #' ############################################################
  #' Compare effect estimates and P values of the ancestry informative SNP 
  #' FP associations
  #' ############################################################
  idx_s=which.min(GWAS[SNPs,3])
  C = SNPs[idx_s]
  print("lm (y ~ geno[,C]  )");
  print(coef(summary(lm (y ~ geno[,C]  ))))
  print("lm (y ~ geno[,C] + as.factor(group) )");
  print(coef(summary(lm (y ~ geno[,C] + as.factor(group) )))[1:2,])
  print("lm (y ~ geno[,C] + PC.svd[,1:3] )");
  print(coef(summary(lm (y ~ geno[,C] + PC.svd[,1:3] )))[1:2,])
  print("lmer (y ~ 1 + geno[,C] + (1|group) )");
  print(coef(summary(lmer (y ~ 1 + geno[,C] + (1|group) )))[1:2,])
  print("lmm.svd(y,  geno[,C, drop = F], sgeno )");
  print(coef(summary(lmm.svd(y,  geno[,C, drop = F], sgeno ) )))
  print("lmm.eig(y,  geno[,C, drop = F], GRM )");
  print(coef(summary(lmm.eig(y,  geno[,C, drop = F], GRM))))
  print("lmm.chl(y,  geno[,C, drop = F], GRM )");
  print(coef(summary(lmm.chl(y,  geno[,C, drop = F], GRM))))
  #' It will work with the standardized genotypes
  print("lmm.eig(y,  geno[,C, drop = F], GRMs )");
  print(coef(summary(lmm.eig(y,  geno[,C, drop = F], GRMs))))
  print("lmm.chl(y,  geno[,C, drop = F], GRMs )");
  print(coef(summary(lmm.chl(y,  geno[,C, drop = F], GRMs))))
  
  
  readline(prompt="Press [enter] to continue")
  
  
  #' ############################################################
  #' Compare effect estimates and P values of the causal SNP 
  #' TP associations
  #' ############################################################
  idx_s=which.min(GWAS[SNPc,3])
  C = SNPc[idx_s]
  print("lm (y ~ geno[,C]  )");
  print(coef(summary(lm (y ~ geno[,C]  ))))
  print("lm (y ~ geno[,C] + as.factor(group) )");
  print(coef(summary(lm (y ~ geno[,C] + as.factor(group) )))[1:2,])
  print("lm (y ~ geno[,C] + PC.svd[,1:3] )");
  print(coef(summary(lm (y ~ geno[,C] + PC.svd[,1:3] )))[1:2,])
  print("lmer (y ~ 1 + geno[,C] + (1|group) )");
  print(coef(summary(lmer (y ~ 1 + geno[,C] + (1|group) )))[1:2,])
  print("lmm.svd(y,  geno[,C, drop = F], sgeno )");
  print(coef(summary(lmm.svd(y,  geno[,C, drop = F], sgeno ) )))
  print("lmm.eig(y,  geno[,C, drop = F], GRM )");
  print(coef(summary(lmm.eig(y,  geno[,C, drop = F], GRM))))
  print("lmm.chl(y,  geno[,C, drop = F], GRM )");
  print(coef(summary(lmm.chl(y,  geno[,C, drop = F], GRM))))
  #' It will work with the standardized genotypes
  print("lmm.eig(y,  geno[,C, drop = F], GRMs )");
  print(coef(summary(lmm.eig(y,  geno[,C, drop = F], GRMs))))
  print("lmm.chl(y,  geno[,C, drop = F], GRMs )");
  print(coef(summary(lmm.chl(y,  geno[,C, drop = F], GRMs))))
  
  #' remove attached variables to environment
  detach(data)
}

  
