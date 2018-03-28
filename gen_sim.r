############################################
#' Simulate genotype (SNP) and phenotype data with 
#' k subpopulations. The phenotype has  genetic and ancestry components. 
#' It is possible to define informative SNPs, causal SNPs, 
#' and non-informative SNPs.
#' 
#' INPUT ARGUMENTS
#'  m  is the number of SNPs
#'  k  is the number of sub-populations
#'  n  is the number of individuals in each sub-population, total = n * k
#'  Ns is the number of ancestry informative SNPs located in the geno matrix between (1:Ns) 
#'  Nc is the number of causal SNPs located in the geno matrix between (Ns+1):(Ns+Nc) 
#'  d  is the proportion of change in allele frequency of ancestry informative SNPs 
#'  Vg is the variance of the allele effects for causal SNPs assuming a random effects model with 0 mean.
#'  Vs is the variance of the allele effects for the ancestry informative SNPs assuming a random effects model with 0 mean.
#'  Ve is the variance of the error term
#'  
#'  MODEL
#'  b    ~ N(0, sqrt(Vg))
#'  u    ~ N(0, sqrt(Vs))
#'  e    ~ N(0, sqrt(Ve))
#'  SNP is a (k*n) by 1 vector with additive allele dosages for a SNP.  
#'  Z   is a (k*n) by k  matrix of ancestry factors (dummy variables)  
#'  y    = SNP * b + Z * u + e
#'  
#'  OUTPUT VALUE IS A LIST OF
#'  geno  is the (k*n) by m genotype matrix 
#'  y     is the (k*n) by 1 phenotype vector  
#'  yc    is the (k*n) by 1 vector with the causal genetic component of the phenotype (usually non-observable)
#'  yz    is the (k*n) by 1 vector with the ancestry-related component of the phenotype (usually non-observable)
#'  Z     is the (k*n) by k matrix with dummy variables for the ancestry subgroups 
#'  SNPc  is the Nc by 1 vector with the geno matrix index of the causal genetic effect SNPs 
#'  SNPs  is the Ns by 1 vector with the geno matrix index of the ancestry informative SNPs 
#'  group is the (k*n) by 1 vector with index of the subpopulations (another way of expressing Z) 
#'  beta  is the Nc by 1 vector with the genetic effect of each causal SNP 
#'  d,m,n,k are the arguments supplied and described above
gen_sim <- function(m=1000, k=3, n=1000, Ns=100, Nc=100, d=c(1.2,1,0.8), Vg = 2, Ve = 4, Vs = 2 ){
  library(MASS)
  
  stopifnot( length(d) == k )
  
  # Generate the frequency of the coding allele for all m SNPs
  f = c(.1,.2,.3,.4,.5)
  p = f[sample(x = 1:5, size = m,replace = T,prob = c(.3,.25,.2,.15,.1))]

  #Ancestry informative SNPs
  SNPs = 1:Ns
  
  #Causal SNPs
  SNPc  = (Ns+1):(Ns+Nc)
  
  # Generate the SNP data
  geno=data.frame()
  for(i in 1:k){
    x = simGeno(m = m, n = n, p = p, Ns = Ns, d = d[i])
    if(i == 1){
      geno = x
    }else{
      geno = rbind(geno,x)
    }
  }
  
  # create a factor to keep track of which "students" are in which "school"
  group <- as.factor(rep(1:k, each = n))
  
  #Generate ancestry effect (confounder of phenotype with ancestry)
  u  <- sqrt(Vs)*scale(rnorm(k,mean = 0,sd = 1))
  Z  <- kronecker( diag(k), rep(1,n) )
  yz <- Z %*% u
  
  #Generate causal effect 
  beta <- sqrt(Vg)*scale(rnorm(n = Nc,mean = 0,1))
  yc   <- geno[, SNPc] %*% beta
  
  #Generate phenotype
  e    <- rnorm(n = Nc,mean = 0,sqrt(Ve)) 
  y    <- yc + yz + e
  
  return(list(geno=geno, y=y, yz=yz, yc=yc, Z=Z, SNPc=SNPc, SNPs=SNPs, group=group, beta=beta, d=d, m=m, n=n, k=k))
  
}

simGeno <- function(m,n,p,Ns,d){
  p[1:Ns] = p[1:Ns] * d 
  geno = matrix(nrow = n,ncol = m)
  for(i in 1:m){
    geno[,i] = rbinom(n, 2, p[i] )
  }
  return(geno)
}


