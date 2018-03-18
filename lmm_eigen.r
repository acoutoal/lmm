#' ##############3
#' Mixed effects model
#' Obtain U by applying SVD to X
#' X is the (genotype) matrix of SNPs  
lmm.svd <- function(X,y,U){
  #' Notice that: X   = USV'
  #'             XX'  = US^2U' 
  #'  Hence: cov(U'X) = S^2
  #'  Demo:   U'XX'U  = U'USV'VS'U'U = S^2             
  
  #Apply rotation
  Yr = t(U) %*% y
  Xr = t(U) %*% as.matrix(data.frame(Intercept = rep(1,nrow(X)),X))
  
  m1 = lm(Yr ~ Xr - 1)
  
  return(m1)
}
#' ##############3
#' Mixed effects model
#' Kinship matrix must be positive semi-definite
#' MZ twins not allowed, exclude them. A pair of MZ twins have
#' roughly the same entries in the K matrix, this leads to problems 
#' in the eigen decomposition  
#' Notice that XX'  = VLV'
#' Hence:     V'X     will diagonalize covariance of X into L,
#'           V'XX'V = V'VLV'V = L
lmm.eig <- function(X,y,K){
  
  #Check if K is symetric
  if ( identical(K,t(K)) == F ) {
    message("Kinship matrix is not symetric, exiting...")
    return(-1)
  }
  
  #Make eigendecomposition of K matrix
  f = eigen(K)
  
  #Check if K is singular
  Lmin = min(f$values)
  if ( Lmin < 0 ) {
    message("K is a singular matrix. Matrix was regularized by adding to its diagonal");
    f = eigen(K + diag(rep(-2*Lmin,nrow(K))))
    message("Smallest eigenvalue of the regularized matrix: ",signif(min(f$values),2)) 
  }
  
  L = f$valu
  V = f$vectors
  
  #Apply rotation
  Yr = t(V) %*% y
  Xr = t(V) %*% as.matrix(data.frame(Intercept = rep(1,nrow(X)),X))
  
  m1 = lm(Yr ~ Xr - 1)
  
  return(m1)
}

#' Runs demo
run_example <- function(){
  K = 50
  S = 5
  N = 1000
  P = 100
  C = 2
  SigmaK = 1 
  classroom_effects  = rnorm(K,mean = 0,sd = SigmaK)
  classroom_clusters = rep(seq(1,K),each=N/K, length.out=N)
  school_effects     = rnorm(S,mean = 0,sd = SigmaK)
  school_clusters    = rep(seq(1,S),each=N/S, length.out=N)  
  X  = matrix(rnorm(P*N), ncol = P)
  E  = rnorm(N)
  B  = c(rep(0,C), rep(c(-1,1),length.out = P-C))
  Z1 = classroom_effects[classroom_clusters] 
  Z2 = school_effects[school_clusters]
  Z  = Z1 + Z2
  X  = X + Z # confounding with outcome for all associations
  X  = scale(X)
  y  = X %*% B + Z + 2*E

  #Generate a matrix with dummies for each school and class
  Ku=model.matrix(object = ~u+v-1,data =  data.frame(u=as.factor(Z1),v=as.factor(Z2))) 
    
  #Factorization of X 
  m.svd = svd(X)
  U     = m.svd$u

  #check decomposition
  X.svd = m.svd$u %*% diag(m.svd$d) %*% t(m.svd$v)
  colMeans( (X.svd - X)^2 )

  #Projections
  PC.svd =  m.svd$u %*% diag(m.svd$d)  
  dim(PC.svd) #N x P

  #Plotting PCs
  layout(matrix(1:2))
  plot(PC.svd[,1],PC.svd[,2])  
  
  #Plotting clustered data
  mcl=hclust(d = dist(X,method = "euclidean"))
  plot(mcl)

  #' Demonstrate 3 methods of correcting multilevel data
  #' Notice that X1 is independent of y conditional on the random effects
  #' y ~ B*X1, B=0
  coef(summary(lm (y ~ X[,1]  ))) #Original model without correction
  coef(summary(lm (y ~ X[,1] + Ku )))[1:2,]
  coef(summary(lm (y ~ X[,1] + PC.svd[,1] )))[1:2,]
  coef(summary(lmm.svd(y, X[,1, drop = F],U)))
  
  #' Demonstrate 3 methods of correcting multilevel data
  #' Notice that X4 is dependent of y conditional on the random effects
  #' y ~ B*X4, B=1
  coef(summary(lm (y ~ X[,4]  )))
  coef(summary(lm (y ~ X[,4] + Ku )))[1:2,]
  coef(summary(lm (y ~ X[,4] + PC.svd[,1] )))[1:2,]
  coef(summary(lmm.svd(y,  X[,4, drop = F],U)))
  
}

