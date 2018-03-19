#' ##############3
#' Mixed effects model
#' Obtain U by applying SVD to G
#' G is the (genotype) matrix of SNPs  
lmm.svd <- function(X,y,G){
  #' Notice that: X   = USV'
  #'             XX'  = US^2U' 
  #'  Hence: cov(U'X) = S^2
  #'  Demo:   U'XX'U  = U'USV'VS'U'U = S^2             
  m.svd = svd(G)
  U     = m.svd$u
  S     = m.svd$d
  D     = diag(1/S) %*% t(U)
  
  #Apply rotation
  Yr = D %*% y
  Xr = D %*% as.matrix(data.frame(Intercept = rep(1,nrow(X)),X))
  
  m1 = lm(Yr ~ Xr - 1)
  
  return(m1)
}
#' ##############3
#' Mixed effects model
#' Kinship matrix must be positive semi-definite
#' MZ twins not allowed, exclude them. A pair of MZ twins have
#' roughly the same entries in the K matrix, this leads to problems 
#' in the eigen decomposition  
#' Notice that XX'  = VLV' = S
#' Hence:     V'X     will diagonalize covariance of X into L,
#'           V'XX'V = V'VLV'V = L
#'        S^{-1/2}  = VL^{-1/2}V'         
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
  
  D = V %*% diag(1/sqrt(L)) %*% t(V)
  
  #Apply rotation
  Yr = D %*% y
  Xr = D %*% as.matrix(data.frame(Intercept = rep(1,nrow(X)),X))
  
  m1 = lm(Yr ~ Xr - 1)
  
  return(m1)
}


#' ##############3
#' Mixed effects model
#' Kinship matrix must be positive semi-definite
#' MZ twins not allowed, exclude them. A pair of MZ twins have
#' roughly the same entries in the K matrix, this leads to problems 
#' in the cholesky decomposition  
#' Notice that XX'  = LL' = S
#'        S^{-1/2}  = L^{-1}         
lmm.chl <- function(X,y,K){
  
  #Check if K is symetric
  if ( identical(K,t(K)) == F ) {
    message("Kinship matrix is not symetric, exiting...")
    return(-1)
  }
  
  #Make choleski decomposition of K matrix and get lower triangular
  L = t(chol(K + diag(rep(4*1E-11,nrow(K)))))

  #S^{-1/2} = L^{-1}
  D = solve(L)
  
  #Apply rotation
  Yr = D %*% y
  Xr = D %*% as.matrix(data.frame(Intercept = rep(1,nrow(X)),X))
  
  m1 = lm(Yr ~ Xr - 1)
  
  return(m1)
}



#' Runs demo
run_example <- function(){
  K = 50
  S = 5
  N = 1000
  P = 100
  C = 10
  SigmaK = 1 
  classroom_effects  = rnorm(K,mean = 0,sd = SigmaK)
  classroom_clusters = rep(seq(1,K),each=N/K, length.out=N)
  school_effects     = rnorm(S,mean = 0,sd = SigmaK)
  school_clusters    = rep(seq(1,S),each=N/S, length.out=N)  
  X  = matrix(rnorm(P*N), ncol = P)
  E  = rnorm(N)
  B  = c(rep(0,C), rep(c(-1,1),length.out = P - C))
  Z1 = classroom_effects[classroom_clusters] 
  Z2 = school_effects[school_clusters]
  Z  = Z1 + Z2
  X[,1:C]  = X[,1:C] + Z # confounding with outcome for C associations
  X  = scale(X)
  y  = X %*% B + Z + 2*E

  #Generate a matrix with dummies for each school and class
  Ku=model.matrix(object = ~u+v-1,data =  data.frame(u=as.factor(Z1),v=as.factor(Z2))) 
  
  #Generate Kinship Matrix using the ancestry informative signals 1:C
  K = X[,1:C] %*% t(X[,1:C])
    
  #Factorization of X 
  m.svd = svd(X[,1:C])
  U     = m.svd$u

  #check decomposition
  X.svd = m.svd$u %*% diag(m.svd$d) %*% t(m.svd$v)
  colMeans( (X.svd - X[,1:C])^2 )

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
  coef(summary(lmm.svd(y, X[,1, drop = F],X[,1:C])))
  coef(summary(lmm.eig(y, X[,1, drop = F],K)))
  coef(summary(lmm.chl(y, X[,1, drop = F],K)))
  
  #' Demonstrate 3 methods of correcting multilevel data
  #' Notice that X4 is dependent of y conditional on the random effects
  #' y ~ B*X4, B=1
  coef(summary(lm (y ~ X[,C+1]  )))
  coef(summary(lm (y ~ X[,C+1] + Ku )))[1:2,]
  coef(summary(lm (y ~ X[,C+1] + PC.svd[,1] )))[1:2,]
  coef(summary(lmm.svd(y,  X[,C+1, drop = F],X[,1:C])))
  coef(summary(lmm.eig(y,  X[,C+1, drop = F],K)))
  coef(summary(lmm.chl(y,  X[,C+1, drop = F],K)))
  
}

