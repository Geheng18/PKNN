
#'  A Lasso-based linear mixed model with generalized method-of-moments estimation for complex phenotype prediction.
#'
#' @param Y The continuous phenotypes of subjects.
#' @param Gen Gen is a list where each element is a n*p genomic matrix with n subjects and p SNPs.
#' @param index index is the index for testing subjects.
#' @param K K is a list where each element is the kernel matrix of the corresponding genomic matrix in the list Gen.
#' @param returnK returnK=T means the kernel matrix will be returned.returnK=F means the kernel matrix will not be returned. The default is returnK=F
#'
#'
#' @return The prediction values and ture values for testing subjects will be returned.
#' @return The effects sizes for each genomic region can also returned.
#' @export
#' @import glmnet
#' @import MASS
#' @importFrom glmnet cv.glmnet
#'
#' @examples
#'OGen=matrix(sample(0:2,500*150,replace = TRUE),500,150)
#'start <- seq(1, by = 3, length = ncol(OGen) / 3)
#'Gen <- lapply(start, function(i, OGen) OGen[,i:(i+2)], OGen = OGen)
#'y = rowSums(scale(Gen[[1]]))+rowSums(scale(Gen[[2]]))+rnorm(500)
#'index= sample(1:length(y),100)
#'fit=GmmLasso::GMMLasso(y=y, Gen=Gen, index=index, K = NULL, returnK = F)


  
MINQUE.selection.numerical = function(y,geno,input.kernel,alpha = 0,KNN.type = 'interaction',lambda = NULL,X = NULL){
  

  K = KNN2LMM(geno, input.kernel, KNN.type)
  list_K=lapply(K, scaleK)
  
  numK=length(list_K)
  
  if(is.null(X)){
    n=length(y);
    V=diag(n)
    E=eigen(V)$vectors
    A=E                       ##get A matrix
  }else{
    n=length(y);
    V=diag(n)-X %*% MASS::ginv(t(X) %*% X) %*% t(X)
    E=eigen(V)$vectors
    nc=ncol(X)
    A=E[,1:(n-nc)]                       ##get A matrix
  }
  
  
  lasso_Pheno = c(t(A)%*%y%*%t(y)%*%A)   ##training Phenotypes used in GMMLasso
  Tlist=lapply(list_K,function(kk) c(t(A)%*%kk%*%A))
  T=matrix(unlist(Tlist),nrow=ncol(A)*ncol(A),ncol=length(list_K))
  
  
  p.fac = c(rep(1,numK))
  p.fac[1] = 0
  fit_cv = cv.glmnet(T, lasso_Pheno, alpha=alpha,lower.limits =0,penalty.factor = p.fac)     ## fit the model
  
  if(is.null(lambda)){
    lambda=MINQUE.selection.numerical.cv(fit_cv,1) 
  }
  ##search the optimal lambda
  
  
  fit_best = glmnet(T, lasso_Pheno, alpha = alpha, lambda = lambda, lower.limits =0, penalty.factor = p.fac)       ## fit the model
  theta= as.numeric(fit_best$beta)   
  return(list(theta=theta))
}

#'  Scaling the Linear kernel matrix for genomic matrix
#'
#' @param K K is the kernel matrix of the corresponding genomic matrix..

#' @return The scaled kernel matrix.
#' @export
scaleK=function(K){
  Ktrace=sum(diag(K))
  Kscale=K/Ktrace*nrow(K)
  Kscale
}









