#' Phenotypes prediction based on VLS
#'
#' @param ytr The continuous phenotypes of training subjects.
#' @param geno A N*p genomic matrix of Ntr training subjects and Nts testing subjects.
#' @param input.kernel The name of input kernels.
#' @param theta Estimated variance components.
#' @param zbeta.tr The fixed effect of phenotypes of training set.
#' @param zbeta.ts The fixed effect of phenotypes of testing set.
#' 
#' @return The prediction of phenotypes of testing subjects will be returned.
#'
#' @examples
#' VLS.predict(ytr,geno,'product',theta)
#' 

MINQUE.predict = function(y,geno,input.kernel,theta,KNN.type='interaction',beta=NULL,Xtr=NULL,Xts=NULL){
  
  N = dim(geno[[1]])[1]
  Ntr = length(y)
  Nts = N-Ntr
  
  # Initialize fix effect
  if (is.null(Xtr) | is.null(Xts)) {
    Xtr = matrix(0,Ntr,1)
    Xts = matrix(0,Nts,1)
    beta = 0
  }
  
  # Get dimension of KNN and transfer to LMM
  variance.component.list = KNN2LMM(geno, input.kernel, KNN.type)
  n.variance.component = length(variance.component.list)
  
  # Predict y.hat
  K.U = Reduce(`+`,mapply(`*`, variance.component.list[2:n.variance.component], theta[2:n.variance.component], SIMPLIFY =FALSE))
  sigma = K.U[(Ntr + 1):N, 1:Ntr]
  Sigma = K.U[1:Ntr, 1:Ntr] + diag(theta[1], Ntr)
  y.hat = as.numeric(Xts %*% beta + sigma %*% ginv(Sigma) %*% (y - Xtr %*%beta))
  
  # Return results
  return(y.hat)
}
