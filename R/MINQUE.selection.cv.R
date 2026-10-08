#' Cross validation for VLS
#'
#' @param y The continuous phenotypes of subjects.
#' @param geno A n*p genomic matrix with n subjects and p SNPs.
#' @param input.kernel The name of input kernel.
#' @param lambda The choice candidate parameter of pernalty term.
#' @param zbeta The fixed effect of phenotypes.
#'
#' @return The estimated variance components values will be returned.
#' @return The matrix of variance components will be returned.
#'
#' @examples
#' VLS.cv(y,geno,'product',lambda=c(0.3,3,10,50,100),criteria='MSE')
#' 

MINQUE.selection.cv = function(y, geno, input.kernel, MINQUE.type='MINQUE1',KNN.type='interaction',lambda = NULL, w=NULL, X=NULL, constrain=TRUE,criteria='COR'){
  # Initialize lambda
  if(is.null(lambda)){
    lambda = c(0.05,0.5,5,50,500,5000,50000,500000)
  }
  
  # Initialize the fix effect
  N = length(y)
  
  if(is.null(X)){
    X = matrix(0,N,N)
  }
  
  # Split training set
  n.lambda = length(lambda)
  n.batch = N%/%n.lambda
  n.region = length(geno)
  ind.temp1 = sample(1:N)
  ind = list()
  for (i in 0:(n.lambda-2)) {
    ind.temp2 = ind.temp1[(i*n.batch+1):((i+1)*n.batch)]
    ind.temp2 = list(ind.temp2)
    ind = c(ind,ind.temp2)
  }
  ind.temp2 = ind.temp1[-c(1:((n.lambda-1)*n.batch))]
  ind.temp2 = list(ind.temp2)
  ind = c(ind,ind.temp2)
  
  # Main selection process
  MSE = c()
  COR = c()
  for (i in 1:n.lambda) {
    ind.train = unlist(ind[-i])
    ind.test = unlist(ind[i])
    ytr = y[ind.train]
    yts = y[ind.test]
    Xfixtr = X[ind.train,]
    Xfixts = X[ind.test,]
    Xtr = list()
    Xts = list()
    for (j in 1:n.region) {
      Xtr.temp1 = geno[[j]][ind.train,]
      Xts.temp1 = geno[[j]][ind.test,]
      Xtr.temp2 = list(Xtr.temp1)
      Xts.temp2 = list(rbind(Xtr.temp1,Xts.temp1))
      Xtr = c(Xtr,Xtr.temp2)
      Xts = c(Xts,Xts.temp2)
    }
  	if(is.null(X)){ 
  	  result.MINQUE = MINQUE.selection(ytr, Xtr, input.kernel, MINQUE.type,KNN.type,lambda[i], w, X=NULL, constrain)
  	  y.hat = MINQUE.predict(ytr,Xts,input.kernel,result.MINQUE$theta,KNN.type,result.MINQUE$beta,Xtr=NULL,Xts=NULL)
    }else{
      result.MINQUE = MINQUE.selection(ytr, Xtr, input.kernel, MINQUE.type,KNN.type,lambda[i], w, Xfixtr, constrain)
  	  y.hat = MINQUE.predict(ytr,Xts,input.kernel,result.MINQUE$theta,KNN.type,result.MINQUE$beta,Xfixtr,Xfixts)
    }
    MSE[i] = mean((yts-y.hat)^2)
    COR[i] = cor(yts,y.hat)
  }
  
  # Return optimal lambda
  if(criteria=='MSE'){
    ind.opt = which.min(MSE)
    lambda.opt = lambda[ind.opt]
  }
  if(criteria=='COR'){
    ind.opt = which.max(COR)
    lambda.opt = lambda[ind.opt]
  }
  return(lambda.opt)
}
