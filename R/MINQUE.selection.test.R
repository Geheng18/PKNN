MINQUE.selection.test = function(y,geno,test.ind,input.kernel,MINQUE.type = 'MINQUE1',KNN.type = 'interaction',lambda = 500,w = NULL,X = NULL,constrain = FALSE){
  result = MINQUE.selection(y,geno,input.kernel,MINQUE.type,KNN.type,lambda,w,X,constrain)
  Cmat = result$Cmat
  beta = result$beta
  v = result$v
  N = length(y)
  sigmaH1 = Reduce(f ="+",x=mapply(FUN="*",v,Cmat[1,],SIMPLIFY=FALSE)) 
  sigmaH2 = Reduce(f ="+",x=mapply(FUN="*",v,Cmat[test.ind,],SIMPLIFY=FALSE)) 
  ratio = beta[test.ind]/beta[1]
  sigmaH = sigmaH2 - ratio*sigmaH1
  eigenvalue  = eigen(sigmaH)$values 
  pvalue = davies(q=0,lambda=eigenvalue)$Qq
  return(list(pvalue=pvalue))
}