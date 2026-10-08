# main function
findKernel = function(kernelName, geno, maf = NULL, weights = 1, constant = 1, deg = 2,sigma = 1,theta = 1)
{
  switch(kernelName,
    ibs = {getIbs(geno,weights = 1)},     
    car = {getCar(geno = geno, maf = maf, weights = weights)},
    identity = {getIdentity(geno = geno)},
    product = {getProduct(geno = geno, weights = weights)},
    absproduct = {getAbsproduct(geno = geno, weights = weights)},
    polynomial = {getPolynomial(geno = geno, weights = weights, constant = constant, deg = deg)},
    gaussian = {getGaussian(geno = geno, weights = weights, sigma = sigma)},
    sigmoid={getSigmoid(geno = geno)},
    relu={getRelu(geno=geno)},
    exponential={getExponential(geno = geno)},
    J = {getJ(geno = geno)},
    sqidentity = {getSqidentity(geno = geno)},
    fstidentity = {getFstidentity(geno = geno)},
    RBF = {getRBF(geno = geno, sigma=sigma)},
    SP = {getSP(geno = geno, theta=theta)},
    stop("Invalied Kernel Name!")
  )
}


# Followings are the definition of base kernels
getIbs = function(geno,weights = 1)
{
  n = nrow(geno)
  m = ncol(geno)
  gtemp1 = geno
  gtemp2 = geno
  gtemp1[geno == 2] = 1
  gtemp2[geno == 1] = 0
  gtemp2[geno == 2]= 1
  gtemp = cbind(gtemp1,gtemp2)
  Inner = gtemp %*% diag(weights,nrow = 2*m,ncol = 2*m) %*% t(gtemp)
  X2 = matrix(diag(Inner),nrow = n,ncol = n)
  Y2 = matrix(diag(Inner),nrow = n,ncol = n,byrow = T)
  Bound = sum(matrix(weights,nrow = 1,ncol = m) * 2)
  Dis = (X2 + Y2 - 2 * Inner)
  IBS = Bound - Dis
  return(IBS)
}

getCar = function(geno, maf, weights)
{
  S=getIbs(geno,weights)/(2*sum(weights))
  if(weights==1){S=getIbs(geno,weights)/(2*ncol(geno))}
  diag(S)=0
  diag=rowSums(S)
  D=diag(diag)
  gamma = mean(cor(geno))
  Va = D - gamma * S
  VaL = chol(Va)
  K = chol2inv(VaL)
  return(K)
}

getIdentity = function(geno)
{
  n = nrow(geno)
  K = diag(rep(1, n))
  return(K)
}

getProduct = function(geno, weights = 1)
{
  m = ncol(geno)
  sqrtW = sqrt(diag(weights, m))
  K = (1/(m * sum(weights))) * tcrossprod(geno %*% sqrtW)
  return(K)
}

getAbsproduct = function(geno, weights = 1)
{
  m = ncol(geno)
  sqrtW = sqrt(diag(weights, m))
  K = (1/(m * sum(weights))) * tcrossprod((geno^2) %*% sqrtW)
  return(K)
}

getPolynomial = function(geno, weights =  1, constant = 1, deg =2)
{
  m = ncol(geno)
  K = (constant + (1/(m * sum(weights))) * tcrossprod(geno)) ^ deg
  return(K)
}

getGaussian = function(geno, weights = 1, sigma = 1)
{
  m = ncol(geno);
  wtGeno = t(t(geno) * sqrt(weights))
  distMat = as.matrix(dist(wtGeno, method = "euclidean", diag = T))
  K = exp(-(1/(2 * m * (sigma) * sum(weights))) * distMat ^ 2)
  return(K)
}

getSigmoid=function(geno,sigma=1)
{
  K=getProduct(geno)
  K=tanh(K+sigma)
  return(K)
}

getRelu=function(geno)
{
  K=getProduct(geno)
  K[K<0]=0 
  return(K)
}

getExponential=function(geno)
{
  K=getProduct(geno)
  K=exp(K)
  return(K)
}

getJ=function(geno)
{
  n = nrow(geno)
  K=matrix(1,n,n)
  return(K)
}

getSqidentity=function(geno){
  K = getProduct(geno)
  K = diag(diag(K)^2)
  return(K)
}

getFstidentity=function(geno){
  K = getProduct(geno)
  K = diag(diag(K))
  return(K)
}

getRBF=function(geno,sigma=1){
  n = nrow(geno)
  XtX = tcrossprod(geno)
  XX = matrix(1, n) %*% diag(XtX)
  D = XX - 2*XtX + t(XX)
  K = exp(-D/(2*sigma))
  return(K)
}

getSP=function(geno,theta=1){
  num = getProduct(geno)
  row.norm = theta + apply(geno,1,L2.norm)
  deno = sqrt(row.norm%*%t(row.norm))
  K = asin(num/deno)/(2*pi)
}

L2.norm = function(X){
  m = length(X)
  y = sum(X^2)/m
  return(y)
}
