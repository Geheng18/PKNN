#' Variance Least Square Estimator
#'
#' @param y The continuous phenotypes of subjects.
#' @param geno A n*p genomic matrix with n subjects and p SNPs.
#' @param input.kernel The name of input kernel.
#' @param lambda The parameter of pernalty term.
#' @param zbeta The fixed effect of phenotypes.
#'
#' @return The estimated variance components values will be returned.
#' @return The matrix of variance components will be returned.
#'
#' @examples
#' VLS(y,geno,'product',lambda=0.3)
#' type = 'interaction' 'second' 'first'

KNN2LMM = function(geno,input.kernel,KNN.type='interaction'){
  
  # Get the dimension of KNN
  N = dim(geno[[1]])[1]
  kernel.matrix = list()
  n.region = length(geno)
  for (i in 1:n.region) {
    n.kernel = length(input.kernel[[i]])
    for (j in 1:n.kernel) {
      km.temp = findKernel(input.kernel[[i]][j],geno[[i]])
      km.temp = list(km.temp)
      kernel.matrix = c(kernel.matrix,km.temp)
    }
  }
  L = length(kernel.matrix)
  
  # List of variance component
  J = list(matrix(1,N,N))
  variance.component = c(J,kernel.matrix)
  
  if (KNN.type=='interaction'){
	  variance.component.list = list()
	  I = list(diag(1,N,N))
	  variance.component.list = c(I,variance.component.list)
	  for (i in 1:(L+1)){
	    for (j in i:(L+1)){
	      variance.component.temp = list(variance.component[[i]]*variance.component[[j]])
	      variance.component.list = c(variance.component.list,variance.component.temp)
	    }
	  }
  }
  
  if (KNN.type=='second'){
	  variance.component.list = list()
	  I = list(diag(1,N,N))
	  variance.component.list = c(I,variance.component)
	  for (i in 1:L){
	    variance.component.temp = list(variance.component[[i]]^2)
	    variance.component.list = c(variance.component.list,variance.component.temp)
	  }
  }
  
  if (KNN.type=='first'){
  	variance.component.list = list()
  	I = list(diag(1,N,N))
	  variance.component.list = c(I,variance.component)
  }
  

  return(variance.component.list)
}
