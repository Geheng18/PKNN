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

