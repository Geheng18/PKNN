#'  Calculate the Linear kernel for genomic matrix
#'
#' @param Gen Gen is a n*p genomic matrix with n subjects and p SNPs.

#' @return The kernel matrix of the corresponding genomic matrix.
#' @export

Fun_Klinear=function(Gen) Gen%*%t(Gen)/ncol(Gen)
