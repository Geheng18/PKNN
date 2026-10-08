#'  How to choose tunning parameters in the process for selecting variables.
#'
#' @param object object is the fitted lasso model.
#' @param d d is to measure the acceptance of standard deviation. For example, d=1 means we can only choose those lambda which are located in the range of minimum error with one standard deviation.d=2 means we can only choose those lambda which are located in the range of minimum error with two standard deviations.

#' @return The lambda value.
#' @export



chooselambda=function(object,d){
  mincv=object$lambda.min
  sdcv=object$cvsd[which(object$lambda == object$lambda.min)]
  range <- c(object$lambda.min - d*sdcv,object$lambda.min + d*sdcv)
  max(object$lambda[object$lambda>=range[1] & object$lambda<=range[2]])
}

