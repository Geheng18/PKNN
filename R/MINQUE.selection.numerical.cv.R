MINQUE.selection.numerical.cv = function(object,d){
  mincv=object$lambda.min
  sdcv=object$cvsd[which(object$lambda == object$lambda.min)]
  range = c(object$lambda.min - d*sdcv,object$lambda.min + d*sdcv)
  return(max(object$lambda[object$lambda>=range[1] & object$lambda<=range[2]]))
}