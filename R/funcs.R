#' mean squared error
#' @export
mse = function(x, y) mean((x - y)^2)

#' l2 norm
#' @export
l2norm = function(x) sum(x^2)

#' mean integrated squared error
#' @export
mise = function(x, y, r) 10000 * mean(apply(x - rep(1, r) %o% y, 1, l2norm)/l2norm(y))

#' data.frame to ggplot data
#' @param df
#' @export
df2gg=function(df){
  n=dim(df)[1]
  p=dim(df)[2]
  name=colnames(df)
  MSE=as.vector(unlist(df))
  method=factor(rep(colnames(df),each=n),levels = name)
  dfgg=data.frame(MSE=MSE,method=method)
  return(dfgg)
}
