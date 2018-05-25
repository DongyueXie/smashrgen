#' Trend filtering allow covariates
#' @param x: covariates, n by p
#' @param y: a vector of data
#' @param ord: order of trend filter
#' @param lambda: reg param, default to be sqrt(n*log(p))
#' @return estimtaed mu and beta
#' @export

tf.gen=function(x,y,ord=1,lambda=NULL){
  n=nrow(x)
  p=ncol(x)
  x.aug=cbind(diag(n),x)
  D=getDtf(n+p,ord=ord)
  idx=which(D[,n]!=0)[1]
  D=D[1:idx,]
  fit=genlasso(y,x.aug,D)
  if(is.null(lambda)){lambda=sqrt(n*log(n+p))}
  coefs=coef.genlasso(fit,lambda=lambda)$beta
  return(list(mu.hat=coefs[1:n],beta.hat=coefs[(n+1):(n+p)]))
}
