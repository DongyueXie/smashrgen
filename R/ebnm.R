#' A function to do eb for normal mean
#' @param x: a vector of data
#' @param s: a vector of sd
#' @param eps: adjustment
#' @return eb estimate of prior mu, sd and posterior mean
#' @export
ebnm_normal=function(x,s,eps=1e-8){
  mu=mean(x)
  sigma=sqrt(max(c(0,mean((x-mu)^2)-s^2)))+eps
  posmean=mu*(1/sigma^2/(1/s^2+1/sigma^2))+x*(1/s^2/(1/s^2+1/sigma^2))
  return(list(mu=mu,sigma=sigma,posmean=posmean))
}
