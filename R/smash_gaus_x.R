#' smash gaus to deal with covariate effect.
#' @param x: covariate, n by p
#' @param y: a vector of data
#' @return a list: mu.hat, beta.hat
#' @export
smash.gaus.x=function(x,y,filter.number=1,family='DaubExPhase'){
  proj=solve((t(x)%*%x))%*%t(x)
  beta.tilde=proj%*%y
  resi=y-x%*%beta.hat
  mu.hat=smash.gaus(resi,filter.number=1,family='DaubExPhase')
  #re-estimate beta
  beta.hat=proj%*%(y-mu.hat)
  return(list(mu.hat=mu.hat,beta.hat=beta.hat))
}
