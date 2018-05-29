#' smash gaus to deal with covariate effect.
#' @param X: covariate, n by p
#' @param y: a vector of data
#' @param homo: wether homoscedastic or heteroscedatic variance.
#' @return a list: mu.hat, beta.hat, se.hat
#' @export
smash.gaus.x=function(X,y,filter.number=1,family='DaubExPhase',homo=T){
  proj=solve((t(X)%*%X))%*%t(X)
  beta.tilde=proj%*%y
  resi=y-X%*%beta.tilde
  smash.fit=smash.gaus(as.numeric(resi),filter.number=filter.number,family=family,
                       v.est=T,homoskedastic=homo,joint = T)
  mu.hat=smash.fit$mu.res
  var.hat=smash.fit$var.res
  #re-estimate beta
  if(homo){
    beta.hat=proj%*%(y-mu.hat)
  }else{
    beta.hat=lm((y-mu.hat)~X-1,weights=1/var.hat)$coefficients
  }
  return(list(mu.hat=mu.hat,beta.hat=beta.hat,se.hat=sqrt(var.hat)))
}
