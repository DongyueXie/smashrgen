#' smash gaus to deal with covariate effect.
#' @param x: covariate, n by p
#' @param y: a vector of data
#' @param homo: wether homoscedastic or heteroscedatic variance.
#' @return a list: mu.hat, beta.hat, se.hat
#' @export
smash.gaus.x=function(x,y,filter.number=1,family='DaubExPhase',homo=T){
  proj=solve((t(x)%*%x))%*%t(x)
  beta.tilde=proj%*%y
  resi=y-x%*%beta.tilde
  smash.fit=smash.gaus(resi,filter.number=filter.number,family=family,
                       v.est=T,homoskedastic=homo,joint = T)
  mu.hat=smash.fit$mu.res
  var.hat=smash.fit$var.res
  #re-estimate beta
  if(homo){
    beta.hat=proj%*%(y-mu.hat)
  }else{
    beta.hat=lm((y-mu.hat)~x-1,weights=1/var.hat)$coefficients
  }
  return(list(mu.hat=mu.hat,beta.hat=beta.hat,se.hat=sqrt(var.hat)))
}
