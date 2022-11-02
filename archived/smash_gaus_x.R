#' smash gaus to deal with covariate effect.
#' @param X: covariate, n by p
#' @param y: a vector of data
#' @param homo: wether homoscedastic or heteroscedatic variance.
#' @param sigma: variance known or not, if not provided, will be estimated.
#' @param niter: how many iterations to do, default is 1
#' @param tol: stop criterion
#' @param verbose: whether print out convergence information.
#' @return a list: mu.hat, beta.hat, se.hat
#' @export
smash.gaus.x=function(X,y,sigma=NULL,niter=1,
                      filter.number=1,family='DaubExPhase',
                      homo=T,tol=1e-2,verbose=F){
  proj=solve((t(X)%*%X))%*%t(X)
  beta.tilde=proj%*%y
  resi=y-X%*%beta.tilde
  for(iter in 1:niter){
    smash.fit=smash.gaus(as.numeric(resi),sigma=sigma,filter.number=filter.number,family=family,
                         v.est=T,homoskedastic=homo,joint = T)
    mu.hat=smash.fit$mu.res
    var.hat=if(is.null(sigma)){smash.fit$var.res}else{sigma^2}

    beta.hat.old=if(iter>1){beta.hat}else{0}

    if(homo|!is.null(sigma)){
      beta.hat=proj%*%(y-mu.hat)
    }else{
      beta.hat=lm((y-mu.hat)~X-1,weights=1/var.hat)$coefficients
    }
    #check if converge, beta

    if(norm(beta.hat-beta.hat.old,'2')<=tol){
      if(verbose&niter!=1){
        message(sprintf('Converge after %i iterations',iter))
      }
      break
    }
    if(verbose&iter==niter&niter!=1){
      message(sprintf('Algorithm does not converge after %i iterations',iter))
    }
    resi=y-X%*%beta.hat
  }
  #re-estimate beta
  return(list(mu.hat=mu.hat,beta.hat=beta.hat,se.hat=sqrt(var.hat)))
}
