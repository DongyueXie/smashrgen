#'@title Smooth Poisson squence, accounting for nugget effect
#'@param x observed Poisson sequence
#'@param nugget nugget effect
#'@param s Scale factor for Poisson observations: y~Pois(scale*lambda), can be a vector.
#'@param transformation transformation of Poisson data, either 'vst' or 'lik_expansion'
#'@param method smoothing method for Gaussian sequence, either 'smash' or 'ti.thresh'. When n is large, ti.thresh is much faster
#'@param lambda If choose lik_expansion, for \tilde{\lambda}, either set it to x('mle') or smash.poiss('smoothing') output
#'@param ash.pm If choose lik_expansion, whehter use ash posterior mean approxiamtion if x=0. If not x = x+eps.
#'@param eps If choose lik_expansion, if x=0, set x = x + eps.
#'@param filter.number,family wavelet basis, see wavethresh pakcage for more details
#'@param maxiter max iterations for estimating nugget effect
#'@param tol tolerance to stop iterations.
#'@return estimated smoothed \lambda, expectation of log(\lambda), estimated nugget effect.
#'@import smashr
#'@import ashr
#'@export

smash.gen.pois = function(x,nugget=NULL,s=1,transformation = 'vst',
                          method='ti.thresh',
                          lambda = 'smoothing',
                          ash.pm=FALSE,eps=0.0001,
                          filter.number = 1,
                          family = "DaubExPhase",
                          maxiter=10,
                          tol=1e-2){

  if(!ispowerof2(length(x))){
    reflect.x = reflect(x)
    x = reflect.x$x
    idx = reflect.x$idx
  }else{
    idx = 1:length(x)
  }

  n = length(x)

  if(transformation == 'vst'){
    y = sqrt(x)/sqrt(s)
    st = rep(sqrt(0.25/s),n)
  }

  if(transformation == 'lik_expansion'){
    if(lambda == 'smoothing'){
      lambda_tilde = smash.poiss(x) / s
    }

    if(lambda == 'mle'){
      lambda_tilde = x/s
      if(min(x)<1){
        if(ash.pm){
          x_pm = ash(rep(0,n),1,lik=lik_pois(x,scale=s),
                     optmethod='mixSQP',pointmass=T)$result$PosteriorMean
          lambda_tilde[x<1] = x_pm[x<1]
        }else{
          lambda_tilde[x<1] = (x[x<1]+eps)/s
        }
      }
    }


    # working data
    st=sqrt(1/(s*lambda_tilde))
    y=log(lambda_tilde)+(x-s*lambda_tilde)/(s*lambda_tilde)

  }


  # estimate nugget effect and estimate mu

  nug.est = nugget

  if(is.null(nugget)){
    fit = NuggetEst(y,st,nug.est,method,filter.number = filter.number,family = family,maxiter,tol)
    nug.est = fit$nug.est
  }else{
    if(method=='smash'){
      fit = smash.gaus(y,sigma=sqrt(st^2+nug.est),
                       filter.number = filter.number,family = family,
                       post.var = TRUE)
    }
    if(method == 'ti.thresh'){
      fit = ti.thresh(y,sigma=sqrt(st^2+nug.est),filter.number = filter.number,family = family)
      fit = list(mu.est = fit, mu.est.var = rep(0,length(fit)))
    }
  }

  mu.est = (fit$mu.est)[idx]
  mu.est.var = (fit$mu.est.var)[idx]

  if(transformation == 'vst'){
    lambda.est = mu.est^2
  }
  if(transformation == 'lik_expansion'){
    lambda.est = exp(mu.est+mu.est.var/2)
  }

  return(list(lambda.est=lambda.est,mu.est=mu.est,nugget.est=nug.est))
}



normaleqn=function(nug,y,mu,st){
  return(sum((y-mu)^2/(nug+st^2)^2)-sum(1/(nug+st^2)))
}

NuggetEst=function(y,st,nug.init=NULL,method,filter.number,family,maxiter,tol){
  #initialize nugget effect sigma^2
  n=length(y)
  if(is.null(nug.init)){
    x.m=c(y[n],y,y[1])
    st.m=c(st[n],st,st[1])
    nug.init = ((x.m[2:n]-x.m[3:(n+1)])^2+(x.m[2:n]-x.m[1:(n-1)])^2-2*st.m[2:n]^2-st.m[1:(n-1)]^2-st.m[3:(n+1)]^2)/4
    nug.init = nug.init[nug.init>0&nug.init<var(y)]
    nug.init = median(nug.init)
  }
  #given st and nug to estimate mean


  nug.est = nug.init
  for(iter in 1:maxiter){
    # update mean
    if(method == 'smash'){
      est = smash.gaus(y,sigma=sqrt(st^2+nug.est),filter.number = filter.number,family = family,post.var = TRUE)
      mu.est = est$mu.est
      mu.est.var = est$mu.est.var
    }
    if(method == 'ti.thresh'){
      mu.est = ti.thresh(y,sigma=sqrt(st^2+nug.est),filter.number = filter.number,family = family)
      mu.est.var = rep(0,n)
    }

    # update nugget effect
    nug.est.new=uniroot(normaleqn,c(-1e6,1e6),y=y,mu=mu.est,st=st)$root

    if(abs(nug.est.new - nug.est)<=tol){
      break
    }else{
      nug.est = nug.est.new
    }
 }

  return(list(mu.est = mu.est,mu.est.var=mu.est.var,nug.est=nug.est))

}


ispowerof2 <- function (x){
  x >= 1 & 2^ceiling(log2(x)) == x
}


