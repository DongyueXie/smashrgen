#'@title Smooth Poisson sequence, accounting for nugget effect
#'@param x observed Poisson sequence
#'@param nugget nugget effect
#'@param s Scale factor for Poisson observations: y~Pois(scale*lambda), can be a vector.
#'@param transformation transformation of Poisson data, either 'vst' or 'lik_expansion'
#'@param robust whether perform robust wavelet regression
#'@param robust.q quantile to determine outliers
#'@param method smoothing method for Gaussian sequence, either 'smash' or 'ti.thresh'. When n is large, ti.thresh is much faster
#'@param nug.init init value of nugget effect, either a scalar or NULL
#'@param ash.pm If choose lik_expansion, whehter use ash posterior mean approxiamtion if x=0. If not x = x+eps.
#'@param eps If choose lik_expansion, if x=0, set x = x + eps. Either input a numerical value or 'estimate'. If estimate, eps = sum(x==1)/sum(x<=1)
#'@param filter.number,family wavelet basis, see wavethresh pakcage for more details
#'@param maxiter max iterations for estimating nugget effect
#'@param tol tolerance to stop iterations.
#'@return estimated smoothed \lambda, estimated nugget effect.
#'@import smashr
#'@import ashr
#'@export

smash.gen.poiss = function(x,nugget=NULL,s=1,
                           transformation = 'vst',
                           method='ti.thresh',
                           robust = F,
                           robust.q = 0.99,
                           est.nug = TRUE,
                           nug.init = NULL,
                           ash.pm=TRUE,ash.link = 'identity',
                           eps=0.01,
                           filter.number = 1,
                           family = "DaubExPhase",
                           maxiter=10,
                           tol=1e-2){

  if(min(x)<1 & ash.pm){
    x_pm = ash(rep(0,length(x)),1,lik=lik_pois(x,scale=s,link=ash.link),
               optmethod='mixSQP',pointmass=T)$result$PosteriorMean
  }

  if(!ispowerof2(length(x))){
    reflect.x = reflect(x)
    x = reflect.x$x
    idx = reflect.x$idx
    if(min(x)<1 & ash.pm){
      x_pm = reflect(x_pm)$x
    }
  }else{
    idx = 1:length(x)
  }

  n = length(x)

  if(transformation == 'vst'){
    y = sqrt(x+3/8)/sqrt(s)
    st = rep(sqrt(0.25),n)
    st[x==0] = 0
    #st[x==0] = 0.1383
    #st[x==1] = 0.2018
    #st[x==2] = 0.2294
    st = st/s
  }

  if(transformation == 'lik_expansion'){

    lambda_tilde = x/s
    if(min(x)<1){
      if(ash.pm){
        #x_pm = ash(rep(0,n),1,lik=lik_pois(x,scale=s),
        #           optmethod='mixSQP',pointmass=T)$result$PosteriorMean
        lambda_tilde[x<1] = x_pm[x<1]
      }else{
        lambda_tilde[x<1] = (x[x<1]+eps)/s
      }
    }
    # working data
    st=sqrt(1/(s*lambda_tilde))
    y=log(lambda_tilde)+(x-s*lambda_tilde)/(s*lambda_tilde)

  }


  # estimate nugget effect and estimate mu

  if(robust){
    win.size = round(sqrt(n)/2)*2+1
    #win.size = round(log(n,2)/2)*2+1
    #browser()
    y.wd = wd(y,filter.number,family,'station')
    y.wd.coefJ = accessD(y.wd,level = log(n,2)-1)
    y.rmed = runmed(y,win.size)

    robust.z = qnorm(0.5+robust.q/2)

    if(is.null(nugget)){
      nug.init = uniroot(normaleqn,c(-1e6,1e6),y=y,mu=y.rmed,st=st)$root
      nug.init = max(c(0,nug.init))
      outlier.idx = which(abs(y-y.rmed)>=(robust.z*sqrt(st^2+nug.init)))
    }else{
      outlier.idx = which(abs(y-y.rmed)>=robust.z*sqrt(st^2+nugget))
    }
    st[outlier.idx] = abs((y.wd.coefJ)[outlier.idx] - median(y.wd.coefJ))
  }


  if(is.null(nugget)){
    if(est.nug){
      fit = NuggetEst(y,st,nug.init,method,filter.number = filter.number,family = family,maxiter,tol)
      nug.est = fit$nug.est
      mu.est = fit$mu.est
    }else{
      mu.est = smash.gaus(y)
      nug.est = uniroot(normaleqn,c(-1e6,1e6),y=y,mu=mu.est,st=st)$root
      nug.est = max(c(nug.est,0))
    }
  }else{
    nug.est = nugget
    if(method=='smash'){
      mu.est = smash.gaus(y,sigma=sqrt(st^2+nug.est),
                       filter.number = filter.number,family = family)
    }
    if(method == 'ti.thresh'){
      mu.est = ti.thresh(y,sigma=sqrt(st^2+nug.est),filter.number = filter.number,family = family)
      #fit = list(mu.est = fit, mu.est.var = rep(0,length(fit)))
    }
  }

  mu.est = (mu.est)[idx]
  #mu.est.var = (fit$mu.est.var)[idx]

  if(transformation == 'vst'){
    lambda.est = mu.est^2+nug.est-3/(8*s)
  }
  if(transformation == 'lik_expansion'){
    lambda.est = exp(mu.est+nug.est/2)
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
    #x.m=c(y[n],y,y[1])
    #st.m=c(st[n],st,st[1])
    #nug.init = ((x.m[2:n]-x.m[3:(n+1)])^2+(x.m[2:n]-x.m[1:(n-1)])^2-2*st.m[2:n]^2-st.m[1:(n-1)]^2-st.m[3:(n+1)]^2)/4
    #nug.init = max(c(median(nug.init),0))
    #nug.init[nug.init<=0] = 0
    #nug.init = nug.init[nug.init<var(y)]
    #nug.init = nug.init[nug.init>0&nug.init<var(y)]
    #nug.init = mean(nug.init)

    mu.init = smash.gaus(y)
    nug.init = uniroot(normaleqn,c(-1e6,1e6),y=y,mu=mu.init,st=st)$root
    nug.init = max(c(nug.init,0))

  }
  #given st and nug to estimate mean


  nug.est = nug.init
  for(iter in 1:maxiter){

    #print(nug.est)

    # update mean
    if(method == 'smash'){
      mu.est = smash.gaus(y,sigma=sqrt(st^2+nug.est),filter.number = filter.number,family = family)
      #mu.est = est$mu.est
      #mu.est.var = est$mu.est.var
    }
    if(method == 'ti.thresh'){
      mu.est = ti.thresh(y,sigma=sqrt(st^2+nug.est),filter.number = filter.number,family = family)
      #mu.est.var = rep(0,n)
    }

    # update nugget effect
    nug.est.new=uniroot(normaleqn,c(-1e6,1e6),y=y,mu=mu.est,st=st)$root
    nug.est.new = max(c(0,nug.est.new))


    if(abs(nug.est.new - nug.est)<=tol){
      break
    }else{
      nug.est = nug.est.new
    }
  }

  #return(list(mu.est = mu.est,mu.est.var=mu.est.var,nug.est=nug.est))
  return(list(mu.est = mu.est,nug.est=nug.est))

}

#'@export
ispowerof2 <- function (x){
  x >= 1 & 2^ceiling(log2(x)) == x
}



