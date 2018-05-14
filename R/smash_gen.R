#' A more general verison of smashgen that works for different data type and variance settings

#' @param x: a vector of observations
#' @param sigma: standard deviations, scalar. If not known, can be estimated.
#' @param wave_family: choice of wavelet basis to be used, as in wavethresh.
#' @param niter: number of iterations for IRLS
#' @param tol: tolerance of the criterion to stop the iterations
#' @param ashp: whether expand around ash posterior mean
#' @param robust: whether set the highest resolution wavelet coeffs to 0
#' @param verbose: whether print out the number of iterations to converge.
#' @param dist_family: family of ditribution: 'poisson','binomial'.
#' @param ntri: Binomial number of trials
#' @param method: method to estimate \sigma or \sigma^2+st^2 if \sigma is unknown.
#' @param reall: whether return all estimates including mean and standard deviation(\sigma, only when method='smashc','ols' or 'moment').
#' @return estimated mean
#' @export

smash_gen=function(x,sigma=NULL,ntri=NULL,wave_family='DaubExPhase',
                   ashp=TRUE,verbose=FALSE, robust=FALSE,
                   niter=1,tol=1e-2,dist_family='poisson',method='smashu',
                   reall=FALSE){
  n=length(x)
  mu=c()
  s=c()
  y=c()
  munorm=c()
  ini=init_smash_gen(x=x,ashp=ashp,ntri=ntri,dist_family=dist_family)
  mu=rbind(mu,ini$m0)
  s=rbind(s,ini$s0)
  y0=ini$y0
  #set wavelet coeffs to 0?
  if(robust){
    wds=wd(y0,family = family,filter.number = filter.number)
    wtd=threshold(wds, levels = wds$nlevels-1,  policy="manual",value = Inf)
    y=rbind(y,wr(wtd))
  }else{
    y=rbind(y,y0)
  }
  for(i in 1:niter){
    vars=ifelse(s[i,]<0,0,s[i,])
    mu.hat=smash.gaus(y[i,],sigma=sqrt(vars))#mu.hat is \mu_t+E(u_t|y)
    mu=rbind(mu,mu.hat)
    munorm[i]=norm(mu.hat-mu[i,],'2')
    if(munorm[i]<tol){
      if(verbose){
        message(sprintf('Converge after %i iterations',i))
      }
      break
    }
    #update m and s_t
    upd=update_smash_gen(mu.hat,x,ntri,dist_family)

    s=rbind(s,upd$st)
    y=rbind(y,upd$yt)
  }
  #give the final estimate
  if(is.null(sigma)){
    if(method=='smashu'){
      mu.hat=smash.gaus(y[i,])
    }
    if(method=='rmad'){
      mu.hat=smash.gaus(y[i,],sigma=sst_est(y[i,],method=method),family=wave_family)
    }
    if(method=='ols'|method=='smashc'|method=='moment'){
      sgg=smash.gaus.gen(y[i,],sqrt(ifelse(s[i,]<0,1e-8,s[i,])),
                         family=wave_family,method=method)
      mu.hat=sgg$mu.hat
      sd.hat=sgg$sd.hat
    }
  }else{
    mu.hat=smash.gaus(y[i,],sigma = sqrt(sigma^2+ifelse(s[i,]<0,1e-8,s[i,])))
  }
  m.est=update_smash_gen(mu.hat,x,ntri,dist_family)$mt
  if(reall&(method=='ols'|method=='smashc'|method=='moment')){
    return(list(m.est=m.est,sd.est=sd.hat))
  }else{
    return(m.est)
  }
}
