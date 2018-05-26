#' A more general verison of smashgen that works for different data type and variance settings

#' @param x: a vector of observations
#' @param X: covariates
#' @param sigma: standard deviations, scalar. If not known, can be estimated.
#' @param wave_family: choice of wavelet basis to be used, as in wavethresh.
#' @param niter: number of iterations for IRLS
#' @param tol: tolerance of the criterion to stop the iterations
#' @param ashp: whether expand around ash posterior mean
#' @param robust: whether set the highest resolution wavelet coeffs to 0
#' @param verbose: whether print out the number of iterations to converge.
#' @param dist_family: family of ditribution: 'poisson','binomial','poi_binom'.
#' @param ntri: Binomial number of trials
#' @param y_var_est: method to estimate \sigma or \sigma^2+st^2 if \sigma is unknown.
#' @param re_all: whether return all estimates including mean and standard deviation(\sigma).
#' @param z_var_est: method to estimate variance of intermidiate variable z: 'rmad', 'smash', 'default'
#' @param k: parameter in huber m estimator
#' @param log_scale: whether return the estimate in log scale(for poisson data only)
#' @return estimated mean or (estimated mean and nuggest effect)
#' @export

smash_gen=function(x,X=NULL,sigma=NULL,ntri=NULL,z_var_est='rmad',wave_family='DaubExPhase',
                   ashp=TRUE,verbose=FALSE,robust=FALSE,filter.number=1,
                   niter=1,tol=1e-2,dist_family='poisson',y_var_est='smashu',
                   re_all=FALSE,log_scale=FALSE){
  n=length(x)
  m=c()
  s=c()
  y=c()
  ini=init_smash_gen(x=x,ashp=ashp,ntri=ntri,dist_family=dist_family)
  mu=rbind(mu,ini$m0)
  s=rbind(s,ini$s0)
  y0=ini$y0
  #set wavelet coeffs to 0?
  if(robust){
    y=rbind(y,fine_coeff_rm(y0,wave_family,filter.number))
  }else{
    y=rbind(y,y0)
  }
  # if niter >1, we do iterations; else, skip it to save time.
  if(niter>1&is.null(X)){
    for(i in 1:niter){
      vars=ifelse(s[i,]<0,0,s[i,])
      mu.hat=smash.gaus(y[i,],sigma=sqrt(vars),family = wave_family,filter.number = filter.number)#mu.hat is \mu_t+E(u_t|y)
      mu=rbind(mu,mu.hat)
      munorm=norm(mu.hat-mu[i,],'2')
      if(munorm<tol){
        if(verbose){
          message(sprintf('Converge after %i iterations',i))
        }
        break
      }
      if(verbose&i==niter){
        message(sprintf('Algorithm does not converge after %i iterations',i))
      }
      #update m and s_t
      upd=update_smash_gen(mu.hat,x,ntri,dist_family)

      s=rbind(s,upd$st)

      if(robust){
        y=rbind(y,fine_coeff_rm(upd$yt,wave_family,filter.number))
      }else{
        y=rbind(y,upd$yt)
      }
    }
  }else{
    i=1
  }

  #give the final estimate
  # if X is not provided, several methods to estimate variance
  # If X isp rovided, for now we directly use smashu
  if(is.null(X)){
    if(is.null(sigma)){
      if(y_var_est=='smashu'){
        mu.hat=smash.gaus(y[i,],family = wave_family,filter.number = filter.number)
      }else if(y_var_est=='rmad'){
        mu.hat=smash.gaus(y[i,],sigma=sst_est(y[i,],y_var_est,filter.number,wave_family),
                          family=wave_family,filter.number=filter.number)
      }else{
        sgg=smash.gaus.gen(y[i,],sqrt(ifelse(s[i,]<0,1e-8,s[i,])),
                           family=wave_family,y_var_est=y_var_est,
                           z_var_est=z_var_est,filter.number=filter.number,k=k)
        mu.hat=sgg$mu.hat
        sd.hat=sgg$sd.hat
      }
    }else{
      mu.hat=smash.gaus(y[i,],sigma = sqrt(sigma^2+ifelse(s[i,]<0,1e-8,s[i,])),family=wave_family,filter.number=filter.number)
    }
  }else{
    fit=smash.gaus.x(X,y[i,],filter.number = filter.number,family = wave_family,homo=F)
    mu.hat=fit$mu.hat
    beta.hat=fit$beta.hat
  }

  #whether return the log mu for poisson data
  if(log_scale&dist_family=='poisson'){
    mu.est=mu.hat
  }else{
    mu.est=update_smash_gen(mu.hat,x,ntri,dist_family)$mt
  }

  #If X is not provided, then return mu.est or mu.est and sd.est(nuggect effect)
  #If X is provided,return mu.est and bea.hat
  if(is.null(X)){
    if(re_all&exists('sd.hat')){
      return(list(mu.est=mu.est,sd.est=sd.hat))
    }else{
      return(mu.est)
    }
  }else{
    return(list(mu.est=mu.est,beta.est=beta.hat))
  }

}
