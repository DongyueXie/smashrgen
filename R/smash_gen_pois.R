#'@title Smooth Poisson sequence, accounting for nugget effect
#'@param x observed Poisson sequence
#'@param s Scale factor for Poisson observations: y~Pois(scale*lambda), can be a vector.
#'@param transformation transformation of Poisson data, either 'vst' or 'lik_expan'; 'vst' for variance stabilizing transformation; 'lik_expansion' for likelihood expansion
#'@param lik_expan_at if transformation='lik_expan', where to expand it? Can be logx, or smash_poi
#'@param robust whether perform robust wavelet regression
#'@param robust.q quantile to determine outliers
#'@param smoother smoothing method for Gaussian sequence, either 'smash' or 'ti.thresh'. When n is large, ti.thresh is much faster
#'@param nug.init init value of nugget effect, either a scalar or NULL
#'@param ash.pm If choose lik_expansion, whether use ash posterior mean approximation if x=0. If not x = x+eps.
#'@param eps If choose lik_expansion, if x=0, set x = x + eps. Either input a numerical value or 'estimate'. If estimate, eps = sum(x==1)/sum(x<=1)
#'@param filter.number,family wavelet basis, see wavethresh package for more details
#'@param maxiter max iterations for estimating nugget effect
#'@param tol tolerance to stop iterations.
#'@return estimated smoothed lambda, estimated nugget effect.
#'@import smashr
#'@import ashr
#'@export

smash_gen_pois = function(x,
                          s=1,
                          #nugget=NULL,
                          nug.init = NULL,
                          est_nugget = TRUE,
                          transformation = 'lik_expan',
                          lik_expan_at ='logx',
                          nug.est.limit = 1,
                          smoother='smash',
                          robust = FALSE,
                          robust.q = 0.99,
                          ash_pm_init_for0=TRUE,
                          eps='estimate',
                          filter.number = 1,
                          family = "DaubExPhase",
                          est_nugget_maxiter=3,
                          est_nugget_tol=1e-3){
  t_start = Sys.time()

  if(!ispowerof2(length(x))){
    reflect.x = reflect(x)
    x = reflect.x$x
    idx = reflect.x$idx
  }else{
    idx = 1:length(x)
  }

  n = length(x)
  if(length(s)==1){
    s = rep(s,n)
  }

  if(transformation == 'vst'){
    y = sqrt(x+3/8)/sqrt(s)
    st = sqrt(0.25/s)
  }

  if(transformation == 'lik_expan'){

    if(lik_expan_at=='logx'){
      lambda_tilde = x/s
      if(min(x)<1){
        if(ash_pm_init_for0){
          x_pm = ash_pois(x,scale=s,link='identity')$result$PosteriorMean
          lambda_tilde[x<1] = x_pm[x<1]
        }else{
          if(eps == 'estimate'){
            eps = sum(round(x)==1)/sum(round(x)<=1)+0.1
          }
          lambda_tilde[x<1] = (x[x<1]+eps)/s
        }
      }
      # working data
      st=sqrt(1/(s*lambda_tilde))
      y=log(lambda_tilde)+(x-s*lambda_tilde)/(s*lambda_tilde)
    }

    if(lik_expan_at=='smash_poi'){
      lambda_tilde = smash.poiss(x)/s
      st=sqrt(1/(s*lambda_tilde))
      y=log(lambda_tilde)+(x-s*lambda_tilde)/(s*lambda_tilde)
    }


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
      nug.init = uniroot(normaleqn_nugget,c(-1e6,1e6),y=y,mu=y.rmed,st=st)$root
      nug.init = max(c(0,nug.init))
      outlier.idx = which(abs(y-y.rmed)>=(robust.z*sqrt(st^2+nug.init)))
    }else{
      outlier.idx = which(abs(y-y.rmed)>=robust.z*sqrt(st^2+nugget))
    }
    st[outlier.idx] = abs((y.wd.coefJ)[outlier.idx] - median(y.wd.coefJ))
  }


  if(est_nugget){
    fit = nugget_est(y,st,nug.init,nug.est.limit,smoother,filter.number = filter.number,family = family,est_nugget_maxiter,est_nugget_tol)
    nug.est = fit$nug.est
    mu.est = (fit$mu.res)[idx]
  }else{
    if(smoother=='smash'){
      mu.est = smash.gaus(y)[idx]
    }
    if(smoother=='ti.thresh'){
      mu.est = ti.thresh(y,method='rmad')[idx]
    }
    nug.est = NULL
  }

  if(transformation == 'vst'){
    lambda.est = mu.est^2-3/(8*s)
  }else{
    lambda.est = exp(mu.est)
  }

  t_end = Sys.time()
  return(list(posterior=list(mean_smooth=lambda.est,
                             mean_latent_smooth = mu.est),
              fitted_g=list(sigma2=nug.est),
              run_time = difftime(t_end,t_start,units='secs')))
  #return(list(lambda.est=lambda.est,mu.est=mu.est,nugget.est=nug.est))
}



normaleqn_nugget=function(nug,y,mu,st){
  return(sum((y-mu)^2/(nug+st^2)^2)-sum(1/(nug+st^2)))
}

nugget_est=function(y,st,nug.init=NULL,nug.est.limit,method,filter.number,family,maxiter,tol){
  #initialize nugget effect sigma^2
  n=length(y)
  if(nug.est.limit<1){
    top.idx = order(st,decreasing = F)[1:round(n*nug.est.limit)]
  }else{
    top.idx = 1:n
  }


  if(is.null(nug.init)){
    x.m=c(y[n],y,y[1])
    st.m=c(st[n],st,st[1])
    nug.init = ((x.m[2:(n+1)]-x.m[3:(n+2)])^2+(x.m[2:(n+1)]-x.m[1:(n)])^2-2*st.m[2:(n+1)]^2-st.m[1:(n)]^2-st.m[3:(n+2)]^2)/4
    #nug.init = ((x.m[2:n]-x.m[3:(n+1)])^2+(x.m[2:n]-x.m[1:(n-1)])^2-2*st.m[2:n]^2-st.m[1:(n-1)]^2-st.m[3:(n+1)]^2)/4
    #print(length(nug.init))
    #print(nug.init)
    nug.init = mean(nug.init[top.idx])
    nug.init = max(0,nug.init)
  }
  #given st and nug to estimate mean
  nug.est = nug.init
  for(iter in 1:maxiter){

    #print(nug.est)

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
    nug.est.new=uniroot(normaleqn_nugget,c(-1e6,1e6),y=y[top.idx],mu=mu.est[top.idx],st=st[top.idx])$root
    nug.est.new = max(c(0,nug.est.new))


    if(abs(nug.est.new - nug.est)<=tol){
      break
    }else{
      nug.est = nug.est.new
    }
  }

  return(list(mu.res = mu.est,mu.est.var=mu.est.var,nug.est=nug.est))

}

ispowerof2 <- function (x){
  x >= 1 & 2^ceiling(log2(x)) == x
}
