#'@title Smooth Poisson sequence, accounting for nugget effect, iterative pseudo-likelihood
#'@param x observed Poisson sequence
#'@param s Scale factor for Poisson observations: y~Pois(scale*lambda), can be a vector.
#'@param lik_expan_init if transformation='lik_expan', where to expand it? Can be mle_pm0, or smash_poi
#'@param robust whether perform robust wavelet regression
#'@param robust.q quantile to determine outliers
#'@param smoother smoothing method for Gaussian sequence, either 'smash' or 'ti.thresh'. When n is large, ti.thresh is much faster
#'@param nug.init init value of nugget effect, either a scalar or NULL
#'@param nugget_type 'homoskedastic' or 'heteroskedastic'
#'@param ash.pm If choose lik_expansion, whether use ash posterior mean approximation if x=0. If not x = x+eps.
#'@param eps If choose lik_expansion, if x=0, set x = x + eps. Either input a numerical value or 'estimate'. If estimate, eps = sum(x==1)/sum(x<=1)
#'@param filter.number,family wavelet basis, see wavethresh package for more details
#'@param maxiter max iterations for estimating nugget effect
#'@param tol tolerance to stop iterations.
#'@return estimated smoothed lambda, estimated nugget effect.
#'@import smashr
#'@import ashr
#'@export

smash_gen_pois_iterative = function(x,
                                    s=1,
                                    nug.init = NULL,
                                    nugget_type = 'heteroskedastic',
                                    lik_expan_init ='mle',
                                    smoother='smash',
                                    #robust = FALSE,
                                    #robust.q = 0.99,
                                    ash_pm_init_for0=TRUE,
                                    eps='estimate',
                                    filter.number = 1,
                                    family = "DaubExPhase",
                                    maxiter=10,
                                    tol=1e-3){

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

  if(lik_expan_init=='mle'){
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
  }

  if(lik_expan_init=='smash_poi'){
    lambda_tilde = smash.poiss(x)/s
  }

  st=sqrt(1/(s*lambda_tilde))
  y=log(lambda_tilde)+(x/(s*lambda_tilde) - 1)

  if(is.null(nug.init)){
    x.m=c(y[n],y,y[1])
    st.m=c(st[n],st,st[1])
    nug.init = ((x.m[2:n]-x.m[3:(n+1)])^2+(x.m[2:n]-x.m[1:(n-1)])^2-2*st.m[2:n]^2-st.m[1:(n-1)]^2-st.m[3:(n+1)]^2)/4
    nug.init = nug.init[nug.init>0&nug.init<var(y)]
    sigma2 = median(nug.init)
  }
  # estimate nugget effect and estimate mu
#
#   if(robust){
#     win.size = round(sqrt(n)/2)*2+1
#     #win.size = round(log(n,2)/2)*2+1
#     #browser()
#     y.wd = wd(y,filter.number,family,'station')
#     y.wd.coefJ = accessD(y.wd,level = log(n,2)-1)
#     y.rmed = runmed(y,win.size)
#
#     robust.z = qnorm(0.5+robust.q/2)
#
#     if(is.null(nugget)){
#       nug.init = uniroot(normaleqn,c(-1e6,1e6),y=y,mu=y.rmed,st=st)$root
#       nug.init = max(c(0,nug.init))
#       outlier.idx = which(abs(y-y.rmed)>=(robust.z*sqrt(st^2+nug.init)))
#     }else{
#       outlier.idx = which(abs(y-y.rmed)>=robust.z*sqrt(st^2+nugget))
#     }
#     st[outlier.idx] = abs((y.wd.coefJ)[outlier.idx] - median(y.wd.coefJ))
#   }



  mu.est.old = rep(Inf,n)
  convergence_trace = c()
  for(iter in 1:maxiter){
    # estimate smooth mean
    if(nugget_type=='heteroskedastic'){
      if(smoother=='smash'){
        fit = smash.gaus(y,v.est = T,joint=T,filter.number = filter.number,family = family)
        mu.est = fit$mu.res
        var.est = fit$var.res
      }
      if(smoother=='ti.thresh'){
        fit = ti.thresh(y,method = 'rmad',filter.number = filter.number,family = family,return_sd=TRUE)
        mu.est = fit$mu.est
        var.est = fit$sigma^2
      }

    }else{
      if(smoother=='smash'){
        mu.est = smash.gaus(y,sigma=sqrt(st^2+sigma2),filter.number = filter.number,family = family)
      }
      if(smoother=='ti.thresh'){
        mu.est = ti.thresh(y,sigma=sqrt(st^2+sigma2),filter.number = filter.number,family = family)
      }
    }

    # check convergence
    convergence_trace[iter] = sqrt(mean((mu.est.old-mu.est)^2))
    if(convergence_trace[iter]<tol){
      break
    }else{
      mu.est.old = mu.est
    }

    # estimate nugget variance
    if(nugget_type=='heteroskedastic'){
      sigma2 = pmax(var.est - st^2,0)
    }else{
      sigma2 = try(uniroot(normaleqn_nugget,c(0,var(y)),y=y,mu=mu.est,st=st)$root,silent=TRUE)
      if(class(sigma2)=='try-error'){
        sigma2 = 0
      }
    }

    # update posterior of nugget
    nugget_pm = (y-mu.est)/(1+st^2/sigma2)
    # new expansion
    lambda_tilde = exp(mu.est + nugget_pm)
    st=sqrt(1/(s*lambda_tilde))
    y = mu.est + nugget_pm + (x/(s*lambda_tilde)-1)

  }

  t_end = Sys.time()
  return(list(posterior=list(mean_smooth=exp(mu.est),
                             mean_latent_smooth = mu.est),
              fitted_g=list(sigma2=sigma2),
              convergence_trace=convergence_trace,
              run_time = difftime(t_end,t_start,units='secs')))
  #return(list(lambda.est=lambda.est,mu.est=mu.est,nugget.est=nug.est))
}


