#'@title Smooth over-dispersed Poisson sequence via splitting method
#'@param x data vector
#'@param maxiter,tol max iteration and tolerance for stopping it.
#'@param Eb_init,sigma2_init initial values of smooth mean and nugget effect.
#'@examples
#' set.seed(12345)
#' n=2^9
#' sigma=0.5
#' mu=c(rep(0.3,n/4), rep(3, n/4), rep(10, n/4), rep(0.3, n/4))
#' x = rpois(n,exp(log(mu)+rnorm(n,sd=sigma)))
#' fit = pois_smooth_split(x,maxiter=30)
#' plot(x,col='grey80')
#' lines(exp(fit$Eb))
#' fit$sigma2
#' plot(fit$obj)
#'@details The problem is
#'\deqn{x_i\sim Poisson(\exp(\mu_i)),}
#'\deqn{\mu_i\sim N(b_i,\sigma^2),}
#'\deqn{\b_i\sim g(.).}

pois_smooth_split = function(x,
                             s = NULL,
                             Eb_init = NULL,
                             sigma2_init = NULL,
                             est_sigma2 = TRUE,
                             maxiter = 100,
                             tol=1e-5,
                             filter.number = 1,
                             family = 'DaubExPhase',
                             verbose=FALSE,
                             printevery = 10,
                             ebnm_params=list(mode=0),
                             optim_method='L-BFGS-B'){

  n = length(x)
  if(is.null(s)){
    s = 1
  }
  if(length(s)==1){
    s = rep(s,n)
  }
  if(is.null(Eb_init)){
    Eb = log(runmed(x/s,1 + 2 * min((n-1)%/% 2, ceiling(0.1*n)))+0.01)
  }else{
    Eb = Eb_init
  }
  if(is.null(sigma2_init)){
    sigma2 = var(log(x/s+0.01)-Eb)
  }else{
    sigma2 = sigma2_init
  }

  W = (t(GenW(n,filter.number,family)))[-1,]

  mu_pm = rep(0,n)
  mu_pv = rep(1/n,n)
  obj = -Inf

  for(iter in 1:maxiter){
    # get m, s^2
    opt = optim(c(mu_pm,log(mu_pv)),
                fn = pois_mean_GG_opt_obj,
                gr = pois_mean_GG_opt_obj_gradient,
                x=x,
                s=s,
                beta=Eb,
                sigma2=sigma2,
                n=n,
                method = optim_method)
    mu_pm = opt$par[1:n]
    mu_pv = exp(opt$par[(n+1):(2*n)])
    qb = smash_dwt(mu_pm,sqrt(sigma2),filter.number=filter.number,family=family,ebnm_params=ebnm_params,W=W)
    Eb = qb$mu.est
    Eb2 = qb$mu.est.var + Eb^2
    # get sigma2
    if(est_sigma2){
      sigma2 = mean(mu_pm^2+mu_pv+Eb2-2*mu_pm*Eb)
    }


    # calc obj
    obj[iter+1] = pois_smooth_split_obj(x,s,mu_pm,mu_pv,Eb,Eb2,sigma2,qb$dKL)

    if(verbose){
      if(iter%%printevery==0){
        print(paste("Done iter",iter,"obj =",obj[iter+1]))
      }

    }

    if((obj[iter+1]-obj[iter])<tol){
      break
    }

  }
  return(list(Emean = exp(mu_pm+mu_pv/2),
              Vmean = exp(mu_pv-1)*exp(2*mu_pm+mu_pv),
              Emu=mu_pm,
              Vmu=mu_pv,
              Eb=Eb,
              Eb2=Eb2,
              sigma2=sigma2,
              obj=obj,
              H = qb$dKL + sum(log(2*pi*mu_pv)/2-log(2*pi*sigma2)/2-(mu_pm^2+mu_pv-2*mu_pm*Eb+Eb2)/2/sigma2)))
}

pois_smooth_split_obj = function(x,s,m,s2,Eb,Eb2,sigma2,KLb){
  return(sum(x*m-s*exp(m+s2/2)+log(s2)/2-log(sigma2)/2-(m^2+s2-2*m*Eb+Eb2)/2/sigma2)+KLb)
}

