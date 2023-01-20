#'@title Smooth over-dispersed Poisson sequence via splitting method
#'@param x data vector
#'@param maxiter,tol max iteration and tolerance for stopping it.
#'@param m_init,sigma2_init initial values of latent variable and nugget effect.
#'@param smooth_init init of smooth function.
#'@param wave_trans dwt or ndwt. If ndwt, stopping criteria cannot be `objabs`
#'@param ndwt_method if wave_trans is ndwt, either use `smash` or `ti.thresh`. When n is large, `ti.thresh` is much faster.
#'@param convergence_criteria 'objabs' for absolute diff in ELBO, 'nugabs' for absolute diff in nugget effect
#'@param warmstart whether warmstart of dwt smoother
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
#'\deqn{x_i\sim Poisson(\lambda_i,}
#'\deqn{\lambda_i = \exp(\mu_i)),}
#'\deqn{\mu_i\sim N(b_i,\sigma^2),}
#'\deqn{\b_i\sim g(.).}
#'@import vebpm
#'@import wavethresh
#'@import smashr
#'@export

pois_smooth_split = function(x,
                             s = NULL,
                             m_init = 'vga',
                             smooth_init = NULL,
                             ash_pm_init_for0 = TRUE,
                             eps_for0 = 'estimate',
                             sigma2_init = NULL,
                             est_sigma2 = TRUE,
                             warmstart=TRUE,
                             maxiter = 100,
                             maxiter_vga = 100,
                             vga_tol = 1e-5,
                             tol=1e-5,
                             filter.number = 1,
                             family = 'DaubExPhase',
                             wave_trans='dwt',
                             ndwt_method='ti.thresh',
                             verbose=FALSE,
                             printevery = 10,
                             ebnm_params=list(mode=0),
                             convergence_criteria = 'objabs',
                             W=NULL,
                             make_power_of_2 = 'reflect',
                             plot_updates = FALSE){

  t_start = Sys.time()
  n = length(x)
  if(is.null(s)){
    s = 1
  }
  if(length(s)==1){
    s = rep(s,n)
  }

  if(!is.null(sigma2_init)){
    if(any(is.na(sigma2_init))){
      sigma2_init = NULL
    }
  }


  J = log(n,2)
  n_orig = n
  if(ceiling(J)!=floor(J)){
    #stop('Length of x must be power of 2')
    # reflect
    if(make_power_of_2=='reflect'){
      x = reflect(x)
      idx = x$idx
      x = x$x
      n = length(x)
      J = log(n,2)
      s = reflect(s)$x
      if(is.numeric(m_init)){
        m_init = reflect(m_init)$x
      }
      if(is.numeric(smooth_init)){
        smooth_init = reflect(smooth_init)$x
      }
    }
    if(make_power_of_2=='extend'){
      x = extend(x)
      idx = x$idx
      x = x$x
      n = length(x)
      J = log(n,2)
      s = extend(s)$x
      if(is.numeric(m_init)){
        m_init = extend(m_init)$x
      }
      if(is.numeric(smooth_init)){
        smooth_init = extend(smooth_init)$x
      }
    }
  }else{
    idx = 1:n
  }

  const = sum(lfactorial(x))
  if(!is.numeric(m_init)|length(m_init)!=n){
    if(m_init == 'smash_poi'){
      m_init = smash.poiss(x,log=TRUE) - log(s)
    }else if(m_init == 'logx'){
      m_init = log(x/s)
      if(min(x)==0){
        idx0 = (x == 0)
        if(ash_pm_init_for0){
          x_pm = ash_pois(x,scale=s,link='identity')$result$PosteriorMean
          m_init[idx0] = log(x_pm[idx0])
        }else{
          if(eps_for0 == 'estimate'){
            eps_for0 = sum(round(x)==1)/sum(round(x)<=1)+0.1
          }
          m_init[idx0] = log((x[idx0]+eps_for0)/s[idx0])
        }
      }
    }else if(m_init == 'vga'){
      if(is.null(sigma2_init)){
        if(is.null(smooth_init)){
          fit_init = pois_mean_GG(x,s,prior_mean = log(sum(x)/sum(s)))
        }else{
          fit_init = pois_mean_GG(x,s,prior_mean = smooth_init)
        }
        m_init = fit_init$posterior$mean_log
        sigma2_init = fit_init$fitted_g$var
      }else{
        if(is.null(smooth_init)){
          fit_init = pois_mean_GG(x,s,prior_mean = log(sum(x)/sum(s)),prior_var = sigma2_init)
        }else{
          fit_init = pois_mean_GG(x,s,prior_mean = smooth_init,prior_var = sigma2_init)
        }
        m_init = fit_init$posterior$mean_log
      }
    }else{
      stop('unknown init method of mu')
    }
  }

  if(is.null(sigma2_init)){
    #sigma2_init = var(m - ti.thresh(m,method='rmad'))
    if(is.null(smooth_init)){
      sigma2_init = pois_mean_GG(x,s,prior_mean = log(sum(x)/sum(s)))$fitted_g$var
      #sigma2_init = var(m - smash.gaus(m))
    }else{
      sigma2_init = var(m - smooth_init)
    }

  }
  sigma2 = sigma2_init
  m = m_init
  v = rep(sigma2/2,n)

  if(wave_trans=='ndwt'){
    convergence_criteria = 'nugabs'
  }

  if(wave_trans=='dwt'&is.null(W)&filter.number != 1&family != 'DaubExPhase'){
    W = (t(GenW(n,filter.number,family)))[-1,]
  }

  if(convergence_criteria=='objabs'){
    obj = -Inf
  }

  #m = rep(0,n)
  #v = rep(1/n,n)

  Eb_old = Inf

  sigma2_trace = c(sigma2)

  if(wave_trans=='dwt' & warmstart){
    qb = list(fitted_g = NULL)
  }

  for(iter in 1:maxiter){

    if(wave_trans=='dwt'){
      if(warmstart){
        qb = suppressWarnings(smash_dwt(m,sqrt(sigma2),filter.number=filter.number,family=family,ebnm_params=list(g_init=qb$fitted_g),W=W))
      }else{
        qb = smash_dwt(m,sqrt(sigma2),filter.number=filter.number,family=family,ebnm_params=ebnm_params,W=W)
      }

      Eb = qb$posterior$mean
      Eb2 = qb$posterior$var + Eb^2
    }
    if(wave_trans=='ndwt'){
      if(ndwt_method=='smash'){
        qb = smash.gaus(m,sqrt(sigma2),filter.number=filter.number,family=family,ebnm_param=ebnm_params,post.var = TRUE)
        Eb = qb$mu.est
        Eb2 = Eb^2+qb$mu.est.var
      }
      if(ndwt_method=='ti.thresh'){
        Eb = ti.thresh(m,sqrt(sigma2),filter.number=filter.number,family=family)
        Eb2 = Eb^2
      }
    }

    if(plot_updates){
      plot(m,col='grey80',ylim=range(m_init))
      lines(Eb)
    }

    # opt = vga_pois_solver(m,x,s,Eb,sigma2,tol=vga_tol)
    if(maxiter_vga==1){
      opt = vga_pois_solver_newton_1iter(m,v,x,s,Eb,sigma2)
    }else{
      opt = vga_pois_solver(m,x,s,Eb,sigma2,tol=vga_tol,maxiter = maxiter_vga)
    }

    m = opt$m
    v = opt$v

    # get sigma2
    if(est_sigma2){
      sigma2_new = mean(m^2+v+Eb2-2*m*Eb)
      sigma2_trace = c(sigma2_trace,sigma2_new)
      if(convergence_criteria=='nugabs'){
        if(abs(sigma2_new-sigma2)<tol){
          break
        }
      }
      #print(sigma2_new)
      sigma2 = sigma2_new
    }else{
      if(convergence_criteria=='nugabs'){
        if(sqrt(mean((Eb-Eb_old)^2))<tol){
          break
        }
        Eb_old = Eb
      }
    }


    # calc obj
    if(convergence_criteria=='objabs'){
      obj[iter+1] = pois_smooth_split_obj(x,s,m,v,Eb,Eb2,sigma2,qb$dKL,const)
      if(verbose){
        if(iter%%printevery==0){
          print(paste("Done iter",iter,"obj =",obj[iter+1]))
        }
      }

      if((obj[iter+1]-obj[iter])/n <tol){
        break
      }
    }

  }
  t_end = Sys.time()
  if(wave_trans=='dwt'){
    return(list(posterior=list(mean=exp(m+v/2)[idx],
                               mean_log = m[idx],
                               mean_smooth = exp(Eb)[idx],
                               mean_log_smooth=Eb[idx]),
                # posterior_full=list(mean_smooth = exp(Eb),
                #                     mean_lambda=exp(m+v/2),
                #                     var_lambda = exp(v-1)*exp(2*m+v),
                #                     mean_mu = m,
                #                     var_mu = v,
                #                     mean_latent_smooth = Eb,
                #                     var_latent_smooth = Eb2-Eb^2),
                fitted_g = list(sigma2=sigma2,sigma2_trace=sigma2_trace),
                elbo=obj[length(obj)]/n*n_orig,
                elbo_trace = obj/n*n_orig,
                H = (qb$dKL + sum(log(2*pi*v)/2-log(2*pi*sigma2)/2-(m^2+v-2*m*Eb+Eb2)/2/sigma2))/n*n_orig,
                log_likelihood = obj[length(obj)]/n*n_orig,
                run_time = difftime(t_end,t_start,units='secs')))
  }else{
    return(list(posterior=list(mean = exp(m+v/2)[idx],
                               mean_log = m[idx],
                               mean_smooth = exp(Eb)[idx],
                               mean_log_smooth=Eb[idx]),
                # posterior_full=list(mean_smooth = exp(Eb),
                #                mean_lambda=exp(m+v/2),
                #                var_lambda = exp(v-1)*exp(2*m+v),
                #                mean_mu = m,
                #                var_mu = v,
                #                mean_latent_smooth = Eb,
                #                var_latent_smooth = Eb2-Eb^2),
                log_likelihood = NULL,
                fitted_g = list(sigma2=sigma2,sigma2_trace=sigma2_trace),
                run_time = difftime(t_end,t_start,units='secs')))
  }
}

pois_smooth_split_obj = function(x,s,m,s2,Eb,Eb2,sigma2,KLb,const){
  return(sum(x*m-s*exp(m+s2/2)+log(s2)/2-log(sigma2)/2-(m^2+s2-2*m*Eb+Eb2)/2/sigma2)+KLb-const)
}

extend = function(x){
  n = length(x)
  J = log2(n)
  if ((J%%1) == 0) {
    return(list(x = x, idx = 1:n))
  }else {
    n.ext = 2^ceiling(J)
    lnum = round((n.ext - n)/2)
    rnum = n.ext - n - lnum
    if (lnum == 0) {
      x.lmir = NULL
    }else {
      x.lmir = x[lnum:1]
    }
    if (rnum == 0) {
      x.rmir = NULL
    }else {
      x.rmir = x[n:(n - rnum + 1)]
    }
    x.ini = c(x.lmir, x, x.rmir)
    return(list(x = x.ini, idx = (lnum + 1):(lnum + n)))
  }

}

