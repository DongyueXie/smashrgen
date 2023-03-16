#'@title empirical Bayes Poisson smooting
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
#' lines(fit$posterior$mean_smooth)
#' fit$sigma2
#' plot(fit$elbo_trace)
#'@details The problem is
#'\deqn{x_i\sim Poisson(\lambda_i,}
#'\deqn{\lambda_i = \exp(\mu_i)),}
#'\deqn{\mu_i\sim N(b_i,\sigma^2),}
#'\deqn{\b_i\sim g(.).}
#'@import vebpm
#'@import wavethresh
#'@import smashr
#'@export




ebps = function(x,
               s = NULL,
               init_control = list(),
               general_control = list(),
               smooth_control = list()
               ){

  t_start = Sys.time()
  n_orig = length(x)
  if(is.null(s)){
    s = 1
  }
  if(length(s)==1){
    s = rep(s,n_orig)
  }

  init_controls = ebps_init_control_default()
  init_controls = modifyList(init_controls,init_control,keep.null = TRUE)

  general_controls = ebps_general_control_default()
  general_controls = modifyList(general_controls,general_control,keep.null = TRUE)

  smooth_controls = ebps_smooth_control_default()
  smooth_controls = modifyList(smooth_controls,smooth_control,keep.null = TRUE)

  init_val = ebps_init(x,s,
                       general_controls$make_power_of_2,
                       init_controls$sigma2_init,
                       init_controls$m_init,
                       init_controls$smooth_init,
                       init_controls$ash_pm_init_for0,
                       init_controls$eps_for0)

  sigma2 = init_val$sigma2_init
  m = init_val$m_init
  x = init_val$x
  s = init_val$s
  idx = init_val$idx
  Eb = init_val$smooth_init
  if(is.null(Eb)){
    Eb = rep(mean(m),length(x))
  }

  const = sum(lfactorial(x))
  v = rep(sigma2/2,length(x))

  if(smooth_controls$wave_trans=='ndwt' | smooth_controls$robust){
    general_controls$convergence_criteria = 'nugabs'
  }
  if(smooth_controls$wave_trans=='dwt'&is.null(smooth_controls$W)&(smooth_controls$filter.number != 1 | smooth_controls$family != 'DaubExPhase')){
    smooth_controls$W = (t(GenW(n,filter.number,family)))[-1,]
  }

  obj = -Inf
  s_update = list(Eb=Eb,
                  qb = list(fitted_g = NULL))

  Eb_old = Inf
  sigma2_trace = c(sigma2)

  for(iter in 1:general_controls$maxiter){

    s_update = ebps_smooth_update(m,sigma2,
                                  smooth_controls$filter.number,
                                  smooth_controls$family,
                                  smooth_controls$wave_trans,
                                  smooth_controls$ndwt_method,
                                  smooth_controls$warmstart,
                                  smooth_controls$ebnm_params,
                                  s_update$qb,
                                  smooth_controls$W,
                                  smooth_controls$robust,
                                  s_update$Eb)

    if(general_controls$plot_updates){
      plot(m,col='grey80')
      lines(s_update$Eb)
    }

    # opt = vga_pois_solver(m,x,s,Eb,sigma2,tol=vga_tol)
    if(general_controls$maxiter_vga==1){
      opt = vga_pois_solver_newton_1iter(m,v,x,s,s_update$Eb+s_update$E_out,sigma2)
    }else{
      opt = vga_pois_solver(m,x,s,s_update$Eb+s_update$E_out,sigma2,tol=general_controls$vga_tol,maxiter = general_controls$maxiter_vga)
    }

    m = opt$m
    v = opt$v

    # get sigma2
    if(general_controls$est_sigma2){
      sigma2_new = mean(m^2+v+s_update$Eb2+s_update$E_out2+2*s_update$Eb*s_update$E_out-2*m*(s_update$Eb+s_update$E_out))
      sigma2_trace = c(sigma2_trace,sigma2_new)
      if(general_controls$convergence_criteria=='nugabs'){
        if(abs(sigma2_new-sigma2)<general_controls$tol){
          break
        }
      }
      #print(sigma2_new)
      sigma2 = sigma2_new
    }else{
      if(general_controls$convergence_criteria=='nugabs'){
        if(sqrt(mean((s_update$Eb-s_update$Eb_old)^2))<general_controls$tol){
          break
        }
        Eb_old = s_update$Eb
      }
    }


    # calc obj
    if(general_controls$convergence_criteria=='objabs'){
      obj[iter+1] = pois_smooth_split_obj(x,s,m,v,s_update$Eb,s_update$Eb2,sigma2,s_update$qb$dKL,const)
      if(general_controls$verbose){
        if(iter%%general_controls$printevery==0){
          print(paste("Done iter",iter,"obj =",obj[iter+1]))
        }
      }

      if((obj[iter+1]-obj[iter])/n <general_controls$tol){
        break
      }
    }

  }
  t_end = Sys.time()
  if(smooth_controls$wave_trans=='dwt'){
    return(list(posterior=list(mean=exp(m+v/2)[idx],
                               mean_log = m[idx],
                               mean_smooth = exp(s_update$Eb)[idx],
                               mean_log_smooth=s_update$Eb[idx],
                               var_log = v[idx],
                               var_log_smooth = (s_update$Eb2-s_update$Eb^2)[idx]),
                fitted_g = list(sigma2=sigma2,sigma2_trace=sigma2_trace,g = s_update$qb$fitted_g),
                elbo=obj[length(obj)]/n*n_orig,
                elbo_trace = obj/n*n_orig,
                H = (s_update$qb$dKL + sum(log(2*pi*v)/2-log(2*pi*sigma2)/2-(m^2+v-2*m*s_update$Eb+s_update$Eb2)/2/sigma2))/n*n_orig,
                log_likelihood = obj[length(obj)]/n*n_orig,
                run_time = difftime(t_end,t_start,units='secs')))
  }else{
    return(list(posterior=list(mean = exp(m+v/2)[idx],
                               mean_log = m[idx],
                               mean_smooth = exp(s_update$Eb)[idx],
                               mean_log_smooth=s_update$Eb[idx]),
                log_likelihood = NULL,
                fitted_g = list(sigma2=sigma2,sigma2_trace=sigma2_trace),
                run_time = difftime(t_end,t_start,units='secs')))
  }
}

ebps_smooth_update = function(m,sigma2,
                              filter.number,
                              family,
                              wave_trans,
                              ndwt_method,
                              warmstart,
                              ebnm_params,
                              qb,
                              W,
                              robust,
                              Eb){

  if(robust){
    #print(range(m-Eb))
    res = ebnm(m-Eb,sqrt(sigma2),
               prior_family = 'point_laplace',
               mode = 0,
               scale =  "estimate",
               g_init = NULL,
               fix_g = FALSE)
    E_out = res$posterior$mean
    E_out2 = res$posterior$sd^2 + E_out^2
  }else{
    E_out = 0
    E_out2 = 0
  }
  m = m - E_out

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
  return(list(Eb=Eb,
              Eb2=Eb2,
              E_out=E_out,
              E_out2=E_out2,
              qb=qb))
}

ebps_init_control_default = function(){
  return(list(m_init = 'vga',
              smooth_init = NULL,
              ash_pm_init_for0 = TRUE,
              eps_for0 = 'estimate',
              sigma2_init = NULL))
}

ebps_general_control_default = function(){
  return(list(est_sigma2 = TRUE,
              maxiter = 100,
              maxiter_vga = 100,
              vga_tol = 1e-5,
              verbose=FALSE,
              tol=1e-5,
              printevery = 10,
              convergence_criteria = 'objabs',
              make_power_of_2 = 'reflect',
              plot_updates = FALSE))
}

ebps_smooth_control_default = function(){
  return(list(filter.number = 1,
              family = 'DaubExPhase',
              wave_trans='dwt',
              ndwt_method='ti.thresh',
              ebnm_params=list(),
              warmstart=TRUE,
              W=NULL,
              robust = FALSE
  ))
}

ebps_init = function(x,s,
                     make_power_of_2,
                     sigma2_init,
                     m_init,
                     smooth_init,
                     ash_pm_init_for0,
                     eps_for0
){
  n = length(x)
  # m_init, sigma2_init, smooth_init
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
          fit_init = ebpm_normal(x,s,g_init = list(mean=log(sum(x)/sum(s)),var=NULL),fix_g = c(T,F))
        }else{
          fit_init = ebpm_normal(x,s,g_init = list(mean=smooth_init,var=NULL),fix_g = c(T,F))
        }
        m_init = fit_init$posterior$mean_log
        sigma2_init = fit_init$fitted_g$var
      }else{
        if(is.null(smooth_init)){
          fit_init = ebpm_normal(x,s,g_init = list(mean=log(sum(x)/sum(s)),var = sigma2_init),fix_g = c(T,T))
        }else{
          fit_init = ebpm_normal(x,s,g_init=list(mean = smooth_init,var = sigma2_init),fix_g = c(T,T))
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
      sigma2_init = ebpm_normal(x,s,g_init = list(mean=log(sum(x)/sum(s)),var=NULL),fix_g=c(T,F))$fitted_g$var
      #sigma2_init = var(m - smash.gaus(m))
    }else{
      sigma2_init = var(m_init - smooth_init)
    }

  }
  #sigma2 = sigma2_init
  #m = m_init
  #v = rep(sigma2/2,n)
  return(list(sigma2_init = sigma2_init,
              smooth_init=smooth_init,
              m_init = m_init,
              idx=idx,
              x=x,
              s=s))
}

# pois_smooth_split_obj = function(x,s,m,s2,Eb,Eb2,sigma2,KLb,const){
#   return(sum(x*m-s*exp(m+s2/2)+log(s2)/2-log(sigma2)/2-(m^2+s2-2*m*Eb+Eb2)/2/sigma2)+KLb-const)
# }
#
# extend = function(x){
#   n = length(x)
#   J = log2(n)
#   if ((J%%1) == 0) {
#     return(list(x = x, idx = 1:n))
#   }else {
#     n.ext = 2^ceiling(J)
#     lnum = round((n.ext - n)/2)
#     rnum = n.ext - n - lnum
#     if (lnum == 0) {
#       x.lmir = NULL
#     }else {
#       x.lmir = x[lnum:1]
#     }
#     if (rnum == 0) {
#       x.rmir = NULL
#     }else {
#       x.rmir = x[n:(n - rnum + 1)]
#     }
#     x.ini = c(x.lmir, x, x.rmir)
#     return(list(x = x.ini, idx = (lnum + 1):(lnum + n)))
#   }
#
# }

