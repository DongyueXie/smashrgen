#'@title Smooth over-dispersed Poisson sequence via splitting method
#'@param x data vector
#'@param maxiter,tol max iteration and tolerance for stopping it.
#'@param Eb_init,sigma2_init initial values of smooth mean and nugget effect.
#'@param wave_trans dwt or ndwt. If ndwt, stopping criteria cannot be `objabs`
#'@param ndwt_method if wave_trans is ndwt, either use `smash` or `ti.thresh`. When n is large, `ti.thresh` is much faster.
#'@param convergence_criteria 'objabs' for absolute diff in ELBO, 'nugabs' for absolute diff in nugget effect
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

pois_smooth_split_init_b = function(x,
                             s = NULL,
                             Eb_init = 'runmed',
                             sigma2_init = NULL,
                             est_sigma2 = TRUE,
                             maxiter = 100,
                             tol=1e-5,
                             filter.number = 1,
                             family = 'DaubExPhase',
                             wave_trans='dwt',
                             ndwt_method='smash',
                             verbose=FALSE,
                             printevery = 10,
                             ebnm_params=list(mode=0),
                             optim_method='L-BFGS-B',
                             convergence_criteria = 'objabs',
                             W=NULL,
                             sigma2_est_top = NULL){

  t_start = Sys.time()
  n = length(x)
  if(is.null(s)){
    s = 1
  }
  if(length(s)==1){
    s = rep(s,n)
  }
  const = sum(lfactorial(x))
  if(!is.numeric(Eb_init)|length(Eb_init)!=n){
    if(Eb_init=='runmed'){
      Eb = log(runmed(x/s,1 + 2 * min((n-1)%/% 2, ceiling(0.1*n)))+0.01)
    }
    if(Eb_init == 'smash_poi'){
      Eb = smash.poiss(x,log=TRUE) - log(s)
    }
    if(Eb_init == 'smooth_gaus'){
      Eb = ti.thresh(log(1/s+x/s),method = 'rmad')
    }
    if(Eb_init == 'log1px'){
      Eb = log(1/s+x/s)
    }
  }else{
    Eb = Eb_init
  }
  if(is.null(sigma2_init)){
    sigma2 = var(log(x/s+0.01)-Eb)
  }else{
    sigma2 = sigma2_init
  }

  if(wave_trans=='ndwt'){
    convergence_criteria = 'nugabs'
  }

  if(wave_trans=='dwt'&is.null(W)){
    W = (t(GenW(n,filter.number,family)))[-1,]
  }

  if(convergence_criteria=='objabs'){
     obj = -Inf
  }

  if(!is.null(sigma2_est_top)&convergence_criteria == 'nugabs'&est_sigma2){
    top_idx = order(x,decreasing = TRUE)[1:round(n*sigma2_est_top)]
  }

  mu_pm = rep(0,n)
  mu_pv = rep(1/n,n)

  Eb_old = Eb

  sigma2_trace = c()

  for(iter in 1:maxiter){
    # get m, s^2
    opt = vga_optimize(c(mu_pm,log(mu_pv)),x,s,Eb,sigma2)
    mu_pm = opt$m
    mu_pv = opt$v

    if(wave_trans=='dwt'){
      qb = smash_dwt(mu_pm,sqrt(sigma2),filter.number=filter.number,family=family,ebnm_params=ebnm_params,W=W)
      Eb = qb$posterior$mean
      Eb2 = qb$posterior$var + Eb^2
    }
    if(wave_trans=='ndwt'){
      if(ndwt_method=='smash'){
        qb = smash.gaus(mu_pm,sqrt(sigma2),filter.number=filter.number,family=family,ebnm_param=ebnm_params,post.var = TRUE)
        Eb = qb$mu.est
        Eb2 = Eb^2+qb$mu.est.var
      }
      if(ndwt_method=='ti.thresh'){
        Eb = ti.thresh(mu_pm,sqrt(sigma2),filter.number=filter.number,family=family)
        Eb2 = Eb^2
      }
    }
    # get sigma2
    if(est_sigma2){
      if(convergence_criteria=='nugabs'&!is.null(sigma2_est_top)){
        sigma2_new = mean((mu_pm^2+mu_pv+Eb2-2*mu_pm*Eb)[top_idx])
      }else{
        sigma2_new = mean(mu_pm^2+mu_pv+Eb2-2*mu_pm*Eb)
      }
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
      obj[iter+1] = pois_smooth_split_obj(x,s,mu_pm,mu_pv,Eb,Eb2,sigma2,qb$dKL,const)
      if(verbose){
        if(iter%%printevery==0){
          print(paste("Done iter",iter,"obj =",obj[iter+1]))
        }
      }

      if((obj[iter+1]-obj[iter])<tol){
        break
      }
    }

  }
  t_end = Sys.time()
  if(wave_trans=='dwt'){
      return(list(posterior=list(mean_smooth = exp(Eb),
                             mean_lambda=exp(mu_pm+mu_pv/2),
                             var_lambda = exp(mu_pv-1)*exp(2*mu_pm+mu_pv),
                             mean_mu = mu_pm,
                             var_mu = mu_pv,
                             mean_latent_smooth = Eb,
                             Var_latent_smooth = Eb2-Eb^2),
              fitted_g = list(sigma2=sigma2,sigma2_trace=sigma2_trace),
              obj_value=obj,
              H = qb$dKL + sum(log(2*pi*mu_pv)/2-log(2*pi*sigma2)/2-(mu_pm^2+mu_pv-2*mu_pm*Eb+Eb2)/2/sigma2),
              run_time = difftime(t_end,t_start,units='secs')))
  }else{
    return(list(posterior=list(mean_smooth = exp(Eb),
                               mean_lambda=exp(mu_pm+mu_pv/2),
                               var_lambda = exp(mu_pv-1)*exp(2*mu_pm+mu_pv),
                               mean_mu = mu_pm,
                               var_mu = mu_pv,
                               mean_latent_smooth = Eb,
                               Var_latent_smooth = Eb2-Eb^2),
                fitted_g = list(sigma2=sigma2,sigma2_trace=sigma2_trace),
                run_time = difftime(t_end,t_start,units='secs')))
  }
}

pois_smooth_split_obj = function(x,s,m,s2,Eb,Eb2,sigma2,KLb,const){
  return(sum(x*m-s*exp(m+s2/2)+log(s2)/2-log(sigma2)/2-(m^2+s2-2*m*Eb+Eb2)/2/sigma2)+KLb-const)
}

