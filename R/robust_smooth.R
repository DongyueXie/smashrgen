#'@title empirical Bayes robust smoothing
#'@import ebnm

robust_smooth = function(y,sigma2,
                         mu_init=NULL,
                         ebnm_control = list(),
                         tol=1e-5,maxiter=100){
  ebnm_controls = robust_smooth_ebnm_control_default()
  ebnm_controls = modifyList(ebnm_controls,ebnm_control,keep.null = TRUE)
  mu = mu_init
  if(is.null(mu)){
    mu = rep(mean(y),length(y))
  }
  mu_old = mu
  for(i in 1:maxiter){
    outlier_res = ebnm(y-mu,sqrt(sigma2),
               prior_family = ebnm_controls$prior_family,
               mode = ebnm_controls$mode,
               scale = ebnm_controls$scale,
               g_init = ebnm_controls$g_init,
               fix_g = ebnm_controls$fix_g)
    smooth_res = smash_dwt(y-outlier_res$posterior$mean,sqrt(sigma2))
    mu = smooth_res$posterior$mean
    if(norm(mu-mu_old,'2') < tol){
      break
    }
    mu_old = mu
  }
  return(list(smooth_res=smooth_res,outlier_res=outlier_res))
}

robust_smooth_ebnm_control_default = function(){
  return(list(prior_family = "point_laplace",
              mode = 0,
              scale = "estimate",
              g_init = NULL,
              fix_g = FALSE))
}
