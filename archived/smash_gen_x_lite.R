#' A light version for smashgen with covariate X, to make it fast.
#' @export
#'
smash_gen_x_lite=function(y,X,sigma=NULL,ntri=NULL,wave_family='DaubExPhase',
                        filter.number=1,
                        dist_family='poisson',
                        log_scale=FALSE,homo=F){
  ini=init_smash_gen(x=y,ashp=T,ntri=ntri,dist_family=dist_family)

  fit=smash.gaus.x(X,ini$y0,family = wave_family,filter.number = filter.number,homo=homo)
  mu.hat=fit$mu.hat
  beta.hat=fit$beta.hat
  if(log_scale&dist_family=='poisson'){
    mu.est=mu.hat
  }else{
    mu.est=update_smash_gen(mu.hat,x=y,ntri,dist_family)$mt
  }
  return(list(mu.est=mu.est,beta.est=beta.hat))
}
