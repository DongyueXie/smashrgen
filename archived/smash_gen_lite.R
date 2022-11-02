#' A light version for smashgen, to make it fast.
#' @export
smash_gen_lite=function(x,sigma=NULL,ntri=NULL,wave_family='DaubExPhase',
                   filter.number=1,
                   dist_family='poisson',
                   log_scale=FALSE){
  ini=init_smash_gen(x=x,ashp=T,ntri=ntri,dist_family=dist_family)
  if(is.null(sigma)){
    mu.hat=smash.gaus(ini$y0,family = wave_family,filter.number = filter.number)
  }else{
    mu.hat=smash.gaus(ini$y0,sigma = sqrt(sigma^2+ifelse(ini$s0<0,1e-8,ini$s0)),
                      family=wave_family,filter.number=filter.number)
  }

  if(log_scale&dist_family=='poisson'){
    mu.est=mu.hat
  }else{
    mu.est=update_smash_gen(mu.hat,x,ntri,dist_family)$mt
  }
  return(mu.est)
}
