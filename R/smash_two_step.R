
#'@title A two step method for smoothing Poisson seq
#'@description This method first fit a smash poisson to the seq, then fit smash.gaus to the log estimation.
#'@param x  data vector
#'@import smashr
#'@export

smash_two_step = function(x,homoskedastic = FALSE){
  t_start = Sys.time()
  lx = smash.poiss(x,log=TRUE)
  fit = smash.gaus(lx,v.est=T,joint=T,homoskedastic=homoskedastic)
  t_end = Sys.time()
  return(list(posterior=list(mean_smooth = exp(fit$mu.res),
                             mean_latent_smooth = fit$mu.res),
              fitted_g = list(sigma2 = fit$var.res),
              run_time = difftime(t_end,t_start,units='secs')))
}
