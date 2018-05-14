#' function to estimate both mu and \sigma

#' @param x:data
#' @param st: known variance
#' @param family: wavelet basis
#' @param niters: number of iterations to estimate \sigma and mu
#' @param method: 'smashc', 'ols' amd 'moment'
#' @return estimated mean and \sigma
#' @export

smash.gaus.gen=function(x,st,family='DaubExPhase',niters=2,method='smashc'){
  #initialize \sigma^2
  sigma0=sigma_est(x,st=st,method = 'moment')
  #sd0=sqrt(sigma0^2+st^2)
  sd.est=c()
  mu.est=c()
  sd.est=c(sd.est,sigma0)
  for(iter in 1:niters){
    #estimate mu given sd
    mu.hat=smash.gaus(x,sigma=sqrt(sd.est[iter]^2+st^2),family=family)
    mu.est=rbind(mu.est,mu.hat)
    #estimate sd given mu
    sd.hat=sigma_est(x,mu.est[iter,],st,method=method)
    sd.est=c(sd.est,sd.hat)
  }
  #mu.hat=smash.gaus(x,sigma=sqrt(sd.hat^2+st^2),family=family)
  return(list(mu.hat=mu.hat,sd.hat=sd.hat))
}
