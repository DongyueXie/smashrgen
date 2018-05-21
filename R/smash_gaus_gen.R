#' function to estimate both mu and \sigma

#' @param x:data
#' @param st: known variance
#' @param family: wavelet basis
#' @param filter.number: as in smash.gaus
#' @param niters: number of iterations to estimate \sigma and mu
#' @param y_var_est: 'mle', 'moment', 'eb', 'huber','wls'
#' @param z_var_est: method to estimate variance: 'rmad', 'smash', 'default'
#' @param k: parameter in huber m estimator
#' @return estimated mean and \sigma
#' @export

smash.gaus.gen=function(x,st,y_var_est='mle',z_var_est='smash',family='DaubExPhase',
                        filter.number=1,niters=2,k=NULL){
  #initialize \sigma^2 using moment method
  sigma0=sigma_est(x,st=st,method = 'moment')
  #sd0=sqrt(sigma0^2+st^2)
  sd.est=c()
  mu.est=c()
  sd.est=c(sd.est,sigma0)
  for(iter in 1:niters){
    #estimate mu given sd
    mu.hat=smash.gaus(x,sigma=sqrt(sd.est[iter]^2+st^2),family=family,filter.number = filter.number)
    mu.est=rbind(mu.est,mu.hat)
    #estimate sd given mu
    sd.hat=sigma_est(x,mu.est[iter,],st,var_est=z_var_est,k=k,method=y_var_est,family=family,filter.number = filter.number)
    sd.est=c(sd.est,sd.hat)
  }
  #mu.hat=smash.gaus(x,sigma=sqrt(sd.hat^2+st^2),family=family)
  return(list(mu.hat=mu.hat,sd.hat=sd.hat))
}
