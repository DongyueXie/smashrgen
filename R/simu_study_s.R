#' Test the sigma_est function(estimate sigma)
#' @return a data frame of sigma hat from each method
#' @export

simu_study_s=function(mu,st,sigma=1,var_est='default',seed=1234){
  wls=c()
  mle=c()
  moment=c()
  eb=c()
  huberh=c()

  n=length(mu)
  set.seed(seed)
  for(i in 1:100){
    y=mu+rnorm(n,0,sigma)+rnorm(n,0,st)
    wls[i]=sigma_est(y,mu,st,'wls',var_est)
    eb[i]=sigma_est(y,mu,st,'eb',var_est)
    moment[i]=sigma_est(y,mu,st,'moment',var_est)
    mle[i]=sigma_est(y,mu,st,'mle',var_est)
    huberh[i]=sigma_est(y,mu,st,'huber',var_est)
  }
  return(data.frame(wls=wls,mle=mle,eb=eb,moment=moment,huberm=huberh))
}
