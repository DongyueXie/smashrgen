#' Function to do comparisons among methods

#' @param est: 's' or 'sst'. 's':sigma only. 'sst': sigma^2+s_t^2
#' @param mu: true mean
#' @param sigma:tue sigma
#' @return a list of values
#' @export

simu_study_var=function(mu,st,sigma,niter=100,seed=12345,var_est='smash'){
  set.seed(seed)
  n=length(mu)

  mle.sst=c()
  moment.sst=c()
  eb.sst=c()
  wls.sst=c()
  huber.sst=c()
  smashu.sst=c()
  rmad.sst=c()

  mle.s=c()
  moment.s=c()
  eb.s=c()
  huber.s=c()
  wls.s=c()


  for(iter in 1:niter){
    y=mu+rnorm(n,0,sigma)+rnorm(n,0,st)
    #mle, moment,eb,huber,wls
    mle=sigma_est(y,mu,st,'mle',var_est)
    mle.s[iter]=mle
    mle.sst[iter]=mse(sqrt(mle^2+st^2),sqrt(sigma^2+st^2))
    moment=sigma_est(y,mu,st,'moment',var_est)
    moment.s[iter]=moment
    moment.sst[iter]=mse(sqrt(moment^2+st^2),sqrt(sigma^2+st^2))
    eb=sigma_est(y,mu,st,'eb',var_est)
    eb.s[iter]=eb
    eb.sst[iter]=mse(sqrt(eb^2+st^2),sqrt(sigma^2+st^2))
    huberh=sigma_est(y,mu,st,'huber',var_est)
    huber.s[iter]=huberh
    huber.sst[iter]=mse(sqrt(huberh^2+st^2),sqrt(sigma^2+st^2))
    wls=sigma_est(y,mu,st,'wls',var_est)
    wls.s[iter]=wls
    wls.sst[iter]=mse(sqrt(wls^2+st^2),sqrt(sigma^2+st^2))
    ###############################
    smashu=sst_est(y,'smashu')
    #smashu.s[iter]= ifelse(smashu^2-st^2)
    smashu.sst[iter]=mse(smashu,sqrt(sigma^2+st^2))
    rmad=sst_est(y,'rmad')
    rmad.sst[iter]=mse(rmad,sqrt(sigma^2+st^2))
  }
  return(list(mle.sst=mle.sst,moment.sst=moment.sst,eb.sst=eb.sst,huber.sst=huber.sst,wls.sst=wls.sst,
              smashu.sst=smashu.sst,rmad.sst=rmad.sst,mle.s=mle.s,moment.s=moment.s,
              eb.s=eb.s,huber.s=huber.s,wls.s=wls.s))
}
