#' Function to do comparisons among methods(total sd)
#' @param mu: true mean
#' @param sigma:tue sigma
#' @param mu.out: whether compare mu.est
#' @return a list of values
#' @export

simu_study_sst=function(mu,st,sigma,niter=100,seed=12345,var_est='rmad',mu.out=T){
  set.seed(seed)
  n=length(mu)

  mle.sst=c()
  moment.sst=c()
  eb.sst=c()
  wls.sst=c()
  huber.sst=c()
  smashu.sst=c()
  rmad.sst=c()
  if(mu.out){
    mle.mu=c()
    moment.mu=c()
    eb.mu=c()
    wls.mu=c()
    huber.mu=c()
    smashu.mu=c()
    rmad.mu=c()
    smashtrue.mu=c()
  }

  for(iter in 1:niter){
    y=mu+rnorm(n,0,sigma)+rnorm(n,0,st)
    #mle, moment,eb,huber,wls
    mle=sigma_est(y,mu,st,'mle',var_est)
    #mle.s[iter]=mle
    mle.sst[iter]=mse(sqrt(mle^2+st^2),sqrt(sigma^2+st^2))
    moment=sigma_est(y,mu,st,'moment',var_est)
    #moment.s[iter]=moment
    moment.sst[iter]=mse(sqrt(moment^2+st^2),sqrt(sigma^2+st^2))
    eb=sigma_est(y,mu,st,'eb',var_est)
    #eb.s[iter]=eb
    eb.sst[iter]=mse(sqrt(eb^2+st^2),sqrt(sigma^2+st^2))
    huberh=sigma_est(y,mu,st,'huber',var_est)
    #huber.s[iter]=huberh
    huber.sst[iter]=mse(sqrt(huberh^2+st^2),sqrt(sigma^2+st^2))
    wls=sigma_est(y,mu,st,'wls',var_est)
    #wls.s[iter]=wls
    wls.sst[iter]=mse(sqrt(wls^2+st^2),sqrt(sigma^2+st^2))
    ###############################
    smashu=sst_est(y,'smashu')
    #smashu.s[iter]= ifelse(smashu^2-st^2)
    smashu.sst[iter]=mse(smashu,sqrt(sigma^2+st^2))
    rmadh=sst_est(y,'rmad')
    rmad.sst[iter]=mse(rmadh,sqrt(sigma^2+st^2))
    if(mu.out){
      mle.mu[iter]=mse(smash.gaus(y,sigma=sqrt(mle^2+st^2)),mu)
      moment.mu[iter]=mse(smash.gaus(y,sigma=sqrt(moment^2+st^2)),mu)
      eb.mu[iter]=mse(smash.gaus(y,sigma=sqrt(eb^2+st^2)),mu)
      wls.mu[iter]=mse(smash.gaus(y,sigma=sqrt(wls^2+st^2)),mu)
      huber.mu[iter]=mse(smash.gaus(y,sigma=sqrt(huberh^2+st^2)),mu)
      smashu.mu[iter]=mse(smash.gaus(y,sigma=smashu),mu)
      rmad.mu[iter]=mse(smash.gaus(y,sigma=rmadh),mu)
      smashtrue.mu[iter]=mse(smash.gaus(y,sigma=sqrt(sigma^2+st^2)),mu)
    }
  }
  if(mu.out){
    return(list(var.est=data.frame(mle=mle.sst,moment=moment.sst,eb=eb.sst,huberm=huber.sst,wls=wls.sst,smashu=smashu.sst,rmad=rmad.sst),mu.est=data.frame(mle=mle.mu,moment=moment.mu,eb=eb.mu,huberm=huber.mu,wls=wls.mu,smashu=smashu.mu,rmad=rmad.mu,smashtrue=smashtrue.mu)))
  }else{
    return(data.frame(mle=mle.sst,moment=moment.sst,eb=eb.sst,huber=huber.sst,wls=wls.sst,smashu=smashu.sst,rmad=rmad.sst))
  }

}
