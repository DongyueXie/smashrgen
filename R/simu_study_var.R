#' Function to do comparisons among methods

#' @param est: 's' or 'sst'. 's':sigma only. 'sst': sigma^2+s_t^2
#' @param mu: true mean
#' @param sigma:tue sigma
#' @return a list of values
#' @export

simu_study_var=function(mu,st,sigma,niter=100,seed=12345){
  set.seed(seed)
  n=length(mu)

  ols.sst=c()
  moment.sst=c()
  smashc.sst=c()
  smashu.sst=c()
  rmad.sst=c()

  ols.s=c()
  moment.s=c()
  smashc.s=c()
  #smashu.s=c()
  #rmad.s=c()

  for(iter in 1:niter){
    y=mu+rnorm(n,0,sigma)+rnorm(n,0,st)
    ols=sigma_est(y,mu,st,'ols')
    ols.s[iter]=ols
    ols.sst[iter]=norm(sqrt(ols^2+st^2)-sqrt(sigma^2+st^2),'2')
    moment=sigma_est(y,mu,st,'moment')
    moment.s[iter]=moment
    moment.sst[iter]=norm(sqrt(moment^2+st^2)-sqrt(sigma^2+st^2),'2')
    smashc=sigma_est(y,mu,st,'smashc')
    smashc.s[iter]=smashc
    smashc.sst[iter]=norm(sqrt(smashc^2+st^2)-sqrt(sigma^2+st^2),'2')
    smashu=sst_est(y,'smashu')
    #smashu.s[iter]= ifelse(smashu^2-st^2)
    smashu.sst[iter]=norm(smashu-sqrt(sigma^2+st^2),'2')
    rmad=sst_est(y,'rmad')
    rmad.sst[iter]=norm(rmad-sqrt(sigma^2+st^2),'2')
  }
  return(list(ols.sst=ols.sst,moment.sst=moment.sst,smashc.sst=smashc.sst,
              smashu.sst=smashu.sst,rmad.sst=rmad.sst,ols.s=ols.s,moment.s=moment.s,
              smashc.s=smashc.s))
}
