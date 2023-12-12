sin_func = function(n){
  x = seq(-0.5,0.5,length.out=n)
  sin(2*pi*x)
}
step_func = function(n){
  x = seq(0,1,length.out=n)
  y = rep(0,n)
  y[x>0.2&x<0.4] = 1
  y[x>0.6&x<0.8] = 3
  y
}
library(smashrgen)
library(smashr)
library(torch)
simu_func = function(n_list,signal_func,snr=3,sigma2=1,n_rep=10,m_list=c(10,30,50,100,150),seed=12345){
  set.seed(seed)
  res = list()
  for(n in n_list){
    res_n = list()
    for(rep in 1:n_rep){
      print(paste('running n=',n,', rep=',rep, sep=''))
      # generate data
      mu = signal_func(n)
      mu = mu*sqrt(snr/var(mu))
      y = mu + rnorm(n,0,sqrt(sigma2))
      # fit smashr
      t0_smash = Sys.time()
      fit_smash = smash.gaus(y,sqrt(sigma2))
      t_smash = Sys.time() - t0_smash
      res_smash = list(posterior=list(mean=fit_smash),run_time = t_smash)
      # fit gp
      res_sgp = list()
      res_sgp_torch=list()
      for(m in m_list){
        fit_sgp = sgp(y,m=m,sigma2=sigma2,fix_sigma2 = T,fix_X_ind = T)
        fit_sgp_torch = sgp_torch(y,m=m,sigma2=sigma2,fix_sigma2 = T,fix_X_ind = T)
        res_sgp[[paste0('m',m,sep='')]] = fit_sgp
        res_sgp_torch[[paste0('m',m,sep='')]] = fit_sgp_torch
      }
      res_n[[paste0('rep_',rep,sep='')]] = list(res_smash=res_smash,res_sgp=res_sgp,res_sgp_torch=res_sgp_torch,data=list(mu=mu,y=y))
    }
    res[[paste0('n_',n,sep='')]] = res_n
  }
  res
}

#test = simu_func(c(200,210),sin_func,n_rep=2,m_list=c(10,20))

output = simu_func(c(300,500,1000,3000,5000,10000),step_func)
