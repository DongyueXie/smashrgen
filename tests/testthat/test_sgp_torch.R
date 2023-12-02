set.seed(1)
library(mvnfast)
library(smashrgen)
n = 1000
#m = 10
X = seq(0,1,length.out = n)
sigma2 = 0.5
#X_ind = seq(0,1,length.out = m)
kernel_param = c(2,-1)
K = Maternkernel(X,X,kernel_param,eps=1e-5,nu=3/2)
K_chol = Rfast::cholesky(K)
mu = 2
f = rmvn(1,rep(mu,n),K_chol,isChol = T)
s2 = runif(n,0.5,2)
#s2 = 1
y = drop(f + rnorm(n,0,sqrt(sigma2*s2)))


mod_torch = sgp_torch(y,fix_X_ind = T,m=100,s2=s2)
mod_reg = sgp(y,fix_X_ind = T,m=100,s2=s2)

mod_torch$fitted_g$kernel_param
mod_reg$fitted_g$kernel_param

mod_torch$fitted_g$mu
mod_reg$fitted_g$mu

mod_torch$fitted_g$sigma2
mod_reg$fitted_g$sigma2

mod_reg$elbo
mod_torch$elbo

mod_reg$fitted_g$X_ind
mod_torch$fitted_g$X_ind

mod_torch$opt_res$counts
mod_reg$opt_res$counts

mod_reg$run_time
mod_torch$run_time

plot(y,pch=20,col='grey70')
lines(c(f))
lines(mod_reg$posterior$mean,col=2,lwd=2)
lines(mod_torch$posterior$mean,col=3,lwd=2)

all_params = torch_tensor(c(mod_reg$fitted_g$kernel_param,log(mod_reg$fitted_g$sigma2)),dtype = torch_float64())
arg = list(y=torch_tensor(y,dtype = torch_float64()),
           X = torch_tensor(X,dtype = torch_float64()),
           X_ind = torch_tensor(mod_reg$fitted_g$X_ind,dtype = torch_float64()),
           kernel_func=Maternkernel_torch,
           mu=NULL,
           s2=torch_unsqueeze(torch_tensor(s2,dtype = torch_float64()),dim=2),
           sigma2 = torch_tensor(mod_reg$fitted_g$sigma2,dtype = torch_float64()),
           fix_sigma2=F,
           Jitter=torch_tensor(1e-5,dtype = torch_float64()),
           d_mat_nm=NULL,
           d_mat_mm=NULL)
sgp_obj_torch(all_params,arg)

X_ind=mod_reg$fitted_g$X_ind
kernel_param = mod_reg$fitted_g$kernel_param
log_sigma2 = log(mod_reg$fitted_g$sigma2)
mu=NULL
s2 = 1
kernel_func = Maternkernel
Jitter=  1e-5
sgp_obj(X_ind,kernel_param,log_sigma2,mu,s2,y,X,kernel_func,Jitter)
d_mat_nm=NULL
d_mat_mm=NULL




X = NULL
X_ind=NULL
m=10
kernel_func=Maternkernel_torch
kernel_param=c(0,0)
mu=NULL
s2=NULL
sigma2=NULL
opt_method='L-BFGS-B'
fix_sigma2=F
fix_X_ind = T
Jitter=1e-5
verbose=FALSE
