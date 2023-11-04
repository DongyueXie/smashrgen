library(mvnfast)
set.seed(2023)
n = 100
m = 10
X = seq(-5,5,length.out = n)
sigma2 = 0.5
X_ind = seq(-5,5,length.out = m)
kernel_param = c(0,0)
#K = RBFkernel(X,X,kernel_param,eps=1e-5)
K = Maternkernel(X,X,kernel_param,eps=1e-5,nu=3/2)
mu = 0
f = rmvn(1,rep(mu,n),K)
s2 = runif(n,0.5,2)
y = drop(f + rnorm(n,0,sqrt(sigma2*s2)))


# True posterior, if X_ind == X
res = sgp_get_posterior(X,kernel_param,log(sigma2),mu,s2=1,y,X,kernel_func=RBFkernel)
plot(X,y,pch=20,col='grey70')
lines(X,f)
lines(X,res$mean,col=4,lwd=2)
legend("topright",c('y',"f","f_hat"),lty=c(NA,1,1),pch=c(20,NA,NA),col=c('grey70',1,4))

# Posterior with X_ind
res = sgp_get_posterior(X_ind,kernel_param,log(sigma2),mu,s2=1,y,X,kernel_func=RBFkernel)
plot(X,y,pch=20,col='grey70')
lines(X,f)
lines(X,res$mean,col=4,lwd=2)
legend("topright",c('y',"f","f_hat"),lty=c(NA,1,1),pch=c(20,NA,NA),col=c('grey70',1,4))

# Try optimization

res = sgp(y,kernel_func = RBFkernel)
plot(y,pch=20,col='grey70')
lines(c(f))
lines(res$posterior$mean+res$fitted_g$mu,col=4,lwd=2)
lines(res$posterior$mean+res$fitted_g$mu+2*res$posterior$sd,col=2,lty=2)
lines(res$posterior$mean+res$fitted_g$mu-2*res$posterior$sd,col=2,lty=2)
legend("topright",c('y',"f","f_hat"),lty=c(NA,1,1),pch=c(20,NA,NA),col=c('grey70',1,4))

res = sgp(y,fix_x_ind = T,fix_sigma2 = T,sigma2=sigma2,s2=s2,m=100)
plot(X,y,pch=20,col='grey70')
lines(X,f)
lines(X,res$posterior$mean+res$fitted_g$mu,col=4,lwd=2)
legend("topright",c('y',"f","f_hat"),lty=c(NA,1,1),pch=c(20,NA,NA),col=c('grey70',1,4))

res = ebnm_sgp(y,sqrt(sigma2*s2))
plot(y,pch=20,col='grey70')
lines(c(f))
lines(res$posterior$mean,col=4,lwd=2)
legend("topright",c('y',"f","f_hat"),lty=c(NA,1,1),pch=c(20,NA,NA),col=c('grey70',1,4))
# s2=NULL
# opt_method='L-BFGS-B'
# fix_x_ind=F
# fix_kernel_param=F
# fix_sigma2=F
# fix_mu=F
#
# m = length(X_ind)
# n_kernel_param = length(kernel_param)
# n = length(y)
# if(is.null(s2)){s2 = rep(1,n)}
# init_val = c(X_ind,kernel_param,log(sigma2))
# kernel_func = RBFkernel
#
# opt_res = optim(init_val,sgp_obj_for_optim,method = opt_method,y=y,X=X,s2=s2,m=m,n_kernel_param=n_kernel_param,kernel_func=kernel_func,
#                 X_ind=X_ind,kernel_param=kernel_param,mu=mu,log_sigma2=log(sigma2),
#                 fix_kernel_param=fix_kernel_param,fix_sigma2=fix_sigma2,fix_x_ind=fix_x_ind,fix_mu=fix_mu)
# params = opt_res$par
# if(!fix_mu){
#   mu = NULL
# }
# if(!fix_x_ind){
#   X_ind = params[1:m]
# }
# if(!fix_kernel_param){
#   kernel_param = params[(m+1):(m+n_kernel_param)]
# }
# if(!fix_sigma2){
#   log_sigma2 = params[length(params)]
# }else{
#   log_sigma2 = log(sigma2)
# }
# elbo = -opt_res$value
# post = sgp_get_posterior(X_ind,kernel_param,log_sigma2,mu,s2,y,X,kernel_func)
# plot(X,y,pch=20,col='grey70')
# lines(X,f)
# lines(X,post$mean+post$mu,col=4,lwd=2)
# legend("topright",c('y',"f","f_hat"),lty=c(NA,1,1),pch=c(20,NA,NA),col=c('grey70',1,4))


###### test speed

library(mvtnorm)
set.seed(2023)
n_list = c(100,500,1000,2000,3000,5000,10000)
#n = 1000
m = 10
sigma2 = 0.5
X_ind = seq(-5,5,length.out = m)
kernel_param = c(0,0)
mu = 0
res_list = list()
for(n in n_list){
  X = seq(-5,5,length.out = n)
  K = RBFkernel(X,X,kernel_param)
  f = rmvnorm(1,rep(mu,n),K)
  rm(K)
  gc()
  y = drop(f + rnorm(n,0,sqrt(sigma2)))
  res = sgp(X,y,X_ind)
  res_list$n = res
}




