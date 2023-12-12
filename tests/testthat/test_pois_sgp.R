library(mvnfast)
library(smashrgen)
set.seed(12)
n = 1000
X = seq(0,1,length.out = n)
# X_ind = seq(-5,5,length.out = m)
kernel_param = c(2, -2)
K = Maternkernel(X,X,kernel_param)
mu = 0
f = drop(rmvn(1,rep(mu,n),K))
f = f - min(f)
f = f*20
# f = c(rep(0.1,n/5),rep(20,n/5),rep(0.1,n/5),rep(100,n/5),rep(0.1,n/5))
# f = f
#s = runif(n,0,5)
s=1
sigma2=0.0
y = rpois(length(f),s*exp(log(f)+rnorm(n,0,sqrt(sigma2))))
plot(y,pch=20,col='grey70')
lines(f)

plot(y/s,pch=20,col='grey70')
lines(f)

res = ebpm_pois_sgp(y,s)

res = pois_sgp(y,sc=s,verbose = T,m=100)
# plot(res$elbo_tace)
plot(y,pch=20,col='grey70')
lines(f,col='grey60')
lines(res$posterior$rate,col=2,lwd=2)
legend("topright",c('y',"f","f_hat"),lty=c(NA,1,1),pch=c(20,NA,NA),col=c('grey70',1,4))
lines(res$posterior$rate+2*res$posterior$rate_sd,col=2,lty=2)
lines(res$posterior$rate-2*res$posterior$rate_sd,col=2,lty=2)
res$fitted_g

#11/23/2023

ly = log(1+y)
plot(ly,pch=20,col='grey70')
fit = sgp(ly,fix_x_ind = T,m=50,verbose = T)
lines(fit$posterior$mean,col=2,lwd=2)
fit$fitted_g$mu
fit$fitted_g$X_ind
fit$fitted_g$kernel_param
fit$fitted_g$sigma2



gp = gp_init(cf_sexp(),lik_poisson(),method = method_fitc())
gp <- gp_optim(gp, X, y,restarts = 10)
pred = gp_pred(gp,X,transform = T)
plot(y,pch=20,col='grey70')
lines(f)
lines(pred$mean,col=4,lwd=2)
lines(res$posterior$rate,col=2,lwd=2)

fit_sgp = sgp(y,m=100,fix_x_ind = T,opt_method = "Nelder-Mead")
lines(fit_sgp$posterior$mean+fit_sgp$fitted_g$mu,col=2)

s = rep(1,n)
kernel_func = Maternkernel
post_mean = log(f)[seq(1,1000,length.out=m)]
V = kernel_func(X[seq(1,1000,length.out=m)],X[seq(1,1000,length.out=m)],kernel_param,eps=0.1)

# fix_kernel_param=F
# fix_X_ind=F
# fix_mu=F
#
# maxiter=100
# tol=1e-5
# maxiter_mean=100
# tol_mean=1e-5
# maxiter_V=100
# tol_V=1e-5

l_b = c(rep(-5,m),-5,-5)
r_b = c(rep(5,m),5,5)

res = pois_sgp(y,X,X_ind,length(X_ind),
               post_mean,V,
               kernel_func=Maternkernel,
               kernel_param=c(0,0),
               mu=0,
               s=NULL,
               opt_method='L-BFGS-B',
               fix_X_ind=F,fix_kernel_param=F,fix_mu=F,
               maxiter=100,tol=1e-5,
               maxiter_mean=100,tol_mean=1e-5,
               maxiter_V=100,tol_V=1e-5,
               l_b=l_b,r_b=r_b,
               verbose=T,printevery=1)

plot(res$elbo_tace)
plot(X,y,pch=20,col='grey70')
lines(X,f)
lines(X,res$posterior$rate,col=4,lwd=2)
legend("topright",c('y',"f","f_hat"),lty=c(NA,1,1),pch=c(20,NA,NA),col=c('grey70',1,4))
lines(X,res$posterior$rate+2*res$posterior$rate_sd,col=2,lty=2)
lines(X,res$posterior$rate-2*res$posterior$rate_sd,col=2,lty=2)

# pois_sgp_elbo(X,y,s,X_ind,kernel_param,kernel_func,post_mean,V,mu)
# grad_hess = grad_post_mean(X,y,s,X_ind,kernel_param,kernel_func,post_mean,V,mu,return_grad=T,return_hess=T)
# grad = grad_hess$grad
# hess =  grad_hess$hess
# direction = c(pcgsolve(hess,post_mean))
#
# solve(-grad_post_mean(X,y,s,X_ind,kernel_param,kernel_func,post_mean,V,mu,return_grad=F,return_hess=T))



pois_sgp_elbo_for_post_mean = function(post_mean,X,y,s,X_ind,kernel_param,kernel_func,V,mu){
  n = length(y)
  m = length(X_ind)
  Knm = kernel_func(X,X_ind,kernel_param)
  Kmm = kernel_func(X_ind,X_ind,kernel_param,eps=1e-5)
  Knn_diag = kernel_func(X,X,kernel_param,diagonal=T)
  L_Kmm = t(chol(Kmm))
  L_Kmm_inv = forwardsolve(L_Kmm,diag(ncol(Kmm)))
  Kmm_inv = crossprod(L_Kmm_inv)
  A = Knm%*%Kmm_inv
  L_V = t(chol(V))
  Ab = drop(A%*%post_mean)
  obj = sum(y*Ab)+mu*sum(y)-exp(mu)*sum(s*exp(Ab+rowSums((A%*%L_V)^2)/2+Knn_diag/2-rowSums((Knm%*%t(L_Kmm_inv))^2)/2))+sum(y*log(s))-sum(lfactorial(y))-sum(log(diag(L_Kmm)))/2 - sum((L_Kmm_inv%*%post_mean)^2)/2 - sum((L_Kmm_inv%*%L_V)^2)/2+sum(log(diag(L_V)))/2-m/2*log(2*pi)+m/2
  return(-drop(obj))
}

pois_sgp_update_post_mean(X,y,s,X_ind,kernel_param,kernel_func,post_mean_new,V,mu,maxiter=maxiter_mean,tol=tol_mean)
optim(post_mean,pois_sgp_elbo_for_post_mean,method="L-BFGS-B",X=X,y=y,s=s,X_ind=X_ind,kernel_param=kernel_param,kernel_func=kernel_func,V=V,mu=mu)
post_mean_new = c(1.939839, 2.283306, 1.318693, 1.264437, 3.122596, 3.243795, 2.343830, 2.412892, 0.867362, 1.626417)
grad_post_mean(X,y,s,X_ind,kernel_param,kernel_func,post_mean_new,V,mu,return_grad=T,return_hess=T)

pois_sgp_elbo(X,y,s,X_ind,kernel_param,kernel_func,post_mean,V,mu)

