#'@import Rfast
grad_kernel_param = function(kernel_param,X,y,sc,X_ind,kernel_func,post_mean,V,mu,Jitter,d_mat_nm,d_mat_mm){
  m=length(X_ind)
  temp = pois_sgp_matrix_helper(X,y,sc,X_ind,kernel_param,kernel_func,post_mean,V,mu,Jitter,get_delta=T,d_mat_nm,d_mat_mm)
  grad_log_c1 = -sum(temp$delta*(temp$Knn_diag/2-Rfast::rowsums(temp$A*temp$Knm)/2))
  grad_log_c2 = - m/2 + sum((temp$L_Kmm_inv%*%post_mean)^2)/2+sum((temp$Kmm_inv*V))/2-mu*sum(temp$Kmm_inv%*%post_mean)+mu^2*sum(temp$Kmm_inv)/2

  g_Knm = kernel_func(X,X_ind,kernel_param,d_mat=d_mat_nm,calc_grad_scales=T)$g_scales
  g_Kmm = kernel_func(X_ind,X_ind,kernel_param,eps=Jitter,d_mat=d_mat_mm,calc_grad_scales=T)$g_scales
  Smmm = temp$Kmm_inv%*%g_Kmm%*%temp$Kmm_inv
  dA = g_Knm%*%temp$Kmm_inv-temp$Knm%*%Smmm
  dAb = dA%*%(post_mean-mu)
  Smmmb = Smmm%*%post_mean
  grad_log_theta1 = sum(y*dAb)-sum(temp$delta*(dAb+(rowsums(mat.mult(dA,V)*temp$A)+rowsums(mat.mult(temp$A,V)*dA)-rowsums(dA*temp$Knm)-rowsums(temp$A*g_Knm))/2))
  grad_log_theta2 = -sum(temp$Kmm_inv*g_Kmm)/2 + sum(post_mean*Smmmb)/2 + sum(Smmm*V)/2 + mu*sum(Smmmb) - mu^2*(sum(Smmm))/2
  return(c(grad_log_c1+grad_log_c2,grad_log_theta1+grad_log_theta2))
}

elbo_kernel_F1 = function(kernel_param2,kernel_param1,mu,X,y,sc,X_ind,kernel_func,post_mean,V,Jitter,d_mat_nm,d_mat_mm){
  m = length(X_ind)
  kernel_param=c(kernel_param1,kernel_param2)
  temp = pois_sgp_matrix_helper(X,y,sc,X_ind,kernel_param,kernel_func,post_mean,V,mu,Jitter,get_delta=T,d_mat_nm,d_mat_mm)
  F1 = sum(y*temp$Ab) +mu*sum(y*temp$tilde1) -sum(temp$delta)+sum(y*log(sc))-sum(lfactorial(y))
  F1
}
elbo_kernel_F2 = function(kernel_param2,kernel_param1,mu,X,y,sc,X_ind,kernel_func,post_mean,V,Jitter,d_mat_nm,d_mat_mm){
  m = length(X_ind)
  kernel_param=c(kernel_param1,kernel_param2)
  temp = pois_sgp_matrix_helper(X,y,sc,X_ind,kernel_param,kernel_func,post_mean,V,mu,Jitter,get_delta=T,d_mat_nm,d_mat_mm)
  F2 = -sum(log(diag(temp$L_Kmm))) - sum((temp$L_Kmm_inv%*%post_mean)^2)/2  - sum((temp$Kmm_inv*V))/2 +sum(log(diag(temp$L_V))) + m*0.5 + mu*sum(temp$Kmm_inv%*%post_mean) - mu^2*sum(temp$Kmm_inv)/2
  F2
}

set.seed(12345)
kernel_param=c(5,3)
n=100
m=10
X=seq(0,1,length.out=n)
y=0:(n-1)
sc=rep(1,n)
X_ind=seq(0,1,length.out=m)
kernel_func=Maternkernel
post_mean=rep(0.1,m)
#V=cov(matrix(rnorm(n*m),nrow=n))
V = diag(m)
mu=0.5
Jitter=1e-5
d_mat_nm=NULL
d_mat_mm=NULL
print(grad_kernel_param(kernel_param,X,y,sc,X_ind,kernel_func,post_mean,V,mu,Jitter,d_mat_nm,d_mat_mm))

# print(c(numDeriv::grad(elbo_kernel_F1,kernel_param[2],kernel_param1=kernel_param[1],mu=mu,X=X,y=y,sc=sc,X_ind=X_ind,kernel_func=kernel_func,post_mean=post_mean,V=V,Jitter=Jitter,d_mat_nm=d_mat_nm,d_mat_mm=d_mat_mm),
# numDeriv::grad(elbo_kernel_F2,kernel_param[2],kernel_param1=kernel_param[1],mu=mu,X=X,y=y,sc=sc,X_ind=X_ind,kernel_func=kernel_func,post_mean=post_mean,V=V,Jitter=Jitter,d_mat_nm=d_mat_nm,d_mat_mm=d_mat_mm),
print(numDeriv::grad(pois_sgp_update_kernel_param_obj,kernel_param,mu=mu,X=X,y=y,sc=sc,X_ind=X_ind,kernel_func=kernel_func,post_mean=post_mean,V=V,Jitter=Jitter,d_mat_nm=d_mat_nm,d_mat_mm=d_mat_mm))

#pois_sgp_update_kernel_param_obj(kernel_param,mu,X,y,sc,X_ind,kernel_func,post_mean,V,Jitter,d_mat_nm,d_mat_mm)


#
# m = length(X_ind)
# temp = pois_sgp_matrix_helper(X,y,sc,X_ind,kernel_param,kernel_func,post_mean,V,mu,Jitter,get_delta=T,d_mat_nm,d_mat_mm)
# F1 = sum(y*temp$Ab)+mu*sum(y*temp$tilde1)-sum(temp$delta)+sum(y*log(sc))-sum(lfactorial(y))
# F2 = -sum(log(diag(temp$L_Kmm)))/2 - sum((temp$L_Kmm_inv%*%post_mean)^2)/2 - sum((temp$L_Kmm_inv%*%temp$L_V)^2)/2+sum(log(diag(temp$L_V)))/2 + m*0.5 + mu*sum(temp$Kmm_inv%*%post_mean) - mu^2*sum(temp$Kmm_inv)/2
#
# F1
# F2
#
#
#
#
# tt = function(kernel_param1,kernel_param2,X,V){
#   kernel_param = c(kernel_param1,kernel_param2)
#   K = Maternkernel(X,X,kernel_param,eps=1e-5)
#   return(-sum(diag(solve(K)%*%V))/2)
# }
# numDeriv::grad(tt,0,kernel_param2=0,X=X_ind,V=V)
#
# tt(0,0,X_ind,V)
#
#
#
#
