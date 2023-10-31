#'@title Sparse Gaussian Process via VFE

sgp = function(X,y,X_ind,
               kernel_func=RBFkernel,
               kernel_param=c(0,0),
               mu=NULL,
               s2=NULL,sigma2=1,
               opt_method='L-BFGS-B',
               fix_x_ind=F,fix_kernel_param=F,fix_sigma2=F,fix_mu=F,l_b=-Inf,r_b=Inf){
  t0 = Sys.time()
  m = length(X_ind)
  n_kernel_param = length(kernel_param)
  n = length(y)
  if(is.null(s2)){s2 = rep(1,n)}
  init_val = c(X_ind,kernel_param,log(sigma2))
  opt_res = optim(par=init_val,fn=sgp_obj_for_optim,y=y,X=X,s2=s2,m=m,n_kernel_param=n_kernel_param,kernel_func=kernel_func,
                  X_ind=X_ind,kernel_param=kernel_param,mu=mu,log_sigma2=log(sigma2),
                  fix_kernel_param=fix_kernel_param,fix_sigma2=fix_sigma2,fix_x_ind=fix_x_ind,fix_mu=fix_mu,method = opt_method,lower=l_b,upper=r_b)
  params = opt_res$par
  if(!fix_mu){
    mu = NULL
  }
  if(!fix_x_ind){
    X_ind = params[1:m]
  }
  if(!fix_kernel_param){
    kernel_param = params[(m+1):(m+n_kernel_param)]
  }
  if(!fix_sigma2){
    log_sigma2 = params[length(params)]
  }else{
    log_sigma2 = log(sigma2)
  }
  elbo = -opt_res$value
  post = sgp_get_posterior(X_ind,kernel_param,log_sigma2,mu,s2,y,X,kernel_func)
  return(list(elbo=elbo,
              posterior=list(mean=post$mean,sd=sqrt(post$v)),
              fitted_g=list(mu=post$mu,X_ind=X_ind,kernel_param=kernel_param,sigma2=exp(log_sigma2)),
              run_time = difftime(Sys.time(),t0)))
  # opt_res = optim(init_val,sgp_obj_opt_kernel,
  #                 X_ind=X_ind,mu=0,
  #                 sigma2=sigma2,
  #                 y=y,X=X,s2=s2,
  #                 kernel_func=kernel_func,
  #                 method = opt_method)
}

# sgp_obj_opt_kernel = function(kernel_param,y,X,X_ind,mu,sigma2,s2,kernel_func){
#   log_sigma2 = log(sigma2)
#   obj = sgp_obj(X_ind,kernel_param,log_sigma2,mu,s2,y,X,kernel_func)
#   print(obj)
#   print(kernel_param)
#   return(obj)
# }

sgp_obj_for_optim = function(params,y,X,s2,m,n_kernel_param,kernel_func,
                             X_ind,kernel_param,mu,log_sigma2,fix_kernel_param,fix_sigma2,fix_x_ind,fix_mu){
  if(!fix_mu){
    mu = NULL
  }
  if(!fix_x_ind){
    X_ind = params[1:m]
  }
  if(!fix_kernel_param){
    kernel_param = params[(m+1):(m+n_kernel_param)]
  }
  if(!fix_sigma2){
    log_sigma2 = params[length(params)]
  }
  obj = sgp_obj(X_ind,kernel_param,log_sigma2,mu,s2,y,X,kernel_func)
  #print(obj)
  #print(params)
  return(obj)
}

sgp_get_posterior = function(X_ind,kernel_param,log_sigma2,mu,s2,y,X,kernel_func){
  n = length(y)
  sigma2 = exp(log_sigma2)
  Knm = kernel_func(X,X_ind,kernel_param)
  Kmm = kernel_func(X_ind,X_ind,kernel_param,eps=1e-5)
  Knn_diag = kernel_func(X,X,kernel_param,diagonal=T)
  L_Kmm = chol(Kmm)
  L_Kmm_inv = backsolve(L_Kmm,diag(ncol(Kmm)))
  Kmm_inv = tcrossprod(L_Kmm_inv)
  A = Knm%*%Kmm_inv
  V = solve(crossprod(A/s2,A)/sigma2 + Kmm_inv)
  L_V = t(chol(V))
  if(is.null(mu)){
    beta = solve(crossprod(A/s2,A)/sigma2 + Kmm_inv- tcrossprod(colSums(A/s2))/sum(1/s2)/sigma2)%*%(t(A/s2)%*%y-sum(y/s2)/sum(1/s2)*colSums(A/s2))/sigma2
    Ab = A%*%beta
    mu = (sum(y/s2)-sum(Ab/s2))/sum(1/s2)
  }else{
    # beta = solve(crossprod(A/s2,A)/sigma2 + Kmm_inv)%*%(t(A/s2)%*%y-mu*colSums(A/s2))/sigma2
    beta = V%*%(t(A/s2)%*%y-mu*colSums(A/s2))/sigma2
    Ab = A%*%beta
  }
  #beta = solve(crossprod(A/s2,A)/sigma2 + Kmm_inv- tcrossprod(colSums(A/s2))/sum(1/s2)/sigma2)%*%(t(A/s2)%*%y-sum(y/s2)/sum(1/s2)*colSums(A/s2))/sigma2
  #beta = Kmm%*%solve(Kmm+crossprod(Knm)/sigma2)%*%t(Knm)%*%y/sigma2
  # beta = solve(crossprod(A/s2,A)/sigma2 + Kmm_inv)%*%(t(A/s2)%*%y-mu*colSums(A/s2))/sigma2
  # V = solve(crossprod(A/s2,A)/sigma2 + Kmm_inv)
  # Ab = A%*%beta
  v = Knn_diag - rowSums((Knm%*%L_Kmm_inv)^2)  + rowSums((A%*%L_V)^2)
  return(list(mean=Ab,var=v,mu=mu))
}


sgp_obj2 = function(X_ind,kernel_param,log_sigma2,mu,s2,y,X,kernel_func){
  n = length(y)
  sigma2 = exp(log_sigma2)
  Knm = kernel_func(X,X_ind,kernel_param)
  Kmm = kernel_func(X_ind,X_ind,kernel_param,eps=1e-5)
  Knn_diag = kernel_func(X,X,kernel_param,diagonal=T)
  Kmm_inv = solve(Kmm)
  temp = Knm%*%Kmm_inv%*%t(Knm)
  dmvnorm(y,rep(0,n),diag(sigma2,n)+temp,log = T) - sum(Knn_diag)/2/sigma2 + sum(diag(temp))/2/sigma2
}

sgp_obj = function(X_ind,kernel_param,log_sigma2,mu,s2,y,X,kernel_func){
  n = length(y)
  sigma2 = exp(log_sigma2)
  Knm = kernel_func(X,X_ind,kernel_param)
  Kmm = kernel_func(X_ind,X_ind,kernel_param,eps=1e-5)
  Knn_diag = kernel_func(X,X,kernel_param,diagonal=T)
  L_Kmm = t(chol(Kmm))
  L_Kmm_inv = forwardsolve(L_Kmm,diag(ncol(Kmm)))
  #L_Kmm_inv = solve(L_Kmm)
  Kmm_inv = crossprod(L_Kmm_inv)
  #Kmm_inv = solve(Kmm)
  A = Knm%*%Kmm_inv
  V = solve(crossprod(A/s2,A)/sigma2 + Kmm_inv)
  L_V = t(chol(V))
  if(is.null(mu)){
    beta = solve(crossprod(A/s2,A)/sigma2 + Kmm_inv- tcrossprod(colSums(A/s2))/sum(1/s2)/sigma2)%*%(t(A/s2)%*%y-sum(y/s2)/sum(1/s2)*colSums(A/s2))/sigma2
    Ab = A%*%beta
    mu = (sum(y/s2)-sum(Ab/s2))/sum(1/s2)
  }else{
    # beta = solve(crossprod(A/s2,A)/sigma2 + Kmm_inv)%*%(t(A/s2)%*%y-mu*colSums(A/s2))/sigma2
    beta = V%*%(t(A/s2)%*%y-mu*colSums(A/s2))/sigma2
    Ab = A%*%beta
  }
  # obj = -n/2*log(2*pi*sigma2) - sum(log(s2))/2 - (sum(y^2/s2)-2*mu*sum(y/s2)-2*sum(y*Ab/s2)+mu^2*sum(1/s2)+2*mu*sum(Ab/s2)+sum(Ab^2/s2)+sum(diag(t(A/s2)%*%A%*%V))+sum(Knn_diag/s2)-sum(diag((A/s2)%*%t(Knm))))/2/sigma2 - as.numeric(determinant(Kmm,logarithm=T)$modulus)/2 - t(beta)%*%Kmm_inv%*%beta/2 - sum(diag(Kmm_inv%*%V)/2)+as.numeric(determinant(V,logarithm = T)$modulus)/2 + n*0.5
  obj = -n/2*log(2*pi*sigma2) - sum(log(s2))/2 - (sum(y^2/s2)-2*mu*sum(y/s2)-2*sum(y*Ab/s2)+mu^2*sum(1/s2)+2*mu*sum(Ab/s2)+sum(Ab^2/s2)+sum(((A/sqrt(s2))%*%L_V)^2)+sum(Knn_diag/s2)-sum(((Knm/sqrt(s2))%*%t(L_Kmm_inv))^2))/2/sigma2 -sum(log(diag(L_Kmm)))/2 - sum((L_Kmm_inv%*%beta)^2)/2 - sum((L_Kmm_inv%*%L_V)^2)/2+sum(log(diag(L_V)))/2 + n*0.5
  return(-drop(obj))
}

RBFkernel = function(X1,X2,kernel_param,eps=NULL,diagonal=F){
  coeff = exp(kernel_param[1])
  scales = exp(kernel_param[2])
  #X1 = cbind(X1)
  #X2 = cbind(X2)
  if(diagonal){
    return(rep(coeff^2,length(X1)))
  }else{
    #n_x <- nrow(X1)
    #n_y <- nrow(X2)
    #sqdist_xy <- (as.matrix(dist(rbind(X1, X2), method = "euclidean"))^2)[1:n_x, (n_x + 1):(n_x + n_y)]
    sqdist_xy <- outer(X1,X2,FUN="-")^2
    if(!is.null(eps)){
      K = coeff^2*exp(-0.5*sqdist_xy/scales^2)+diag(eps,length(X1))
    }else{
      K = coeff^2*exp(-0.5*sqdist_xy/scales^2)
    }
    return(K)
  }
}

