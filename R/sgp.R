#'@title Sparse Gaussian process in EBNM syntax
#'@param x A vector of observations.
#'@param s A vector of standard errors
#'@param g_init a list of X_ind,kernel_param,kernel_func,mu
#'@param fix_g If TRUE, the prior is fixed at g_init instead of estimated
#'@param output ignore this.
#'@return a list of posterior, log-likelihood, fitted_g
#'@export
ebnm_sgp = function(x,s,g_init=NULL,fix_g=FALSE,output){
  n = length(x)
  X = seq(0,1,length.out=n)
  if(is.null(g_init)){
    m=min(30,n)
    kernel_param = c(0,0)
    kernel_func = Maternkernel
    mu = NULL
    X_ind = seq(0,1,length.out=m)
  }else{
    m = length(g_init$X_ind)
    kernel_param = g_init$kernel_param
    kernel_func = g_init$kernel_func
    mu = g_init$mu
    X_ind = g_init$X_ind
  }
  if(!is.null(g_init) & fix_g){
    fix_x_ind=T
    fix_kernel_param=T
    fix_mu=T
  }else{
    fix_x_ind=T
    fix_kernel_param=F
    fix_mu=F
  }
  res = sgp(x,X,X_ind,m,
           kernel_func=kernel_func,
           kernel_param=kernel_param,
           mu=mu,
           s2=s^2,sigma2=1,
           opt_method='L-BFGS-B',
           fix_x_ind=fix_x_ind,fix_kernel_param=fix_kernel_param,
            fix_sigma2=T,fix_mu=fix_mu,l_b=-Inf,r_b=Inf)
  return(list(posterior=data.frame(mean=res$posterior$mean,second_moment=res$posterior$mean^2+res$posterior$sd^2),
              fitted_g=list(X_ind=res$fitted_g$X_ind,kernel_param=res$fitted_g$kernel_param,kernel_func=res$fitted_g$kernel_func,mu=res$fitted_g$mu),
              log_likelihood=res$elbo))
}


#'@title Fit Sparse Gaussian Process via variational inference
#'@param X,y training data X and response y.
#'@param X_ind inducing point locations.
#'@param m default number of X_ind is 10, if X_ind is not given as an input
#'@param kernel_func,kernel_param Kernel functions to use, and their parameters
#'@param mu prior mean
#'@param s2 known variances
#'@param sigma2 noise variance
#'@param opt_method optimization method for estimating prior parameters in `optim`
#'@param fix_ whether fix those parameters
#'@param l_b,r_b lower and upper bound of prior parameters
#'@param Jitter added to diagonal of the Kernel matrix for numerical stability
#'@param n_restart number of re-start of different kernel params values
#'@export
sgp = function(y,
               X=NULL,X_ind=NULL,m=10,
               kernel_func=Maternkernel,
               kernel_param=c(0,0),
               mu=NULL,
               s2=NULL,sigma2=NULL,
               opt_method='L-BFGS-B',
               fix_x_ind=F,fix_kernel_param=F,fix_sigma2=F,fix_mu=F,
               l_b=-Inf,r_b=Inf,
               Jitter=1e-5,
               n_restart=5,
               verbose=FALSE){
  t0 = Sys.time()
  n = length(y)
  if(is.null(X)){
    X = seq(0,1,length.out=n)
  }
  X_range = range(X)
  if(is.null(X_ind)){
    m = min(m,n)
    X_ind = seq(X_range[1],X_range[2],length.out=m)
  }
  m = length(X_ind)
  n_kernel_param = length(kernel_param)
  if(is.null(s2)){s2 = 1}

  if(fix_x_ind){
    d_mat_nm = Rfast::Outer(X_ind,X,"-")
    d_mat_mm = Rfast::Outer(X_ind,X_ind,"-")
  }else{
    d_mat_nm=NULL
    d_mat_mm=NULL
  }
  if(is.null(sigma2)){
    sigma2 = sd_est_diff2(y)^2
  }

  if(!(fix_x_ind&fix_kernel_param&fix_sigma2&fix_mu)){
    if(kernel_param[1]==0&!fix_kernel_param){
      kernel_param = c(log(max(var(y)-sigma2,1)),0)
    }
    param_hat = try(sgp_workhorse(y,
                              X,X_ind,m,
                              kernel_func,
                              kernel_param,
                              mu,
                              s2,sigma2,
                              opt_method,
                              fix_x_ind,fix_kernel_param,fix_sigma2,fix_mu,
                              l_b,r_b,
                              Jitter,
                              n_kernel_param,
                              X_range,
                              d_mat_nm,
                              d_mat_mm),silent = !verbose)
    if(class(param_hat)=="try-error"){
      param_hat = try(sgp_workhorse(y,
                                    X,X_ind,m,
                                    kernel_func,
                                    kernel_param,
                                    mu,
                                    s2,sigma2,
                                    "Nelder-Mead",
                                    fix_x_ind,fix_kernel_param,fix_sigma2,fix_mu,
                                    l_b,r_b,
                                    Jitter,
                                    n_kernel_param,
                                    X_range,
                                    d_mat_nm,
                                    d_mat_mm),silent = !verbose)
    }
    num_try=0
    while(class(param_hat)=='try-error'&num_try<n_restart){
      num_try = num_try + 1
      if(verbose){
        cat(paste("re-initializing",num_try))
        cat('\n')
      }
      if(!fix_x_ind){
        X_ind=runif(m,min=X_range[1],max=X_range[2])
      }
      if(!fix_kernel_param){
        kernel_param = rnorm(n_kernel_param)
      }
      param_hat = try(sgp_workhorse(y,
                                    X,X_ind,m,
                                    kernel_func,
                                    kernel_param,
                                    mu,
                                    s2,sigma2,
                                    opt_method,
                                    fix_x_ind,fix_kernel_param,fix_sigma2,fix_mu,
                                    l_b,r_b,
                                    Jitter,
                                    n_kernel_param,
                                    X_range,
                                    d_mat_nm,
                                    d_mat_mm),silent = !verbose)
    }
    if(class(param_hat)=='try-error'){
      stop('See error messages.')
    }
    X_ind=param_hat$X_ind
    mu=param_hat$mu
    kernel_param=param_hat$kernel_param
    sigma2=param_hat$sigma2
  }
  elbo = sgp_obj(X_ind,kernel_param,log(sigma2),mu,s2,y,X,kernel_func,Jitter,d_mat_nm,d_mat_mm)
  post = sgp_get_posterior(X_ind,kernel_param,log(sigma2),mu,s2,y,X,kernel_func,Jitter,d_mat_nm,d_mat_mm)
  return(list(elbo=-elbo,
              posterior=list(mean=post$mean,sd=sqrt(post$v),mean_ind=post$mean_ind,V_ind = post$V_ind),
              fitted_g=list(mu=post$mu,X_ind=X_ind,kernel_param=kernel_param,sigma2=sigma2,kernel_func=kernel_func),
              run_time = difftime(Sys.time(),t0),
              opt_res=param_hat$opt_res))
}
sgp_workhorse = function(y,
                         X,X_ind,m,
                         kernel_func,
                         kernel_param,
                         mu,
                         s2,sigma2,
                         opt_method,
                         fix_x_ind,fix_kernel_param,fix_sigma2,fix_mu,
                         l_b,r_b,
                         Jitter,
                         n_kernel_param,
                         X_range,
                         d_mat_nm,
                         d_mat_mm){

  if(!fix_mu){
    mu = NULL
  }
  if(!fix_kernel_param & (fix_x_ind&fix_sigma2)){
    opt_res = optim(par=kernel_param,fn=sgp_obj_kernel_param_only,
                    y=y,X=X,s2=s2,kernel_func=kernel_func,
                    X_ind=X_ind,mu=mu,log_sigma2=log(sigma2),
                    Jitter=Jitter,
                    d_mat_nm=d_mat_nm,
                    d_mat_mm=d_mat_mm,
                    method = opt_method,
                    lower=l_b,upper=r_b)
    kernel_param = opt_res$par
  }else if(!fix_kernel_param & !fix_sigma2 & fix_x_ind){
    opt_res = optim(par=c(kernel_param,log(sigma2)),fn=sgp_obj_kernel_param_and_sigma2,
                    y=y,X=X,s2=s2,kernel_func=kernel_func,
                    X_ind=X_ind,mu=mu,
                    Jitter=Jitter,
                    n_kernel_param=n_kernel_param,
                    d_mat_nm=d_mat_nm,
                    d_mat_mm=d_mat_mm,
                    method = opt_method,
                    lower=l_b,upper=r_b)
    kernel_param = opt_res$par[1:n_kernel_param]
    sigma2 = exp(opt_res$par[length(opt_res$par)])
  }else{
    if(l_b==-Inf){
      l_b = c(rep(X_range[1],m),rep(-5,n_kernel_param),-5)
    }
    if(r_b==Inf){
      r_b=c(rep(X_range[2],m),rep(5,n_kernel_param),5)
    }
    init_val = c(X_ind,kernel_param,log(sigma2))
    opt_res = optim(par=init_val,fn=sgp_obj_for_optim,y=y,X=X,s2=s2,m=m,n_kernel_param=n_kernel_param,kernel_func=kernel_func,
                    X_ind=X_ind,kernel_param=kernel_param,mu=mu,log_sigma2=log(sigma2),
                    fix_kernel_param=fix_kernel_param,fix_sigma2=fix_sigma2,fix_x_ind=fix_x_ind,fix_mu=fix_mu,method = opt_method,Jitter=Jitter,
                    lower=l_b,upper=r_b)
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
      sigma2 = exp(params[length(params)])
    }
  }
  return(list(X_ind=X_ind,mu=mu,kernel_param=kernel_param,sigma2=sigma2,opt_res=opt_res))
}
# sgp_obj_opt_kernel = function(kernel_param,y,X,X_ind,mu,sigma2,s2,kernel_func){
#   log_sigma2 = log(sigma2)
#   obj = sgp_obj(X_ind,kernel_param,log_sigma2,mu,s2,y,X,kernel_func)
#   print(obj)
#   print(kernel_param)
#   return(obj)
# }

sgp_obj_kernel_param_only = function(params,y,X,s2,kernel_func,
                             X_ind,mu,log_sigma2,Jitter,d_mat_nm,d_mat_mm){
  obj = sgp_obj(X_ind,params,log_sigma2,mu,s2,y,X,kernel_func,Jitter,d_mat_nm,d_mat_mm)
  return(obj)
}

sgp_obj_kernel_param_and_sigma2 = function(params,y,X,s2,kernel_func,
                                     X_ind,mu,Jitter,n_kernel_param,d_mat_nm,d_mat_mm){
  obj = sgp_obj(X_ind,params[1:n_kernel_param],params[length(params)],mu,s2,y,X,kernel_func,Jitter,d_mat_nm,d_mat_mm)
  return(obj)
}

sgp_obj_for_optim = function(params,y,X,s2,m,n_kernel_param,kernel_func,
                             X_ind,kernel_param,mu,log_sigma2,fix_kernel_param,fix_sigma2,fix_x_ind,fix_mu,Jitter){
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
  obj = sgp_obj(X_ind,kernel_param,log_sigma2,mu,s2,y,X,kernel_func,Jitter)
  return(obj)
}

sgp_matrix_helper = function(X_ind,kernel_param,log_sigma2,mu,s2,y,X,kernel_func,Jitter,d_mat_nm=NULL,d_mat_mm=NULL){
  sigma2 = exp(log_sigma2)
  Knm = kernel_func(X,X_ind,kernel_param,d_mat=d_mat_nm)
  Kmm = kernel_func(X_ind,X_ind,kernel_param,eps=Jitter,d_mat=d_mat_mm)
  Knn_diag = kernel_func(X,X,kernel_param,diagonal=T)
  L_Kmm = t(Rfast::cholesky(Kmm))
  L_Kmm_inv = forwardsolve(L_Kmm,diag(ncol(Kmm)))
  Kmm_inv = crossprod(L_Kmm_inv)
  A = Rfast::mat.mult(Knm,Kmm_inv)
  ASA = Rfast::Crossprod(A/s2,A)/sigma2
  V = solve(ASA + Kmm_inv)
  L_V = t(Rfast::cholesky(V))
  colsumAS = Rfast::colsums(A/s2)
  return(list(Knm=Knm,
              Kmm=Kmm,
              Knn_diag=Knn_diag,
              L_Kmm=L_Kmm,
              L_Kmm_inv=L_Kmm_inv,
              Kmm_inv=Kmm_inv,
              A=A,
              ASA=ASA,
              V=V,
              L_V=L_V,
              colsumAS=colsumAS))
}


sgp_get_posterior = function(X_ind,kernel_param,log_sigma2,mu,s2,y,X,kernel_func,Jitter,d_mat_nm,d_mat_mm){
  n = length(y)
  sigma2 = exp(log_sigma2)
  temp = sgp_matrix_helper(X_ind,kernel_param,log_sigma2,mu,s2,y,X,kernel_func,Jitter,d_mat_nm,d_mat_mm)
  tilde1 = rep(1,n) - c(Rfast::rowsums(temp$A))
  a = (sum(y/s2*tilde1))/sigma2/(sum(tilde1/s2*tilde1)/sigma2+sum(temp$Kmm_inv))
  b = (rowSums(temp$Kmm_inv)-crossprod(temp$A/s2,tilde1)/sigma2)/(sum(tilde1/s2*tilde1)/sigma2+sum(temp$Kmm_inv))
  if(is.null(mu)){
    beta = solve(temp$ASA + temp$Kmm_inv + tcrossprod(crossprod(temp$A/s2,tilde1),b)/sigma2- tcrossprod(rowSums(temp$Kmm_inv),b), crossprod(temp$A/s2,y)/sigma2-crossprod(temp$A/s2,tilde1)*a/sigma2+rowSums(temp$Kmm_inv)*a)
    Ab = temp$A%*%beta
    mu = a + sum(b*beta)
  }else{
    beta = temp$V%*%(crossprod(temp$A/s2,y)/sigma2+mu*(rowSums(temp$Kmm_inv)-crossprod(temp$A/s2,tilde1)/sigma2))
    Ab = temp$A%*%beta
  }
  v = temp$Knn_diag - Rfast::rowsums(Rfast::mat.mult(temp$Knm,t(temp$L_Kmm_inv))^2)  + Rfast::rowsums(Rfast::mat.mult(temp$A,temp$L_V)^2)
  return(list(mean=Ab+mu*tilde1,var=v,mu=mu,mean_ind=beta,V_ind=temp$V))
}

sgp_obj = function(X_ind,kernel_param,log_sigma2,mu,s2,y,X,kernel_func,Jitter,d_mat_nm=NULL,d_mat_mm=NULL){

  n = length(y)
  m = length(X_ind)
  sigma2 = exp(log_sigma2)
  temp = sgp_matrix_helper(X_ind,kernel_param,log_sigma2,mu,s2,y,X,kernel_func,Jitter,d_mat_nm,d_mat_mm)
  tilde1 = rep(1,n) - c(Rfast::rowsums(temp$A))
  Kmm_inv1 = Rfast::rowsums(temp$Kmm_inv)
  ADtilde1 = crossprod(temp$A/s2,tilde1)
  tilde1Dtilde1 = sum(tilde1^2/s2)
  sum_Kmminv = sum(temp$Kmm_inv)
  if(is.null(mu)){
    a = (sum(y/s2*tilde1))/sigma2/(tilde1Dtilde1/sigma2+sum_Kmminv)
    b = (Kmm_inv1-ADtilde1/sigma2)/(tilde1Dtilde1/sigma2+sum_Kmminv)
    beta = solve(temp$ASA + temp$Kmm_inv + tcrossprod(ADtilde1,b)/sigma2- tcrossprod(Kmm_inv1,b), crossprod(temp$A/s2,y)/sigma2-ADtilde1*a/sigma2+Kmm_inv1*a)
    Ab = temp$A%*%beta
    mu = a + sum(b*beta)
  }else{
    beta = temp$V%*%(crossprod(temp$A/s2,y)/sigma2+mu*(Kmm_inv1-ADtilde1/sigma2))
    Ab = temp$A%*%beta
  }
  F1=-n/2*log(2*pi*sigma2)-sum(log(s2))/2-(sum(y^2/s2)-2*sum(y/s2*Ab)-2*mu*sum(y/s2*tilde1)+sum(Ab^2/s2)+2*mu*sum(Ab/s2*tilde1)+mu^2*tilde1Dtilde1+sum(Rfast::mat.mult(temp$A/sqrt(s2),temp$L_V)^2)+sum(temp$Knn_diag/s2)-sum(Rfast::Tcrossprod(temp$Knm/sqrt(s2),temp$L_Kmm_inv)^2))/2/sigma2
  F2 = -sum(log(diag(temp$L_Kmm))) - sum((temp$L_Kmm_inv%*%beta)^2)/2 - sum((temp$L_Kmm_inv%*%temp$L_V)^2)/2+sum(log(diag(temp$L_V))) + m*0.5 + mu*sum(temp$Kmm_inv%*%beta) - mu^2*sum_Kmminv/2
  elbo=drop(F1+F2)
  return(-elbo)
}

#'@title RBF Kernel
#'@param d_mat output of outer(X1,X2,"-")
#'@export
RBFkernel = function(X1,X2,kernel_param,eps=NULL,diagonal=F,d_mat=NULL){
  coeff = exp(kernel_param[1])
  scales = exp(kernel_param[2])
  #X1 = cbind(X1)
  #X2 = cbind(X2)
  if(diagonal){
    return(rep(coeff,length(X1)))
  }else{
    if(!is.null(d_mat)){
      sqdist_xy = d_mat^2
    }else{
      sqdist_xy <- Rfast::Outer(X2,X1,"-")^2
    }
    if(!is.null(eps)){
      K = coeff*exp(-0.5*sqdist_xy/scales)+diag(eps,length(X1))
    }else{
      K = coeff*exp(-0.5*sqdist_xy/scales)
    }
    return(K)
  }
}

#'@title Matern Kernel
#'@export
Maternkernel = function(X1,X2,kernel_param,eps=NULL,diagonal=F,nu=3/2,d_mat=NULL,calc_grad_scales=FALSE){
  coeff = exp(kernel_param[1])
  scales = exp(kernel_param[2])
  if(diagonal){
    return(rep(coeff,length(X1)))
  }else{
    if(!is.null(d_mat)){
      abs_dist = abs(d_mat)
    }else{
      abs_dist <- abs(Rfast::Outer(X2,X1,"-"))
    }

    if(nu==1/2){
      R = exp(-abs_dist/scales)
      if(calc_grad_scales){
        g_scales = coeff*abs_dist*exp(-abs_dist/scales)/scales
      }
    }else if(nu==3/2){
      R = (1+sqrt(3)*abs_dist/scales)*exp(-sqrt(3)*abs_dist/scales)
      if(calc_grad_scales){
        g_scales = coeff*(3*abs_dist^2/scales^2*exp(-sqrt(3)*abs_dist/scales))
      }
    }else if(nu==5/2){
      R = (1+sqrt(5)*abs_dist/scales+5*abs_dist^2/3/scales^2)*exp(-sqrt(5)*abs_dist/scales)
      if(calc_grad_scales){
        g_scales = coeff*5*abs_dist^2*(sqrt(5)*abs_dist+scales)*exp(-sqrt(5)*abs_dist/scales)/(3*scales^3)
      }
    }else{
      stop('nu can only be 1/3, 3/2, 5/2')
    }
    if(!is.null(eps)){
      R = R+diag(eps,length(X1))
    }
    if(calc_grad_scales){
      return(list(K=coeff*R,g_scales=g_scales))
    }else{
      return(coeff*R)
    }
  }
}

# Maternkernel_grad = function(X1,X2,kernel_param,nu=3/2,d_mat=NULL){
#   coeff = exp(kernel_param[1])
#   scales = exp(kernel_param[2])
#   if(!is.null(d_mat)){
#     abs_dist = abs(d_mat)
#   }else{
#     abs_dist <- abs(Rfast::Outer(X2,X1,"-"))
#   }
#   if(nu==1/2){
#     g_scales = coeff*abs_dist*exp(-abs_dist/scales)/scales
#     #g_coeff = coeff*exp(-abs_dist/scales)
#   }else if(nu==3/2){
#     g_scales = coeff*(3*abs_dist^2/scales^2*exp(-sqrt(3)*abs_dist/scales))
#     #g_coeff = coeff*(1+sqrt(3)*abs_dist/scales)*exp(-sqrt(3)*abs_dist/scales)
#   }else if(nu==5/2){
#     g_scales = coeff*5*abs_dist^2*(sqrt(5)*abs_dist+scales)*exp(-sqrt(5)*abs_dist/scales)/(3*scales^3)
#     #g_coeff = coeff*(1+sqrt(5)*abs_dist/scales+5*abs_dist^2/3/scales^2)*exp(-sqrt(5)*abs_dist/scales)
#   }else{
#     stop('nu can only be 1/3, 3/2, 5/2')
#   }
#   return(g_scales)
# }

# sgp_obj2 = function(X_ind,kernel_param,log_sigma2,mu,s2,y,X,kernel_func){
#   n = length(y)
#   sigma2 = exp(log_sigma2)
#   Knm = kernel_func(X,X_ind,kernel_param)
#   Kmm = kernel_func(X_ind,X_ind,kernel_param,eps=1e-5)
#   Knn_diag = kernel_func(X,X,kernel_param,diagonal=T)
#   Kmm_inv = solve(Kmm)
#   temp = Knm%*%Kmm_inv%*%t(Knm)
#   dmvnorm(y,rep(0,n),diag(sigma2,n)+temp,log = T) - sum(Knn_diag)/2/sigma2 + sum(diag(temp))/2/sigma2
# }


################ Old obj where the model is y~N(mu+f,sigma2^2D), f is mean 0. ################
################################################################################################
# sgp_obj = function(X_ind,kernel_param,log_sigma2,mu,s2,y,X,kernel_func,Jitter,d_mat_nm=NULL,d_mat_mm=NULL){
#   n = length(y)
#   m = length(X_ind)
#   sigma2 = exp(log_sigma2)
#   Knm = kernel_func(X,X_ind,kernel_param,d_mat=d_mat_nm)
#   Kmm = kernel_func(X_ind,X_ind,kernel_param,eps=Jitter,d_mat=d_mat_mm)
#   Knn_diag = kernel_func(X,X,kernel_param,diagonal=T)
#   L_Kmm = t(Rfast::cholesky(Kmm))
#   L_Kmm_inv = forwardsolve(L_Kmm,diag(ncol(Kmm)))
#   #L_Kmm_inv = solve(L_Kmm)
#   Kmm_inv = crossprod(L_Kmm_inv)
#   #Kmm_inv = solve(Kmm)
#   #A = Knm%*%Kmm_inv
#   A = Rfast::mat.mult(Knm,Kmm_inv)
#   ASA = Rfast::Crossprod(A/s2,A)/sigma2
#   V = solve(ASA + Kmm_inv)
#   L_V = t(Rfast::cholesky(V))
#   colsumAS = Rfast::colsums(A/s2)
#   if(is.null(mu)){
#     # beta = solve(ASA + Kmm_inv- tcrossprod(colsumAS)/sum(1/s2)/sigma2)%*%(t(A/s2)%*%y-sum(y/s2)/sum(1/s2)*colsumAS)/sigma2
#     beta = solve(ASA + Kmm_inv- tcrossprod(colsumAS)/sum(1/s2)/sigma2, crossprod(A/s2,y)-sum(y/s2)/sum(1/s2)*colsumAS)/sigma2
#     Ab = A%*%beta
#     mu = (sum(y/s2)-sum(Ab/s2))/sum(1/s2)
#   }else{
#     # beta = solve(crossprod(A/s2,A)/sigma2 + Kmm_inv)%*%(t(A/s2)%*%y-mu*colSums(A/s2))/sigma2
#     beta = V%*%(crossprod(A/s2,y)-mu*colsumAS)/sigma2
#     Ab = A%*%beta
#   }
#   # obj = -n/2*log(2*pi*sigma2) - sum(log(s2))/2 - (sum(y^2/s2)-2*mu*sum(y/s2)-2*sum(y*Ab/s2)+mu^2*sum(1/s2)+2*mu*sum(Ab/s2)+sum(Ab^2/s2)+sum(diag(t(A/s2)%*%A%*%V))+sum(Knn_diag/s2)-sum(diag((A/s2)%*%t(Knm))))/2/sigma2 - as.numeric(determinant(Kmm,logarithm=T)$modulus)/2 - t(beta)%*%Kmm_inv%*%beta/2 - sum(diag(Kmm_inv%*%V)/2)+as.numeric(determinant(V,logarithm = T)$modulus)/2 + n*0.5
#   obj = -n/2*log(2*pi*sigma2) - sum(log(s2))/2 - (sum(y^2/s2)-2*mu*sum(y/s2)-2*sum(y*Ab/s2)+mu^2*sum(1/s2)+2*mu*sum(Ab/s2)+sum(Ab^2/s2)+sum(Rfast::mat.mult(A/sqrt(s2),L_V)^2)+sum(Knn_diag/s2)-sum(Rfast::Tcrossprod(Knm/sqrt(s2),L_Kmm_inv)^2))/2/sigma2 -sum(log(diag(L_Kmm)))/2 - sum((L_Kmm_inv%*%beta)^2)/2 - sum((L_Kmm_inv%*%L_V)^2)/2+sum(log(diag(L_V)))/2 + m*0.5
#   print(kernel_param)
#   print(mu)
#   print(drop(obj))
#   return(-drop(obj))
# }

# sgp_get_posterior = function(X_ind,kernel_param,log_sigma2,mu,s2,y,X,kernel_func,Jitter,d_mat_nm,d_mat_mm){
#   n = length(y)
#   sigma2 = exp(log_sigma2)
#   Knm = kernel_func(X,X_ind,kernel_param,d_mat=d_mat_nm)
#   Kmm = kernel_func(X_ind,X_ind,kernel_param,eps=Jitter,d_mat=d_mat_mm)
#   Knn_diag = kernel_func(X,X,kernel_param,diagonal=T)
#   L_Kmm = Rfast::cholesky(Kmm)
#   L_Kmm_inv = backsolve(L_Kmm,diag(ncol(Kmm)))
#   Kmm_inv = tcrossprod(L_Kmm_inv)
#   A = Rfast::mat.mult(Knm,Kmm_inv)
#   ASA = Rfast::Crossprod(A/s2,A)/sigma2
#   V = solve(ASA + Kmm_inv)
#   L_V = t(Rfast::cholesky(V))
#   colsumAS = Rfast::colsums(A/s2)
#   if(is.null(mu)){
#     beta = solve(ASA + Kmm_inv- tcrossprod(colsumAS)/sum(1/s2)/sigma2, crossprod(A/s2,y)-sum(y/s2)/sum(1/s2)*colsumAS)/sigma2
#     Ab = A%*%beta
#     mu = (sum(y/s2)-sum(Ab/s2))/sum(1/s2)
#   }else{
#     beta = V%*%(crossprod(A/s2,y)-mu*colsumAS)/sigma2
#     Ab = A%*%beta
#   }
#   v = Knn_diag - Rfast::rowsums(Rfast::mat.mult(Knm,L_Kmm_inv)^2)  + Rfast::rowsums(Rfast::mat.mult(A,L_V)^2)
#   return(list(mean=Ab,var=v,mu=mu,mean_ind=beta,V_ind=V))
# }
