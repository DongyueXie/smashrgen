#'@title Sparse Gaussian process for Poisson data via variational inference
#'@param y count vector
pois_sgp = function(y,X=NULL,X_ind=NULL,m=30,
                    post_mean=NULL,V=NULL,
                    kernel_func=Maternkernel,
                    kernel_param=NULL,
                    mu=NULL,
                    s=NULL,
                    opt_method='L-BFGS-B',
                    fix_X_ind=F,fix_kernel_param=F,fix_mu=F,
                    maxiter=100,tol=1e-5,
                    maxiter_mean=100,tol_mean=1e-5,
                    maxiter_V=100,tol_V=1e-5,
                    l_b=-Inf,r_b=Inf,
                    Jitter=1e-5,
                    verbose=T,printevery=1){
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
  if(is.null(s)){s = rep(1,n)}
  if(is.null(post_mean)|is.null(V)|is.null(kernel_param)){
    temp = sgp(log(1+y/s),X,X_ind,mu=mean(log(1+y/s)),fix_x_ind = T,fix_mu = T,kernel_func = kernel_func)
    if(is.null(post_mean)){
      post_mean = temp$posterior$mean_ind
    }
    if(is.null(V)){
      V = temp$posterior$V_ind
    }
    if(is.null(kernel_param)){
      kernel_param = temp$fitted_g$kernel_param
    }
  }
  if(is.null(mu)){
    mu = mean(log(1+y/s))
  }
  if(l_b[1]==-Inf){
    l_b = c(rep(X_range[1],m),rep(-6,n_kernel_param))
  }
  if(r_b[1]==Inf){
    r_b=c(rep(X_range[2],m),rep(6,n_kernel_param))
  }

  elbo_tace = -Inf
  for(iter in 1:maxiter){
    # update post mean
    post_mean = pois_sgp_update_post_mean(X,y,s,X_ind,kernel_param,kernel_func,post_mean,V,mu,maxiter=maxiter_mean,tol=tol_mean,Jitter=Jitter)
    #print(paste("elbo after post_mean update:",round(pois_sgp_elbo(X,y,s,X_ind,kernel_param,kernel_func,post_mean,V,mu,Jitter=Jitter),3)))
    # update post var
    V = pois_sgp_update_V(X,y,s,X_ind,kernel_param,kernel_func,post_mean,V,mu,maxiter=maxiter_V,tol=tol_V,Jitter=Jitter)
    #print(paste("elbo after V update:",round(pois_sgp_elbo(X,y,s,X_ind,kernel_param,kernel_func,post_mean,V,mu,Jitter=Jitter),3)))
    # update kernel and x_ind
    prior = pois_sgp_update_prior(X,y,s,X_ind,kernel_param,kernel_func,post_mean,V,mu,fix_kernel_param,fix_X_ind,m,n_kernel_param,opt_method,l_b,r_b,Jitter=Jitter)
    X_ind = prior$X_ind
    kernel_param = prior$kernel_param
    #print(paste("elbo after prior update:",round(pois_sgp_elbo(X,y,s,X_ind,kernel_param,kernel_func,post_mean,V,mu,Jitter=Jitter),3)))
    # update mu
    mu = pois_sgp_update_mu(X,y,s,X_ind,kernel_param,kernel_func,post_mean,V,mu,fix_mu,Jitter)
    #print(paste("elbo after mu update:",round(pois_sgp_elbo(X,y,s,X_ind,kernel_param,kernel_func,post_mean,V,mu,Jitter=Jitter),3)))

    elbo = pois_sgp_elbo(X,y,s,X_ind,kernel_param,kernel_func,post_mean,V,mu,Jitter=Jitter)
    elbo_tace[iter+1] = elbo
    if(verbose){
      if(iter%%printevery==0){
        print(paste("At iter ",iter, ", elbo = ",elbo,sep=''))
      }
    }
    if((elbo - elbo_tace[iter])/n < tol){
      break
    }


  }
  post = pois_sgp_get_posterior(X,y,s,X_ind,kernel_param,kernel_func,post_mean,V,mu,Jitter)
  return(list(elbo_tace=elbo_tace,elbo=elbo,
              fitted_g = list(X_ind=X_ind,kernel_param=kernel_param,mu=mu,kernel_func=kernel_func),
              posterior=list(mean=post$mean,sd=sqrt(post$v),rate=post$rate,rate_sd = sqrt(post$rate_v),mean_ind=post_mean,V_ind=V),
              run_time = difftime(Sys.time(),t0)))
}

pois_sgp_update_prior = function(X,y,s,X_ind,kernel_param,kernel_func,post_mean,V,mu,
                                 fix_kernel_param,fix_X_ind,m,n_kernel_param,opt_method,l_b,r_b,Jitter){
  if(!(fix_kernel_param&fix_X_ind)){
    if(fix_X_ind&(!fix_kernel_param)){
      res = try(optim(kernel_param,pois_sgp_obj_for_optim_kernel_only,
                  X=X,y=y,s=s,X_ind=X_ind,
                  kernel_func=kernel_func,post_mean=post_mean,V=V,mu=mu,
                  lower = rep(-6,n_kernel_param),upper = rep(6,n_kernel_param),Jitter=Jitter,
                  method = opt_method))
      if(class(res)!='try-error'){
        kernel_param = res$par
      }
    }else{
      res = try(optim(c(X_ind,kernel_param),pois_sgp_obj_for_optim,
                  m=m,n_kernel_param=n_kernel_param,X=X,y=y,s=s,X_ind=X_ind,kernel_param=kernel_param,
                  kernel_func=kernel_func,post_mean=post_mean,V=V,mu=mu,fix_X_ind=fix_X_ind,fix_kernel_param=fix_kernel_param,
                  lower = l_b,upper = r_b,Jitter=Jitter,
                  method = opt_method))
      if(class(res)!="try-error"){
        if(!fix_X_ind){
          X_ind = res$par[1:m]
        }
        if(!fix_kernel_param){
          kernel_param = res$par[(m+1):(m+n_kernel_param)]
        }
      }
    }
  }

  #print(pois_sgp_elbo(X,y,s,X_ind,kernel_param,kernel_func,post_mean,V,mu))
  return(list(X_ind = X_ind,kernel_param=kernel_param))
}

pois_sgp_update_mu = function(X,y,s,X_ind,kernel_param,kernel_func,post_mean,V,mu,fix_mu,Jitter){
  if(!fix_mu){
    n = length(y)
    m = length(X_ind)
    Knm = kernel_func(X,X_ind,kernel_param)
    Kmm = kernel_func(X_ind,X_ind,kernel_param,eps=Jitter)
    Knn_diag = kernel_func(X,X,kernel_param,diagonal=T)
    L_Kmm = t(chol(Kmm))
    L_Kmm_inv = forwardsolve(L_Kmm,diag(ncol(Kmm)))
    Kmm_inv = crossprod(L_Kmm_inv)
    A = Knm%*%Kmm_inv
    L_V = t(chol(V))
    Ab = drop(A%*%post_mean)
    delta = s*exp(Ab+rowSums((A%*%L_V)^2)/2+Knn_diag/2-rowSums((Knm%*%t(L_Kmm_inv))^2)/2)
    mu = log(sum(y)/sum(delta))
  }
  #print(pois_sgp_elbo(X,y,s,X_ind,kernel_param,kernel_func,post_mean,V,mu))
  return(mu)
}

pois_sgp_update_V = function(X,y,s,X_ind,kernel_param,kernel_func,post_mean,V,mu,maxiter=100,tol=1e-5,Jitter){
  V_old = V
  elbo_old = pois_sgp_elbo(X,y,s,X_ind,kernel_param,kernel_func,post_mean,V,mu,Jitter)
  for(i in 1:maxiter){
    V = solve(-grad_post_mean(X,y,s,X_ind,kernel_param,kernel_func,post_mean,V_old,mu,return_grad=F,return_hess=T,Jitter))
    #print(pois_sgp_elbo(X,y,s,X_ind,kernel_param,kernel_func,post_mean,V,mu,Jitter))
    #print(norm(V-V_old,type ='2'))
    elbo = pois_sgp_elbo(X,y,s,X_ind,kernel_param,kernel_func,post_mean,V,mu,Jitter)
    # if(norm(V-V_old,type ='2')<tol){
    #   break
    # }
    if((elbo-elbo_old)<tol){
      break
    }else{
      V_old = V
    }
  }
  if(i==maxiter){
    print("max iteration reached when updating posterior variance")
  }
  if((elbo-elbo_old)<0){
    return(V_old)
  }else{
    return(V)
  }
}

pois_sgp_update_post_mean = function(X,y,s,X_ind,kernel_param,kernel_func,post_mean,V,mu,maxiter=100,tol=1e-5,Jitter){
  post_mean_old = post_mean
  elbo_old = pois_sgp_elbo(X,y,s,X_ind,kernel_param,kernel_func,post_mean,V,mu,Jitter)
  for(iter in 1:maxiter){
    grad_hess = grad_post_mean(X,y,s,X_ind,kernel_param,kernel_func,post_mean_old,V,mu,return_grad=T,return_hess=T,Jitter)
    grad = -grad_hess$grad
    hess = -grad_hess$hess
    # direction = c(pcgsolve(hess,post_mean))
    direction = c(solve(hess)%*%grad)
    post_mean = post_mean_old - direction
    elbo = pois_sgp_elbo(X,y,s,X_ind,kernel_param,kernel_func,post_mean,V,mu,Jitter)
    #print(pois_sgp_elbo(X,y,s,X_ind,kernel_param,kernel_func,post_mean,V,mu))
    # if(norm(post_mean-post_mean_old,type='2')<tol){
    #   break
    # }
    if((elbo-elbo_old)<tol){
      break
    }else{
      post_mean_old = post_mean
    }
  }
  if(iter==maxiter){
    print("max iteration reached when updating posterior mean")
  }
  if((elbo-elbo_old)<0){
    return(post_mean_old)
  }else{
    return(post_mean)
  }
}

grad_post_mean = function(X,y,s,X_ind,kernel_param,kernel_func,post_mean,V,mu,return_grad=T,return_hess=T,Jitter){
  n = length(y)
  m = length(X_ind)
  Knm = kernel_func(X,X_ind,kernel_param)
  Kmm = kernel_func(X_ind,X_ind,kernel_param,eps=Jitter)
  Knn_diag = kernel_func(X,X,kernel_param,diagonal=T)
  L_Kmm = t(chol(Kmm))
  L_Kmm_inv = forwardsolve(L_Kmm,diag(ncol(Kmm)))
  Kmm_inv = crossprod(L_Kmm_inv)
  A = Knm%*%Kmm_inv
  L_V = t(chol(V))
  Ab = drop(A%*%post_mean)

  delta = s*exp(Ab+rowSums((A%*%L_V)^2)/2+Knn_diag/2-rowSums((Knm%*%t(L_Kmm_inv))^2)/2)
  if(return_grad&return_hess){
    grad = crossprod(A,y-exp(mu)*delta)-Kmm_inv%*%post_mean
    hess = -exp(mu)*crossprod(A*sqrt(delta))-Kmm_inv
    return(list(grad=grad,hess=hess))
  }else if(return_grad){
    grad = crossprod(A,y-exp(mu)*delta)-Kmm_inv%*%post_mean
    return(grad)
  }else if(return_hess){
    hess = -exp(mu)*crossprod(A*sqrt(delta))-Kmm_inv
    return(hess)
  }else{
    return(NULL)
  }

}


pois_sgp_elbo = function(X,y,s,X_ind,kernel_param,kernel_func,post_mean,V,mu,Jitter){
  n = length(y)
  m = length(X_ind)
  Knm = kernel_func(X,X_ind,kernel_param)
  Kmm = kernel_func(X_ind,X_ind,kernel_param,eps=Jitter)
  Knn_diag = kernel_func(X,X,kernel_param,diagonal=T)
  L_Kmm = t(chol(Kmm))
  L_Kmm_inv = forwardsolve(L_Kmm,diag(ncol(Kmm)))
  Kmm_inv = crossprod(L_Kmm_inv)
  A = Knm%*%Kmm_inv
  L_V = t(chol(V))
  Ab = drop(A%*%post_mean)
  obj = sum(y*Ab)+mu*sum(y)-exp(mu)*sum(s*exp(Ab+rowSums((A%*%L_V)^2)/2+Knn_diag/2-rowSums((Knm%*%t(L_Kmm_inv))^2)/2))+sum(y*log(s))-sum(lfactorial(y))-sum(log(diag(L_Kmm)))/2 - sum((L_Kmm_inv%*%post_mean)^2)/2 - sum((L_Kmm_inv%*%L_V)^2)/2+sum(log(diag(L_V)))/2-m/2*log(2*pi)+m/2
  return(drop(obj))
}

pois_sgp_obj_for_optim = function(params,m,n_kernel_param,X,y,s,X_ind,kernel_param,kernel_func,post_mean,V,mu,fix_X_ind,fix_kernel_param,Jitter){
  #n = length(y)
  if(!fix_X_ind){
    X_ind = params[1:m]
  }
  if(!fix_kernel_param){
    kernel_param = params[(m+1):(m+n_kernel_param)]
  }
  obj=-pois_sgp_elbo(X,y,s,X_ind,kernel_param,kernel_func,post_mean,V,mu,Jitter)
  return(obj)
}

pois_sgp_obj_for_optim_kernel_only = function(kernel_param,X,y,s,X_ind,kernel_func,post_mean,V,mu,Jitter){
  obj = -pois_sgp_elbo(X,y,s,X_ind,kernel_param,kernel_func,post_mean,V,mu,Jitter)
  #print(kernel_param)
  #print(obj)
  return(obj)
}

pois_sgp_get_posterior = function(X,y,s,X_ind,kernel_param,kernel_func,post_mean,V,mu,Jitter){
  n = length(y)
  Knm = kernel_func(X,X_ind,kernel_param)
  Kmm = kernel_func(X_ind,X_ind,kernel_param,eps=Jitter)
  Knn_diag = kernel_func(X,X,kernel_param,diagonal=T)
  L_Kmm = chol(Kmm)
  L_Kmm_inv = backsolve(L_Kmm,diag(ncol(Kmm)))
  Kmm_inv = tcrossprod(L_Kmm_inv)
  A = Knm%*%Kmm_inv
  L_V = t(chol(V))
  Ab = drop(A%*%post_mean)
  v = Knn_diag - rowSums((Knm%*%L_Kmm_inv)^2)  + rowSums((A%*%L_V)^2)
  return(list(mean=Ab,var=v,rate=exp(mu)*exp(Ab+v/2),rate_v = (exp(v)-1)*exp(2*Ab+v)))
}

# RBFkernel = function(X1,X2,kernel_param,eps=NULL,diagonal=F){
#   coeff = exp(kernel_param[1])
#   scales = exp(kernel_param[2])
#   if(diagonal){
#     return(rep(coeff^2,length(X1)))
#   }else{
#     sqdist_xy <- outer(X1,X2,FUN="-")^2
#     if(!is.null(eps)){
#       K = coeff^2*exp(-0.5*sqdist_xy/scales^2)+diag(eps,length(X1))
#     }else{
#       K = coeff^2*exp(-0.5*sqdist_xy/scales^2)
#     }
#     return(K)
#   }
# }
