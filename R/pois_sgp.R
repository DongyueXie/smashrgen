
#'@title Poisson Sparse Gaussian Process
#'@param x A vector of Poisson observations.
#'@param s A vector of scaling factors for Poisson observations
#'@param g_init a list of prior parameters: X_ind, kernel_param, mu, kernel_func
#'@param fix_g a boolean or a vector of boolean of length 3. If TRUE, set the corresponding prior at g_init.
#'@param q_init a list of init values of posterior: mean_log_ind, v_log_ind. These are the posterior mean and var of the f at the inducing points.
#'@param control a list of other parameters passed to the `pois_sgp` function.
#'@return a list of posterior, log_likelihood, fitted_g
#'@details The model is
#'\deqn{x_{i}\sim \text{Poisson}(s_i\exp(\mu+f_{i})),}
#'\deqn{f\sim N(0,\Sigma(\theta)).}
#'@export
ebpm_pois_sgp = function(x,s=1,g_init=NULL,fix_g=FALSE,q_init=NULL,control=list()){
  if(!is.null(g_init)){
    X_ind = g_init$X_ind
    kernel_param = g_init$kernel_param
    mu=g_init$mu
    kernel_func=g_init$kernel_func
  }else{
    X_ind = NULL
    kernel_param = NULL
    mu = NULL
    kernel_func = Maternkernel
  }
  if(!is.null(q_init)){
    post_mean = q_init$mean_log_ind
    V = q_init$v_log_ind
  }else{
    post_mean = NULL
    V = NULL
  }
  if(length(fix_g)>1){
    fix_X_ind = fix_g[1]
    fix_kernel_param = fix_g[2]
    fix_mu = fix_g[3]
  }else{
    if(fix_g){
      fix_X_ind = T
      fix_kernel_param = T
      fix_mu = T
    }else{
      fix_X_ind = T
      fix_kernel_param = F
      fix_mu = F
    }
  }
  controls = modifyList(ebpm_pois_sgp_control_default(),control,keep.null = TRUE)
  res = pois_sgp(x,X=controls$X,X_ind=X_ind,m=controls$m,
                            post_mean=post_mean,V=V,
                            kernel_func=kernel_func,
                            kernel_param=kernel_param,
                            mu=mu,
                            sc=s,
                            opt_method=controls$opt_method,
                            fix_X_ind=fix_X_ind,fix_kernel_param=fix_kernel_param,fix_mu=fix_mu,
                            maxiter=controls$maxiter,tol=controls$tol,
                            maxiter_mean=controls$maxiter_mean,tol_mean=controls$tol_mean,
                            maxiter_V=controls$maxiter_V,tol_V=controls$tol_V,
                            l_b=controls$l_b,r_b=controls$r_b,
                            Jitter=controls$Jitter,
                            verbose=controls$verbose,printevery=controls$printevery)
  return(list(posterior=list(mean=res$posterior$rate,
                             mean_log=res$posterior$mean,
                             sd=res$posterior$rate_sd,
                             sd_log = res$posterior$sd,
                             mean_log_ind=res$posterior$mean_ind,
                             v_log_ind=res$posterior$V_ind),
              log_likelihood = res$elbo,
              fitted_g=res$fitted_g))
}

ebpm_pois_sgp_control_default = function(){
  list(X=NULL,
       m=30,
       opt_method='L-BFGS-B',maxiter=100,tol=1e-5,
       maxiter_mean=100,tol_mean=1e-5,
       maxiter_V=100,tol_V=1e-5,
       l_b=-Inf,r_b=Inf,
       Jitter=1e-5,
       verbose=F,printevery=1)
}

#'@title Sparse Gaussian process for Poisson data via variational inference
#'@param y count vector
#'@param X,X_ind grids, and inducing points
#'@param m number of inducing points
#'@param post_mean,V posterior mean and var of f at inducing point
#'@param sc scaling scalar or vectors
#'@param kernel_func,kernel_param functions and their parameters
#'@export
pois_sgp = function(y,X=NULL,X_ind=NULL,m=30,
                    post_mean=NULL,V=NULL,
                    kernel_func=Maternkernel,
                    kernel_param=NULL,
                    mu=NULL,
                    sc=NULL,
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
  if(is.null(sc)){sc = rep(1,n)}
  if(length(sc)==1){sc = rep(sc,n)}
  if(is.null(post_mean)|is.null(V)|is.null(kernel_param)){
    #temp = sgp(log(1+y/s),X,X_ind,mu=mean(log(1+y/s)),fix_x_ind = T,fix_mu = T,kernel_func = kernel_func)
    #ly = smash.poiss(y/s,log=T)
    temp = sgp(log(1+y/sc),X,X_ind,mu=NULL,fix_x_ind = T,kernel_func = kernel_func)
    mu = temp$fitted_g$mu
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
    mu = log(sum(y)/sum(sc))
  }
  if(l_b[1]==-Inf){
    l_b = c(rep(X_range[1],m),rep(-6,n_kernel_param))
  }
  if(r_b[1]==Inf){
    r_b=c(rep(X_range[2],m),rep(6,n_kernel_param))
  }

  elbo_tace = pois_sgp_elbo(X,y,sc,X_ind,kernel_param,kernel_func,post_mean,V,mu,Jitter=Jitter)
  #print(paste("elbo init:",round(elbo_tace,3)))
  for(iter in 1:maxiter){
    # update post mean
    post_mean = pois_sgp_update_post_mean(X,y,sc,X_ind,kernel_param,kernel_func,post_mean,V,mu,maxiter=maxiter_mean,tol=tol_mean,Jitter=Jitter)
    #print(paste("elbo after post_mean update:",round(pois_sgp_elbo(X,y,sc,X_ind,kernel_param,kernel_func,post_mean,V,mu,Jitter=Jitter),3)))
    # update post var
    V = pois_sgp_update_V(X,y,sc,X_ind,kernel_param,kernel_func,post_mean,V,mu,maxiter=maxiter_V,tol=tol_V,Jitter=Jitter)
    #print(paste("elbo after V update:",round(pois_sgp_elbo(X,y,sc,X_ind,kernel_param,kernel_func,post_mean,V,mu,Jitter=Jitter),3)))
    # update kernel and x_ind
    prior = pois_sgp_update_prior(X,y,sc,X_ind,kernel_param,kernel_func,post_mean,V,mu,fix_kernel_param,fix_X_ind,m,n_kernel_param,opt_method,l_b,r_b,Jitter=Jitter)
    X_ind = prior$X_ind
    kernel_param = prior$kernel_param
    #print(paste("elbo after prior update:",round(pois_sgp_elbo(X,y,sc,X_ind,kernel_param,kernel_func,post_mean,V,mu,Jitter=Jitter),3)))
    # update mu
    mu = pois_sgp_update_mu(X,y,sc,X_ind,kernel_param,kernel_func,post_mean,V,mu,fix_mu,Jitter)
    #print(paste("elbo after mu update:",round(pois_sgp_elbo(X,y,sc,X_ind,kernel_param,kernel_func,post_mean,V,mu,Jitter=Jitter),3)))

    elbo = pois_sgp_elbo(X,y,sc,X_ind,kernel_param,kernel_func,post_mean,V,mu,Jitter=Jitter)
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
  post = pois_sgp_get_posterior(X,y,sc,X_ind,kernel_param,kernel_func,post_mean,V,mu,Jitter)
  return(list(elbo_tace=elbo_tace,elbo=elbo,
              fitted_g = list(X_ind=X_ind,kernel_param=kernel_param,mu=mu,kernel_func=kernel_func),
              posterior=list(mean=post$mean,sd=sqrt(post$v),rate=post$rate,rate_sd = sqrt(post$rate_v),mean_ind=post_mean,V_ind=V),
              run_time = difftime(Sys.time(),t0)))
}

#'@title Function for optimize kernel parameters in the prior
pois_sgp_update_prior_kernel = function(X,y,sc,X_ind,kernel_param,kernel_func,post_mean,V,mu,Jitter){
  n_kernel_param = length(kernel_param)
  elbo_old = pois_sgp_elbo(X,y,sc,X_ind,kernel_param,kernel_func,post_mean,V,mu,Jitter)
  res = try(optim(kernel_param,pois_sgp_obj_for_optim_kernel_only,
                  X=X,y=y,sc=sc,X_ind=X_ind,
                  kernel_func=kernel_func,post_mean=post_mean,V=V,mu=mu,
                  lower = rep(-6,n_kernel_param),upper = rep(6,n_kernel_param),Jitter=Jitter,
                  method = opt_method),silent = TRUE)
  if(class(res)!='try-error'){
    if((-res$value)>elbo_old){
      kernel_param = res$par
    }else{
      refit = TRUE
    }
  }else{
    refit=TRUE
  }
  if(refit){
    res = try(pois_sgp_optimize_kernel_newton(kernel_param,X,y,sc,X_ind,kernel_func,post_mean,V,mu,Jitter))
    if(class(res)!='try-error'){
      kernel_param = res
    }
  }
  kernel_param
}

#'@title Master function for optimizing priors (X_ind, kernel_param)
pois_sgp_update_prior = function(X,y,sc,X_ind,kernel_param,kernel_func,post_mean,V,mu,
                                 fix_kernel_param,fix_X_ind,m,n_kernel_param,opt_method,l_b,r_b,Jitter){

  if(!(fix_kernel_param&fix_X_ind)){
    # if only optimize kernel_param
    if(fix_X_ind&(!fix_kernel_param)){
      kernel_param = pois_sgp_update_prior_kernel(X,y,sc,X_ind,kernel_param,kernel_func,post_mean,V,mu,Jitter)
    }else{
      elbo_old = pois_sgp_elbo(X,y,sc,X_ind,kernel_param,kernel_func,post_mean,V,mu,Jitter)
      # optimize both X_ind and kernel_param
      res = try(optim(c(X_ind,kernel_param),pois_sgp_obj_for_optim,
                  m=m,n_kernel_param=n_kernel_param,X=X,y=y,sc=sc,X_ind=X_ind,kernel_param=kernel_param,
                  kernel_func=kernel_func,post_mean=post_mean,V=V,mu=mu,fix_X_ind=fix_X_ind,fix_kernel_param=fix_kernel_param,
                  lower = l_b,upper = r_b,Jitter=Jitter,
                  method = opt_method),silent = TRUE)
      if(class(res)!="try-error"){
        if((-res$value)>elbo_old){
          if(!fix_X_ind){
            X_ind = res$par[1:m]
          }
          if(!fix_kernel_param){
            kernel_param = res$par[(m+1):(m+n_kernel_param)]
          }
        }else{
          refit_kernel_only=TRUE
        }
      }else{
        refit_kernel_only=TRUE
      }
      if(refit_kernel_only){
        # if optim not working well when optimizing both X_ind, and kernel param, only optimize kernel_param
        if(!fix_kernel_param&!fix_X_ind){
          kernel_param = pois_sgp_update_prior_kernel(X,y,sc,X_ind,kernel_param,kernel_func,post_mean,V,mu,Jitter)
        }
      }
    }
  }

  #print(pois_sgp_elbo(X,y,s,X_ind,kernel_param,kernel_func,post_mean,V,mu))
  return(list(X_ind = X_ind,kernel_param=kernel_param))
}

#'@title Function for optimize mu
pois_sgp_update_mu = function(X,y,sc,X_ind,kernel_param,kernel_func,post_mean,V,mu,fix_mu,Jitter){
  if(!fix_mu){
    delta = pois_sgp_matrix_helper(X,y,sc,X_ind,kernel_param,kernel_func,post_mean,V,mu,Jitter,get_delta=T)$delta
    mu = log(sum(y)/sum(delta))
  }
  #print(pois_sgp_elbo(X,y,s,X_ind,kernel_param,kernel_func,post_mean,V,mu))
  return(mu)
}

#'@title Function for updating posterior var
pois_sgp_update_V = function(X,y,sc,X_ind,kernel_param,kernel_func,post_mean,V,mu,maxiter=100,tol=1e-5,Jitter){
  V_old = V
  n = length(y)
  elbo_old = pois_sgp_elbo(X,y,sc,X_ind,kernel_param,kernel_func,post_mean,V,mu,Jitter)
  for(i in 1:maxiter){
    V = solve(-grad_post_mean(X,y,sc,X_ind,kernel_param,kernel_func,post_mean,V_old,mu,return_grad=F,return_hess=T,Jitter))
    #print(pois_sgp_elbo(X,y,s,X_ind,kernel_param,kernel_func,post_mean,V,mu,Jitter))
    #print(norm(V-V_old,type ='2'))
    elbo = pois_sgp_elbo(X,y,sc,X_ind,kernel_param,kernel_func,post_mean,V,mu,Jitter)
    # if(norm(V-V_old,type ='2')<tol){
    #   break
    # }
    if((elbo-elbo_old)/n<tol){
      break
    }else{
      V_old = V
      elbo_old = elbo
    }
  }
  # if(i==maxiter){
  #   print("max iteration reached when updating posterior variance")
  # }
  if((elbo-elbo_old)<0){
    return(V_old)
  }else{
    return(V)
  }
}

#'@title Function for optimizing posterior mean
pois_sgp_update_post_mean = function(X,y,sc,X_ind,kernel_param,kernel_func,post_mean,V,mu,maxiter=100,tol=1e-5,Jitter){
  post_mean_old = post_mean
  elbo_old = pois_sgp_elbo(X,y,sc,X_ind,kernel_param,kernel_func,post_mean,V,mu,Jitter)
  n = length(y)
  for(iter in 1:maxiter){
    grad_hess = grad_post_mean(X,y,sc,X_ind,kernel_param,kernel_func,post_mean_old,V,mu,return_grad=T,return_hess=T,Jitter)
    grad = -grad_hess$grad
    hess = -grad_hess$hess
    # direction = c(pcgsolve(hess,post_mean))
    direction = c(solve(hess)%*%grad)
    post_mean = post_mean_old - direction
    elbo = pois_sgp_elbo(X,y,sc,X_ind,kernel_param,kernel_func,post_mean,V,mu,Jitter)
    #print(pois_sgp_elbo(X,y,s,X_ind,kernel_param,kernel_func,post_mean,V,mu))
    # if(norm(post_mean-post_mean_old,type='2')<tol){
    #   break
    # }
    if((elbo-elbo_old)/n<tol){
      break
    }else{
      post_mean_old = post_mean
      elbo_old = elbo
    }
  }
  # if(iter==maxiter){
  #   print("max iteration reached when updating posterior mean")
  # }
  if((elbo-elbo_old)<0){
    return(post_mean_old)
  }else{
    return(post_mean)
  }
}

#'@title Gradient and hessian of eblo w.r.t post_mean
grad_post_mean = function(X,y,sc,X_ind,kernel_param,kernel_func,post_mean,V,mu,return_grad=T,return_hess=T,Jitter){
  temp = pois_sgp_matrix_helper(X,y,sc,X_ind,kernel_param,kernel_func,post_mean,V,mu,Jitter,get_delta=T)
  if(return_grad&return_hess){
    grad = crossprod(temp$A,y-exp(mu)*temp$delta)-temp$Kmm_inv%*%post_mean
    hess = -exp(mu)*crossprod(temp$A*sqrt(temp$delta))-temp$Kmm_inv
    return(list(grad=grad,hess=hess))
  }else if(return_grad){
    grad = crossprod(temp$A,y-exp(mu)*temp$delta)-temp$Kmm_inv%*%post_mean
    return(grad)
  }else if(return_hess){
    hess = -exp(mu)*crossprod(temp$A*sqrt(temp$delta))-temp$Kmm_inv
    return(hess)
  }else{
    return(NULL)
  }

}

#'@title calc elbo
pois_sgp_elbo = function(X,y,sc,X_ind,kernel_param,kernel_func,post_mean,V,mu,Jitter){
  m = length(X_ind)
  temp = pois_sgp_matrix_helper(X,y,sc,X_ind,kernel_param,kernel_func,post_mean,V,mu,Jitter,get_delta=T)
  obj = sum(y*temp$Ab)+mu*sum(y)-exp(mu)*sum(temp$delta)+sum(y*log(sc))-sum(lfactorial(y))-sum(log(diag(temp$L_Kmm)))/2 - sum((temp$L_Kmm_inv%*%post_mean)^2)/2 - sum((temp$L_Kmm_inv%*%temp$L_V)^2)/2+sum(log(diag(temp$L_V)))/2-m/2*log(2*pi)+m/2
  return(drop(obj))
}

#'@title objective function for optimizing elbo w.r.t X-ind and kernel_param for `optim`
pois_sgp_obj_for_optim = function(params,m,n_kernel_param,X,y,sc,X_ind,kernel_param,kernel_func,post_mean,V,mu,fix_X_ind,fix_kernel_param,Jitter){
  #n = length(y)
  if(!fix_X_ind){
    X_ind = params[1:m]
  }
  if(!fix_kernel_param){
    kernel_param = params[(m+1):(m+n_kernel_param)]
  }
  obj=-pois_sgp_elbo(X,y,sc,X_ind,kernel_param,kernel_func,post_mean,V,mu,Jitter)
  return(obj)
}

#'@title objective function for optimizing elbo w.r.t kernel_param for `optim`
pois_sgp_obj_for_optim_kernel_only = function(kernel_param,X,y,sc,X_ind,kernel_func,post_mean,V,mu,Jitter){
  obj = -pois_sgp_elbo(X,y,sc,X_ind,kernel_param,kernel_func,post_mean,V,mu,Jitter)
  #print(kernel_param)
  #print(obj)
  return(obj)
}

#'@title Newton's method for optimizing elbo w.r.t kernel_param, using numerical grad and hess
pois_sgp_optimize_kernel_newton = function(kernel_param,X,y,sc,X_ind,kernel_func,post_mean,V,mu,Jitter,tol=1e-5,maxiter=100){
  kernel_param_old = kernel_param
  elbo_old = pois_sgp_elbo(X,y,sc,X_ind,kernel_param,kernel_func,post_mean,V,mu,Jitter)
  elbo_init = elbo_old
  kernel_param_init = kernel_param_old
  n = length(y)
  for(iter in 1:maxiter){
    hess_mat = numDeriv::hessian(pois_sgp_obj_for_optim_kernel_only,kernel_param_old,X=X,y=y,sc=sc,X_ind=X_ind,
                             kernel_func=kernel_func,post_mean=post_mean,V=V,mu=mu,Jitter=Jitter)
    grad_vec = numDeriv::grad(pois_sgp_obj_for_optim_kernel_only,kernel_param_old,X=X,y=y,sc=sc,X_ind=X_ind,
                             kernel_func=kernel_func,post_mean=post_mean,V=V,mu=mu,Jitter=Jitter)
    if(length(c(grad_vec))==2){
      direction = c(hess_mat[2,2]*grad_vec[1]-hess_mat[1,2]*grad_vec[2], hess_mat[1,1]*grad_vec[2]-hess_mat[2,1]*grad_vec[1])/(prod(diag(hess_mat))-hess_mat[1,2]*hess_mat[2,1])
    }else{
      direction = c(solve(hess_mat)%*%grad_vec)
    }

    # print(direction)
    kernel_param = kernel_param_old - direction
    elbo = pois_sgp_elbo(X,y,sc,X_ind,kernel_param,kernel_func,post_mean,V,mu,Jitter)
    # print(elbo)
    if(abs(elbo-elbo_old)/n<tol){
      break
    }else{
      kernel_param_old = kernel_param
      elbo_old = elbo
    }
  }
  if((elbo-elbo_init)<0){
    return(kernel_param_init)
  }else{
    return(kernel_param)
  }
}

#'@title Obtain full posterior distribution
pois_sgp_get_posterior = function(X,y,sc,X_ind,kernel_param,kernel_func,post_mean,V,mu,Jitter){
  temp = pois_sgp_matrix_helper(X,y,sc,X_ind,kernel_param,kernel_func,post_mean,V,mu,Jitter,get_delta=F)
  v = temp$Knn_diag - rowSums((temp$Knm%*%t(temp$L_Kmm_inv))^2)  + rowSums((temp$A%*%temp$L_V)^2)
  return(list(mean=temp$Ab+mu,var=v,rate=exp(mu+temp$Ab+v/2),rate_v = (exp(v)-1)*exp(2*(temp$Ab+mu)+v)))
}

#'@title Helper function for matrix calculations
pois_sgp_matrix_helper = function(X,y,sc,X_ind,kernel_param,kernel_func,post_mean,V,mu,Jitter,get_delta=T){
  Knm = kernel_func(X,X_ind,kernel_param)
  Kmm = kernel_func(X_ind,X_ind,kernel_param,eps=Jitter)
  Knn_diag = kernel_func(X,X,kernel_param,diagonal=T)
  L_Kmm = t(chol(Kmm))
  L_Kmm_inv = forwardsolve(L_Kmm,diag(ncol(Kmm)))
  Kmm_inv = crossprod(L_Kmm_inv)
  A = Knm%*%Kmm_inv
  L_V = t(chol(V))
  Ab = drop(A%*%post_mean)
  if(get_delta){
    delta=sc*exp(Ab+rowSums((A%*%L_V)^2)/2+Knn_diag/2-rowSums((Knm%*%t(L_Kmm_inv))^2)/2)
  }else{
    delta=NULL
  }
  return(list(delta=delta,
              Knm=Knm,
              Kmm=Kmm,
              Knn_diag=Knn_diag,
              L_Kmm=L_Kmm,
              L_Kmm_inv=L_Kmm_inv,
              Kmm_inv=Kmm_inv,
              A=A,
              L_V=L_V,
              Ab=Ab))
}

