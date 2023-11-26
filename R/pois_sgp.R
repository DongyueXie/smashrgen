
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
                 fix_X_ind=T,fix_kernel_param=fix_kernel_param,fix_mu=fix_mu,
                 maxiter=controls$maxiter,tol=controls$tol,
                 maxiter_mean=controls$maxiter_mean,tol_mean=controls$tol_mean,
                 maxiter_V=controls$maxiter_V,tol_V=controls$tol_V,
                 #l_b=controls$l_b,r_b=controls$r_b,
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
pois_sgp = function(y,sc=NULL,
                    X=NULL,X_ind=NULL,m=30,
                    kernel_func=Maternkernel,
                    kernel_param=NULL,
                    mu=NULL,
                    post_mean=NULL,V=NULL,
                    opt_method='L-BFGS-B',
                    fix_X_ind=T,fix_kernel_param=F,fix_mu=F,
                    maxiter=100,tol=1e-5,
                    maxiter_mean=100,tol_mean=1e-5,
                    maxiter_V=100,tol_V=1e-5,
                    Jitter=1e-5,
                    verbose=T,printevery=1){
  t0 = Sys.time()
  if(verbose){cat("initializing...");cat("\n")}
  n = length(y)
  if(is.null(X)){
    if(!is.null(X_ind)){
      X = seq(range(X_ind)[1],range(X_ind)[2],length.out=n)
    }else{
      X = seq(0,1,length.out=n)
    }
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
    ly = log(1+y/sc)
    temp = sgp(ly,X,X_ind,mu=NULL,fix_x_ind = T,kernel_func = kernel_func,
               sigma2 = sd_est_diff2(ly)^2,fix_sigma2 = TRUE)
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
  init_val = list(mu=mu,post_mean=post_mean,V=V,kernel_param=kernel_param)
  # if(l_b[1]==-Inf & opt_method=="L-BFGS-B"){
  #   l_b = c(rep(X_range[1],m),rep(-5,n_kernel_param))
  # }
  # if(r_b[1]==Inf & opt_method=="L-BFGS-B"){
  #   r_b=c(rep(X_range[2],m),rep(5,n_kernel_param))
  # }

  if(fix_X_ind){
    d_mat_nm = Rfast::Outer(X_ind,X,"-")
    d_mat_mm = Rfast::Outer(X_ind,X_ind,"-")
  }else{
    d_mat_nm=NULL
    d_mat_mm=NULL
  }

  elbo_tace = pois_sgp_elbo(X,y,sc,X_ind,kernel_param,kernel_func,post_mean,V,mu,Jitter,d_mat_nm,d_mat_mm)
  #print(paste("elbo init:",round(elbo_tace,3)))
  if(verbose){cat("running iterations...");cat("\n")}
  for(iter in 1:maxiter){
    # update post mean
    post_mean = pois_sgp_update_post_mean(X,y,sc,X_ind,kernel_param,kernel_func,post_mean,V,mu,maxiter=maxiter_mean,tol=tol_mean,Jitter,d_mat_nm,d_mat_mm,verbose)
    #print(paste("elbo after post_mean update:",round(pois_sgp_elbo(X,y,sc,X_ind,kernel_param,kernel_func,post_mean,V,mu,Jitter,d_mat_nm,d_mat_mm),3)))
    # update post var
    V = pois_sgp_update_V(X,y,sc,X_ind,kernel_param,kernel_func,post_mean,V,mu,maxiter=maxiter_V,tol=tol_V,Jitter,d_mat_nm,d_mat_mm,verbose)
    #print(paste("elbo after V update:",round(pois_sgp_elbo(X,y,sc,X_ind,kernel_param,kernel_func,post_mean,V,mu,Jitter,d_mat_nm,d_mat_mm),3)))
    # update kernel and x_ind
    # prior = pois_sgp_update_prior(X,y,sc,X_ind,kernel_param,kernel_func,post_mean,V,mu,fix_kernel_param,fix_X_ind,m,n_kernel_param,opt_method,l_b,r_b,Jitter,d_mat_nm,d_mat_mm,verbose)
    # X_ind = prior$X_ind
    # kernel_param = prior$kernel_param
    vari = pois_sgp_update_variational_param(X,y,sc,X_ind,kernel_param,kernel_func,post_mean,V,mu,
                                             fix_kernel_param,fix_mu,m,opt_method,Jitter,d_mat_nm,d_mat_mm,verbose)
    kernel_param = vari$kernel_param
    mu = vari$mu
    #print(paste("elbo after prior update:",round(pois_sgp_elbo(X,y,sc,X_ind,kernel_param,kernel_func,post_mean,V,mu,Jitter,d_mat_nm,d_mat_mm),3)))
    # update mu
    # mu = pois_sgp_update_mu(X,y,sc,X_ind,kernel_param,kernel_func,post_mean,V,mu,fix_mu,Jitter,d_mat_nm,d_mat_mm,verbose)
    #print(paste("elbo after mu update:",round(pois_sgp_elbo(X,y,sc,X_ind,kernel_param,kernel_func,post_mean,V,mu,Jitter,d_mat_nm,d_mat_mm),3)))

    elbo = pois_sgp_elbo(X,y,sc,X_ind,kernel_param,kernel_func,post_mean,V,mu,Jitter,d_mat_nm,d_mat_mm)
    elbo_tace[iter+1] = elbo
    if(verbose){
      if(iter%%printevery==0){
        cat(paste("At iter ",iter, ", elbo = ",elbo,sep=''))
        cat("\n")
      }
    }
    if((elbo - elbo_tace[iter])/n < tol){
      break
    }


  }
  post = pois_sgp_get_posterior(X,y,sc,X_ind,kernel_param,kernel_func,post_mean,V,mu,Jitter,d_mat_nm,d_mat_mm)
  return(list(elbo_tace=elbo_tace,elbo=elbo,init_val=init_val,
              fitted_g = list(X_ind=X_ind,kernel_param=kernel_param,mu=mu,kernel_func=kernel_func),
              posterior=list(mean=post$mean,sd=sqrt(post$v),rate=post$rate,rate_sd = sqrt(post$rate_v),mean_ind=post_mean,V_ind=V),
              run_time = difftime(Sys.time(),t0)))
}



#'@title Master function for optimizing variational parameters (Kernel_param, mu)
#'
pois_sgp_update_variational_param = function(X,y,sc,X_ind,kernel_param,kernel_func,post_mean,V,mu,
                                             fix_kernel_param,fix_mu,m,opt_method,Jitter,d_mat_nm,d_mat_mm,verbose){
  if(!(fix_kernel_param & fix_mu)){
    elbo_old = pois_sgp_elbo(X,y,sc,X_ind,kernel_param,kernel_func,post_mean,V,mu,Jitter,d_mat_nm,d_mat_mm)
    res = try(optim(c(kernel_param,mu),pois_sgp_update_variational_param_obj,
                    gr = pois_sgp_update_variational_param_obj_grad,
                    X=X,y=y,sc=sc,X_ind=X_ind,
                    kernel_func=kernel_func,post_mean=post_mean,V=V,
                    Jitter=Jitter,
                    d_mat_nm=d_mat_nm,d_mat_mm=d_mat_mm,
                    method = opt_method),silent = !verbose)

    if(class(res)!='try-error'){
      if((-res$value)>elbo_old){
        kernel_param = res$par[1:(length(res$par)-1)]
        mu = res$par[length(res$par)]
        refit = FALSE
      }else{
        refit = TRUE
      }
    }else{
      refit=TRUE
    }

    if(refit & opt_method!="Nelder-Mead"){
      if(verbose){cat("re-optimizing using Nelder-Mead");cat("\n")}
      res = try(optim(c(kernel_param,mu),pois_sgp_update_variational_param_obj,
                      X=X,y=y,sc=sc,X_ind=X_ind,
                      kernel_func=kernel_func,post_mean=post_mean,V=V,
                      Jitter=Jitter,
                      d_mat_nm=d_mat_nm,d_mat_mm=d_mat_mm,
                      method = "Nelder-Mead"),silent = !verbose)
      #################
      # if this also fails, may consider trust region method from trust or trust.Optim packages
      ################
      if(class(res)!='try-error'){
        if((-res$value)>elbo_old){
          kernel_param = res$par[1:(length(res$par)-1)]
          mu = res$par[length(res$par)]
          refit = FALSE
        }else{
          if(verbose){cat("Failed to increase ELBO")}
        }
      }
    }

  }
  return(list(kernel_param = kernel_param,mu=mu))
}

pois_sgp_update_variational_param_obj = function(vari_param,X,y,sc,X_ind,kernel_func,post_mean,V,Jitter,d_mat_nm,d_mat_mm){
  kernel_param = vari_param[1:(length(vari_param)-1)]
  mu = vari_param[length(vari_param)]
  obj = -pois_sgp_elbo(X,y,sc,X_ind,kernel_param,kernel_func,post_mean,V,mu,Jitter,d_mat_nm,d_mat_mm)
  return(obj)
}
pois_sgp_update_kernel_param_obj = function(kernel_param,mu,X,y,sc,X_ind,kernel_func,post_mean,V,Jitter,d_mat_nm,d_mat_mm){
  obj = -pois_sgp_elbo(X,y,sc,X_ind,kernel_param,kernel_func,post_mean,V,mu,Jitter,d_mat_nm,d_mat_mm)
  return(obj)
}
pois_sgp_update_variational_param_obj_grad = function(vari_param,X,y,sc,X_ind,kernel_func,post_mean,V,Jitter,d_mat_nm,d_mat_mm){
  kernel_param = vari_param[1:(length(vari_param)-1)]
  mu = vari_param[length(vari_param)]
  # temp = pois_sgp_matrix_helper(X,y,sc,X_ind,kernel_param,kernel_func,post_mean,V,mu,Jitter,get_delta=T,d_mat_nm,d_mat_mm)
  # d_mu = sum(y*temp$tilde1) - sum(temp$delta*temp$tilde1)+sum(temp$Kmm_inv%*%post_mean)-mu*sum(temp$Kmm_inv)
  # d_kernel_param = numDeriv::grad(pois_sgp_update_kernel_param_obj,kernel_param,
  #                                 mu=mu,X=X,y=y,sc=sc,X_ind=X_ind,kernel_func=kernel_func,
  #                                 post_mean=post_mean,V=V,Jitter=Jitter,d_mat_nm=d_mat_nm,d_mat_mm=d_mat_mm)
  d_vari_param = grad_kernel_param(kernel_param,X,y,sc,X_ind,kernel_func,post_mean,V,mu,Jitter,d_mat_nm,d_mat_mm)
  return(-d_vari_param)

}

#'@title Get derivative of log kernel parameters and mean
#'@importFrom Rfast rowsums
#'@importFrom Rfast mat.mult
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

  d_mu = sum(y*temp$tilde1) - sum(temp$delta*temp$tilde1)+sum(temp$Kmm_inv%*%post_mean)-mu*sum(temp$Kmm_inv)
  return(c(grad_log_c1+grad_log_c2,grad_log_theta1+grad_log_theta2,d_mu))
}



#'@title Function for updating posterior var
pois_sgp_update_V = function(X,y,sc,X_ind,kernel_param,kernel_func,post_mean,V,mu,maxiter=100,tol=1e-5,Jitter,d_mat_nm,d_mat_mm,verbose){
  V_old = V
  n = length(y)
  elbo_old = pois_sgp_elbo(X,y,sc,X_ind,kernel_param,kernel_func,post_mean,V,mu,Jitter,d_mat_nm,d_mat_mm)
  for(i in 1:maxiter){
    V = solve(-grad_post_mean(X,y,sc,X_ind,kernel_param,kernel_func,post_mean,V_old,mu,return_grad=F,return_hess=T,Jitter,d_mat_nm,d_mat_mm))
    #print(pois_sgp_elbo(X,y,s,X_ind,kernel_param,kernel_func,post_mean,V,mu,Jitter))
    #print(norm(V-V_old,type ='2'))
    elbo = pois_sgp_elbo(X,y,sc,X_ind,kernel_param,kernel_func,post_mean,V,mu,Jitter,d_mat_nm,d_mat_mm)
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
  if(verbose){
    if(i==maxiter){
      print("max iteration reached when updating posterior variance")
    }
  }

  if((elbo-elbo_old)<0){
    return(V_old)
  }else{
    return(V)
  }
}

#'@title Function for optimizing posterior mean
pois_sgp_update_post_mean = function(X,y,sc,X_ind,kernel_param,kernel_func,post_mean,V,mu,maxiter=100,tol=1e-5,Jitter,d_mat_nm,d_mat_mm,verbose){
  post_mean_old = post_mean
  elbo_old = pois_sgp_elbo(X,y,sc,X_ind,kernel_param,kernel_func,post_mean,V,mu,Jitter,d_mat_nm,d_mat_mm)
  n = length(y)
  for(iter in 1:maxiter){
    grad_hess = grad_post_mean(X,y,sc,X_ind,kernel_param,kernel_func,post_mean_old,V,mu,return_grad=T,return_hess=T,Jitter,d_mat_nm,d_mat_mm)
    grad = -grad_hess$grad
    hess = -grad_hess$hess
    # direction = c(pcgsolve(hess,post_mean))
    direction = c(solve(hess,grad))
    post_mean = post_mean_old - direction
    elbo = pois_sgp_elbo(X,y,sc,X_ind,kernel_param,kernel_func,post_mean,V,mu,Jitter,d_mat_nm,d_mat_mm)
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
  if(verbose){
    if(iter==maxiter){
      print("max iteration reached when updating posterior mean")
    }
  }

  if((elbo-elbo_old)<0){
    return(post_mean_old)
  }else{
    return(post_mean)
  }
}

#'@title Gradient and hessian of eblo w.r.t post_mean
grad_post_mean = function(X,y,sc,X_ind,kernel_param,kernel_func,post_mean,V,mu,return_grad=T,return_hess=T,Jitter,d_mat_nm,d_mat_mm){
  temp = pois_sgp_matrix_helper(X,y,sc,X_ind,kernel_param,kernel_func,post_mean,V,mu,Jitter,get_delta=T,d_mat_nm,d_mat_mm)
  if(return_grad&return_hess){
    grad = crossprod(temp$A,y-temp$delta)-temp$Kmm_inv%*%post_mean + mu*sum(temp$Kmm_inv)
    hess = -crossprod(temp$A*sqrt(temp$delta))-temp$Kmm_inv
    return(list(grad=grad,hess=hess))
  }else if(return_grad){
    grad = crossprod(temp$A,y-temp$delta)-temp$Kmm_inv%*%post_mean + mu*sum(temp$Kmm_inv)
    return(grad)
  }else if(return_hess){
    hess = -crossprod(temp$A*sqrt(temp$delta))-temp$Kmm_inv
    return(hess)
  }else{
    return(NULL)
  }

}

#'@title calc elbo
pois_sgp_elbo = function(X,y,sc,X_ind,kernel_param,kernel_func,post_mean,V,mu,Jitter,d_mat_nm,d_mat_mm){
  m = length(X_ind)
  temp = pois_sgp_matrix_helper(X,y,sc,X_ind,kernel_param,kernel_func,post_mean,V,mu,Jitter,get_delta=T,d_mat_nm,d_mat_mm)
  F1 = sum(y*temp$Ab)+mu*sum(y*temp$tilde1)-sum(temp$delta)+sum(y*log(sc))-sum(lfactorial(y))
  F2 = -sum(log(diag(temp$L_Kmm))) - sum((temp$L_Kmm_inv%*%post_mean)^2)/2 - sum((temp$L_Kmm_inv%*%temp$L_V)^2)/2+sum(log(diag(temp$L_V))) + m*0.5 + mu*sum(temp$Kmm_inv%*%post_mean) - mu^2*sum(temp$Kmm_inv)/2
  elbo=drop(F1+F2)
  # if(elbo>1e6){
  #   browser()
  # }
  return(elbo)
}

#'@title Obtain full posterior distribution
pois_sgp_get_posterior = function(X,y,sc,X_ind,kernel_param,kernel_func,post_mean,V,mu,Jitter,d_mat_nm,d_mat_mm){
  temp = pois_sgp_matrix_helper(X,y,sc,X_ind,kernel_param,kernel_func,post_mean,V,mu,Jitter,get_delta=F,d_mat_nm,d_mat_mm)
  v = temp$Knn_diag - Rfast::rowsums(Rfast::Tcrossprod(temp$Knm,temp$L_Kmm_inv)^2)  + Rfast::rowsums(Rfast::mat.mult(temp$A,temp$L_V)^2)
  #v = temp$Knn_diag - rowSums((temp$Knm%*%t(temp$L_Kmm_inv))^2)  + rowSums((temp$A%*%temp$L_V)^2)
  return(list(mean=temp$Ab+mu*temp$tilde1,var=v,rate=exp(temp$Ab+mu*temp$tilde1+v/2),rate_v = (exp(v)-1)*exp(2*(temp$Ab+mu*temp$tilde1)+v)))
}

#'@title Helper function for matrix calculations
pois_sgp_matrix_helper = function(X,y,sc,X_ind,kernel_param,kernel_func,post_mean,V,mu,Jitter,get_delta=T,d_mat_nm,d_mat_mm){
  Knm = kernel_func(X,X_ind,kernel_param,d_mat=d_mat_nm)
  Kmm = kernel_func(X_ind,X_ind,kernel_param,eps=Jitter,d_mat=d_mat_mm)
  Knn_diag = kernel_func(X,X,kernel_param,diagonal=T)
  L_Kmm = t(Rfast::cholesky(Kmm))
  L_Kmm_inv = forwardsolve(L_Kmm,diag(ncol(Kmm)))
  Kmm_inv = crossprod(L_Kmm_inv)
  A = Rfast::mat.mult(Knm,Kmm_inv)
  L_V = t(Rfast::cholesky(V))
  Ab = drop(A%*%post_mean)
  tilde1 = 1-Rfast::rowsums(A)
  if(get_delta){
    delta=sc*exp(Ab+mu*tilde1+Rfast::rowsums(A*(Rfast::mat.mult(A,V)-Knm))/2+Knn_diag/2)
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
              Ab=Ab,
              tilde1=tilde1))
}
