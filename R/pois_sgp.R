pois_sgp = function(X,y,X_ind,
                    post_mean,V,
                    kernel_func=RBFkernel,
                    kernel_param=c(0,0),
                    mu=NULL,
                    s=NULL,
                    opt_method='L-BFGS-B',
                    fix_x_ind=F,fix_kernel_param=F,fix_mu=F,
                    maxiter=100,tol=1e-5,
                    l_b=-Inf,r_b=Inf){
  n = length(y)
  m = length(X_ind)
  n_kernel_param = length(kernel_param)
  elbo_tace = -Inf
  for(iter in 1:maxiter){
    # update post mean

    # update post var

    # update kernel and x_ind

    # update mu

    elbo = pois_sgp_elbo()
    elbo_tace[iter+1] = elbo
    if(elbo - elbo_tace[iter] < tol){
      break
    }

  }
  return(list())
}


pois_sgp_update_beta = function(){

  for(iter in 1:maxiter){
    grad = grad_beta()
    hessian = hessian_beta()
    direction = solve(hessian)%*%grad
    beta = beta - direction
  }
  return(beta)
}

grad_beta = function(){

}

hessian_beta = function(){

}


