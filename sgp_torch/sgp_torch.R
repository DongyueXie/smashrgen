#'@title sparse GP via torch R
#'@import torch
#'@importFrom Rfast Outer
#'@export
sgp_torch = function(y,
                     X=NULL,X_ind=NULL,m=10,
                     kernel_func=Maternkernel_torch,
                     kernel_param=c(0,0),
                     mu=NULL,
                     s2=NULL,sigma2=NULL,
                     opt_method='L-BFGS-B',
                     fix_X_ind=T,
                     fix_sigma2=F,
                     Jitter=1e-5){
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
  if(is.null(s2)){s2 = rep(1,n)}
  if(length(s2)==1){s2 = rep(s2,n)}

  if(fix_X_ind){
    d_mat_nm = torch_tensor(Outer(X_ind,X,"-"),dtype = torch_float64())
    d_mat_mm = torch_tensor(Outer(X_ind,X_ind,"-"),dtype = torch_float64())
  }else{
    d_mat_nm=NULL
    d_mat_mm=NULL
  }
  if(is.null(sigma2)){
    sigma2 = sd_est_diff2(y)^2
    fix_sigma2 = F
  }
  if(is.null(kernel_param)){
    kernel_param=c(0,0)
  }
  if(kernel_param[1]==0){
    kernel_param = c(log(max(var(y)-sigma2,0)),0)
  }
  arg = list(y=torch_tensor(y,dtype = torch_float64()),
              X = torch_tensor(X,dtype = torch_float64()),
              X_ind = torch_tensor(X_ind,dtype = torch_float64()),
              kernel_func=kernel_func,
              mu=NULL,
              s2=torch_unsqueeze(torch_tensor(s2,dtype = torch_float64()),dim=2),
              sigma2 = torch_tensor(sigma2,dtype = torch_float64()),
              fix_sigma2=fix_sigma2,
              fix_X_ind=fix_X_ind,
              Jitter=torch_tensor(Jitter,dtype = torch_float64()),
             d_mat_nm=d_mat_nm,
             d_mat_mm=d_mat_mm)
  sgp_obj_torch_wrapper = function(all_params){
    loss = sgp_obj_torch(torch_tensor(all_params,dtype = torch_float64()),arg)
    return(as.numeric(loss))
  }
  sgp_obj_grad_torch = function(all_params){
    all_params = torch_tensor(all_params,requires_grad = T,dtype = torch_float64())
    loss = sgp_obj_torch(all_params,arg)
    grad = autograd_grad(loss,all_params)[[1]]
    return(as.numeric(grad))
  }
  pars = optim(c(kernel_param,log(sigma2),X_ind),fn=sgp_obj_torch_wrapper,gr = sgp_obj_grad_torch,method = opt_method)
  post = sgp_get_posterior_torch(torch_tensor(pars$par,dtype = torch_float64()),arg)
  return(list(elbo=-pars$value*n,
              posterior=list(mean=as.numeric(post$mean),sd=as.numeric(sqrt(post$v)),mean_ind=as.numeric(post$mean_ind),V_ind = as.matrix(post$V_ind)),
              fitted_g=list(mu=as.numeric(post$mu),X_ind=as.numeric(pars$par[4:length(pars$par)]),kernel_param=as.numeric(pars$par[1:2]),sigma2=as.numeric(exp(pars$par[3])),kernel_func=kernel_func),
              run_time = difftime(Sys.time(),t0),
              opt_res=pars))
}



sgp_matrix_helper_torch = function(X_ind,kernel_param,log_sigma2,mu,s2,y,X,kernel_func,Jitter,d_mat_nm=NULL,d_mat_mm=NULL){
  sigma2 = torch_exp(log_sigma2)
  Knm = kernel_func(X,X_ind,kernel_param,d_mat=d_mat_nm)
  Kmm = kernel_func(X_ind,X_ind,kernel_param,eps=Jitter,d_mat=d_mat_mm)
  Knn_diag = kernel_func(X,X,kernel_param,diagonal=T)
  L_Kmm = torch_cholesky(Kmm)
  L_Kmm_inv = torch_triangular_solve(torch_eye(ncol(Kmm),dtype = torch_float64()),L_Kmm,upper = F)[[1]]
  Kmm_inv = torch_matmul(torch_t(L_Kmm_inv),L_Kmm_inv)
  # Kmm_inv = torch_inverse(Kmm)
  # Kmm_inv = torch_cholesky_inverse(L_Kmm)
  A = torch_matmul(Knm,Kmm_inv)
  ASA = torch_matmul(torch_t(A/s2),A)/sigma2
  V = torch_inverse(ASA + Kmm_inv)
  L_V = torch_cholesky(V)
  colsumAS = torch_sum(A/s2,dim=1)
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


sgp_get_posterior_torch = function(all_params,arg){
  n = length(arg$y)
  kernel_param = all_params[1:2]
  if(arg$fix_sigma2){
    log_sigma2 = torch_log(torch_tensor(arg$sigma2,dtype = torch_float64()))
  }else{
    log_sigma2 = all_params[3]
  }
  if(arg$fix_X_ind){
    X_ind = arg$X_ind
  }else{
    X_ind = all_params[4:length(all_params)]
  }
  sigma2 = torch_exp(log_sigma2)
  temp = sgp_matrix_helper_torch(X_ind,kernel_param,log_sigma2,arg$mu,arg$s2,arg$y,arg$X,arg$kernel_func,arg$Jitter,arg$d_mat_nm,arg$d_mat_mm)
  tilde1 = torch_ones(n,dtype = torch_float64()) - torch_sum(temp$A,dim=2)
  a = (torch_sum(arg$y/arg$s2$squeeze()*tilde1))/sigma2/(torch_sum(tilde1/arg$s2$squeeze()*tilde1)/sigma2+torch_sum(temp$Kmm_inv))
  b = (torch_sum(temp$Kmm_inv,dim=2)-torch_matmul(torch_t(temp$A/arg$s2),tilde1)/sigma2)/(torch_sum(tilde1/arg$s2$squeeze()*tilde1)/sigma2+torch_sum(temp$Kmm_inv))
  Kmm_inv1 = torch_sum(temp$Kmm_inv,dim=2)
  ADtilde1 = torch_matmul(torch_t(temp$A/arg$s2),tilde1)
  # tilde1Dtilde1 = torch_sum(tilde1^2/arg$s2)
  if(is.null(arg$mu)){
    # beta = solve(temp$ASA + temp$Kmm_inv + tcrossprod(crossprod(temp$A/s2,tilde1),b)/sigma2- tcrossprod(rowSums(temp$Kmm_inv),b), crossprod(temp$A/s2,y)/sigma2-crossprod(temp$A/s2,tilde1)*a/sigma2+rowSums(temp$Kmm_inv)*a)
    beta = torch_matmul(torch_inverse(temp$ASA + temp$Kmm_inv + torch_matmul(ADtilde1$view(c(-1,1)),b$view(c(1,-1)))/sigma2- torch_matmul(Kmm_inv1$view(c(-1,1)),b$view(c(1,-1)))), torch_matmul(torch_t(temp$A/arg$s2),arg$y)/sigma2-ADtilde1*a/sigma2+Kmm_inv1*a)
    Ab = torch_matmul(temp$A,beta)
    mu = a + torch_sum(b*beta)
  }else{
    # beta = temp$V%*%(torch_cross(temp$A/arg$s2,arg$y)/sigma2+arg$mu*(Kmm_inv1-ADtilde1/sigma2))
    # Ab = temp$A%*%beta
    # mu = arg$mu
    stop('not implemented yet')
  }
  v = temp$Knn_diag - torch_sum(torch_matmul(temp$Knm,torch_t(temp$L_Kmm_inv))^2,dim=2)  + torch_sum(torch_matmul(temp$A,temp$L_V)^2,dim=2)
  return(list(mean=Ab+mu*tilde1,var=v,mu=mu,mean_ind=beta,V_ind=temp$V))
}



#'@param all_params kernel params,log_sigma2,X_ind,args
#'@param arg other args input as a list
sgp_obj_torch = function(all_params,arg){

  #X_ind,kernel_param,log_sigma2,mu,s2,y,X,kernel_func,Jitter,d_mat_nm=NULL,d_mat_mm=NULL
  kernel_param = all_params[1:2]
  if(arg$fix_sigma2){
    log_sigma2 = torch_log(torch_tensor(arg$sigma2,dtype = torch_float64()))
  }else{
    log_sigma2 = all_params[3]
  }
  if(arg$fix_X_ind){
    X_ind = arg$X_ind
  }else{
    X_ind = all_params[4:length(all_params)]
  }
  n = length(arg$y)
  m = length(arg$X_ind)
  sigma2 = torch_exp(log_sigma2)
  temp = sgp_matrix_helper_torch(X_ind,kernel_param,log_sigma2,arg$mu,arg$s2,arg$y,arg$X,arg$kernel_func,arg$Jitter,arg$d_mat_nm,arg$d_mat_mm)
  tilde1 = torch_ones(n,dtype = torch_float64()) - torch_sum(temp$A,dim=2)
  Kmm_inv1 = torch_sum(temp$Kmm_inv,dim=2)
  ADtilde1 = torch_matmul(torch_t(temp$A/arg$s2),tilde1)
  tilde1Dtilde1 = torch_sum(tilde1^2/arg$s2$squeeze())
  sum_Kmminv = torch_sum(temp$Kmm_inv)
  if(is.null(arg$mu)){
    a = (sum(arg$y/arg$s2$squeeze()*tilde1))/sigma2/(tilde1Dtilde1/sigma2+sum_Kmminv)
    b = (Kmm_inv1-ADtilde1/sigma2)/(tilde1Dtilde1/sigma2+sum_Kmminv)
    # beta = solve(temp$ASA + temp$Kmm_inv + tcrossprod(ADtilde1,b)/sigma2- tcrossprod(Kmm_inv1,b), crossprod(temp$A/s2,y)/sigma2-ADtilde1*a/sigma2+Kmm_inv1*a)
    beta = torch_matmul(torch_inverse(temp$ASA + temp$Kmm_inv + torch_matmul(ADtilde1$view(c(-1,1)),b$view(c(1,-1)))/sigma2- torch_matmul(Kmm_inv1$view(c(-1,1)),b$view(c(1,-1)))), torch_matmul(torch_t(temp$A/arg$s2),arg$y)/sigma2-ADtilde1*a/sigma2+Kmm_inv1*a)
    Ab = torch_matmul(temp$A,beta)
    mu = a + torch_sum(b*beta)
  }else{
    # beta = temp$V%*%(torch_cross(temp$A/arg$s2,arg$y)/sigma2+arg$mu*(Kmm_inv1-ADtilde1/sigma2))
    # Ab = temp$A%*%beta
    # mu = arg$mu
    stop('Not implemented yet')
  }
  F1=-n/2*torch_log(2*pi*sigma2)-torch_sum(torch_log(arg$s2))/2-(torch_sum(arg$y^2/arg$s2$squeeze())-2*torch_sum(arg$y/arg$s2$squeeze()*Ab)-2*mu*torch_sum(arg$y/arg$s2$squeeze()*tilde1)+torch_sum(Ab^2/arg$s2$squeeze())+2*mu*torch_sum(Ab/arg$s2$squeeze()*tilde1)+mu^2*tilde1Dtilde1+torch_sum(torch_matmul(temp$A/torch_sqrt(arg$s2),temp$L_V)^2)+torch_sum(temp$Knn_diag/arg$s2$squeeze())-torch_sum(torch_matmul(temp$Knm/torch_sqrt(arg$s2),torch_t(temp$L_Kmm_inv))^2))/2/sigma2
  F2 = -torch_sum(torch_log(torch_diag(temp$L_Kmm))) - torch_sum(torch_matmul(temp$L_Kmm_inv,beta)^2)/2 - torch_sum(torch_matmul(temp$L_Kmm_inv,temp$L_V)^2)/2+torch_sum(torch_log(torch_diag(temp$L_V))) + m*0.5 + mu*torch_sum(torch_matmul(temp$Kmm_inv,beta)) - mu^2*sum_Kmminv/2
  elbo=F1+F2
  return(-elbo/n)
}


#'@title Matern Kernel
#'@param d_mat this is x_i - x_j
#'@import torch
#'@export

Maternkernel_torch = function(X1,X2,kernel_param,eps=NULL,diagonal=F,nu=3/2,d_mat=NULL){
  coeff = torch_exp(kernel_param[1])
  scales = torch_exp(kernel_param[2])
  if(length(dim(X1))==1){X1=X1$view(c(-1,1))}
  if(length(dim(X2))==1){X2=X2$view(c(-1,1))}
  n1 = nrow(X1)
  if(diagonal){
    return(coeff*torch_ones(n1,dtype = torch_float64()))
  }else{
    if(!is.null(d_mat)){
      abs_dist = torch_abs(d_mat)
    }else{
      X1 = torch_unsqueeze(X1,dim=2)
      X2 = torch_unsqueeze(X2,dim=1)
      abs_dist <- torch_sum(torch_abs(X1-X2),dim=3)
    }

    if(nu==1/2){
      R = torch_exp(-abs_dist/scales)
    }else if(nu==3/2){
      R = (1+torch_sqrt(3.0)*abs_dist/scales)*torch_exp(-torch_sqrt(3.0)*abs_dist/scales)
    }else if(nu==5/2){
      R = (1+torch_sqrt(5.0)*abs_dist/scales+5*abs_dist^2/3/scales^2)*torch_exp(-torch_sqrt(5.0)*abs_dist/scales)
    }else{
      stop('nu can only be 1/3, 3/2, 5/2')
    }
    if(!is.null(eps)){
      R = R+eps*torch_eye(n1,dtype = torch_float64())
    }
    return(coeff*R)
  }
}
