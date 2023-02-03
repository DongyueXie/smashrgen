#'@title this function fits empirical Bayes Sparse Poisson regression
#'@importFrom mr.ash fit_mr_ash
pois_reg_splitting = function(X,y,
                              fix_g=FALSE,
                              g_init=NULL,
                              maxiter=100,
                              maxiter_mrash=100,
                              printevery=10,
                              tol=1e-5){
  t_start = Sys.time()
  n = nrow(X)
  p = ncol(X)
  XtX = crossprod(X)

  m = log(1+y)
  fit_mrash = fit_mr_ash(X,m,intercept = FALSE,verbose = 'none')
  if(fix_g){
    fit_mrash = mr_ash(X,m,se=(fit_mrash$resid.sd)^2,s0=g_init$prior.sd^2,w0=g_init$prior.weights,b=fit_mrash$b,maxiter=maxiter_mrash,verbose=F,update.w0=F)
  }else{
    fit_mrash$prior$sd[1] = 0.01
    fit_mrash = mr_ash(X,m,se=(fit_mrash$resid.sd)^2,s0=fit_mrash$prior$sd^2,w0=fit_mrash$prior$weights,b=fit_mrash$b,maxiter=maxiter_mrash,verbose=F)
  }

  sigma2 = fit_mrash$se
  fitted_mrash = c(X%*%fit_mrash$b)

  obj = -Inf
  for(iter in 1:maxiter){
    opt=vga_pois_solver(m,y,s=1,beta=fitted_mrash,sigma2=sigma2)
    m = opt$m
    v = opt$v

    if(fix_g){
      fit_mrash = mr_ash(X,m,se=(sigma2),s0=g_init$prior.sd^2,w0=g_init$prior.weights,
                         b=fit_mrash$b,maxiter=maxiter_mrash,verbose=F,update.w0=F,update.se = F)
      #fit0 = init_mr_ash(X,m,b=fit_mrash$b,prior.sd = g_init$prior.sd,prior.weights = g_init$prior.weights,resid.sd=sqrt(sigma2))
    }else{
      fit_mrash = mr_ash(X,m,se=(sigma2),s0=fit_mrash$s0,w0=fit_mrash$w0,
                         b=fit_mrash$b,maxiter=maxiter_mrash,verbose=F,update.w0=T,update.se = F)
      #fit0 = init_mr_ash(X,m,b=fit_mrash$b,prior.sd = fit.mr.ash$prior$sd,prior.weights = fit.mr.ash$prior$weights,resid.sd=sqrt(sigma2))
    }


    fitted_mrash = c(X%*%fit_mrash$b)
    vxb = c(t(fit_mrash$b)%*%XtX%*%fit_mrash$b+sum(diag(XtX)*fit_mrash$vars))

    H = fit_mrash$elbo - (-n/2*log(sigma2)-(sum(m^2)-sum(2*crossprod(m,X)*fit_mrash$b) + vxb)/2/sigma2)
    sigma2 = (sum(m^2+v-2*m*fitted_mrash)+vxb)/n


    obj[iter+1] = pois_reg_splitting_obj(y,m,v,sigma2,fitted_mrash,XtX,fit_mrash$b,fit_mrash$vars,H,n)
    if((obj[iter+1]-obj[iter])<tol){
      break
    }
    if(iter%%printevery==0){
      print(paste('iter ',iter,'elbo =',obj[iter+1],sep=' '))
    }

  }
  return(list(fit_mrash=fit_mrash,obj=obj))

}

pois_reg_splitting_obj = function(y,m,v,sigma2,beta,XtX,Eb,Vb,H,n){
  sum(y*m-exp(m+v/2)+log(v)/2) - n/2*log(sigma2) - sum(m^2+v-2*m*beta)/2/sigma2 - t(Eb)%*%XtX%*%Eb/2/sigma2 - sum(diag(XtX)*Vb)/2/sigma2+H
}

# get_mrash_pvar = function(fit.mr.ash){
#   rowSums(fit.mr.ash$phi*(fit.mr.ash$m^2+fit.mr.ash$s2))-fit.mr.ash$b^2
# }
