
#'@import mr.ash
pois_reg_splitting = function(X,y,
                              fix_g=FALSE,
                              g_init=NULL,
                              maxiter=100,
                              tol=1e-5){
  t_start = Sys.time()
  n = nrow(X)
  p = ncol(X)
  XtX = crossprod(X)

  m = log(1+y)
  sigma2 = var(m)
  fit_mrash = fit_mr_ash(X,m,intercept = FALSE)
  sigma2 = fit.mrash$resid.sd^2
  fitted_mrash = fit_mrash$fitted

  obj = -Inf
  for(iter in 1:maxiter){
    opt=vga_pois_solver(m,y,s=1,beta=fitted_mrash,sigma2=sigma2)
    m = opt$m
    v = opt$v

    if(fix_g){
      fit0 = init_mr_ash(X,m,b=fit_mrash$b,prior.sd = g_init$prior.sd,prior.weights = g_init$prior.weights,resid.sd=sqrt(sigma2))
    }else{
      fit0 = init_mr_ash(X,m,b=fit_mrash$b,prior.sd = fit.mr.ash$prior$sd,prior.weights = fit.mr.ash$prior$weights,resid.sd=sqrt(sigma2))
    }

    fit_mrash = fit_mr_ash(X,m,fit0,standardize = FALSE,intercept = FALSE,
                         control = list(update.prior.weights=fix_g,update.resid.sd=FALSE))
    fitted_mrash = fit_mrash$fitted
    vxb = t(fit_mrash$b)%*%XtX%*%fit_mrash$b+sum(diag(XtX)*(get_mrash_pvar(fit_mrash)))

    H = fit_mrash$elbo - (-n/2*log(sigma2)-(sum(m^2)-sum(2*crossprod(m,X)*fit_mrash$b) + vxb)/2/sigma2)
    sigma2 = (sum(m^2+v-2*m*fitted_mrash+vxb))


    obj[iter+1] = pois_reg_splitting_obj(y,m,v,sigma2,fit_mrash$fitted,XtX,fit_mrash$b,get_mrash_pvar(fit_mrash),H,n)
    if((obj[iter+1]-obj[iter])>tol){
      break
    }

  }
  return(list(fit_mrash=fit_mrash,obj=obj))

}

pois_reg_splitting_obj = function(y,m,v,sigma2,beta,XtX,Eb,Vb,H,n){
  sum(y*m-exp(m+v/2)+log(v)/2) - n/2*log(sigma2) - sum(m^2+v-2*m*beta)/2/sigma2 - t(Eb)%*%XtX%*%Eb/2/sigma2 - sum(diag(XtX)*Vb)/2/sigma2+H
}

get_mrash_pvar = function(fit.mr.ash){
  rowSums(fit.mr.ash$phi*(fit.mr.ash$m^2+fit.mr.ash$s2))-fit.mr.ash$b^2
}
