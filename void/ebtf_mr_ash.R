#'@title Empirical bayes Trend filtering via mr ash
#'@param x observed Poisson sequence
#'@return
#'@import mr.ash

ebtf_mr_ash = function(x,k){
  n = length(x)
  H = getH(n,k,remove_col1 = T)
  fit = fit_mr_ash(H,y,fit0 = init_mr_ash(H, y,init.method = 'glmnet'),intercept = T)
  y_hat_var = diag(H%*%diag(c(rowSums(fit$phi*(fit$m^2+fit$s2)) - rowSums(fit$phi*fit$m)^2))%*%t(H))
  return(list(posterior = list(mean = fit$b,var = y_hat_var),fit_mrash = fit))
}


getH = function(n,k,sparse=FALSE,remove_col1 =FALSE){

  if(k==0){
    H = matrix(0,n,n)
    H[lower.tri(H,diag = T)] = 1
  }else{
    H = matrix(0,n,n)
    A = list()
    for(i in 1:(k+1)){
      A[[i]] = rep(1,n)
    }
    for(i in 2:(k+1)){
      A[[i]] = cumsum(A[[i-1]])
    }
    for(j in 1:k){
      H[j:n,j] = A[[j]][1:(n-j+1)]
    }
    for(j in (k+1):n){
      H[j:n,j] = A[[k+1]][1:(n-j+1)]
    }
  }
  if(remove_col1){
    H = H[,-1]
  }
  if(sparse){
    H = Matrix::Matrix(H,sparse = T)
  }else{
    H
  }
}
