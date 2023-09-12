#'@title Empirical Bayes wavelet smoothing via DWT
#'@description Smooth homogeneous Gaussian data.
#'@param x data
#'@param sigma known standard error
#'@param filter.number,family wavelet family and filter number as in `wavethresh` package
#'@param ebnm_params a list of `ebnm` parameters
#'@param W the dwt matrix for calc posterior variance. Remove the first row which is all 1/sqrt(n).
#'@return a list of
#'  \item{mu.est:}{posterior mean}
#'  \item{mu.var:}{posterior variance}
#'  \item{loglik:}{log likelihood}
#'  \item{dKL:}{KL divergence between g(the prior) and q(the posterior)}
#'@import wavethresh
#'@import ebnm
#'@export
smash_dwt = function(x,sigma,filter.number=1,
                     family="DaubExPhase",
                     ebnm_params=list(),W=NULL){

  n = length(x)
  J = log(n,2)
  if(ceiling(J)!=floor(J)){
    stop('Length of x must be power of 2')
  }
  if(filter.number==1&family=='DaubExPhase'){
    wavelet_name = "haar"
  }else{
    wavelet_name = 'non-haar'
  }
  ebnm_params = modifyList(ebnm_params_default_plr(),ebnm_params,keep.null =  TRUE)
  tsum = sum(x)/sqrt(n)
  x.w = wd(x, filter.number = filter.number,
           family = family, type = "wavelet")

  data.var = sigma^2
  if(length(data.var==1)){
    data.var = rep(data.var,n)
  }

  if(wavelet_name!='haar'){
    if(is.null(W)){
      W = (t(GenW(n,filter.number,family)))[-1,]
    }
  }


  if(length(sigma)==1){
    x.w.v = rep(sigma^2,n-1)
    tsum.var = sigma^2
  }else{
    x.w.v =  data.var
    tsum.var = x.w.v[1]
    x.w.v = x.w.v[-1]
  }

  dKL = 0
  loglik.scale = c()
  fitted_g = list()
  x.w.v.s = rep(0, 2^J-1)
  for (j in 0:(J - 1)) {
    x.pm = rep(0, 2^j)
    #index = (((J - 1) - j) * n + 1):((J - j) * n)
    index = (n-2^(j+1)+1):(n-2^j)
    x.w.j = accessD(x.w, j)
    x.w.v.j = x.w.v[index]
    ind.nnull = (x.w.v.j != 0)

    a = ebnm(x.w.j[ind.nnull],sqrt(x.w.v.j[ind.nnull]),
             mode=ebnm_params$mode,
             prior_family=ebnm_params$prior_family,
             scale = ebnm_params$scale,
             g_init = ebnm_params$g_init[[j+1]],
             fix_g = ebnm_params$fix_g,
             output = ebnm_params$output,
             optmethod = ebnm_params$optmethod,
             control = ebnm_params$control)

    dKL = dKL + a$log_likelihood - Eloglik(x.w.j[ind.nnull], sqrt(x.w.v.j[ind.nnull]),a$posterior$mean, a$posterior$mean^2+a$posterior$sd^2)
    x.pm[ind.nnull] = a$posterior$mean
    x.pm[!ind.nnull] = 0
    x.w = putD(x.w, j, x.pm)
    loglik.scale[j + 1] = a$log_likelihood
    fitted_g[[j+1]] = a$fitted_g
    x.w.v.s[index[ind.nnull]] = a$posterior$sd^2
    x.w.v.s[index[!ind.nnull]] = 0
  }
  mu.est = wr(x.w)
  loglik = sum(loglik.scale)
  #x.w.v.s = c(tsum.var,x.w.v.s)
  if(wavelet_name=='haar'){
    mu.est.var = haar_inv_var(c(x.w.v.s,0))
  }else{
    mu.est.var = colSums(W^2*x.w.v.s)
  }

  return(list(posterior=list(mean=mu.est,var=mu.est.var),
              fitted_g = fitted_g,
              loglik = loglik,
              dKL = dKL,
              x.w.v.s=x.w.v.s))
}



Eloglik = function(x, s, Et, Et2) {
  # Deal with infinite SEs:
  idx = is.finite(s)
  x = x[idx]
  s = s[idx]
  Et = Et[idx]
  Et2 = Et2[idx]
  return(-0.5 * sum(log(2*pi*s^2) + (1/s^2) * (Et2 - 2*x*Et + x^2)))
}

#'@importFrom ebnm ebnm_output_default
ebnm_params_default_plr = function(){
  return(list(prior_family='normal_scale_mixture',
              mode=0,
              scale = "estimate",
              g_init = NULL,
              fix_g = FALSE,
              output = ebnm_output_default(),
              optmethod = NULL,
              control = NULL))
}


haar = function(x,scale= sqrt(2)){
  if(length(x)==1){
    return(x)
  }
  else{
    x = matrix(x,nrow=2)
    diff = (x[1,]-x[2,])/scale
    sum = (x[1,]+x[2,])/scale
    return(c(diff, haar(sum)))
  }
}

haar_inv = function(x,scale=sqrt(2)){
  n=length(x)
  if(n==1){
    return(x)
  }
  x = matrix(scale*x,nrow=2,byrow=TRUE)
  smoothed = haar_inv(x[2,])
  return(as.vector(rbind(smoothed+x[1,], smoothed-x[1,]))/2)
}

haar_inv_var = function(v,scale=sqrt(2)){
  n=length(v)
  if(n==1){
    return(v)
  }
  v = matrix(scale^2*v,nrow=2,byrow=TRUE)
  smoothed = haar_inv_var(v[2,])
  return(as.vector(rbind(smoothed+v[1,], smoothed+v[1,]))/4)
  #return(rep((smoothed+v[1,])/4,each=2))
}
