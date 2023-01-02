#'@title Empirical Bayes wavelet smoothing for flashier
#'@description Smooth homogeneous Gaussian data.
#'@param x data
#'@param s known standard error
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
#'@importFrom smashr reflect
#'@export
dwt.fn = function(x, s, g_init, fix_g, output){

  filter.number=1
  family="DaubExPhase"
  ebnm_params=list()
  W=NULL

  n = length(x)
  J = log(n,2)
  n_orig = n
  if(ceiling(J)!=floor(J)){
    #stop('Length of x must be power of 2')
    # reflect
    x = reflect(x)
    idx = x$idx
    x = x$x
    n = length(x)
    J = log(n,2)
  }else{
    idx = 1:n
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

  data.var = s^2
  if(length(data.var==1)){
    data.var = rep(data.var,n)
  }else if(length(unique(data.var))==1){
    data.var = rep(data.var[1],n)
  }else{
    stop('sigma must be constant for all observations.')
  }

  if(wavelet_name!='haar'){
    if(is.null(W)){
      W = (t(GenW(n,filter.number,family)))[-1,]
    }
  }



  x.w.v =  data.var
  tsum.var = x.w.v[1]
  x.w.v = x.w.v[-1]


  dKL = 0
  loglik.scale = c()
  x.w.v.s = rep(0, 2^J-1)
  fitted_g = list()
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
             g_init = g_init[[j+1]],
             fix_g = fix_g,
             output = ebnm_params$output,
             optmethod = ebnm_params$optmethod,
             control = ebnm_params$control)
    fitted_g[[j+1]] = a$fitted_g
    dKL = dKL + a$log_likelihood - Eloglik(x.w.j[ind.nnull], sqrt(x.w.v.j[ind.nnull]),a$posterior$mean, a$posterior$mean^2+a$posterior$sd^2)
    x.pm[ind.nnull] = a$posterior$mean
    x.pm[!ind.nnull] = 0
    x.w = putD(x.w, j, x.pm)
    loglik.scale[j + 1] = a$log_likelihood
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

  return(list(posterior=data.frame(mean = mu.est[idx],second_moment=(mu.est^2+mu.est.var)[idx]),
              fitted_g=fitted_g,
              log_likelihood = loglik/n*n_orig))
}

