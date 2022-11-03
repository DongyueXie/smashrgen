#'@title generate Poisson sequence data
#'@param n_simu number of simulation reps
#'@param n length of Poisson seq
#'@param snr signal  to noise ratio
#'@param count_size max of exp(b)
#'@param smooth_func blocks, bumps, heavi,doppler
#'@return a list of
#'  \item{X:}{data}
#'  \item{L:}{mean}
#'  \item{other imputs:}{sigma2,snr, ...}
#'@import wavethresh
#'@export
sim_data_smooth = function(n_simu,n=2^9,snr=3,count_size,smooth_func='blocks',seed=12345){
  set.seed(seed)
  b = DJ.EX(n=n,signal=1,noisy=FALSE,plotfn = FALSE)[[smooth_func]]
  b = b - min(b)
  b = b/(max(b)/log(count_size))
  sigma2 = var(b)/snr
  X = matrix(nrow=n_simu,ncol=n)
  L = matrix(nrow=n_simu,ncol=n)
  for(i in 1:n_simu){
    l = exp(b+rnorm(n,0,sd=sqrt(sigma2)))
    L[i,] = l
    X[i,] = rpois(n,l)
  }
  return(list(X=X,latentX=L,b=b,snr=snr,sigma2=sigma2,count_size=count_size,smooth_func=smooth_func,seed=seed))
}

#'@title compare methods for smoothing Poisson sequence
#'@import parallel
#'@import vebpm
#'@import smashr
#'@export
simu_study_poisson_smooth = function(simdata,save_data=TRUE,
                                     method_list=c('vst','lik_exp','split_dwt',
                                                   'split_ndwt','smash','smash_two_step'),
                                     n_cores = 10,
                                     filter.number = 1,
                                     family='DaubExPhase',
                                     maxiter=100){
  n_simu = nrow(simdata$X)
  n = ncol(simdata$X)
  n_method = length(method_list)
  res = mclapply(1:n_simu,function(i){
    fitted_model = vector('list',n_method)
    names(fitted_model) = method_list

    if('vst'%in%method_list){
      res_vst = try(smash_gen_pois(simdata$X[i,],transformation='vst',method='smash',
                                   filter.number = filter.number,family = family,maxiter=1))
      fitted_model$vst = res_vst
    }
    if('lik_exp'%in%method_list){
      res_lik = try(smash_gen_pois(simdata$X[i,],transformation='lik_expansion',method='smash',
                                   filter.number = filter.number,family = family,maxiter=1))
      fitted_model$lik_exp = res_lik
    }
    if('split_dwt'%in%method_list){
      res_split = try(pois_smooth_split(simdata$X[i,],wave_trans='dwt',
                                        filter.number = filter.number,family = family,maxiter=maxiter))
      fitted_model$split_dwt = res_split
    }
    if('split_ndwt'%in%method_list){
      res_ndwt = try(pois_smooth_split(simdata$X[i,],wave_trans='ndwt',
                                       filter.number = filter.number,family = family,maxiter=maxiter))
      fitted_model$split_ndwt = res_ndwt
    }
    if('smash'%in%method_list){
      res_smash = try(smash.poiss(simdata$X[i,]))
      res_smash = list(posterior=list(mean_smooth=res_smash,mean_latent_smooth = log(res_smash)))
      fitted_model$smash = res_smash
    }
    if('smash_two_step'%in%method_list){
      res_smash2 = try(smash_two_step(simdata$X[i,]))
      fitted_model$smash_two_step = res_smash2
    }
    mse_smooth = simplify2array(lapply(fitted_model, function(x) {
      mse(x$posterior$mean_smooth, exp(simdata$b))}))
    mae_smooth = simplify2array(lapply(fitted_model, function(x) {
      mae(x$posterior$mean_smooth, exp(simdata$b))}))
    mse_latent_smooth = simplify2array(lapply(fitted_model, function(x) {
      mse(x$posterior$mean_latent_smooth,simdata$b)}))
    mae_latent_smooth = simplify2array(lapply(fitted_model, function(x) {
      mae(x$posterior$mean_latent_smooth,simdata$b)}))
    return(list(fitted_model=fitted_model,mse_smooth=mse_smooth,mae_smooth=mae_smooth,
                mse_latent_smooth=mse_latent_smooth,mae_latent_smooth=mae_latent_smooth))
  },mc.cores = n_cores
  )
  if(save_data){
    return(list(sim_data = simdata, output = res))
  }else{
    return(res)
  }
}
