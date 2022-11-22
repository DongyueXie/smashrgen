#'@title generate Poisson sequence data
#'@param n_simu number of simulation reps
#'@param n length of Poisson seq
#'@param snr signal  to noise ratio
#'@param count_size max of exp(b)
#'@param smooth_func cblocks, sblocks, heavi, angles, bursts, spike
#'@return a list of
#'  \item{X:}{data}
#'  \item{L:}{mean}
#'  \item{other imputs:}{sigma2,snr, ...}
#'@import wavethresh
#'@export
sim_data_smooth = function(n_simu,n=2^9,snr=3,count_size,smooth_func='blocks',seed=12345){
  set.seed(seed)
  #b = DJ.EX(n=n,signal=1,noisy=FALSE,plotfn = FALSE)[[smooth_func]]
  t = seq(0,1,length.out = n)
  if(smooth_func=='sblocks'){
    b = sblocks.f(t)
  }
  if(smooth_func=='cblocks'){
    b = cblocks.f(t)
  }
  if(smooth_func=='heavi'){
    b = heavi.f(t)
  }
  if(smooth_func=='angles'){
    b = angles.f(t)
  }
  if(smooth_func=='bursts'){
    b = bursts.f(t)
  }
  if(smooth_func=='spike'){
    b = spike.f(t)
  }

  #b = b - min(b)
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
#'
#'@export

simu_study_poisson_smooth = function(simdata,save_data=TRUE,
                                     method_list=c('vst_nug',
                                                   'vst_nug_top',
                                                   'vst_smooth',
                                                   'lik_exp_mle_nug',
                                                   'lik_exp_mle_nug_top',
                                                   'lik_exp_mle_smooth',
                                                   'lik_exp_smashpoi_nug',
                                                   'lik_exp_smashpoi_nug_top',
                                                   'lik_exp_smashpoi_smooth',
                                                   'split_runmed_dwt',
                                                   'split_runmed_ndwt',
                                                   'split_runmed_ndwt_top',
                                                   'split_log1px_ndwt',
                                                   'split_runmed_ndwt_fix_nug',
                                                   'split_log1px_ndwt_fix_nug',
                                                   'split_true',
                                                   'smash_two_step_homo',
                                                   'smash_two_step_hetero',
                                                   'smash'),
                                     smoother = 'smash',
                                     n_cores = 20,
                                     filter.number = 1,
                                     family='DaubExPhase',
                                     nug.est.limit = 0.2,
                                     maxiter=100){
  n_simu = nrow(simdata$X)
  n = ncol(simdata$X)
  n_method = length(method_list)
  res = mclapply(1:n_simu,function(i){
    fitted_model = vector('list',n_method)
    names(fitted_model) = method_list

    if('vst_nug'%in%method_list){
      res_vst_nug = try(smash_gen_pois(simdata$X[i,],transformation='vst',smoother=smoother,
                                   filter.number = filter.number,family = family,est_nugget_maxiter = maxiter,est_nugget = TRUE))
      fitted_model$vst_nug = res_vst_nug
    }

    if('vst_nug_top'%in%method_list){
      fitted_model$vst_nug_top = try(smash_gen_pois(simdata$X[i,],transformation='vst',smoother=smoother,
                                       filter.number = filter.number,family = family,
                                       est_nugget_maxiter = maxiter,est_nugget = TRUE,nug.est.limit=nug.est.limit))
    }

    if('vst_smooth'%in%method_list){
      res_vst_smooth = try(smash_gen_pois(simdata$X[i,],transformation='vst',smoother=smoother,
                                       filter.number = filter.number,family = family,est_nugget_maxiter = maxiter,est_nugget = FALSE))
      fitted_model$vst_smooth = res_vst_smooth
    }

    if('lik_exp_mle_nug'%in%method_list){
      res_lik_exp_mle_nug = try(smash_gen_pois(simdata$X[i,],transformation='lik_expan',smoother=smoother,
                                   filter.number = filter.number,family = family,est_nugget_maxiter=maxiter,lik_expan_at = 'mle',est_nugget = TRUE))
      fitted_model$lik_exp_mle_nug = res_lik_exp_mle_nug
    }
    if('lik_exp_mle_nug_top'%in%method_list){
      fitted_model$lik_exp_mle_nug_top = try(smash_gen_pois(simdata$X[i,],transformation='lik_expan',smoother=smoother,
                                                            filter.number = filter.number,family = family,
                                                            est_nugget_maxiter=maxiter,
                                                            lik_expan_at = 'mle',est_nugget = TRUE,nug.est.limit=nug.est.limit))
    }
    if('lik_exp_mle_smooth'%in%method_list){
      res_lik_exp_mle_smooth = try(smash_gen_pois(simdata$X[i,],transformation='lik_expan',smoother=smoother,
                                               filter.number = filter.number,family = family,est_nugget_maxiter=maxiter,lik_expan_at = 'mle',est_nugget = FALSE))
      fitted_model$lik_exp_mle_smooth = res_lik_exp_mle_smooth
    }
    if('lik_exp_smashpoi_nug'%in%method_list){
      res_lik_exp_smashpoi_nug = try(smash_gen_pois(simdata$X[i,],transformation='lik_expan',smoother=smoother,
                                               filter.number = filter.number,family = family,est_nugget_maxiter=maxiter,lik_expan_at = 'smash_poi',est_nugget = TRUE))
      fitted_model$lik_exp_smashpoi_nug = res_lik_exp_smashpoi_nug
    }
    if('lik_exp_smashpoi_nug_top'%in%method_list){
      fitted_model$lik_exp_smashpoi_nug_top = try(smash_gen_pois(simdata$X[i,],transformation='lik_expan',smoother=smoother,
                                                                 filter.number = filter.number,
                                                                 family = family,est_nugget_maxiter=maxiter,
                                                                 lik_expan_at = 'smash_poi',est_nugget = TRUE,
                                                                 nug.est.limit=nug.est.limit))
    }
    if('lik_exp_smashpoi_smooth'%in%method_list){
      res_lik_exp_smashpoi_smooth = try(smash_gen_pois(simdata$X[i,],transformation='lik_expan',smoother=smoother,
                                                    filter.number = filter.number,family = family,est_nugget_maxiter=maxiter,lik_expan_at = 'smash_poi',est_nugget = FALSE))
      fitted_model$lik_exp_smashpoi_smooth = res_lik_exp_smashpoi_smooth
    }
    if('lik_exp_iter_homo'%in%method_list){
      res_lik_exp_iter_homo = try(smash_gen_pois_iterative(simdata$X[i,],nugget_type = 'homoskedastic',
                                                           lik_expan_init ='mle',smoother = smoother,filter.number = filter.number,family = family,
                                                           maxiter=maxiter))
      fitted_model$lik_exp_iter_homo = res_lik_exp_iter_homo
    }
    if('lik_exp_iter_hetero'%in%method_list){
      res_lik_exp_iter_hetero = try(smash_gen_pois_iterative(simdata$X[i,],nugget_type = 'heteroskedastic',
                                                           lik_expan_init ='mle',smoother = smoother,filter.number = filter.number,family = family,
                                                           maxiter=maxiter))
      fitted_model$lik_exp_iter_hetero = res_lik_exp_iter_hetero
    }


    if('split_runmed_dwt'%in%method_list){
      fitted_model$split_runmed_dwt = try(pois_smooth_split(simdata$X[i,],wave_trans='dwt',Emu_init='runmed',
                                                            filter.number = filter.number,family = family,maxiter=maxiter))
    }
    if('split_smoothgaus_dwt'%in%method_list){
      fitted_model$split_smoothgaus_dwt = try(pois_smooth_split(simdata$X[i,],wave_trans='dwt',Emu_init='smooth_gaus',
                                                            filter.number = filter.number,family = family,maxiter=maxiter))
    }
    if('split_runmed_ndwt'%in%method_list){
      fitted_model$split_runmed_ndwt = try(pois_smooth_split(simdata$X[i,],wave_trans='ndwt',Emu_init='runmed',
                                                             filter.number = filter.number,family = family,maxiter=maxiter))
    }
    if('split_runmed_ndwt_top'%in%method_list){
      fitted_model$split_runmed_ndwt_top = try(pois_smooth_split(simdata$X[i,],wave_trans='ndwt',Emu_init='runmed',
                                                             filter.number = filter.number,family = family,maxiter=maxiter,
                                                             sigma2_est_top = nug.est.limit))
    }
    if('split_log1px_ndwt'%in%method_list){
      fitted_model$split_log1px_ndwt = try(pois_smooth_split(simdata$X[i,],wave_trans='ndwt',Emu_init='log1px',
                                                                 filter.number = filter.number,family = family,maxiter=maxiter))
    }
    if('split_log1px_ndwt_fix_nug'%in%method_list){
      fitted_model$split_log1px_ndwt_fix_nug = try(pois_smooth_split(simdata$X[i,],wave_trans='ndwt',Emu_init='log1px',
                                                             filter.number = filter.number,family = family,maxiter=maxiter,
                                                             est_sigma2 = FALSE, sigma2_init = simdata$sigma2))
    }
    if('split_true'%in%method_list){
      fitted_model$split_true = try(pois_smooth_split(simdata$X[i,],wave_trans='ndwt',Emu_init=simdata$b,
                                                      filter.number = filter.number,family = family,maxiter=maxiter,
                                                      est_sigma2 = FALSE, sigma2_init = simdata$sigma2))
    }
    if('split_runmed_ndwt_fix_nug'%in%method_list){
      fitted_model$split_runmed_ndwt_fix_nug = try(pois_smooth_split(simdata$X[i,],wave_trans='ndwt',Emu_init='runmed',
                                                                     filter.number = filter.number,family = family,maxiter=maxiter,
                                                                     est_sigma2 = FALSE, sigma2_init = simdata$sigma2))
    }
    if('split_smoothgaus_ndwt'%in%method_list){
      fitted_model$split_smoothgaus_ndwt = try(pois_smooth_split(simdata$X[i,],wave_trans='ndwt',Emu_init='smooth_gaus',
                                                             filter.number = filter.number,family = family,maxiter=maxiter))
    }
    if('smash'%in%method_list){
      t1 = Sys.time()
      res_smash = try(smash.poiss(simdata$X[i,]))
      res_smash = list(posterior=list(mean_smooth=res_smash,mean_latent_smooth = log(res_smash)),run_time = difftime(Sys.time(),t1,unit='secs'))
      fitted_model$smash = res_smash
    }

    if('smash_two_step_homo'%in%method_list){
      res_smash2 = try(smash_two_step(simdata$X[i,],homoskedastic=TRUE))
      fitted_model$smash_two_step_homo = res_smash2
    }
    if('smash_two_step_hetero'%in%method_list){
      res_smash3 = try(smash_two_step(simdata$X[i,],homoskedastic=FALSE))
      fitted_model$smash_two_step_hetero = res_smash3
    }
    mse_smooth = simplify2array(lapply(fitted_model, function(x) {
      if(class(x)!='try-error'){
        mse(x$posterior$mean_smooth, exp(simdata$b))
      }else{
        NA
      }
      }))
    mae_smooth = simplify2array(lapply(fitted_model, function(x) {
      if(class(x)!='try-error'){
        mae(x$posterior$mean_smooth, exp(simdata$b))
      }else{
        NA
      }
      }))
    mse_latent_smooth = simplify2array(lapply(fitted_model, function(x) {
      if(class(x)!='try-error'){
        mse(x$posterior$mean_latent_smooth,simdata$b)
      }else{
        NA
      }
      }))
    mae_latent_smooth = simplify2array(lapply(fitted_model, function(x) {
      if(class(x)!='try-error'){
        mae(x$posterior$mean_latent_smooth,simdata$b)
      }else{
        NA
      }
      }))

    run_times = simplify2array(lapply(fitted_model,function(x){
      if(class(x)!='try-error'){
        x$run_time
      }else{
        NA
      }
    }))

    return(list(fitted_model=fitted_model,mse_smooth=mse_smooth,mae_smooth=mae_smooth,run_times=run_times,
                mse_latent_smooth=mse_latent_smooth,mae_latent_smooth=mae_latent_smooth))
  },mc.cores = n_cores
  )
  if(save_data){
    return(list(sim_data = simdata, output = res))
  }else{
    return(res)
  }
}





spike.f <- function (t)
  0.75 * exp(-500   * (t - 0.23)^2) +
  1.5  * exp(-2000  * (t - 0.33)^2) +
  3    * exp(-8000  * (t - 0.47)^2) +
  2.25 * exp(-16000 * (t - 0.69)^2) +
  0.5  * exp(-32000 * (t - 0.83)^2)


angles.f <- function (t) {
  s <- ((2 * t + 0.5) * (t <= 0.15)) +
    ((-12 * (t - 0.15) + 0.8) * (t > 0.15 & t <= 0.2)) +
    0.2 * (t > 0.2 & t <= 0.5) +
    ((6 * (t - 0.5) + 0.2) * (t > 0.5 & t <= 0.6)) +
    ((-10 * (t - 0.6) + 0.8) * (t > 0.6 & t <= 0.65)) +
    ((-0.5 * (t - 0.65) + 0.3) * (t > 0.65 & t <= 0.85)) +
    ((2 * (t - 0.85) + 0.2) * (t > 0.85))
  f <- 3/5 * ((5/(max(s) - min(s))) * s - 1.6) - 0.0419569
  return(f)
}


cblocks.f <- function (t) {
  pos <- c(0.1,0.13,0.15,0.23,0.25,0.4,0.44,0.65,0.76,0.78,0.81)
  hgt <- 2.88/5 * c(4,-5,3,-4,5,-4.2,2.1,4.3,-3.1,2.1,-4.2)
  f   <- rep(0,length(t))
  for (i in 1:length(pos))
    f <- f + (1 + sign(t - pos[i])) * (hgt[i]/2)
  f[f < 0] <- 0
  return(f)
}

sblocks.f <- function (t) {

  n = length(t)
  return(c(rep(0,n/8),rep(1.5,n/4),rep(0,n/4),rep(3,n/4),rep(0,n/8)))
}


# This defines the "Heavisine" signal.
heavi.f <- function (t) {
  heavi <- 4 * sin(4 * pi * t) - sign(t - 0.3) - sign(0.72 - t)
  f <- heavi/sqrt(var(heavi)) * 1 * 2.99/3.366185
  f <- f - min(f)
  return(f)
}

# This defines the "Bursts" signal.
bursts.f <- function (t) {
  I1 <- exp(-(abs(t - 0.2)/0.01)^1.2) * (t <= 0.2) +
    exp(-(abs(t - 0.2)/0.03)^1.2) * (t > 0.2)
  I2 <- exp(-(abs(t - 0.3)/0.01)^1.2) * (t <= 0.3) +
    exp(-(abs(t - 0.3)/0.03)^1.2) * (t > 0.3)
  I3 <- exp(-(abs(t - 0.4)/0.01)^1.2) * (t <= 0.4) +
    exp(-(abs(t - 0.4)/0.03)^1.2) * (t > 0.4)
  f  <- 2.99/4.51804 * (4*I1 + 3*I2 + 4.5*I3)
  return(f)
}


