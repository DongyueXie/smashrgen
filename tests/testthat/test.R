set.seed(12345)
n=2^9
sigma=0.5
mu=c(rep(0.3,n/4), rep(3, n/4), rep(10, n/4), rep(0.3, n/4))
x = rpois(n,exp(log(mu)+rnorm(n,sd=sigma)))
fit = pois_smooth_split(x,maxiter=30,wave_trans = 'dwt')
plot(x,col='grey80')
lines(mu,col='grey50')
lines(fit$posterior$mean_smooth)

fit = pois_smooth_split(x,maxiter=10,wave_trans = 'ndwt',ndwt_method = 'smash')
plot(x,col='grey80')
lines(mu,col='grey50')
lines(fit$posterior$mean_smooth)

fit = pois_smooth_split(x,maxiter=30,wave_trans = 'ndwt',ndwt_method = 'ti.thresh')
plot(x,col='grey80')
lines(mu,col='grey50')
lines(fit$posterior$mean_smooth)

fit = smash_gen_pois(x,smoother = 'ti.thresh',transformation='vst',est_nugget = F)
plot(x,col='grey80')
lines(mu,col='grey50')
lines(fit$posterior$mean_smooth)

fit = smash_gen_pois(x,est_nugget_maxiter = 3,smoother = 'smash',
                     transformation='lik_expan',
                     lik_expan_at = 'smash_poi',est_nugget = T)
plot(x,col='grey80')
lines(mu,col='grey50')
lines(fit$posterior$mean_smooth)

fit = smash_gen_pois_iterative(x,maxiter = 10)


simdata = sim_data_smooth(2,count_size=10,n=128)
out = simu_study_poisson_smooth(simdata,n_cores=1)


