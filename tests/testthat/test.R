set.seed(12345)
n=2^9+100
sigma=0.5
mu=c(rep(0.3,n/4), rep(3, n/4), rep(10, n/4), rep(0.3, n/4))
x = rpois(length(mu),exp(log(mu)+rnorm(length(mu),sd=sigma)))

res_vst_nug = try(smash_gen_pois(x,transformation='vst',smoother='smash',F))
plot(x,col='grey80')
lines(res_vst_nug$posterior$mean_smooth)

res_lik_exp_logx_nug = try(smash_gen_pois(x,transformation='lik_expan',smoother='smash',
                                          est_nugget_maxiter=30,
                                          lik_expan_at = 'logx',
                                          est_nugget = TRUE))
plot(x,col='grey80')
lines(res_lik_exp_logx_nug$posterior$mean_smooth)

split_smashpoi_dwt = try(ebps(x,
                              init_control = list(m_init_method='smash_poi'),
                              smooth_control = list(wave_trans='dwt'),
                              general_control = list(maxiter=100)))

plot(x,col='grey80')
lines(mu,col='grey50')
lines(split_smashpoi_dwt$posterior$mean_smooth)
split_smashpoi_dwt$elbo

split_vga_dwt = try(ebps(x,
                         init_control = list(m_init_method='vga'),
                         smooth_control = list(wave_trans='dwt'),
                         general_control = list(maxiter=100)))
plot(x,col='grey80',pch=20)
lines(mu,col='grey50')
lines(split_vga_dwt$posterior$mean_smooth)
split_vga_dwt$log_likelihood


fit = ebps(x,general_control = list(verbose=T,tol=1e-5,maxiter=1),
           g_init = list(sigma2=fit$fitted_g$sigma2),
           q_init = list(smooth=fit$posterior$mean_log_smooth),
                       smooth_control = list(wave_trans = 'dwt'))
plot(x,col='grey80')
lines(mu,col='grey50')
lines(fit$posterior$mean_smooth)
fit$elbo

split_smashpoi_ndwt = try(ebps(x,
                               init_control = list(m_init_method='smash_poi'),
                               smooth_control = list(wave_trans='ndwt',
                                                     ndwt_method = 'smash'),
                               general_control = list(maxiter=100)))
plot(x,col='grey80')
lines(mu,col='grey50')
lines(split_smashpoi_ndwt$posterior$mean_smooth)

split_vga_ndwt = try(ebps(x,
                          init_control = list(m_init_method='vga'),
                          smooth_control = list(wave_trans='ndwt',
                                                ndwt_method = 'smash'),
                          general_control = list(maxiter=100,tol=0.01)))
plot(x,col='grey80')
lines(mu,col='grey50')
lines(split_vga_ndwt$posterior$mean_smooth)

res_smash2 = try(smash_two_step(x,homoskedastic=TRUE))
plot(x,col='grey80')
lines(mu,col='grey50')
lines(res_smash2$posterior$mean_smooth)

res_smash2 = try(smash_two_step(x,homoskedastic=F))
plot(x,col='grey80')
lines(mu,col='grey50')
lines(res_smash2$posterior$mean_smooth)

fit = smash_gen_pois(x,est_nugget_maxiter = 3,smoother = 'smash',
                     transformation='lik_expan',
                     lik_expan_at = 'smash_poi',est_nugget = T)
plot(x,col='grey80')
lines(mu,col='grey50')
lines(fit$posterior$mean_smooth)

fit = smash_gen_pois_iterative(x,maxiter = 10)


simdata = sim_data_smooth(1,count_size=10,n=128)
out = simu_study_poisson_smooth(simdata,n_cores=1)

set.seed(12345)
n=2^12
count_size = 10
sigma=0.5
t = seq(0,1,length.out = n)
b = spike.f(t)
b = b/(max(b)/count_size)
x = rpois(n,exp(log(b)+rnorm(n,sd=sigma)))
plot(x,col='grey80')
lines(b,col='grey50')

fit = pois_smooth_split(x,maxiter=100,wave_trans = 'dwt',verbose = T,Emu_init = 'vga',maxiter_vga = 1,ndwt_method='ti.thresh')
lines(fit$posterior$mean_smooth)


##############
x = TPM3[1,]

fit = ebps(x,
           init_control = list(),
           general_control = list(verbose=T,tol=1e-2,maxiter=30),
           smooth_control = list(wave_trans = 'ndwt',
                                 ndwt_method = 'ti.thresh',
                                 robust=T))
plot(x,col='grey80',pch='.',cex=2)
lines(fit$posterior$mean_smooth)
plot(fit$fitted_g$sigma2_trace)





