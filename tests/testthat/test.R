set.seed(12345)
n=2^11
sigma=0
mu=c(rep(0.3,n/4), rep(3, n/4), rep(10, n/4), rep(0.3, n/4))
x = rpois(length(mu),exp(log(mu)+rnorm(length(mu),sd=sigma)))
split_vga_dwt = try(ebps(x,
                         init_control = list(m_init_method='vga'),
                         smooth_control = list(wave_trans='dwt'),
                         general_control = list(maxiter=1000,tol=1e-10,printevery=10,verbose=T)))
plot(x,col='grey80',pch=20)
lines(mu,col='grey50')
lines(split_vga_dwt$posterior$mean_smooth)
split_vga_dwt$log_likelihood
split_vga_dwt$fitted_g$sigma2_trace

plot(split_vga_dwt$posterior$mean_log_smooth,split_vga_dwt$posterior$mean_log,col='grey80')
abline(a=0,b=1)

plot(split_vga_dwt$posterior$var_log)
plot(split_vga_dwt$posterior$var_log_smooth)

split_vga_dwt0 = try(ebps(x,
                         init_control = list(m_init_method='vga'),
                         smooth_control = list(wave_trans='dwt'),
                         general_control = list(maxiter=1,tol=1e-10,printevery=10,verbose=T)))
plot(x,col='grey80',pch=20)
lines(mu,col='grey50')
lines(split_vga_dwt0$posterior$mean_smooth)
split_vga_dwt0$log_likelihood
split_vga_dwt0$fitted_g$sigma2_trace

plot(split_vga_dwt0$posterior$mean_log_smooth,split_vga_dwt0$posterior$mean_log,col='grey80')
abline(a=0,b=1)

plot(split_vga_dwt0$posterior$var_log)
plot(split_vga_dwt0$posterior$var_log_smooth)

plot_densities <- function(mu0, sigma0, mu1, sigma1, x_low, x_high) {
  # Create a sequence of x values
  x <- seq(from = x_low, to = x_high, by = 0.01)

  # Compute density values for the two normal distributions
  density0 <- dnorm(x, mean = mu0, sd = sigma0)
  density1 <- dnorm(x, mean = mu1, sd = sigma1)

  # Plot the densities
  plot(x, density0, type = "l", col = "blue",
       xlab = "x", ylab = "Density",
       main = "Normal Densities", ylim = c(0, max(density0, density1)))
  lines(x, density1, type = "l", col = "red")

  # Add a legend
  legend("topright", legend = c("Density 0", "Density 1"),
         lty = c(1, 1), col = c("blue", "red"))
}


# Test the function
ii = 1330
plot_densities(mu0 = split_vga_dwt0$posterior$mean_log_smooth[ii],
               sigma0 = sqrt(split_vga_dwt0$posterior$var_log_smooth[ii]),
               mu1 = split_vga_dwt0$posterior$mean_log[ii],
               sigma1 = sqrt(split_vga_dwt0$posterior$var_log[ii]),
               1,5)


plot_densities(mu0 = split_vga_dwt$posterior$mean_log_smooth[ii],
               sigma0 = sqrt(split_vga_dwt$posterior$var_log_smooth[ii]),
               mu1 = split_vga_dwt$posterior$mean_log[ii],
               sigma1 = sqrt(split_vga_dwt$posterior$var_log[ii]),
               1,5)




#############################
#############################
#############################

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




###########
m = smash.poiss(y,log = T)
yhat = smash.gaus(m,sigma=sd_est_diff2(m))
yhat = smash.gaus(m,v.est = T,joint = T)
plot(y,col='grey80',pch=19)
lines(exp(yhat$mu.res))
lines(fit_vst$posterior$mean_smooth)

set.seed(12345)
n = 512
mu = c(rep(-2,n/4),rep(0,n/4),rep(-2,n/2))
ss = runmed(yhat$mu.res,k=31)
mu = c(rep(-3,40),ss,rep(-3,47))
sigma2 = mean(yhat$var.res)
lambda = exp(mu+rnorm(n,0,sqrt(sigma2)))
y = rpois(n,lambda)
plot(y,col='grey80',pch=19)
lines(exp(mu))

lines(exp(smash.gaus(smash.poiss(y,log=T))),col=2)

fit_split = ebps(y, init_control = list(m_init_method='vga'),
                             smooth_control = list(wave_trans='ndwt',
                                                   ndwt_method = 'smash'),
                             general_control = list(maxiter=20,printevery=1,verbose=T))

plot(y,col='grey80',pch=19)
lines(exp(mu),col='grey60')
lines(fit_split$posterior$mean_smooth)
fit_split$elbo_trace

############ chip ###########
plot(y,col='grey80',pch=19)
y2 = extend(y)
y2 = y2$x
plot(y2,col='grey80',pch=19)
fit_split = ebps(y2, init_control = list(m_init_method='vga'),
                 smooth_control = list(wave_trans='ndwt',
                                       ndwt_method = 'smash'),
                 general_control = list(maxiter=30,printevery=1,verbose=T))
lines(fit_split$posterior$mean_smooth)
plot(fit_split$fitted_g$sigma2_trace)
plot(fit_split$posterior$mean_log_smooth)



# generate data from split results
y_fake = rpois(length(y2),exp(fit_split$posterior$mean_log_smooth+rnorm(length(y2),sd=sqrt(fit_split$fitted_g$sigma2))))
plot(y_fake)
cor(y2,y_fake)



mu_hat_two_step = (smash.gaus((smash.poiss(y2,log = T)),v.est = T,joint=T))
y_fake = rpois(length(y2),exp(mu_hat_two_step$mu.res+rnorm(length(y2),sd=sqrt(mu_hat_two_step$var.res))))
cor(y2,y_fake)
plot(y_fake,col='grey80',pch=19)
lines(exp(mu_hat_two_step$mu.res),col='grey60')
lines(exp(smash.gaus((smash.poiss(y_fake,log=T)))))
# init
init = ebpm_normal(y_fake,1,g_init = list(mean=log(mean(y)),var=NULL),fix_g = c(T,F))

plot(init$posterior$mean_log)
lines(smash.gaus(init$posterior$mean_log))
lines(smash.gaus(init$posterior$mean_log,sigma=sqrt(init$fitted_g$var)))

plot(y_fake)
lines(exp(smash.gaus(init$posterior$mean_log,sigma=sqrt(init$fitted_g$var))))
lines(exp(smash.gaus(init$posterior$mean_log)))

fit_split = ebps(y_fake, init_control = list(m_init_method='smash_poi'),
                 smooth_control = list(wave_trans='dwt',
                                       ndwt_method = 'smash'),
                 general_control = list(maxiter=300,printevery=1,verbose=T))
plot(y_fake,col='grey80',pch=19)
lines(exp(mu_hat_two_step$mu.res),col='grey60')
lines(fit_split$posterior$mean_smooth,col='grey80',pch=19)
plot(fit_split$fitted_g$sigma2_trace)
fit_split$fitted_g$sigma2_trace
