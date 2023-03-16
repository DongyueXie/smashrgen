set.seed(12345)
n=2^9
sigma=1
mu=c(rep(0.3,n/4), rep(3, n/4), rep(10, n/4), rep(0.3, n/4))
x = mu + rnorm(n,0,sigma)
x[128] = x[128] + 10
x[256] = x[256] + 10
plot(x,col='grey80',pch=20)
lines(mu,col='grey60')
fit0 = smash_dwt(x,1)
lines(fit0$posterior$mean,col=2,lwd=2)
fit1 = robust_smooth(x,1,mu_init = rep(mean(x,n)))
plot(x,col='grey80',pch=20)
lines(mu,col='grey60')
lines(fit1$smooth_res$posterior$mean,col=2,lwd=2)

fit1 = robust_smooth(x,1,mu_init = runmed(x,k=7))
plot(x,col='grey80',pch=20)
lines(mu,col='grey60')
lines(fit1$smooth_res$posterior$mean,col=2,lwd=2)
