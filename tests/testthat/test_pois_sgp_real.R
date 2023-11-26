library(smashrgen)
datax = read.csv('/home/dxie/Downloads/SRSF3.Counts.csv.gz')
rownames(datax)=datax[,1]
datax[1:5,1:5]
datax = datax[,-1]
dim(datax)
datax = as.matrix(datax)
sum(datax==0)/prod(dim(datax))

m=100
fit3 = pois_sgp(datax[1,],m=m,maxiter=10)
fit3$run_time
plot(datax[1,],col='grey80',pch=20)
lines(fit3$posterior$rate,col=2,lwd=2)
fit3$fitted_g

ly = log(1+datax[1,])
temp = sgp(ly,m=100,fix_x_ind = T,
           sigma2 = smashrgen:::sd_est_diff2(ly)^2,fix_sigma2 = TRUE)
plot(ly,col='grey80',pch=20)
lines(temp$posterior$mean)
temp$fitted_g$mu
fit3_5 = pois_sgp(datax[1,],fix_X_ind = T,m=100,kernel_param = temp$fitted_g$kernel_param,mu=temp$fitted_g$mu,post_mean = temp$posterior$mean_ind,V=temp$posterior$V_ind)
plot(datax[1,],col='grey80',pch=20)
lines(fit3_5$posterior$rate,col=4,lwd=2)

fit4 = pois_sgp(datax[1,],fix_X_ind = T,m=50)
fit4$run_time
plot(datax[2,],col='grey80',pch=20)
lines(fit4$posterior$rate,col=4,lwd=2)

fit5 = pois_sgp(datax[1,],fix_X_ind = T,m=100)
fit5$run_time
plot(datax[1,],col='grey80',pch=20)
lines(fit5$posterior$rate,col=4,lwd=2)


dd <- readRDS("~/myRpackages/ebnmf/iter39_pois_sgp.rds")
y = dd$y
sc = dd$sc
X_ind = dd$X_ind
post_mean = dd$post_mean
V = dd$V
kernel_param = dd$kernel_param
mu = dd$mu
res = pois_sgp(y,sc=sc,fix_X_ind = T,X_ind=X_ind,opt_method = 'L-BFGS-B',V=V,kernel_param=kernel_param,mu=mu,post_mean=post_mean,verbose = T)
res = pois_sgp(y,sc=sc,X_ind=X_ind,fix_X_ind = T,V=diag(100),kernel_param=c(0,0),mu=0,post_mean=rep(0,100),verbose = T)
res = pois_sgp(y,sc=sc,m=100)

plot(y/sc,col='grey80',pch=20)
lines(res$posterior$rate)
res$fitted_g$kernel_param
res$fitted_g$mu

plot(res$posterior$mean)










plot(debug$y,col='grey80',pch=20)
lines(smashr::smash.poiss(debug$y))
fit1 = pois_sgp(debug$y,X_ind = debug$g_init$X_ind,post_mean = debug$q_init$mean_log_ind,V=debug$q_init$v_log_ind,
                kernel_func = debug$g_init$kernel_func,kernel_param = debug$g_init$kernel_param,mu=debug$g_init$mu,s=debug$s,fix_X_ind = T)

fit_init = sgp(log(1+debug$y/debug$s),X=NULL,X_ind=debug$g_init$X_ind,mu=NULL,fix_x_ind = T,kernel_func = debug$g_init$kernel_func)

fit_init$fitted_g$mu
fit_init$fitted_g$kernel_param
fit_init$fitted_g$sigma2

debug$g_init$mu
debug$g_init$kernel_param

fit2 = pois_sgp(debug$y,X_ind = debug$g_init$X_ind,post_mean = debug$q_init$mean_log_ind,V=debug$q_init$v_log_ind,
                kernel_func = debug$g_init$kernel_func,kernel_param = fit_init$fitted_g$kernel_param,mu=debug$g_init$mu,s=debug$s,fix_X_ind = T)

fit3 = pois_sgp(debug$y,s=debug$s,fix_X_ind = T)

plot(debug$y/debug$s,col='grey80',pch=20)
lines(fit3$posterior$rate)

fit3$fitted_g$kernel_param
fit3$fitted_g$mu

plot(debug$q_init$mean_log_ind,fit3$posterior$mean_ind)
abline(a=0,b=1)

fit4 = pois_sgp(debug$y,X_ind = debug$g_init$X_ind,post_mean = fit3$posterior$mean_ind,V=fit3$posterior$V_ind,
                kernel_func = fit3$fitted_g$kernel_func,kernel_param = debug$g_init$kernel_param,mu=fit3$fitted_g$mu,s=debug$s,fix_X_ind = T)



