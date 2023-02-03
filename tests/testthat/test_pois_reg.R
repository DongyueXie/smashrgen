set.seed(12345)
n = 100
p = 100
X = matrix(runif(n*p),nrow=n,ncol=p)
b = c(rep(2,5),rep(-2,5),rep(0,p-10))
lambda = exp(X%*%b)
y = rpois(n,lambda)

# fit glmnet
library(glmnet)
fit_glmnet = glmnet(X,y,family='poisson',intercept = F)
fit_cv = cv.glmnet(X,y,family='poisson',intercept = F)
plot(fit_cv)
fit_cv$lambda.min
plot(b,coef.glmnet(fit_glmnet,s=fit_cv$lambda.1se)[-1])

# fit
fit_split = pois_reg_splitting(X,y,printevery = 1)
plot(b,fit_split$fit_mrash$b)

# Simulate a data set.
library(mr.ash)
set.seed(1)
n <- 200
p <- 300

data <- simulate_regression_data(n, p, s = 250)
fit.mr.ash <- fit_mr_ash(data$X,data$y)

# X = matrix(runif(n*p),nrow=n,ncol=p)
# b = c(rep(1,5),rep(-1,5),rep(0,p-10))
# y = X%*%b + rnorm(n)
# ### fit Mr.ASH
# fit.mr.ash <- fit_mr_ash(X,y,intercept = F)

plot(fit.mr.ash$b,rowSums(fit.mr.ash$m*fit.mr.ash$phi),type='p',col='grey80')
abline(a=0,b=1,lty=2)

