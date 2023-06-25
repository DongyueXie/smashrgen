

#'@title Define the squared exponential kernel
kernel <- function(grid, lengthscale) {
  matrix <- outer(grid, grid, function(x, y) exp(-(x-y)^2/2/lengthscale^2))
  return(matrix)
}

#'@title Define the negative marginal log-likelihood
#'@importFrom mvtnorm dmvnorm
neg_log_likelihood <- function(params, x, s2, grid) {
  sigma2 <- params[1]
  theta <- params[2]
  l = params[3]
  mu <- rep(theta, length(x))
  cov <- sigma2 * kernel(grid, l) + diag(s2, length(x))
  -dmvnorm(x, mean = mu, sigma = cov, log = TRUE)
}

#'@title Define the gradient of the negative log-likelihood
neg_log_likelihood_gradient <- function(params, x, s2, grid) {
  sigma2 <- params[1]
  theta <- params[2]
  l = params[3]
  n <- length(x)

  mu <- rep(theta, n)
  K <- kernel(grid, l)
  cov <- sigma2 * K + diag(s2, n)
  inv_cov <- solve(cov)

  # Compute the gradient with respect to sigma2
  d_sigma2 <- 0.5 * sum(diag(inv_cov %*% K)) - 0.5 * t(x - mu) %*% inv_cov %*% K %*% inv_cov %*% (x - mu)

  # Compute the gradient with respect to theta
  d_theta <- -sum(inv_cov %*% (x - mu))

  #
  #d_l <- sigma2 * sum(sapply(1:n, function(i) sum(sapply(1:n, function(j) ((grid[i] - grid[j])^2 / l^3) * K[i, j] * (x[i] - mu[i]) * (x[j] - mu[j])))))

  # Compute the gradient for l
  # d_l <- 0
  # for (i in 1:(n - 1)) {
  #   for (j in (i + 1):n) {
  #     d_kern_ij <- sigma2 * ((grid[i] - grid[j])^2 / l^3) * exp(-(grid[i] - grid[j])^2 / (2 * l^2))
  #     d_l <- d_l + (1 / 2) * ((x[i] - theta) * (x[j] - theta) * K_inv[i, j]^2 - K_inv[i, j]) * d_kern_ij
  #   }
  # }
  # d_l <- d_l * 2  # account for symmetry
  c(d_sigma2, d_theta,d_l)
}

#'@title Define the empirical Bayes function
GP <- function(x, s2,init_val=NULL) {
  grid <- seq(1, length(x)) / length(x)

  # Optimization
  if(is.null(init_val)){
    start <- c(1, mean(x),1) # Initialize sigma2 and theta
  }else{
    start <- init_val # Initialize sigma2 and theta
  }

  opt_res <- optim(start, neg_log_likelihood, x = x, s2 = s2, grid = grid, method = "BFGS")

  # Extract the estimates
  sigma2_hat <- opt_res$par[1]
  theta_hat <- opt_res$par[2]
  l_hat = opt_res$par[3]

  # Compute posterior mean and covariance
  K <- kernel(grid, l_hat)
  cov_prior <- sigma2_hat * K
  cov_likelihood <- diag(s2, length(x))
  #post_cov <- solve(solve(cov_prior) + solve(cov_likelihood))
  post_cov = cov_prior%*%solve(cov_prior+cov_likelihood)%*%cov_likelihood
  #post_mean <- post_cov %*% (solve(cov_prior) %*% rep(theta_hat, length(x)) + solve(cov_likelihood) %*% x)
  post_mean = c(cov_likelihood%*%solve(cov_prior+cov_likelihood)%*% rep(theta_hat, length(x)) + cov_prior%*%solve(cov_prior+cov_likelihood)%*%x)

  list(sigma2_hat = sigma2_hat, theta_hat = theta_hat, l_hat=l_hat,post_mean = post_mean, post_cov = post_cov)
}




# # Test the function
# set.seed(123)
# n <- 100
# mu <- rep(2, n)
# sigma <- 1
# x <- rnorm(n, mu, sigma)
# s2 <- rep(sigma^2, n)
#
# # Optimization with gradient
# start <- c(1, mean(x))  # Initialize sigma2 and theta
# grid <- seq(1, length(x)) / length(x)
# opt_res <- optim(start, neg_log_likelihood, gr = neg_log_likelihood_gradient, x = x, s2 = s2, grid = grid, method = "BFGS")
#
# print(opt_res$par)

# set.seed(123)
# n <- 200
# mu <- sin(2 * pi * (1:n)/n)  # Mu follows a sinusoidal curve
# sigma <- 1
# x <- rnorm(n, mu, sigma)
# s2 <- rep(sigma^2, n)

# # Test the function
# set.seed(123)
# n <- 100
# mu <- rep(2, n)
# sigma <- 1
# x <- rnorm(n, mu, sigma)
# s2 <- rep(sigma^2, n)

# result <- GP(x, s2,init_val = c(1,1,1))
# print(result$sigma2_hat)
# print(result$theta_hat)
# print(result$l_hat)
#
# plot(x)
# lines(mu)
# lines(result$post_mean,col=4)



#'@title Fit Poisson + Gaussian process
#'@export
Pois_GP = function(x,maxiter=100,tol=1e-8,printevery=1,verbose=TRUE){
  fit_init = ebpm_normal(x, NULL, g_init = list(mean = log(mean(x)),var = NULL), fix_g = c(T, F))
  m = fit_init$posterior$mean_log
  sigma2 = fit_init$fitted_g$var
  sigma2_trace = sigma2
  obj_trace = -Inf
  init_val_gp = NULL
  for(iter in 1:maxiter){
    fit.gp = GP(m,sigma2,init_val=init_val_gp)
    init_val_gp = c(fit.gp$sigma2_hat, fit.gp$theta_hat, fit.gp$l_hat)
    opt = vga_pois_solver(m, x, rep(1,length(x)), fit.gp$post_mean,
                          sigma2, tol = 1e-10, maxiter = 1000)
    m = opt$m
    v = opt$v
    sigma2 = mean(m^2+v+fit.gp$post_mean^2+diag(fit.gp$post_cov)-2*m*fit.gp$post_mean)
    sigma2_trace = c(sigma2_trace,sigma2)
    obj_trace[iter+1] = Pois_GP_obj(x,m,v,sigma2,fit.gp$post_mean,fit.gp$post_cov)
    if(abs(obj_trace[iter+1] - obj_trace[iter]) < tol){
      break
    }
    if(verbose){
      if(iter%%printevery==0){
        print(paste("Done iter", iter, "obj =", obj_trace[iter + 1]))
      }
    }
  }
  return(list(obj_trace=obj_trace,
              fitted_g = list(sigma2=sigma2,sigma2_trace=sigma2_trace,gp_sigma2 = fit.gp$sigma2_hat,gp_lengthscale = fit.gp$l_hat,gp_mean = fit.gp$theta_hat),
              posterior=list(mean_log=m,
                             var_log=v,
                             mean_b = fit.gp$post_mean,
                             var_b=fit.gp$post_cov)))
}

Pois_GP_obj = function(y,m,v,sigma2,b,Vb){
  sum(y*m-exp(m + v/2)) - n/2*log(sigma2) - 1/2/sigma2*sum(m^2+v+b^2+diag(Vb)-2*m*b) + sum(log(v))/2 + as.numeric(determinant(Vb,log=T)$modulus)/2
}


