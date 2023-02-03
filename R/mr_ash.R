# Perform EM updates for the mr-ash model, in which the (approximate)
# E-step is implemented via coordinate-wise updates (method = "cd") or
# using optim (method = "nm" or method = "bfgs").
#
# This implementation is meant to be "instructive"---that is, I've
# tried to make the code as simple as possible, with an emphasis on
# clarity. Very little effort has been devoted to making the
# implementation efficient, or the code concise.

# s0, se are variances
mr_ash <- function (X, y, se, s0, w0, b = rep(0,ncol(X)),
                    method = c("cd", "nm", "bfgs"),
                    maxiter = 100, tol = 1e-8, update.se = TRUE,
                    update.w0 = TRUE, verbose = TRUE) {
  method <- match.arg(method)
#
#   # Center X and y.
#   X <- scale(X,scale = FALSE)
#   y <- y - mean(y)

  # These two variables are used to keep track of the algorithm's
  # progress: "elbo" stores the value of the objective (the
  # variational lower bound, or "ELBO") at each iteration; "maxd"
  # stores the largest difference in the posterior mean coefficients
  # between two successive iterations.
  elbo <- rep(as.numeric(NA),maxiter)
  maxd <- rep(as.numeric(NA),maxiter)

  # Iterate the EM updates.
  if (verbose)
    cat("iter                elbo max|b-b'|\n")
  for (i in 1:maxiter) {

    # Save the current estimates of the posterior means.
    b0 <- b

    # E STEP
    # ------
    # Update the posterior means of the regression coefficients.
    if (method == "cd") {
      out  <- mr_ash_update_cd(X,y,b,se,s0,w0)
      vars <- out$post_var
    } else {
      out  <- mr_ash_update_optim(X,y,b,se,s0,w0,method)
      vars <- out$vars
    }
    b     <- out$b
    w0.em <- out$w0.em

    # Record the algorithm's progress.
    elbo[i] <- out$elbo
    maxd[i] <- max(abs(b - b0))

    # M STEP
    # ------
    # Update the residual variance, if requested.
    if (update.se)
      se <- out$erss/n

    # Update the mixture weights, if requested.
    if (update.w0)
      w0 <- w0.em

    # Report progress, if requested.
    if (verbose)
      cat(sprintf("%4d %0.12e %0.3e\n",i,elbo[i],maxd[i]))

    # Stop if the largest change in the posterior means of the
    # regression coefficients is small.
    if (maxd[i] <= tol)
      break
  }

  # Return the updated posterior means of the regression coefficicents
  # (b); the prior variances of the mixture components (s0); the
  # updated mixture weights (w0); the EM update for the mixture
  # weights (w0.em); residual variance (se); the value of the
  # objective after the (approximate) E-step at each iteration
  # (elbo); and the maximum change in the regression coefficients at
  # each iteration (maxd).
  #
  # Note that w0 and w0.em will be the same whenever update.w0 = TRUE.
  return(list(b     = b,
              vars  = vars,
              s0    = s0,
              w0    = w0,
              w0.em = w0.em,
              se    = se,
              elbo_trace = elbo[1:i],
              elbo  = elbo[i],
              maxd  = maxd[1:i]))
}

# Perform a single pass of the coordinate-wise updates for the
# mr-ash model.
#
# This implementation is meant to be "instructive"---that is, I've
# tried to make the code as simple as possible, with an emphasis on
# clarity. Very little effort has been devoted to making the
# implementation efficient, or the code concise.
mr_ash_update_cd <- function (X, y, b, se, s0, w0) {

  # Get the number of samples (n), the number of predictors (p), and
  # the number of mixture components (k).
  n <- nrow(X)
  p <- ncol(X)
  k <- length(w0)

  # This will be the M-step update for the weights in the
  # mixture-of-normals.
  w0.em <- rep(0,k)

  # Initialize two variables which will contain the expected residual
  # sum of squares (erss) and the variational lower bound (elbo) upon
  # completion of the loop. Also, store the term in the varaiationaal
  # lower bound that does not change as the variational approximation
  # is updated during the loop (elbo_const).
  elbo_const <- -n*log(2*pi*se)/2
  elbo <- elbo_const
  erss <- 0

  # Compute the expected residuals.
  r  <- drop(y - X %*% b)
  xx <- diag(crossprod(X))

  post_var = c()
  # Repeat for each predictor.
  for (i in 1:p) {
    x <- X[,i]

    # Disregard the ith predictor in the expected residuals.
    r <- r + x*b[i]

    # Compute the least-squares estimate of the coefficients (bhat)
    # and its standard error (shat) in the siimple linear regression
    # r = x*b + e, e ~ N(0,se).
    bhat <- dot(x,r)/xx[i]
    shat <- se/xx[i]

    # Update the posterior distribution of the regression coefficients
    # for the ith predictor.
    out   <- bayes_lr_mix_ss(bhat,shat,se,s0,w0)
    b[i]  <- out$mu1
    post_var[i] <- out$s1
    w0.em <- w0.em + out$w1
    llik0 <- ldnorm(r,se)

    # Calculate the ith term in the "sum of variances" for the ELBO.
    v <- norm2(x)^2 * out$s1

    # Update the expected residuals.
    r <- r - x*b[i]

    # Compute the expected residual sum of of squares (erssi) for the
    # linear regression r = xi*bi + e, and the KL-divergence between
    # the posterior and prior distributions of the regression
    # coefficient for the ith predictor.
    erssi <- norm2(r)^2 + v
    d     <- elbo_const - (out$lbf + llik0 + erssi/(2*se))
    elbo  <- elbo - d
    erss  <- erss + v
  }

  # Complete the computation of the expected residual sum of squares
  # (erss) and the variational lower bound (elbo).
  erss <- erss + norm2(r)^2
  elbo <- elbo - erss/(2*se)

  # Output the updated posterior mean coefficients (b), the M-step
  # update for the mixture weights (w0.em), the updated variational
  # lower bound (elbo), and the updated expected residual sum of
  # squares (erss).
  return(list(b = b,post_var=post_var,w0.em = w0.em/p,elbo = elbo,erss = erss))
}

# Compute approximate mr-ash posteriors by maximizing the variational
# lower bound using optim.
#
# This implementation is meant to be "instructive"---that is, I've
# tried to make the code as simple as possible, with an emphasis on
# clarity. Very little effort has been devoted to making the
# implementation efficient, or the code concise.
mr_ash_update_optim <- function (X, y, b, se, s0, w0,
                                 method = c("nm", "bfgs")) {
  method <- match.arg(method)

  # Recover the free parameters (vars) given the posterior mean
  # estimates (b).
  vars <- b2vars(X,y,b)

  # Find the free parameters (locally) maximizing the variational
  # lower bound.
  if (method == "nm")
    out <- optim(vars,
                 function (vars) -compute_elbo_ss(vars,X,y,se,s0,w0)$elbo,
                 control = list(maxit = 200,reltol = 1e-12),
                 method = "Nelder-Mead")
  else if (method == "bfgs") {
    out <- optim(vars,
                 function (vars) -compute_elbo_ss(vars,X,y,se,s0,w0)$elbo,
                 function (vars) -compute_elbo_ss(vars,X,y,se,s0,w0)$elbograd,
                 control = list(maxit = 200,reltol = 1e-12),
                 method = "BFGS")
  }

  # Output the updated free parameters (vars) and the compute_elbo_ss
  # outputs.
  return(c(list(vars = out$par),compute_elbo_ss(out$par,X,y,se,s0,w0)))
}

# Compute the variational lower bound using the "summary statistics"
# parameterization of the variational distribution.
#
# This implementation is meant to be "instructive"---that is, I've
# tried to make the code as simple as possible, with an emphasis on
# clarity. Very little effort has been devoted to making the
# implementation efficient, or the code concise.
compute_elbo_ss <- function (vars, X, y, se, s0, w0) {

  # Get the number of rows (n) and columns (p) of X, and get the
  # number of mixture components in the mixture prior (k).
  n <- nrow(X)
  p <- ncol(X)
  k <- length(w0)

  # Compute the standard errors of the least-squares estimates.
  xx   <- diag(crossprod(X))
  shat <- se/xx

  # For each variable i, compute the posterior mean coefficient (b),
  # the posterior variance (s), the log-Bayes factor (lbf) and the
  # posterior mixture assignment probabilities (w1) for
  # BSR-mix-sum-stat(bhat,shat,se,w0,s0), where bhat = vars[i].
  b     <- rep(0,p)
  s     <- rep(0,p)
  lbf   <- rep(0,p)
  w0.em <- rep(0,k)
  for (i in 1:p) {
    out    <- bayes_lr_mix_ss(vars[i],shat[i],se,s0,w0)
    b[i]   <- out$mu1
    s[i]   <- out$s1
    lbf[i] <- out$lbf
    w0.em  <- w0.em + out$w1
  }

  # Compute the variational lower bound (the ELBO).
  erss <- norm2(y - X %*% b)^2
  elbo <- -n*log(2*pi*se)/2 - erss/(2*se) +
    sum(lbf + ((vars - b)^2 - vars^2)/(2*shat))

  # Compute the gradient of the ELBO.
  yhat     <- drop(X %*% b)
  elbograd <- s/shat * (drop((y - yhat) %*% X)/se + (b - vars)/shat)

  # Output the updated posterior mean coefficients (b), the M-step
  # update for the mixture weights (w0.em), the updated variational
  # lower bound (elbo), and the updated expected residual sum of
  # squares (erss).
  return(list(b = b,w0.em = w0.em/p,erss = erss,
              elbo = elbo,elbograd = elbograd))
}

# Recover the free parameters (vars) for the variational distribution
# given posterior mean estimates (b) maximizing the ELBO.
b2vars <- function (X, y, b) {
  xx   <- diag(crossprod(X))
  yhat <- drop(X %*% b)
  return(b + drop((y - yhat) %*% X)/xx)
}

# Fit a univariate linear regression model in which the regression
# coefficient is assigned a normal prior with mean zero and variance
# s0. Here, bhat = sum(x*y)/sum(x*x) and shat = se/sum(x*x).
#
# This implementation is meant to be "instructive"---that is, I've
# tried to make the code as simple as possible, with an emphasis on
# clarity. Very little effort has been devoted to making the
# implementation efficient, or the code concise.
bayes_lr_ridge_ss <- function (bhat, shat, se, s0) {

  # Compute the posterior mean (mu1) and variance (s1) assuming a
  # normal prior with zero mean and variance s0.
  s1  <- s0/(1 + s0/shat)
  mu1 <- s1/shat * bhat

  # Compute the log-Bayes factor.
  lbf <- log(shat/(s0 + shat))/2 + mu1^2/(2*s1)

  # Return the posterior mean (mu1) and variance (s1), the
  # log-likelihood for the "null" model (when b = 0), and the
  # log-Bayes factor (lbf).
  return(list(mu1 = mu1,s1 = s1,lbf = lbf))
}

# Fit a univariate linear regression model in which the regression
# coefficient is assigned a mixture-of-normals prior. Here, bhat =
# sum(x*y)/sum(x*x) and shat = se/sum(x*x).
#
# This implementation is meant to be "instructive"---that is, I've
# tried to make the code as simple as possible, with an emphasis on
# clarity. Very little effort has been devoted to making the
# implementation efficient, or the code concise.
bayes_lr_mix_ss <- function (bhat, shat, se, s0, w0) {

  # Get the number of mixture components.
  k <- length(w0)

  # Compute the quantities separately for each mixture component.
  out <- vector("list",k)
  for (i in 1:k)
    out[[i]] <- bayes_lr_ridge_ss(bhat,shat,se,s0[i])

  # Compute the posterior assignment probabilities for the latent
  # indicator variable.
  lbf <- sapply(out,"[[","lbf")
  w1  <- softmax(lbf + log(w0))

  # Compute the log-Bayes factor as a linear combination of the
  # individual Bayes factors for each mixture component.
  z   <- log(w0) + lbf
  u   <- max(z)
  lbf <- u + log(sum(exp(z - u)))

  # Compute the posterior mean (mu1) and variance (s1) of the
  # regression coefficients.
  s1  <- 0
  mu1 <- 0
  for (i in 1:k) {
    w   <- w1[i]
    mu  <- out[[i]]$mu1
    s   <- out[[i]]$s1
    mu1 <- mu1 + w*mu
    s1  <- s1 + w*(s + mu^2)
  }
  s1 <- s1 - mu1^2

  # Return the posterior mean (mu1) and variance (s1), the
  # posterior assignment probabilities (w1) and the log-Bayes factor
  # (lbf).
  return(list(mu1 = mu1,s1 = drop(s1),w1 = w1,lbf = lbf))
}
