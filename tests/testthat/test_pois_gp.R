set.seed(1234)
n = 100
mu = sin(2 * pi * (1:n)/n)
mu = mu - min(mu)
mu = mu * 3+ 0.1
x = rpois(length(mu),exp(log(mu)+rnorm(n,0,0.5)))
fit = Pois_GP(x,maxiter=30)

plot(fit$elbo_trace,type='l')
plot(fit$fitted_g$sigma2_trace,type='l',ylim=c(-0.1,max(fit$fitted_g$sigma2_trace)))
fit$fitted_g$sigma2

n <- length(x)  # assuming x is your data vector
df <- data.frame(
  x = (1:n) / n,
  y = x,
  mu = mu,  # assuming mu is your mean vector
  y_fit = exp(fit$posterior$mean_b + diag(fit$posterior$var_b) / 2)  # fit data
)

ggplot(df, aes(x = x)) +
  geom_point(aes(y = y), color = 'grey80') +
  geom_line(aes(y = mu), color = 'grey50') +
  geom_line(aes(y = y_fit), color = 'blue', size = 1)

#lines((1:n)/n,smashr::smash.poiss(x),col=2)

plot(fit$posterior$mean_log,col='grey80')
lines(fit$posterior$mean_b)

plot(fit$posterior$mean_log,fit$posterior$mean_b,col='grey80')
abline(a = 0,b=1)


fit0 = Pois_GP(x,maxiter=1)

# ii = 100
# xl = 3
# xh = 6
#
# plot_densities <- function(mu0, sigma0, mu1, sigma1, x_low, x_high) {
#   # Create a sequence of x values
#   x <- seq(from = x_low, to = x_high, by = 0.01)
#
#   # Compute density values for the two normal distributions
#   density0 <- dnorm(x, mean = mu0, sd = sigma0)
#   density1 <- dnorm(x, mean = mu1, sd = sigma1)
#
#   # Plot the densities
#   plot(x, density0, type = "l", col = "blue",
#        xlab = "x", ylab = "Density",
#        main = "Normal Densities", ylim = c(0, max(density0, density1)))
#   lines(x, density1, type = "l", col = "red")
#
#   # Add a legend
#   legend("topright", legend = c("Density 0", "Density 1"),
#          lty = c(1, 1), col = c("blue", "red"))
# }
#
# plot_densities(mu0 = fit$posterior$mean_log[ii], sigma0 = sqrt(fit$posterior$var_log[ii]),
#                mu1 = fit$posterior$mean_b[ii], sigma1 = sqrt(fit$posterior$var_b[ii,ii]),
#                x_low = xl, x_high = xh)
# #abline(v = log(mu)[ii],lty=2)
#
# plot_densities(mu0 = fit0$posterior$mean_log[ii], sigma0 = sqrt(fit0$posterior$var_log[ii]),
#                mu1 = fit0$posterior$mean_b[ii], sigma1 = sqrt(fit0$posterior$var_b[ii,ii]),
#                x_low = xl, x_high = xh)
# #abline(v = log(mu)[ii],lty=2)


library(ggplot2)


plot_conf_bands <- function(mean1, sd1, mean2, sd2, true_mean, ylim_range) {
  n <- length(mean1)
  df <- data.frame(
    x = 1:n,
    lower1 = mean1 - 2*sd1,
    upper1 = mean1 + 2*sd1,
    mean1 = mean1,
    lower2 = mean2 - 2*sd2,
    upper2 = mean2 + 2*sd2,
    mean2 = mean2,
    true_mean = true_mean
  )

  ggplot(df, aes(x = x)) +
    geom_ribbon(aes(ymin = lower1, ymax = upper1), fill = "red", alpha = 0.2) +
    geom_line(aes(y = mean1), color = "red") +
    geom_ribbon(aes(ymin = lower2, ymax = upper2), fill = "blue", alpha = 0.2) +
    geom_line(aes(y = mean2), color = "blue") +
    geom_line(aes(y = true_mean), color = "grey80", lty = 2) +
    coord_cartesian(ylim = ylim_range) +
    #theme_minimal() +
    theme(axis.title = element_blank())
}


plot_conf_bands(fit0$posterior$mean_log,sqrt(fit0$posterior$var_log),
                fit0$posterior$mean_b,sqrt(diag(fit0$posterior$var_b)),
                log(mu),
                c(-2,7))
plot_conf_bands(fit$posterior$mean_log,sqrt(fit$posterior$var_log),
                fit$posterior$mean_b,sqrt(diag(fit$posterior$var_b)),
                log(mu),
                c(-2,7))



