% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ebps.R
\name{ebps}
\alias{ebps}
\title{empirical Bayes Poisson smoothing}
\usage{
ebps(
  x,
  s = NULL,
  g_init = NULL,
  q_init = NULL,
  init_control = list(),
  general_control = list(),
  smooth_control = list()
)
}
\arguments{
\item{x, s}{data vector and scaling factor. s can be a vector of the same length as x, or a scalar.}

\item{g_init}{a list of initial value of sigma2, and g_smooth. g_smooth is the initial prior g of the smoothing method. Can be NULL.}

\item{q_init}{a list of initial values of m, smooth. m is the posterior mean of mu, smooth the posterior mean of b. See the details below.}

\item{init_control}{See function ebps_init_control_default}

\item{general_control}{See function ebps_general_control_default}

\item{smooth_control}{See function ebps_smooth_control_default}
}
\description{
empirical Bayes Poisson smoothing
}
\details{
The problem is
\deqn{x_i\sim Poisson(\lambda_i,}
\deqn{\lambda_i = \exp(\mu_i)),}
\deqn{\mu_i\sim N(b_i,\sigma^2),}
\deqn{\b_i\sim g(.).}

The \code{init_control} argument is a list in which any of the following
named components will override the default algorithm settings (as
defined by \code{ebps_init_control_default}):

\describe{
\item{\code{m_init_method}}{'vga' or 'smash_poi'}
}

The \code{general_control} argument is a list in which any of the following
named components will override the default algorithm settings (as
defined by \code{ebps_general_control_default}):

\describe{
\item{\code{est_sigma2}}{whether estiamte sigma2 or fix it}
\item{\code{maxiter}}{max iteration of the main algorithm, default is 100}
\item{\code{maxiter_vga}}{max iteration of the vga step}
\item{\code{vga_tol}}{tolerance for vga step stopping}
\item{\code{verbose}}{print progress?}
\item{\code{tol}}{tolerance for stopping the main algorithm}
\item{\code{convergence_criteria}}{'objabs' or 'nugabs'}
\item{\code{make_power_of_2}}{'reflect' or 'extend'}
\item{\code{plot_updates}}{internal use only}
}

The \code{smooth_control} argument is a list in which any of the following
named components will override the default algorithm settings (as
defined by \code{ebps_smooth_control_default}):

\describe{
\item{\code{wave_trans}}{'dwt' or 'ndwt'}
\item{\code{ndwt_method}}{'smash' or 'ti.thresh'}
\item{\code{ebnm_params}}{parameters for ebnm used in wavelet smoothing}
\item{\code{warmstart}}{init posterior using last iteration's results}
\item{\code{W}}{DWT matrix for non-haar wavelet basis}
}
}
\examples{
set.seed(12345)
n=2^9
sigma=0.5
mu=c(rep(0.3,n/4), rep(3, n/4), rep(10, n/4), rep(0.3, n/4))
x = rpois(n,exp(log(mu)+rnorm(n,sd=sigma)))
fit = ebps(x)
plot(x,col='grey80')
lines(fit$posterior$mean_smooth)
fit$sigma2
plot(fit$elbo_trace)
}
