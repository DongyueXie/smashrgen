% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/smash_gen_pois.R
\name{smash_gen_pois}
\alias{smash_gen_pois}
\title{Smooth Poisson sequence, accounting for nugget effect}
\usage{
smash_gen_pois(
  x,
  s = 1,
  nug.init = NULL,
  est_nugget = TRUE,
  transformation = "lik_expan",
  lik_expan_at = "logx",
  nug.est.limit = 1,
  smoother = "smash",
  robust = FALSE,
  robust.q = 0.99,
  ash_pm_init_for0 = TRUE,
  eps = "estimate",
  filter.number = 1,
  family = "DaubExPhase",
  homoskedastic = FALSE,
  est_nugget_maxiter = 2,
  est_nugget_tol = 0.01
)
}
\arguments{
\item{x}{observed Poisson sequence}

\item{s}{Scale factor for Poisson observations: y~Pois(scale*lambda), can be a vector.}

\item{nug.init}{init value of nugget effect, either a scalar or NULL}

\item{transformation}{transformation of Poisson data, either 'vst' or 'lik_expan'; 'vst' for variance stabilizing transformation; 'lik_expansion' for likelihood expansion}

\item{lik_expan_at}{if transformation='lik_expan', where to expand it? Can be logx, or smash_poi}

\item{smoother}{smoothing method for Gaussian sequence, either 'smash' or 'ti.thresh'. When n is large, ti.thresh is much faster}

\item{robust}{whether perform robust wavelet regression}

\item{robust.q}{quantile to determine outliers}

\item{eps}{If choose lik_expansion, if x=0, set x = x + eps. Either input a numerical value or 'estimate'. If estimate, eps = sum(x==1)/sum(x<=1)}

\item{filter.number, family}{wavelet basis, see wavethresh package for more details}

\item{ash.pm}{If choose lik_expansion, whether use ash posterior mean approximation if x=0. If not x = x+eps.}

\item{maxiter}{max iterations for estimating nugget effect}

\item{tol}{tolerance to stop iterations.}
}
\value{
estimated smoothed lambda, estimated nugget effect.
}
\description{
Smooth Poisson sequence, accounting for nugget effect
}
