% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sigma_est.R
\name{sigma_est}
\alias{sigma_est}
\title{A function estimate \sigma^2}
\usage{
sigma_est(x, mu = NULL, st, method = "mle", var_est = "rmad",
  k = NULL, family = "DaubExPhase", filter.number = 1)
}
\arguments{
\item{x:}{data}

\item{method:}{'eb', 'moment', 'wls', 'mle','huber','smash'}

\item{var_est:}{method to estimate variance: 'rmad', 'smash', 'default'}

\item{k:}{parameter in huber m estimator}
}
\value{
estimated sd
}
\description{
A function estimate \sigma^2
}
