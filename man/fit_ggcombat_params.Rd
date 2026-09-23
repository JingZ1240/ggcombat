\name{fit_ggcombat_params}
\alias{fit_ggcombat_params}
\alias{estimate_ggcombat_hyperparams}
\alias{update_ggcombat_params}
\title{Estimate GG-ComBat empirical-Bayes parameters}
\usage{
fit_ggcombat_params(Z_list, L, singleton_group = ncol(L), tol = 1e-05,
  max_iter = 30, verbose = FALSE)
}
\arguments{
\item{Z_list}{List of site-specific standardized data matrices.}
\item{L}{Feature by community membership matrix.}
\item{singleton_group}{Column of \code{L} containing singleton features.}
\item{tol}{Convergence tolerance.}
\item{max_iter}{Maximum number of iterations.}
\item{verbose}{Print convergence messages.}
}
\value{A list of fitted parameter objects, one per site.}
\description{
Estimates groupwise empirical-Bayes mean, variance, and latent covariance
parameters used by GG-ComBat.
}
