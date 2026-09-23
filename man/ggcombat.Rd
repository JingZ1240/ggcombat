\name{ggcombat}
\alias{ggcombat}
\title{Fit and apply Graph-Guided ComBat harmonization}
\usage{
ggcombat(Y, batch, covariates = NULL, L, target = c("average", "reference"),
  reference = NULL, singleton_group = ncol(L),
  variance_calibration = c("empirical_target", "model_target"),
  singleton_calibration = c("standardize", "target"), tol = 1e-05,
  max_iter = 30, verbose = FALSE)
}
\arguments{
\item{Y}{Numeric subjects by features matrix.}
\item{batch}{Site or scanner labels, one per row of \code{Y}.}
\item{covariates}{Optional data frame of biological covariates to preserve during standardization.}
\item{L}{Feature by community membership matrix. The last column is treated as singleton features by default.}
\item{target}{Use a sample-size weighted average target or a reference-site target.}
\item{reference}{Reference site index or label when \code{target = "reference"}.}
\item{singleton_group}{Column of \code{L} containing singleton or weakly connected features.}
\item{variance_calibration}{Whether to use the direct model-based target
marginal standard deviations or the historical empirical variance calibration.}
\item{singleton_calibration}{Whether singleton features use submitted
ComBat-style standardization or target residual variance scaling.}
\item{tol}{Convergence tolerance for empirical-Bayes updates.}
\item{max_iter}{Maximum empirical-Bayes iterations.}
\item{verbose}{Print convergence messages.}
}
\value{A list containing harmonized standardized data, harmonized data on the input scale, fitted site parameters, target covariance, and metadata.}
\description{
High-level wrapper for the experimental GG-ComBat implementation. By default,
harmonization preserves the submitted design: correlated features use empirical
marginal calibration after correlation alignment, and singleton features use
univariate ComBat-style standardization.
}
