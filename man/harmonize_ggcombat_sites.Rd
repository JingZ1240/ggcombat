\name{harmonize_ggcombat_sites}
\alias{harmonize_ggcombat_sites}
\title{Apply unified GG-ComBat covariance harmonization}
\usage{
harmonize_ggcombat_sites(Z_list, fits = NULL, gamma_hats = NULL,
  delta2_hats = NULL, Sigma_fhats = NULL, L, target = c("average",
  "reference"), reference = 1, singleton_group = ncol(L),
  variance_calibration = c("empirical_target", "model_target"),
  singleton_calibration = c("standardize", "target"))
}
\arguments{
\item{Z_list}{List of site-specific standardized data matrices.}
\item{fits}{Optional list returned by \code{fit_ggcombat_params}.}
\item{gamma_hats, delta2_hats, Sigma_fhats}{Optional low-level fitted parameters.}
\item{L}{Feature by community membership matrix.}
\item{target}{Average or reference target distribution.}
\item{reference}{Reference site index for reference-target harmonization.}
\item{singleton_group}{Column of \code{L} containing singleton features.}
\item{variance_calibration}{Whether to use the direct model-based target
marginal standard deviations or the historical empirical correction based on
the observed post-correlation-alignment standard deviations.}
\item{singleton_calibration}{Whether singleton features use the submitted
ComBat-style standardization by their own residual variance or target residual
variance scaling.}
}
\value{A list of harmonized site-specific matrices. The
\code{cur_sd_after_correlation_alignment} attribute stores the empirical
post-correlation-alignment standard deviations used to diagnose whether the
historical empirical correction is negligible.}
\description{
Applies the GG-ComBat correlation-decomposition transform after removing
estimated additive site effects. For correlated features, the default preserves
the submitted empirical marginal calibration after the row-vector mapping
\code{D_i^{-1} R_i^{-1/2} R_target^{1/2}}. Singleton features use submitted
univariate ComBat-style standardization by default.
}
