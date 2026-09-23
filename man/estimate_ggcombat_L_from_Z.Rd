\name{estimate_ggcombat_L_from_Z}
\alias{estimate_ggcombat_L_from_Z}
\title{Estimate a graph-guided membership matrix from standardized data}
\usage{
estimate_ggcombat_L_from_Z(Z_list, n_blocks = 4,
  prctile_vec = seq(80, 90, by = 0.5),
  lambda_vec = seq(0.4, 0.8, length.out = 5), use_abs = TRUE,
  ncores = 1, use_parallel = FALSE)
}
\arguments{
\item{Z_list}{List of site-specific standardized subject by feature matrices.}
\item{n_blocks}{Number of leading dense ICONS groups to retain as correlated graph-guided blocks. Remaining features are assigned to the final singleton group.}
\item{prctile_vec}{Candidate percentile thresholds passed to \code{ICONS::param_tuning_sigmau()}.}
\item{lambda_vec}{Candidate lambda values passed to \code{ICONS::param_tuning_sigmau()}.}
\item{use_abs}{Whether to run ICONS on the absolute weighted correlation matrix.}
\item{ncores}{Number of cores passed to ICONS tuning.}
\item{use_parallel}{Whether ICONS tuning should use parallel computation.}
}
\value{A list containing \code{L}, dense extraction results, tuning parameters,
the weighted correlation matrix, site-wise correlations, retained block sizes,
and original-order indices of correlated and singleton features. Rows of
\code{L} are in the ICONS graph order \code{results$Clist}.}
\description{
Constructs the GG-ComBat graph scaffold from the standardized matrix
\eqn{Z_{ij}} used in parameter estimation. The function computes site-wise
feature correlations, forms the sample-size weighted pooled correlation matrix,
tunes ICONS parameters, extracts dense groups, retains the first
\code{n_blocks} groups, and assigns all remaining features to a final singleton
or weak-feature group.
}
\details{
The returned \code{L} is in graph order, not original feature order. Reorder
site-specific standardized matrices by \code{results$Clist} before passing them
to \code{fit_ggcombat_params()} or \code{harmonize_ggcombat_sites()}.
}
