\name{simulate_ggcombat_data}
\alias{simulate_ggcombat_data}
\title{Simulate block-structured multi-feature data}
\usage{
simulate_ggcombat_data(n, block_sizes, latent_cov, batch_means, batch_sds,
  lambdas, thetas, singleton_last = TRUE)
}
\arguments{
\item{n}{Number of subjects.}
\item{block_sizes}{Number of features in each community.}
\item{latent_cov}{Covariance matrix for latent community factors.}
\item{batch_means}{Community-level additive batch-effect means.}
\item{batch_sds}{Community-level additive batch-effect standard deviations.}
\item{lambdas, thetas}{Inverse-gamma hyperparameters for residual variances.}
\item{singleton_last}{Whether the last community represents singleton features.}
}
\value{A list containing simulated data \code{Z}, membership matrix \code{L}, and true parameters.}
\description{
Generates toy data for examples, tests, and simulation design.
}
