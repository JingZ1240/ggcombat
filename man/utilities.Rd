\name{ggcombat-utilities}
\alias{ggcombat_standardize}
\alias{restore_standardized}
\alias{spd_project}
\alias{evaluate_site_distances}
\title{GG-ComBat utility functions}
\usage{
ggcombat_standardize(Y, covariates = NULL, intercept = TRUE)
restore_standardized(Z, standardization)
spd_project(sigma, eps = NULL)
evaluate_site_distances(Z_list)
}
\arguments{
\item{Y}{Numeric subjects by features matrix.}
\item{covariates}{Optional covariate data frame.}
\item{intercept}{Whether to include an intercept in standardization.}
\item{Z}{Standardized data matrix.}
\item{standardization}{Object returned by \code{ggcombat_standardize}.}
\item{sigma}{Covariance matrix.}
\item{eps}{Minimum eigenvalue for positive-definite projection.}
\item{Z_list}{List of site-specific matrices.}
}
\value{Utility-specific return values.}
\description{
Small helpers for standardization, positive-definite projection, and evaluation.
}
