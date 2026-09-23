ggcombat_standardize <- function(Y, covariates = NULL, intercept = TRUE) {
  Y <- as.matrix(Y)
  n <- nrow(Y)

  if (is.null(covariates)) {
    design <- if (intercept) matrix(1, nrow = n, ncol = 1) else matrix(nrow = n, ncol = 0)
  } else {
    covariates <- as.data.frame(covariates)
    design <- stats::model.matrix(~ ., data = covariates)
    if (!intercept) {
      design <- design[, colnames(design) != "(Intercept)", drop = FALSE]
    }
  }

  if (ncol(design) == 0) {
    fitted <- matrix(0, nrow = n, ncol = ncol(Y))
    beta <- matrix(numeric(0), nrow = 0, ncol = ncol(Y))
  } else {
    fit <- stats::lm.fit(design, Y)
    beta <- fit$coefficients
    fitted <- design %*% beta
  }

  residual <- Y - fitted
  scale <- sqrt(pmax(colMeans(residual^2), 1e-12))
  Z <- sweep(residual, 2, scale, "/")

  list(
    Z = Z,
    fitted = fitted,
    beta = beta,
    scale = scale,
    design = design,
    feature_names = colnames(Y)
  )
}

restore_standardized <- function(Z, standardization) {
  sweep(Z, 2, standardization$scale, "*") + standardization$fitted
}
