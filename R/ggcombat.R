ggcombat <- function(Y, batch, covariates = NULL, L,
                     target = c("average", "reference"),
                     reference = NULL,
                     singleton_group = ncol(L),
                     variance_calibration = c("empirical_target", "model_target"),
                     singleton_calibration = c("standardize", "target"),
                     tol = 1e-5, max_iter = 30, verbose = FALSE) {
  target <- match.arg(target)
  variance_calibration <- match.arg(variance_calibration)
  singleton_calibration <- match.arg(singleton_calibration)

  Y <- as.matrix(Y)
  batch <- as.factor(batch)
  if (nrow(Y) != length(batch)) {
    stop("Y must have one row per batch label.")
  }
  if (nrow(L) != ncol(Y)) {
    stop("L must have one row per feature in Y.")
  }

  sites <- levels(batch)
  if (is.null(reference)) {
    reference <- 1
  } else if (is.character(reference)) {
    reference <- match(reference, sites)
    if (is.na(reference)) {
      stop("reference was not found among batch levels.")
    }
  }

  std <- ggcombat_standardize(Y, covariates = covariates)
  split_index <- split(seq_len(nrow(Y)), batch)
  Z_list <- lapply(split_index, function(idx) std$Z[idx, , drop = FALSE])

  fits <- fit_ggcombat_params(
    Z_list = Z_list,
    L = L,
    singleton_group = singleton_group,
    tol = tol,
    max_iter = max_iter,
    verbose = verbose
  )

  Z_harmonized_list <- harmonize_ggcombat_sites(
    Z_list = Z_list,
    fits = fits,
    L = L,
    target = target,
    reference = reference,
    singleton_group = singleton_group,
    variance_calibration = variance_calibration,
    singleton_calibration = singleton_calibration
  )

  Z_harmonized <- matrix(NA_real_, nrow = nrow(Y), ncol = ncol(Y))
  for (site in seq_along(split_index)) {
    Z_harmonized[split_index[[site]], ] <- Z_harmonized_list[[site]]
  }
  colnames(Z_harmonized) <- colnames(Y)

  list(
    Z_harmonized = Z_harmonized,
    Y_harmonized = restore_standardized(Z_harmonized, std),
    fits = fits,
    L = L,
    sites = sites,
    target_covariance = attr(Z_harmonized_list, "target_covariance"),
    target_correlation = attr(Z_harmonized_list, "target_correlation"),
    target_delta2 = attr(Z_harmonized_list, "target_delta2"),
    correlated_idx = attr(Z_harmonized_list, "correlated_idx"),
    singleton_idx = attr(Z_harmonized_list, "singleton_idx"),
    cur_sd_after_correlation_alignment = attr(Z_harmonized_list, "cur_sd_after_correlation_alignment"),
    variance_calibration = attr(Z_harmonized_list, "variance_calibration"),
    singleton_calibration = attr(Z_harmonized_list, "singleton_calibration"),
    standardization = std,
    call = match.call()
  )
}
