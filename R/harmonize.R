model_covariance <- function(L, sigma_fhat, delta2) {
  symmetrize_matrix(L %*% sigma_fhat %*% t(L) + diag(delta2, nrow = length(delta2)))
}

harmonize_ggcombat_sites <- function(Z_list, fits = NULL, gamma_hats = NULL,
                                     delta2_hats = NULL, Sigma_fhats = NULL,
                                     L, target = c("average", "reference"),
                                     reference = 1,
                                     singleton_group = ncol(L),
                                     variance_calibration = c("empirical_target", "model_target"),
                                     singleton_calibration = c("standardize", "target")) {
  target <- match.arg(target)
  variance_calibration <- match.arg(variance_calibration)
  singleton_calibration <- match.arg(singleton_calibration)

  n_sites <- length(Z_list)
  if (!is.null(fits)) {
    gamma_hats <- lapply(fits, `[[`, "gamma")
    delta2_hats <- lapply(fits, `[[`, "delta2")
    Sigma_fhats <- lapply(fits, `[[`, "Sigma_fhat")
  }

  stopifnot(
    length(gamma_hats) == n_sites,
    length(delta2_hats) == n_sites,
    length(Sigma_fhats) == n_sites
  )

  n_features <- ncol(Z_list[[1]])
  stopifnot(nrow(L) == n_features)
  singleton_idx <- if (is.null(singleton_group)) integer(0) else which(L[, singleton_group] != 0)
  correlated_idx <- setdiff(seq_len(n_features), singleton_idx)
  if (length(correlated_idx) == 0) {
    stop("At least one non-singleton feature is required for graph-guided covariance alignment.")
  }

  sigma_list <- vector("list", n_sites)
  for (site in seq_len(n_sites)) {
    sigma_full <- model_covariance(L, Sigma_fhats[[site]], delta2_hats[[site]])
    sigma_corr <- sigma_full[correlated_idx, correlated_idx, drop = FALSE]
    sigma_list[[site]] <- spd_project(sigma_corr)$mat
  }

  if (target == "average") {
    weights <- vapply(Z_list, nrow, numeric(1))
    weights <- weights / sum(weights)
    sigma_target <- Reduce(`+`, Map(function(sigma, weight) weight * sigma, sigma_list, weights))
    sigma_target <- spd_project(sigma_target)$mat
    delta2_target <- Reduce(`+`, Map(function(delta2, weight) weight * delta2, delta2_hats, weights))
  } else {
    sigma_target <- sigma_list[[reference]]
    delta2_target <- delta2_hats[[reference]]
  }

  site_sd <- lapply(sigma_list, function(sigma) sqrt(pmax(diag(sigma), 1e-12)))
  site_cor <- Map(function(sigma, sd) spd_project(sigma / (sd %o% sd))$mat, sigma_list, site_sd)
  target_sd <- sqrt(pmax(diag(sigma_target), 1e-12))
  target_cor <- spd_project(sigma_target / (target_sd %o% target_sd))$mat
  target_cor_sqrt <- matrix_sqrt(target_cor)$sqrt
  site_cor_inv_sqrt <- lapply(site_cor, function(cor_mat) matrix_sqrt(cor_mat)$inv_sqrt)

  z_corrected <- vector("list", n_sites)
  cur_sd_after_correlation <- vector("list", n_sites)
  for (site in seq_len(n_sites)) {
    z_site <- Z_list[[site]]
    centered <- sweep(z_site, 2, gamma_hats[[site]], "-")
    z_new <- matrix(0, nrow = nrow(z_site), ncol = n_features)

    z_block <- centered[, correlated_idx, drop = FALSE]
    z_std <- sweep(z_block, 2, site_sd[[site]], "/")
    z_cor <- z_std %*% site_cor_inv_sqrt[[site]] %*% target_cor_sqrt
    cur_sd_after_correlation[[site]] <- sqrt(pmax(apply(z_cor, 2, stats::var), 1e-12))
    if (variance_calibration == "model_target") {
      z_new[, correlated_idx] <- sweep(z_cor, 2, target_sd, "*")
    } else {
      scale_vec <- target_sd / cur_sd_after_correlation[[site]]
      z_new[, correlated_idx] <- sweep(z_cor, 2, scale_vec, "*")
    }

    if (length(singleton_idx) > 0) {
      if (singleton_calibration == "standardize") {
        singleton_scale <- sqrt(pmax(delta2_hats[[site]][singleton_idx], 1e-12))
        z_new[, singleton_idx] <- sweep(centered[, singleton_idx, drop = FALSE], 2, singleton_scale, "/")
      } else {
        singleton_scale <- sqrt(pmax(delta2_target[singleton_idx], 1e-12) /
                                  pmax(delta2_hats[[site]][singleton_idx], 1e-12))
        z_new[, singleton_idx] <- sweep(centered[, singleton_idx, drop = FALSE], 2, singleton_scale, "*")
      }
    }

    z_corrected[[site]] <- z_new
  }

  attr(z_corrected, "target_covariance") <- sigma_target
  attr(z_corrected, "target_correlation") <- target_cor
  attr(z_corrected, "target_delta2") <- delta2_target
  attr(z_corrected, "correlated_idx") <- correlated_idx
  attr(z_corrected, "singleton_idx") <- singleton_idx
  attr(z_corrected, "cur_sd_after_correlation_alignment") <- cur_sd_after_correlation
  attr(z_corrected, "variance_calibration") <- variance_calibration
  attr(z_corrected, "singleton_calibration") <- singleton_calibration
  z_corrected
}
