harmonize_sites_correlation_target <- function(Z_list, gamma_hats, delta2_hats,
                                               Sigma_fhats, L,
                                               reference = "average",
                                               variance_calibration = c("empirical_target", "model_target"),
                                               singleton_calibration = c("standardize", "target")) {
  variance_calibration <- match.arg(variance_calibration)
  singleton_calibration <- match.arg(singleton_calibration)
  n_sites <- length(Z_list)
  stopifnot(
    length(gamma_hats) == n_sites,
    length(delta2_hats) == n_sites,
    length(Sigma_fhats) == n_sites
  )

  n_features <- ncol(Z_list[[1]])
  n_singleton <- sum(L[, ncol(L)] != 0)
  correlated_idx <- seq_len(n_features - n_singleton)
  singleton_idx <- if (n_singleton > 0) (n_features - n_singleton + 1):n_features else integer(0)

  eig_sqrt <- function(S, eps = 1e-10) {
    S <- (S + t(S)) / 2
    ee <- eigen(S, symmetric = TRUE)
    val <- pmax(ee$values, eps)
    Q <- ee$vectors
    list(
      S12 = Q %*% diag(sqrt(val), nrow = length(val)) %*% t(Q),
      Sm12 = Q %*% diag(1 / sqrt(val), nrow = length(val)) %*% t(Q)
    )
  }

  sigma_list <- vector("list", n_sites)
  for (site in seq_len(n_sites)) {
    sigma_full <- L %*% Sigma_fhats[[site]] %*% t(L) + diag(delta2_hats[[site]])
    sigma_corr <- sigma_full[correlated_idx, correlated_idx, drop = FALSE]
    sigma_list[[site]] <- Spd.proj(SigA.hat = sigma_corr)$mat
  }

  if (identical(reference, "average")) {
    weights <- vapply(Z_list, nrow, numeric(1))
    weights <- weights / sum(weights)
    sigma_target <- Reduce(`+`, Map(function(S, w) w * S, sigma_list, weights))
    sigma_target <- Spd.proj(SigA.hat = sigma_target)$mat
    delta2_target <- Reduce(`+`, Map(function(delta2, w) w * delta2, delta2_hats, weights))
  } else {
    sigma_target <- sigma_list[[reference]]
    delta2_target <- delta2_hats[[reference]]
  }

  site_sd <- lapply(sigma_list, function(S) sqrt(pmax(diag(S), 1e-12)))
  site_cor <- Map(function(S, sd) Spd.proj(SigA.hat = S / (sd %o% sd))$mat, sigma_list, site_sd)
  target_sd <- sqrt(pmax(diag(sigma_target), 1e-12))
  target_cor <- Spd.proj(SigA.hat = sigma_target / (target_sd %o% target_sd))$mat
  target_cor_sqrt <- eig_sqrt(target_cor)$S12
  site_cor_inv_sqrt <- lapply(site_cor, function(R) eig_sqrt(R)$Sm12)

  Z_corrected <- vector("list", n_sites)
  for (site in seq_len(n_sites)) {
    Z_centered <- sweep(Z_list[[site]], 2, gamma_hats[[site]], "-")
    Z_new <- matrix(0, nrow = nrow(Z_list[[site]]), ncol = n_features)

    Z_std <- sweep(Z_centered[, correlated_idx, drop = FALSE], 2, site_sd[[site]], "/")
    Z_cor <- Z_std %*% site_cor_inv_sqrt[[site]] %*% target_cor_sqrt
    if (variance_calibration == "empirical_target") {
      cur_sd <- sqrt(pmax(apply(Z_cor, 2, stats::var), 1e-12))
      Z_new[, correlated_idx] <- sweep(Z_cor, 2, target_sd / cur_sd, "*")
    } else {
      Z_new[, correlated_idx] <- sweep(Z_cor, 2, target_sd, "*")
    }

    if (length(singleton_idx) > 0) {
      if (singleton_calibration == "standardize") {
        Z_new[, singleton_idx] <- sweep(
          Z_centered[, singleton_idx, drop = FALSE],
          2,
          sqrt(pmax(delta2_hats[[site]][singleton_idx], 1e-12)),
          "/"
        )
      } else {
        singleton_scale <- sqrt(pmax(delta2_target[singleton_idx], 1e-12) /
                                  pmax(delta2_hats[[site]][singleton_idx], 1e-12))
        Z_new[, singleton_idx] <- sweep(
          Z_centered[, singleton_idx, drop = FALSE],
          2,
          singleton_scale,
          "*"
        )
      }
    }

    Z_corrected[[site]] <- Z_new
  }

  Z_corrected
}

harmonize_sites_sgt_2 <- function(Z_list, gamma_hats, delta2_hats, Sigma_fhats,
                                  L, reference = "average") {
  harmonize_sites_correlation_target(
    Z_list = Z_list,
    gamma_hats = gamma_hats,
    delta2_hats = delta2_hats,
    Sigma_fhats = Sigma_fhats,
    L = L,
    reference = reference,
    variance_calibration = "empirical_target",
    singleton_calibration = "standardize"
  )
}

harmonize_sites_sgt_3 <- function(Z_list, gamma_hats, delta2_hats, Sigma_fhats,
                                  L, reference = 1) {
  harmonize_sites_correlation_target(
    Z_list = Z_list,
    gamma_hats = gamma_hats,
    delta2_hats = delta2_hats,
    Sigma_fhats = Sigma_fhats,
    L = L,
    reference = reference,
    variance_calibration = "empirical_target",
    singleton_calibration = "target"
  )
}
