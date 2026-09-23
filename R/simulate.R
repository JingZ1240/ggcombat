rinvgamma_base <- function(n, shape, scale) {
  1 / stats::rgamma(n, shape = shape, rate = scale)
}

simulate_block <- function(n, latent, n_features, batch_mean, batch_sd,
                           lambda, theta) {
  gamma <- stats::rnorm(n_features, mean = batch_mean, sd = batch_sd)
  delta2 <- rinvgamma_base(n_features, shape = lambda, scale = theta)
  noise <- vapply(delta2, function(v) stats::rnorm(n, sd = sqrt(v)), numeric(n))
  latent_mat <- matrix(latent, nrow = n, ncol = n_features)
  batch_mat <- matrix(gamma, nrow = n, ncol = n_features, byrow = TRUE)

  list(
    data = batch_mat + latent_mat + noise,
    gamma = gamma,
    delta2 = delta2
  )
}

simulate_ggcombat_data <- function(n, block_sizes, latent_cov,
                                   batch_means, batch_sds,
                                   lambdas, thetas,
                                   singleton_last = TRUE) {
  n_blocks <- length(block_sizes)
  n_latent <- if (singleton_last) n_blocks - 1 else n_blocks

  stopifnot(
    is.matrix(latent_cov),
    nrow(latent_cov) == n_latent,
    ncol(latent_cov) == n_latent,
    length(batch_means) == n_blocks,
    length(batch_sds) == n_blocks,
    length(lambdas) == n_blocks,
    length(thetas) == n_blocks
  )

  latent <- matrix(stats::rnorm(n * n_latent), nrow = n)
  chol_cov <- chol(spd_project(latent_cov)$mat)
  latent <- latent %*% chol_cov
  if (singleton_last) {
    latent <- cbind(latent, rep(0, n))
  }

  blocks <- vector("list", n_blocks)
  for (k in seq_len(n_blocks)) {
    blocks[[k]] <- simulate_block(
      n = n,
      latent = latent[, k],
      n_features = block_sizes[k],
      batch_mean = batch_means[k],
      batch_sd = batch_sds[k],
      lambda = lambdas[k],
      theta = thetas[k]
    )
  }

  y <- do.call(cbind, lapply(blocks, `[[`, "data"))
  gamma <- unlist(lapply(blocks, `[[`, "gamma"), use.names = FALSE)
  delta2 <- unlist(lapply(blocks, `[[`, "delta2"), use.names = FALSE)

  l_mat <- matrix(0, nrow = sum(block_sizes), ncol = n_blocks)
  offset <- 0
  for (k in seq_len(n_blocks)) {
    idx <- seq.int(offset + 1, offset + block_sizes[k])
    l_mat[idx, k] <- 1
    offset <- offset + block_sizes[k]
  }

  list(
    Z = y,
    L = l_mat,
    gamma_true = gamma,
    delta2_true = delta2,
    F_true = latent
  )
}
