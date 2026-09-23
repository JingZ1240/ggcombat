estimate_ggcombat_hyperparams <- function(Z, L, singleton_group = ncol(L)) {
  n <- nrow(Z)
  n_features <- ncol(Z)
  n_groups <- ncol(L)
  stopifnot(nrow(L) == n_features)

  mu_hat <- colMeans(Z)
  block_features <- lapply(seq_len(n_groups), function(k) which(L[, k] != 0))
  block_sizes <- vapply(block_features, length, integer(1))
  if (any(block_sizes == 0)) {
    stop("Each group in L must contain at least one feature.")
  }

  mu_bar <- vapply(block_features, function(idx) mean(mu_hat[idx]), numeric(1))
  tau2_bar <- vapply(block_features, function(idx) stats::var(mu_hat[idx]), numeric(1))
  tau2_bar <- pmax(tau2_bar, 1e-8)

  z_centered <- sweep(Z, 2, mu_hat, "-")
  v_tilde <- stats::cov(z_centered)

  L_factor <- L
  if (!is.null(singleton_group)) {
    L_factor[, singleton_group] <- 0
  }
  ltl_inv <- diag(1 / block_sizes, nrow = n_groups)
  f_hat <- ltl_inv %*% t(L_factor) %*% t(z_centered)
  sigma_fhat <- ltl_inv %*% t(L_factor) %*% v_tilde %*% L_factor %*% ltl_inv

  sigma_hat <- v_tilde - L_factor %*% sigma_fhat %*% t(L_factor)
  delta2_hat <- pmax(diag(sigma_hat), 1e-8)

  v_bar <- s2_bar <- numeric(n_groups)
  for (k in seq_along(block_features)) {
    idx <- block_features[[k]]
    vals <- delta2_hat[idx]
    vals_w <- pmin(vals, stats::quantile(vals, 0.95, names = FALSE))
    v_bar[k] <- mean(vals_w)
    s2_bar[k] <- max(stats::var(vals_w), 1e-8)
  }

  lambda_bar <- v_bar^2 / s2_bar + 2
  theta_bar <- v_bar^3 / s2_bar + v_bar

  if (!is.null(singleton_group)) {
    vals_w <- pmin(delta2_hat, stats::quantile(delta2_hat, 0.95, names = FALSE))
    mu_bar[singleton_group] <- mean(mu_hat)
    tau2_bar[singleton_group] <- max(stats::var(mu_hat), 1e-8)
    v_global <- mean(vals_w)
    s2_global <- max(stats::var(vals_w), 1e-8)
    lambda_bar[singleton_group] <- v_global^2 / s2_global + 2
    theta_bar[singleton_group] <- v_global^3 / s2_global + v_global
  }

  list(
    mu_bar = mu_bar,
    tau2_bar = tau2_bar,
    lambda_bar = lambda_bar,
    theta_bar = theta_bar,
    V_tilde = v_tilde,
    f_hat = f_hat,
    Sigma_fhat = sigma_fhat,
    Sigma_hat = sigma_hat,
    block_features = block_features,
    mu_hat = mu_hat,
    singleton_group = singleton_group
  )
}

update_ggcombat_params <- function(Z, L, hyper, delta_prev = NULL) {
  n <- nrow(Z)
  n_features <- ncol(Z)

  psi2 <- diag(L %*% hyper$Sigma_fhat %*% t(L))
  v2 <- if (is.null(delta_prev)) {
    apply(Z, 2, stats::var)
  } else {
    psi2 + delta_prev
  }
  v2 <- pmax(v2, 1e-8)

  gamma_star <- numeric(n_features)
  for (k in seq_along(hyper$block_features)) {
    idx <- hyper$block_features[[k]]
    gamma_star[idx] <- (
      n * hyper$tau2_bar[k] * hyper$mu_hat[idx] + v2[idx] * hyper$mu_bar[k]
    ) / (n * hyper$tau2_bar[k] + v2[idx])
  }

  residual <- Z - matrix(gamma_star, n, n_features, byrow = TRUE) -
    t(L %*% hyper$f_hat)
  ssq <- colSums(residual^2)

  delta_star <- numeric(n_features)
  for (k in seq_along(hyper$block_features)) {
    idx <- hyper$block_features[[k]]
    delta_star[idx] <- (
      hyper$theta_bar[k] + 0.5 * ssq[idx]
    ) / (n / 2 + hyper$lambda_bar[k] - 1)
  }

  list(gamma = gamma_star, delta2 = pmax(delta_star, 1e-8))
}

fit_ggcombat_site <- function(Z, L, singleton_group = ncol(L),
                              tol = 1e-5, max_iter = 30, verbose = FALSE) {
  hyper <- estimate_ggcombat_hyperparams(Z, L, singleton_group = singleton_group)
  delta_prev <- NULL
  gamma_old <- delta_old <- NULL
  converged <- FALSE

  for (iter in seq_len(max_iter)) {
    update <- update_ggcombat_params(Z, L, hyper, delta_prev)

    if (!is.null(gamma_old)) {
      d_gamma <- sum(abs(update$gamma - gamma_old))
      d_delta <- sum(abs(update$delta2 - delta_old))
      converged <- max(d_gamma, d_delta) < tol
      if (verbose) {
        message(sprintf("[Iter %d] d_gamma=%0.4g d_delta=%0.4g", iter, d_gamma, d_delta))
      }
      if (converged) {
        break
      }
    }

    gamma_old <- update$gamma
    delta_old <- update$delta2
    delta_prev <- update$delta2
  }

  list(
    gamma = update$gamma,
    delta2 = update$delta2,
    Sigma_fhat = hyper$Sigma_fhat,
    hyper = hyper,
    iterations = iter,
    converged = converged
  )
}

fit_ggcombat_params <- function(Z_list, L, singleton_group = ncol(L),
                                tol = 1e-5, max_iter = 30, verbose = FALSE) {
  lapply(
    Z_list,
    fit_ggcombat_site,
    L = L,
    singleton_group = singleton_group,
    tol = tol,
    max_iter = max_iter,
    verbose = verbose
  )
}
