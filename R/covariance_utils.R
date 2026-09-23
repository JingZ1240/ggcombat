symmetrize_matrix <- function(x) {
  (x + t(x)) / 2
}

matrix_sqrt <- function(sigma, eps = 1e-10) {
  sigma <- symmetrize_matrix(sigma)
  eig <- eigen(sigma, symmetric = TRUE)
  values <- pmax(eig$values, eps)
  vectors <- eig$vectors

  list(
    sqrt = vectors %*% diag(sqrt(values), nrow = length(values)) %*% t(vectors),
    inv_sqrt = vectors %*% diag(1 / sqrt(values), nrow = length(values)) %*% t(vectors)
  )
}

spd_project <- function(sigma, eps = NULL) {
  p <- ncol(sigma)
  if (is.null(eps)) {
    eps <- 5 / p
  }

  sigma <- symmetrize_matrix(sigma)
  eig <- eigen(sigma, symmetric = TRUE)
  if (min(eig$values) > eps) {
    return(list(mat = sigma, conv = 1))
  }

  vectors <- eig$vectors
  values <- pmax(eig$values, eps)
  projected <- vectors %*% diag(values, nrow = length(values)) %*% t(vectors)
  list(mat = symmetrize_matrix(projected), conv = 2)
}

l1_project <- function(v, radius) {
  stopifnot(radius > 0)
  if (sum(abs(v)) <= radius) {
    return(v)
  }

  u <- sort(abs(v), decreasing = TRUE)
  sv <- cumsum(u)
  rho <- max(which(u > (sv - radius) / seq_along(u)))
  theta <- max(0, (sv[rho] - radius) / rho)
  sign(v) * pmax(abs(v) - theta, 0)
}

admm_project_spd <- function(mat, epsilon = 1e-4, mu = 10, it_max = 1e3,
                             etol = 1e-4, etol_distance = 1e-4) {
  p <- nrow(mat)
  mat <- symmetrize_matrix(mat)
  r_mat <- diag(diag(mat), p)
  s_mat <- matrix(0, p, p)
  lambda <- matrix(0, p, p)

  diagnostics <- vector("list", it_max)
  for (iter in seq_len(it_max)) {
    r_prev <- r_mat
    s_prev <- s_mat

    w_mat <- mat + s_mat + mu * lambda
    eig <- eigen(symmetrize_matrix(w_mat), symmetric = TRUE)
    r_mat <- eig$vectors %*% diag(pmax(eig$values, epsilon), p) %*% t(eig$vectors)

    m_mat <- r_mat - mat - mu * lambda
    lower <- lower.tri(m_mat, diag = TRUE)
    s_mat[lower] <- m_mat[lower] - l1_project(m_mat[lower], radius = mu / 2)
    s_mat[upper.tri(s_mat)] <- t(s_mat)[upper.tri(s_mat)]

    lambda <- lambda - (r_mat - s_mat - mat) / mu

    eps_r <- max(abs(r_mat - r_prev))
    eps_s <- max(abs(s_mat - s_prev))
    eps_primal <- max(abs(r_mat - s_mat - mat))
    distance <- max(abs(r_mat - mat))
    diagnostics[[iter]] <- data.frame(
      iteration = iter,
      eps_R = eps_r,
      eps_S = eps_s,
      eps_primal = eps_primal,
      distance = distance
    )

    if ((eps_r < etol && eps_s < etol && eps_primal < etol) ||
        abs(max(abs(r_prev - mat)) - distance) < etol_distance) {
      diagnostics <- diagnostics[seq_len(iter)]
      break
    }

    if (iter %% 20 == 0) {
      mu <- mu / 2
    }
  }

  list(mat = symmetrize_matrix(r_mat), df_ADMM = do.call(rbind, diagnostics))
}
