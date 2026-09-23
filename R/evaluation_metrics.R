evaluate_site_distances <- function(Z_list) {
  n_sites <- length(Z_list)
  pairs <- utils::combn(seq_len(n_sites), 2, simplify = FALSE)

  out <- lapply(pairs, function(pair) {
    a <- Z_list[[pair[1]]]
    b <- Z_list[[pair[2]]]
    cov_a <- stats::cov(a)
    cov_b <- stats::cov(b)
    eig <- eigen(cov_a - cov_b, symmetric = TRUE, only.values = TRUE)$values

    data.frame(
      site_a = pair[1],
      site_b = pair[2],
      mean_distance = sqrt(sum((colMeans(a) - colMeans(b))^2)),
      variance_distance = sqrt(sum((apply(a, 2, stats::var) - apply(b, 2, stats::var))^2)),
      frobenius_covariance_distance = sqrt(sum((cov_a - cov_b)^2)),
      spectral_covariance_distance = max(abs(eig))
    )
  })

  do.call(rbind, out)
}
