load_study_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(
      "Missing required packages: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(lapply(packages, library, character.only = TRUE))
}

source_legacy_method <- function(root = ".") {
  source(file.path(root, "study/provenance/source_snapshots/str_combat_source_20250527.R"))
  source(file.path(root, "study/analysis/00_common/legacy_harmonize_sites_sgt_3.R"))
}

make_L_from_block_sizes <- function(G_vec, n_groups = length(G_vec)) {
  p <- sum(G_vec)
  L <- matrix(0, nrow = p, ncol = n_groups)
  start <- 1
  for (b in seq_len(min(length(G_vec), n_groups))) {
    end <- start + G_vec[b] - 1
    L[start:end, b] <- 1
    start <- end + 1
  }
  L
}

make_main_simulation_design <- function() {
  list(
    methods = c("Raw", "ComBat", "CovBat", "GGCombat"),
    n = 300,
    K_site = 3,
    G_vec = c(50, 50, 50, 50),
    batch_means = list(
      c(1, 3, 1, 2),
      c(2, 3, 3, 1),
      c(3, 2, 3, 2)
    ),
    batch_sds = list(
      rep(0.5, 4),
      rep(0.5, 4),
      rep(0.5, 4)
    ),
    lambdas = list(
      c(3, 5, 3, 3),
      c(5, 3, 5, 5),
      c(5, 3, 3, 5)
    ),
    thetas = list(
      c(0.4, 0.5, 0.3, 0.2),
      c(0.3, 0.5, 0.5, 0.5),
      c(0.4, 0.5, 0.3, 0.3)
    )
  )
}

gen_sigma_f_list <- function(rep_seed,
                             min.block = 0.8,
                             max.block = 1.1,
                             min.bg = 0.2,
                             max.bg = 0.6) {
  lapply(seq_len(3), function(j) {
    set.seed(rep_seed + j)
    diag_vals <- runif(3, min.block, max.block)
    rho <- runif(3, min.bg, max.bg)

    Sigma_f <- matrix(c(
      diag_vals[1], -rho[1], rho[2],
      -rho[1], diag_vals[2], rho[3] / 2,
      rho[2], rho[3] / 2, diag_vals[3]
    ), nrow = 3, byrow = TRUE)

    Spd.proj(SigA.hat = Sigma_f, eps = 0.001)$mat
  })
}

compute_distances <- function(Z1, Z2) {
  S1 <- cov(Z1)
  S2 <- cov(Z2)

  c(
    Frobenius = sqrt(sum((S1 - S2)^2)),
    Spectral = max(abs(eigen(S1 - S2, symmetric = TRUE)$values)),
    Var = sqrt(sum((diag(S1) - diag(S2))^2)),
    Mean = norm(colMeans(Z1) - colMeans(Z2), type = "2")
  )
}

compute_correlation_distance <- function(Z1, Z2) {
  R1 <- cor(Z1)
  R2 <- cor(Z2)
  sqrt(sum((R1 - R2)^2))
}

compute_total_avg_corr_dist <- function(Z, L) {
  R <- cor(Z)
  K <- ncol(L)
  total <- 0
  for (b in seq_len(K)) {
    idx <- which(L[, b] == 1)
    vals <- R[idx, idx][upper.tri(R[idx, idx])]
    total <- total + mean(1 - vals)
  }
  total
}

split_by_site <- function(x, n_per_site, n_site) {
  lapply(seq_len(n_site), function(j) {
    idx1 <- (j - 1) * n_per_site + 1
    idx2 <- j * n_per_site
    x[idx1:idx2, , drop = FALSE]
  })
}
