make_L_from_dense_results <- function(results, n_blocks) {
  if (n_blocks < 1) {
    stop("n_blocks must be at least 1.", call. = FALSE)
  }
  if (n_blocks >= length(results$CID)) {
    stop("n_blocks must be smaller than the number of detected dense groups.", call. = FALSE)
  }

  cid <- results$CID[seq_len(n_blocks)]
  n_features <- length(results$Clist)
  L <- matrix(0, nrow = n_features, ncol = n_blocks + 1)
  start <- 1
  for (block in seq_len(n_blocks + 1)) {
    end <- if (block <= n_blocks) start + cid[[block]] - 1 else n_features
    L[start:end, block] <- 1
    start <- end + 1
  }

  list(
    L = L,
    CID = cid,
    block_index = results$Clist[seq_len(sum(cid))],
    singleton_index = results$Clist[(sum(cid) + 1):n_features]
  )
}

estimate_ggcombat_L_from_Z <- function(Z_list,
                                       n_blocks = 4,
                                       prctile_vec = seq(80, 90, by = 0.5),
                                       lambda_vec = seq(0.4, 0.8, length.out = 5),
                                       use_abs = TRUE,
                                       ncores = 1,
                                       use_parallel = FALSE) {
  if (!requireNamespace("ICONS", quietly = TRUE)) {
    stop("The ICONS package is required to estimate L from Z.", call. = FALSE)
  }
  if (!is.list(Z_list) || length(Z_list) < 2) {
    stop("Z_list must be a list of at least two site-specific matrices.", call. = FALSE)
  }

  Z_list <- lapply(Z_list, as.matrix)
  n_features <- vapply(Z_list, ncol, integer(1))
  if (length(unique(n_features)) != 1) {
    stop("All matrices in Z_list must have the same number of columns.", call. = FALSE)
  }

  z_cor_data <- lapply(Z_list, stats::cor)
  site_weights <- vapply(Z_list, nrow, numeric(1))
  site_weights <- site_weights / sum(site_weights)
  sigma_comb <- Reduce(`+`, Map(function(R, w) w * R, z_cor_data, site_weights))
  graph_matrix <- if (use_abs) abs(sigma_comb) else sigma_comb
  z_all <- do.call(rbind, Z_list)

  param <- ICONS::param_tuning_sigmau(
    graph_matrix,
    z_all,
    prctile_vec = prctile_vec,
    lam_vec = lambda_vec,
    ncores = ncores,
    use_parallel = use_parallel
  )
  results <- ICONS::dense(
    W_original = graph_matrix,
    threshold = param$cut_out,
    lambda = param$lambda_out
  )
  l_fit <- make_L_from_dense_results(results, n_blocks)

  list(
    input = "Z",
    Z_cor_data = z_cor_data,
    sigma_comb = sigma_comb,
    graph_matrix = graph_matrix,
    param = param,
    results = results,
    L = l_fit$L,
    CID = l_fit$CID,
    block_index = l_fit$block_index,
    singleton_index = l_fit$singleton_index,
    n_blocks = n_blocks
  )
}

reorder_Z_list_by_graph <- function(Z_list, graph_fit) {
  lapply(Z_list, function(Z) Z[, graph_fit$results$Clist, drop = FALSE])
}
