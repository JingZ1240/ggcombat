# Formal operating-range simulation for AOAS revision.
#
# Reviewer 1 Major Comment 1 asks whether GG-ComBat's advantage is limited to
# moderate ABCD-like site effects or persists as additive, variance, and
# covariance-level site heterogeneity become stronger and as block coherence
# weakens. This script starts from the submitted simulation design and varies
# those axes systematically.

root <- normalizePath(file.path(getwd(), "../../../.."))
source(file.path(root, "study/analysis/00_common/simulation_common.R"))
source(file.path(root, "study/analysis/00_common/ml_evaluation.R"))
source(file.path(root, "R/evaluation_metrics.R"))

load_study_packages(c(
  "dplyr",
  "tibble",
  "purrr",
  "tidyr",
  "readr",
  "MASS",
  "MCMCpack",
  "neuroCombat",
  "CovBat",
  "parsnip",
  "recipes",
  "workflows",
  "ranger",
  "pROC",
  "yardstick"
))
source_legacy_method(root)

design <- make_main_simulation_design()

seed_global <- as.integer(Sys.getenv("GGCOMBAT_OP_SEED", "20260909"))
n_reps <- as.integer(Sys.getenv("GGCOMBAT_OP_N_REPS", "100"))
run_prediction <- as.logical(as.integer(Sys.getenv("GGCOMBAT_OP_RUN_PREDICTION", "0")))
n_splits_site <- as.integer(Sys.getenv("GGCOMBAT_OP_SITE_SPLITS", "5"))
trees_site <- as.integer(Sys.getenv("GGCOMBAT_OP_TREES_SITE", "100"))
trees_loso <- as.integer(Sys.getenv("GGCOMBAT_OP_TREES_LOSO", "200"))

additive_scale_grid <- as.numeric(strsplit(Sys.getenv("GGCOMBAT_OP_ADDITIVE", "0.5,1,2"), ",")[[1]])
variance_scale_grid <- as.numeric(strsplit(Sys.getenv("GGCOMBAT_OP_VARIANCE", "0.5,1,2"), ",")[[1]])
covariance_scale_grid <- as.numeric(strsplit(Sys.getenv("GGCOMBAT_OP_COVARIANCE", "0.5,1,2"), ",")[[1]])
coherence_grid <- as.numeric(strsplit(Sys.getenv("GGCOMBAT_OP_COHERENCE", "1,0.75,0.5,0.25"), ",")[[1]])

n <- design$n
K_site <- design$K_site
G_vec <- design$G_vec
p <- sum(G_vec)
methods <- c("Raw", "ComBat", "CovBat", "GG-ComBat")
# The submitted SGT simulator uses the final block as singleton/no-latent
# features. The legacy fitting code expects that singleton block to be present
# as the final column of L and then zeros it internally when estimating latent
# covariance.
L_true <- make_L_from_block_sizes(G_vec, n_groups = length(G_vec))

condition_grid <- expand.grid(
  additive_scale = additive_scale_grid,
  variance_scale = variance_scale_grid,
  covariance_scale = covariance_scale_grid,
  block_coherence = coherence_grid,
  KEEP.OUT.ATTRS = FALSE
)

scale_deviations_by_block <- function(site_vectors, scale) {
  mat <- do.call(rbind, site_vectors)
  center <- colMeans(mat)
  out <- sweep(mat, 2, center, "-") * scale
  out <- sweep(out, 2, center, "+")
  lapply(seq_len(nrow(out)), function(i) as.numeric(out[i, ]))
}

scale_variance_thetas <- function(lambdas, thetas, scale) {
  lambda_mat <- do.call(rbind, lambdas)
  theta_mat <- do.call(rbind, thetas)
  mean_mat <- theta_mat / (lambda_mat - 1)
  center <- colMeans(mean_mat)
  scaled_mean <- sweep(sweep(mean_mat, 2, center, "-") * scale, 2, center, "+")
  scaled_mean <- pmax(scaled_mean, 0.02)
  scaled_theta <- scaled_mean * (lambda_mat - 1)
  lapply(seq_len(nrow(scaled_theta)), function(i) as.numeric(scaled_theta[i, ]))
}

scale_covariance_deviations <- function(sigma_list, scale) {
  avg <- Reduce(`+`, sigma_list) / length(sigma_list)
  lapply(sigma_list, function(S) {
    S_scaled <- avg + scale * (S - avg)
    Spd.proj(SigA.hat = S_scaled, eps = 0.001)$mat
  })
}

mix_block_coherence <- function(values, block_sizes, rho) {
  out <- values
  start <- 1
  for (block in seq_along(block_sizes)) {
    end <- start + block_sizes[block] - 1
    idx <- start:end
    block_mean <- mean(values[idx])
    if (rho >= 1) {
      out[idx] <- values[idx]
    } else {
      centered <- values[idx] - block_mean
      shuffled <- sample(values[idx]) - block_mean
      out[idx] <- block_mean + sqrt(rho) * centered + sqrt(1 - rho) * shuffled
    }
    start <- end + 1
  }
  out
}

compute_pairwise_distances <- function(data_list, batch) {
  split_idx <- split(seq_len(length(batch)), batch)
  dplyr::bind_rows(lapply(names(data_list), function(method) {
    Zm <- lapply(split_idx, function(idx) data_list[[method]][idx, , drop = FALSE])
    pair_df <- evaluate_site_distances(Zm)
    pair_df$Method <- method
    pair_df
  }))
}

long_distance_metrics <- function(pair_df) {
  metric_cols <- c(
    mean_distance = "Mean",
    variance_distance = "Var",
    frobenius_covariance_distance = "Frobenius",
    spectral_covariance_distance = "Spectral"
  )
  dplyr::bind_rows(lapply(names(metric_cols), function(col) {
    tibble::tibble(
      site_a = pair_df$site_a,
      site_b = pair_df$site_b,
      Pair = paste(pair_df$site_a, pair_df$site_b, sep = "-"),
      Method = pair_df$Method,
      Metric = unname(metric_cols[col]),
      Value = pair_df[[col]]
    )
  }))
}

summarize_metric <- function(df, group_cols) {
  df |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
    dplyr::summarise(
      median = stats::median(Value, na.rm = TRUE),
      mean = mean(Value, na.rm = TRUE),
      sd = stats::sd(Value, na.rm = TRUE),
      q25 = stats::quantile(Value, 0.25, na.rm = TRUE),
      q75 = stats::quantile(Value, 0.75, na.rm = TRUE),
      .groups = "drop"
    )
}

pairwise_method_tests <- function(df) {
  df |>
    dplyr::select(
      replicate, condition_id, additive_scale, variance_scale,
      covariance_scale, block_coherence, Pair, Method, Metric, Value
    ) |>
    tidyr::pivot_wider(names_from = Method, values_from = Value) |>
    dplyr::group_by(condition_id, additive_scale, variance_scale, covariance_scale, block_coherence, Metric) |>
    dplyr::summarise(
      p_GG_vs_ComBat = stats::wilcox.test(`GG-ComBat`, ComBat, paired = TRUE, exact = FALSE)$p.value,
      p_GG_vs_CovBat = stats::wilcox.test(`GG-ComBat`, CovBat, paired = TRUE, exact = FALSE)$p.value,
      p_GG_vs_Raw = stats::wilcox.test(`GG-ComBat`, Raw, paired = TRUE, exact = FALSE)$p.value,
      .groups = "drop"
    )
}

set.seed(seed_global)
distance_results <- tibble::tibble()
site_pred_results <- tibble::tibble()
sex_pred_results <- tibble::tibble()

for (row_id in seq_len(nrow(condition_grid))) {
  cond <- condition_grid[row_id, ]
  cond_id <- sprintf(
    "a%s_v%s_c%s_r%s",
    cond$additive_scale,
    cond$variance_scale,
    cond$covariance_scale,
    cond$block_coherence
  )

  batch_means_scaled <- scale_deviations_by_block(design$batch_means, cond$additive_scale)
  thetas_scaled <- scale_variance_thetas(design$lambdas, design$thetas, cond$variance_scale)

  for (r in seq_len(n_reps)) {
    rep_seed <- seed_global + row_id * 100000 + r
    set.seed(rep_seed)

    beta_sex <- rep(0, p)
    beta_sex[1:50] <- stats::runif(50, 0.2, 0.4)

    sex_list <- lapply(seq_len(K_site), function(k) {
      probs <- stats::runif(1, 0.2, 0.8)
      factor(stats::rbinom(n, 1, probs), labels = c("F", "M"))
    })

    sigma_f_base <- gen_sigma_f_list(rep_seed = 2026 + r)
    sigma_f_list <- scale_covariance_deviations(sigma_f_base, cond$covariance_scale)

    Z_list <- gamma_truth <- delta2_truth <- vector("list", K_site)
    for (k in seq_len(K_site)) {
      out <- simulate_multiblock_sgt(
        n = n,
        G_vec = G_vec,
        Sigma_f = sigma_f_list[[k]],
        batch_means = batch_means_scaled[[k]],
        batch_sds = design$batch_sds[[k]],
        lambdas = design$lambdas[[k]],
        thetas = thetas_scaled[[k]]
      )

      gamma_coh <- mix_block_coherence(out$gamma_true, G_vec, cond$block_coherence)
      delta2_coh <- exp(mix_block_coherence(log(out$delta2_true), G_vec, cond$block_coherence))
      centered_noise <- out$Z -
        matrix(out$gamma_true, nrow = n, ncol = p, byrow = TRUE) -
        out$F_true[, rep(seq_along(G_vec), G_vec), drop = FALSE]
      rescaled_noise <- sweep(centered_noise, 2, sqrt(delta2_coh / out$delta2_true), "*")
      Zk <- matrix(gamma_coh, nrow = n, ncol = p, byrow = TRUE) +
        out$F_true[, rep(seq_along(G_vec), G_vec), drop = FALSE] +
        rescaled_noise

      sex_num <- as.numeric(sex_list[[k]]) - 1
      Z_list[[k]] <- Zk + (sex_num %o% beta_sex)
      gamma_truth[[k]] <- gamma_coh
      delta2_truth[[k]] <- delta2_coh
    }

    combined <- do.call(rbind, Z_list)
    batch <- factor(rep(seq_len(K_site), each = n), labels = paste0("Site", seq_len(K_site)))
    sex <- factor(unlist(lapply(sex_list, as.character)), levels = c("F", "M"))
    mod <- stats::model.matrix(~ sex)

    cb <- neuroCombat::neuroCombat(
      dat = t(combined),
      batch = batch,
      mod = mod,
      verbose = FALSE
    )
    comb_std <- t(cb$dat.standardized)
    delta2_g <- cb$estimates$var.pooled
    delta_g <- sqrt(delta2_g)
    alpha_g <- cb$estimates$stand.mean[, 1]
    beta_hat <- cb$estimates$beta.hat["sexM", ]
    Xbeta <- mod[, 2] %o% beta_hat

    out_cov <- CovBat::covbat(dat = t(combined), bat = batch, mod = mod, n.pc = 10)
    Z_list_std <- split_by_site(comb_std, n_per_site = n, n_site = K_site)

    hyper_list <- lapply(Z_list_std, estimate_hyperparams_sgt, L = L_true, sgt = TRUE)
    out_list <- lapply(
      Z_list_std,
      run_combat_iterative,
      L = L_true,
      tol = 1e-5,
      max_iter = 30,
      verbose = FALSE,
      sgt = TRUE
    )

    gamma_hats <- lapply(out_list, `[[`, "gamma")
    delta2_hats <- lapply(out_list, `[[`, "delta2")
    Sigma_fhats <- lapply(hyper_list, `[[`, "Sigma_fhat")

    Z_gg_list <- harmonize_sites_sgt_2(
      Z_list = Z_list_std,
      gamma_hats = gamma_hats,
      delta2_hats = delta2_hats,
      Sigma_fhats = Sigma_fhats,
      L = L_true,
      reference = "average"
    )
    Z_gg_std <- do.call(rbind, Z_gg_list)
    Y_gg <- sweep(sweep(Z_gg_std, 2, delta_g, "*"), 2, alpha_g, "+") + Xbeta

    data_list <- list(
      Raw = combined,
      ComBat = t(cb$dat.combat),
      CovBat = t(out_cov$dat.covbat),
      `GG-ComBat` = Y_gg
    )

    std_data_list <- list(
      Raw = sweep(sweep(combined - Xbeta, 2, alpha_g, "-"), 2, delta_g, "/"),
      ComBat = sweep(sweep(t(cb$dat.combat) - Xbeta, 2, alpha_g, "-"), 2, delta_g, "/"),
      CovBat = sweep(sweep(t(out_cov$dat.covbat) - Xbeta, 2, alpha_g, "-"), 2, delta_g, "/"),
      `GG-ComBat` = Z_gg_std
    )

    dist_long <- long_distance_metrics(compute_pairwise_distances(std_data_list, batch))
    dist_long$replicate <- r
    dist_long$condition_id <- cond_id
    dist_long$additive_scale <- cond$additive_scale
    dist_long$variance_scale <- cond$variance_scale
    dist_long$covariance_scale <- cond$covariance_scale
    dist_long$block_coherence <- cond$block_coherence
    distance_results <- dplyr::bind_rows(distance_results, dist_long)

    if (run_prediction) {
      site_pred_results <- dplyr::bind_rows(
        site_pred_results,
        evaluate_methods_site_multiclass(
          data_list = data_list,
          batch = batch,
          p = p,
          r = r,
          n_splits = n_splits_site,
          trees = trees_site,
          mtry = floor(sqrt(p)),
          base_seed = rep_seed + 5000
        ) |>
          dplyr::mutate(
            condition_id = cond_id,
            additive_scale = cond$additive_scale,
            variance_scale = cond$variance_scale,
            covariance_scale = cond$covariance_scale,
            block_coherence = cond$block_coherence
          )
      )
      sex_pred_results <- dplyr::bind_rows(
        sex_pred_results,
        evaluate_methods_loso_sex(
          data_list = data_list,
          batch = batch,
          sex = sex,
          p = p,
          r = r,
          trees = trees_loso,
          mtry = floor(sqrt(p))
        ) |>
          dplyr::mutate(
            condition_id = cond_id,
            additive_scale = cond$additive_scale,
            variance_scale = cond$variance_scale,
            covariance_scale = cond$covariance_scale,
            block_coherence = cond$block_coherence
          )
      )
    }

    message(sprintf(
      "Finished condition %s (%d/%d), replicate %d/%d",
      cond_id, row_id, nrow(condition_grid), r, n_reps
    ))
  }
}

out_dir <- file.path(root, "study/analysis/04_revision_robustness/01_operating_range/results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

distance_summary <- summarize_metric(
  distance_results,
  c("condition_id", "additive_scale", "variance_scale", "covariance_scale",
    "block_coherence", "Method", "Metric")
)
distance_tests <- pairwise_method_tests(distance_results)

readr::write_csv(distance_results, file.path(out_dir, "operating_range_distances.csv"))
readr::write_csv(distance_summary, file.path(out_dir, "operating_range_distance_summary.csv"))
readr::write_csv(distance_tests, file.path(out_dir, "operating_range_distance_pairwise_tests.csv"))

if (run_prediction) {
  readr::write_csv(site_pred_results, file.path(out_dir, "operating_range_site_prediction.csv"))
  readr::write_csv(sex_pred_results, file.path(out_dir, "operating_range_sex_prediction_loso.csv"))
}

meta <- tibble::tibble(
  seed_global = seed_global,
  n_reps = n_reps,
  run_prediction = run_prediction,
  n_per_site = n,
  K_site = K_site,
  p_features = p,
  G_vec = paste(G_vec, collapse = ","),
  additive_scale = paste(additive_scale_grid, collapse = ","),
  variance_scale = paste(variance_scale_grid, collapse = ","),
  covariance_scale = paste(covariance_scale_grid, collapse = ","),
  block_coherence = paste(coherence_grid, collapse = ","),
  n_conditions = nrow(condition_grid),
  output_dir = out_dir
)
readr::write_csv(meta, file.path(out_dir, "operating_range_meta.csv"))

print(distance_summary)
