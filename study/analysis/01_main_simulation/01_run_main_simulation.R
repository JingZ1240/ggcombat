# Reproducible main simulation workflow.
#
# This script is a cleaned version of the submitted `harm_ml.Rmd` analysis. It
# preserves the original estimands and writes regenerated outputs to
# `results/reproduced/` so submitted files remain unchanged.

root <- normalizePath(file.path(getwd(), "../../.."))
source(file.path(root, "study/analysis/00_common/simulation_common.R"))
source(file.path(root, "study/analysis/00_common/ml_evaluation.R"))

load_study_packages(c(
  "dplyr",
  "tibble",
  "purrr",
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

seed_global <- 927
n_reps <- as.integer(Sys.getenv("GGCOMBAT_N_REPS", "50"))
n_splits_site <- as.integer(Sys.getenv("GGCOMBAT_SITE_SPLITS", "10"))
trees_site <- as.integer(Sys.getenv("GGCOMBAT_TREES_SITE", "100"))
trees_loso <- as.integer(Sys.getenv("GGCOMBAT_TREES_LOSO", "400"))

set.seed(seed_global)

n <- design$n
K_site <- design$K_site
G_vec <- design$G_vec
p <- sum(G_vec)
methods <- design$methods

# Legacy representation: first three 50-feature communities are explicit; the
# final 50 features are singleton/weakly connected by position.
L_true <- make_L_from_block_sizes(G_vec, n_groups = 3)

gamma_long <- tibble::tibble()
delta2_long <- tibble::tibble()
distances_df <- tibble::tibble()
site_pred <- tibble::tibble()
sex_pred_loso <- tibble::tibble()

for (r in seq_len(n_reps)) {
  beta_sex <- rep(0, p)
  beta_sex[1:50] <- runif(50, 0.2, 0.4)

  sex_list <- vector("list", K_site)
  for (k in seq_len(K_site)) {
    probs <- runif(1, 0.2, 0.8)
    sex_bin <- rbinom(n, 1, probs)
    sex_list[[k]] <- factor(sex_bin, labels = c("F", "M"))
  }

  Z_list <- vector("list", K_site)
  gamma_truth <- vector("list", K_site)
  delta2_truth <- vector("list", K_site)
  Sigma_f_list <- gen_sigma_f_list(rep_seed = 2026 + r)

  for (k in seq_len(K_site)) {
    out <- simulate_multiblock_sgt(
      n = n,
      G_vec = G_vec,
      Sigma_f = Sigma_f_list[[k]],
      batch_means = design$batch_means[[k]],
      batch_sds = design$batch_sds[[k]],
      lambdas = design$lambdas[[k]],
      thetas = design$thetas[[k]]
    )

    sex_num <- as.numeric(sex_list[[k]]) - 1
    Zk <- out$Z + (sex_num %o% beta_sex)

    Z_list[[k]] <- Zk
    gamma_truth[[k]] <- out$gamma_true
    delta2_truth[[k]] <- out$delta2_true
  }

  comb <- do.call(rbind, Z_list)
  batch <- factor(rep(seq_len(K_site), each = n), labels = paste0("Site", seq_len(K_site)))
  sex <- factor(unlist(lapply(sex_list, as.character)), levels = c("F", "M"))
  mod <- model.matrix(~ sex)

  cb <- neuroCombat::neuroCombat(
    dat = t(comb),
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

  out_cov <- CovBat::covbat(dat = t(comb), bat = batch, mod = mod, n.pc = 10)
  Z_list_std <- split_by_site(comb_std, n_per_site = n, n_site = K_site)

  avg_dists <- sapply(Z_list_std, compute_total_avg_corr_dist, L = L_true)
  ref_site <- which.min(avg_dists)

  hyper_list <- lapply(
    Z_list_std,
    estimate_hyperparams_sgt,
    L = L_true,
    sgt = TRUE
  )
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

  Z_corr_list <- harmonize_sites_sgt_3(
    Z_list = Z_list_std,
    gamma_hats = gamma_hats,
    delta2_hats = delta2_hats,
    Sigma_fhats = Sigma_fhats,
    L = L_true,
    reference = ref_site
  )

  Z_corr_std <- do.call(rbind, Z_corr_list)
  Z_final <- sweep(sweep(Z_corr_std, 2, delta_g, "*"), 2, alpha_g, "+") + Xbeta

  data_list <- list(
    Raw = comb,
    ComBat = t(cb$dat.combat),
    CovBat = t(out_cov$dat.covbat),
    GGCombat = Z_final
  )

  for (m in methods) {
    dat_std <- if (m == "GGCombat") {
      do.call(rbind, Z_corr_list)
    } else {
      sweep(sweep(data_list[[m]] - Xbeta, 2, alpha_g, "-"), 2, delta_g, "/")
    }
    Zm <- split_by_site(dat_std, n_per_site = n, n_site = K_site)
    site_pairs <- utils::combn(seq_len(K_site), 2, simplify = FALSE)

    pair_dist <- purrr::map_dfr(site_pairs, function(pair) {
      d <- compute_distances(Zm[[pair[1]]], Zm[[pair[2]]])
      tibble::tibble(
        Rep = r,
        Method = m,
        Pair = paste(pair, collapse = "-"),
        Metric = names(d),
        Value = as.numeric(d)
      )
    })
    distances_df <- dplyr::bind_rows(distances_df, pair_dist)
  }

  for (k in seq_len(K_site)) {
    inv_order <- seq_len(p)
    site_lab <- paste0("Site", k)
    gamma_true_k <- as.numeric(gamma_truth[[k]])
    delta2_true_k <- as.numeric(delta2_truth[[k]])
    gamma_cb_k <- as.numeric(alpha_g + delta_g * cb$estimates$gamma.star[k, ])
    delta2_cb_k <- as.numeric(delta2_g * cb$estimates$delta.star[k, ])
    gamma_sg_k <- as.numeric(alpha_g + delta_g * gamma_hats[[k]][inv_order])
    delta2_sg_k <- as.numeric(delta2_g * delta2_hats[[k]][inv_order])

    gamma_long <- dplyr::bind_rows(
      gamma_long,
      tibble::tibble(
        rep = r,
        site = site_lab,
        feature = seq_along(gamma_true_k),
        truth = gamma_true_k,
        combat = gamma_cb_k,
        ggcombat = gamma_sg_k
      )
    )

    delta2_long <- dplyr::bind_rows(
      delta2_long,
      tibble::tibble(
        rep = r,
        site = site_lab,
        feature = seq_along(delta2_true_k),
        truth = delta2_true_k,
        combat = delta2_cb_k,
        ggcombat = delta2_sg_k
      )
    )
  }

  sex_pred_loso <- dplyr::bind_rows(
    sex_pred_loso,
    evaluate_methods_loso_sex(
      data_list = data_list,
      batch = batch,
      sex = sex,
      p = p,
      r = r,
      trees = trees_loso,
      mtry = floor(sqrt(p))
    )
  )

  site_pred <- dplyr::bind_rows(
    site_pred,
    evaluate_methods_site_multiclass(
      data_list = data_list,
      batch = batch,
      p = p,
      r = r,
      n_splits = n_splits_site,
      trees = trees_site,
      mtry = floor(sqrt(p))
    )
  )

  message(sprintf("Finished replicate %d/%d; reference site = %d", r, n_reps, ref_site))
}

out_dir <- file.path(root, "study/analysis/01_main_simulation/results/reproduced")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

readr::write_csv(gamma_long, file.path(out_dir, "gamma_long.csv"))
readr::write_csv(delta2_long, file.path(out_dir, "delta2_long.csv"))
readr::write_csv(distances_df, file.path(out_dir, "distances.csv"))
readr::write_csv(site_pred, file.path(out_dir, "site_pred.csv"))
readr::write_csv(sex_pred_loso, file.path(out_dir, "sex_pred_loso.csv"))

meta <- tibble::tibble(
  seed_global = seed_global,
  n_per_site = n,
  p_features = p,
  K_site = K_site,
  G_vec = paste(G_vec, collapse = ","),
  n_reps = n_reps,
  n_splits_site = n_splits_site,
  trees_site = trees_site,
  trees_loso = trees_loso,
  output_dir = out_dir
)
readr::write_csv(meta, file.path(out_dir, "run_meta.csv"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "session_info.txt"))
