# Graph-scaffold sensitivity analysis for AOAS revision.
#
# Reviewer 2 asks how GG-ComBat behaves when the graph/community scaffold used
# by the method is incorrect or uncertain, and how it differs from Block-ComBat.
# This script uses a submitted-scale simulation: 3 sites, 300 subjects per site,
# 200 features, three correlated 50-feature blocks, and 50 weak/singleton
# features. The final column of L is always the singleton/weak-feature group.

find_project_root <- function(start_dir) {
  current <- normalizePath(start_dir, mustWork = TRUE)
  repeat {
    if (file.exists(file.path(current, "DESCRIPTION")) &&
        dir.exists(file.path(current, "R"))) {
      return(current)
    }
    parent <- dirname(current)
    if (identical(parent, current)) {
      stop("Could not locate ggcombat project root.", call. = FALSE)
    }
    current <- parent
  }
}

args_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(args_file) > 0) {
  sub("^--file=", "", args_file[1])
} else {
  tryCatch(sys.frame(1)$ofile, error = function(e) NA_character_)
}
script_dir <- if (!is.na(script_path) && file.exists(script_path)) dirname(normalizePath(script_path)) else getwd()
root <- find_project_root(script_dir)

out_dir <- file.path(root, "study/analysis/04_revision_robustness/03_graph_scaffold_sensitivity/results")
fig_dir <- file.path(root, "study/analysis/04_revision_robustness/03_graph_scaffold_sensitivity/figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

load_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop("Missing required packages: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(lapply(packages, library, character.only = TRUE))
}

load_packages(c("dplyr", "tidyr", "tibble", "readr", "MASS", "neuroCombat", "CovBat", "ggplot2"))

source(file.path(root, "R/covariance_utils.R"))
source(file.path(root, "R/standardize.R"))
source(file.path(root, "R/fit.R"))
source(file.path(root, "R/harmonize.R"))
source(file.path(root, "R/ggcombat.R"))
source(file.path(root, "R/evaluation_metrics.R"))

seed_global <- as.integer(Sys.getenv("GGCOMBAT_GS_SEED", "20260916"))
n_reps <- as.integer(Sys.getenv("GGCOMBAT_GS_N_REPS", "100"))
n_site <- 3
n_per_site <- as.integer(Sys.getenv("GGCOMBAT_GS_N_PER_SITE", "300"))
block_sizes <- c(50, 50, 50, 50)
n_correlated_blocks <- 3
p <- sum(block_sizes)
n_pc_covbat <- as.integer(Sys.getenv("GGCOMBAT_GS_COVBAT_N_PC", "10"))
max_iter_gg <- as.integer(Sys.getenv("GGCOMBAT_GS_GG_MAX_ITER", "20"))

make_L <- function(labels, n_groups = max(labels)) {
  L <- matrix(0, nrow = length(labels), ncol = n_groups)
  L[cbind(seq_along(labels), labels)] <- 1
  L
}

true_labels <- rep(seq_along(block_sizes), block_sizes)
L_true <- make_L(true_labels, n_groups = length(block_sizes))
singleton_group <- ncol(L_true)

split_by_batch <- function(Y, batch) {
  split_idx <- split(seq_len(nrow(Y)), batch)
  lapply(split_idx, function(idx) Y[idx, , drop = FALSE])
}

make_site_factor_covariances <- function(coherence = 1) {
  base <- list(
    matrix(c(1.00, 0.30, 0.10,
             0.30, 0.90, 0.15,
             0.10, 0.15, 1.10), 3, 3, byrow = TRUE),
    matrix(c(1.20, 0.05, 0.30,
             0.05, 1.10, 0.25,
             0.30, 0.25, 0.85), 3, 3, byrow = TRUE),
    matrix(c(0.85, 0.25, 0.00,
             0.25, 1.20, 0.35,
             0.00, 0.35, 0.95), 3, 3, byrow = TRUE)
  )
  avg <- Reduce(`+`, base) / length(base)
  lapply(base, function(S) spd_project(avg + coherence * (S - avg))$mat)
}

simulate_dataset <- function(rep_id, structure_strength = 1) {
  set.seed(seed_global + rep_id)
  batch <- factor(rep(seq_len(n_site), each = n_per_site), labels = paste0("Site", seq_len(n_site)))
  sex <- unlist(lapply(c(0.45, 0.55, 0.50), function(prob) {
    stats::rbinom(n_per_site, 1, prob)
  }))
  sex_factor <- factor(sex, levels = c(0, 1), labels = c("F", "M"))

  beta_sex <- rep(0, p)
  beta_sex[1:70] <- stats::runif(70, 0.25, 0.55)

  gamma_blocks <- list(
    c(0.55, -0.25, 0.40, 0.05),
    c(-0.20, 0.45, -0.35, -0.05),
    c(-0.45, -0.10, 0.25, 0.00)
  )
  delta_blocks <- list(
    c(0.75, 1.05, 0.85, 1.00),
    c(1.20, 0.90, 1.10, 1.00),
    c(0.90, 1.25, 0.95, 1.00)
  )
  sigmas <- make_site_factor_covariances(coherence = structure_strength)

  Y_list <- vector("list", n_site)
  for (site in seq_len(n_site)) {
    idx <- ((site - 1) * n_per_site + 1):(site * n_per_site)
    factors <- MASS::mvrnorm(n_per_site, mu = rep(0, n_correlated_blocks), Sigma = sigmas[[site]])
    latent <- matrix(0, nrow = n_per_site, ncol = p)
    latent[, 1:50] <- structure_strength * factors[, 1]
    latent[, 51:100] <- structure_strength * factors[, 2]
    latent[, 101:150] <- structure_strength * factors[, 3]

    gamma <- rep(gamma_blocks[[site]], block_sizes) + stats::rnorm(p, sd = 0.06)
    resid_sd <- rep(delta_blocks[[site]], block_sizes) * stats::runif(p, 0.85, 1.15)
    noise <- matrix(stats::rnorm(n_per_site * p), n_per_site, p)
    noise <- sweep(noise, 2, resid_sd, "*")
    bio <- sex[idx] %o% beta_sex
    Y_list[[site]] <- bio + matrix(gamma, n_per_site, p, byrow = TRUE) + latent + noise
  }

  Y <- do.call(rbind, Y_list)
  colnames(Y) <- paste0("feature_", seq_len(p))
  list(Y = Y, batch = batch, sex = sex_factor, beta_sex = beta_sex)
}

safe_combat <- function(Y, batch, mod) {
  out <- neuroCombat::neuroCombat(dat = t(Y), batch = batch, mod = mod, verbose = FALSE)
  t(out$dat.combat)
}

safe_covbat <- function(Y, batch, mod) {
  out <- CovBat::covbat(dat = t(Y), bat = batch, mod = mod, n.pc = n_pc_covbat)
  if (!is.null(out$dat.covbat)) {
    t(out$dat.covbat)
  } else if (!is.null(out$dat.combat)) {
    t(out$dat.combat)
  } else {
    stop("Unexpected CovBat return object.", call. = FALSE)
  }
}

blockwise_combat <- function(Y, batch, mod, L, singleton_group = ncol(L)) {
  out <- matrix(NA_real_, nrow = nrow(Y), ncol = ncol(Y))
  for (group in seq_len(ncol(L))) {
    idx <- which(L[, group] != 0)
    if (length(idx) == 0) next
    if (length(idx) == 1) {
      out[, idx] <- safe_combat(Y[, idx, drop = FALSE], batch, mod)
    } else {
      out[, idx] <- safe_combat(Y[, idx, drop = FALSE], batch, mod)
    }
  }
  colnames(out) <- colnames(Y)
  out
}

run_gg <- function(Y, batch, mod, L) {
  cb <- neuroCombat::neuroCombat(dat = t(Y), batch = batch, mod = mod, verbose = FALSE)
  Z_std <- t(cb$dat.standardized)
  beta_hat <- stats::lm.fit(mod, Y)$coefficients
  bio_hat <- mod[, -1, drop = FALSE] %*% beta_hat[-1, , drop = FALSE]
  Z_list <- split_by_batch(Z_std, batch)
  fits <- fit_ggcombat_params(Z_list, L, singleton_group = ncol(L), max_iter = max_iter_gg)
  gg_list <- harmonize_ggcombat_sites(
    Z_list,
    fits = fits,
    L = L,
    target = "average",
    singleton_group = ncol(L),
    variance_calibration = "empirical_target",
    singleton_calibration = "standardize"
  )
  Y_gg <- do.call(rbind, gg_list)
  Y_gg <- Y_gg + bio_hat
  rownames(Y_gg) <- rownames(Y)
  colnames(Y_gg) <- colnames(Y)
  Y_gg
}

labels_from_L <- function(L) max.col(L, ties.method = "first")

perturb_swap_L <- function(L, prop_swap) {
  labels <- labels_from_L(L)
  candidate <- which(labels != singleton_group)
  n_move <- floor(length(candidate) * prop_swap)
  if (n_move == 0) return(L)
  moved <- sample(candidate, n_move)
  for (idx in moved) {
    labels[idx] <- sample(setdiff(seq_len(n_correlated_blocks), labels[idx]), 1)
  }
  make_L(labels, n_groups = ncol(L))
}

split_first_block_L <- function(L) {
  original <- labels_from_L(L)
  labels <- original
  first <- which(original == 1)
  half <- first[seq_len(floor(length(first) / 2))]
  labels[original == 2] <- 3
  labels[original == 3] <- 4
  labels[original == singleton_group] <- 5
  labels[half] <- 2
  labels[setdiff(first, half)] <- 1
  make_L(labels, n_groups = 5)
}

merge_first_two_blocks_L <- function(L) {
  labels <- labels_from_L(L)
  labels[labels == 2] <- 1
  labels[labels == 3] <- 2
  labels[labels == singleton_group] <- 3
  make_L(labels, n_groups = 3)
}

one_block_L <- function(L) {
  labels <- labels_from_L(L)
  out_labels <- ifelse(labels == singleton_group, 2, 1)
  make_L(out_labels, n_groups = 2)
}

random_L <- function(L) {
  labels <- labels_from_L(L)
  non_singleton <- which(labels != singleton_group)
  labels[non_singleton] <- sample(seq_len(n_correlated_blocks), length(non_singleton), replace = TRUE)
  make_L(labels, n_groups = ncol(L))
}

estimated_L_from_correlation <- function(Y, batch, mod) {
  residual <- stats::lm.fit(mod, Y)$residuals
  corr <- stats::cor(residual)
  corr[is.na(corr)] <- 0
  dist_mat <- stats::as.dist(1 - abs(corr[seq_len(150), seq_len(150), drop = FALSE]))
  hc <- stats::hclust(dist_mat, method = "average")
  labels_cor <- stats::cutree(hc, k = n_correlated_blocks)
  labels <- integer(p)
  labels[seq_len(150)] <- labels_cor
  labels[151:200] <- singleton_group
  make_L(labels, n_groups = singleton_group)
}

residualize_by_mod <- function(Y, mod) {
  stats::lm.fit(mod, Y)$residuals
}

distance_long <- function(Y, batch, mod, method, scaffold = NA_character_, scenario = "block") {
  Y_resid <- residualize_by_mod(Y, mod)
  pair_df <- evaluate_site_distances(split_by_batch(Y_resid, batch))
  dplyr::bind_rows(
    tibble::tibble(method = method, scaffold = scaffold, scenario = scenario, metric = "Mean", value = pair_df$mean_distance),
    tibble::tibble(method = method, scaffold = scaffold, scenario = scenario, metric = "Var", value = pair_df$variance_distance),
    tibble::tibble(method = method, scaffold = scaffold, scenario = scenario, metric = "Frobenius", value = pair_df$frobenius_covariance_distance),
    tibble::tibble(method = method, scaffold = scaffold, scenario = scenario, metric = "Spectral", value = pair_df$spectral_covariance_distance)
  )
}

biology_metrics <- function(Y, dat, method, scaffold = NA_character_, scenario = "block") {
  X <- stats::model.matrix(~ sex, data = data.frame(sex = dat$sex))
  coef <- stats::lm.fit(X, Y)$coefficients
  beta_hat <- coef["sexM", ]
  tibble::tibble(
    method = method,
    scaffold = scaffold,
    scenario = scenario,
    beta_sex_rmse = sqrt(mean((beta_hat - dat$beta_sex)^2)),
    beta_sex_cor = suppressWarnings(stats::cor(beta_hat, dat$beta_sex))
  )
}

summarize_distance <- function(df) {
  df |>
    dplyr::group_by(scenario, method, scaffold, metric) |>
    dplyr::summarise(
      median = stats::median(value, na.rm = TRUE),
      mean = mean(value, na.rm = TRUE),
      q25 = stats::quantile(value, 0.25, na.rm = TRUE),
      q75 = stats::quantile(value, 0.75, na.rm = TRUE),
      .groups = "drop"
    )
}

summarize_biology <- function(df) {
  df |>
    dplyr::group_by(scenario, method, scaffold) |>
    dplyr::summarise(
      beta_sex_rmse_median = stats::median(beta_sex_rmse, na.rm = TRUE),
      beta_sex_cor_median = stats::median(beta_sex_cor, na.rm = TRUE),
      .groups = "drop"
    )
}

set.seed(seed_global)
distance_results <- tibble::tibble()
biology_results <- tibble::tibble()
scaffold_meta <- tibble::tibble()

scenario_grid <- tibble::tibble(
  scenario = c("block_coherent", "weak_block"),
  structure_strength = c(1, 0.25)
)

for (rep_id in seq_len(n_reps)) {
  for (scenario_id in seq_len(nrow(scenario_grid))) {
    scenario_name <- scenario_grid$scenario[[scenario_id]]
    dat <- simulate_dataset(rep_id + scenario_id * 100000, scenario_grid$structure_strength[[scenario_id]])
    mod <- stats::model.matrix(~ sex, data = data.frame(sex = dat$sex))

    L_est <- estimated_L_from_correlation(dat$Y, dat$batch, mod)
    scaffolds <- list(
      oracle_L = L_true,
      estimated_L = L_est,
      swap_10pct = perturb_swap_L(L_true, 0.10),
      swap_30pct = perturb_swap_L(L_true, 0.30),
      split_block = split_first_block_L(L_true),
      merge_blocks = merge_first_two_blocks_L(L_true),
      one_block = one_block_L(L_true),
      random_L = random_L(L_true)
    )

    raw_Y <- dat$Y
    combat_Y <- safe_combat(dat$Y, dat$batch, mod)
    covbat_Y <- safe_covbat(dat$Y, dat$batch, mod)

    method_outputs <- list(
      Raw = list(Y = raw_Y, scaffold = NA_character_),
      ComBat = list(Y = combat_Y, scaffold = NA_character_),
      CovBat = list(Y = covbat_Y, scaffold = NA_character_)
    )

    for (nm in names(method_outputs)) {
      distance_results <- dplyr::bind_rows(
        distance_results,
        distance_long(method_outputs[[nm]]$Y, dat$batch, mod, nm, method_outputs[[nm]]$scaffold, scenario_name) |>
          dplyr::mutate(replicate = rep_id)
      )
      biology_results <- dplyr::bind_rows(
        biology_results,
        biology_metrics(method_outputs[[nm]]$Y, dat, nm, method_outputs[[nm]]$scaffold, scenario_name) |>
          dplyr::mutate(replicate = rep_id)
      )
    }

    for (scaffold_name in names(scaffolds)) {
      L_work <- scaffolds[[scaffold_name]]
      gg_Y <- run_gg(dat$Y, dat$batch, mod, L_work)
      block_Y <- blockwise_combat(dat$Y, dat$batch, mod, L_work)

      for (method_name in c("GG-ComBat", "Blockwise ComBat")) {
        Y_method <- if (method_name == "GG-ComBat") gg_Y else block_Y
        distance_results <- dplyr::bind_rows(
          distance_results,
          distance_long(Y_method, dat$batch, mod, method_name, scaffold_name, scenario_name) |>
            dplyr::mutate(replicate = rep_id)
        )
        biology_results <- dplyr::bind_rows(
          biology_results,
          biology_metrics(Y_method, dat, method_name, scaffold_name, scenario_name) |>
            dplyr::mutate(replicate = rep_id)
        )
      }

      true_members <- labels_from_L(L_true)
      work_members <- labels_from_L(L_work)
      comparable <- seq_len(p)
      scaffold_meta <- dplyr::bind_rows(
        scaffold_meta,
        tibble::tibble(
          replicate = rep_id,
          scenario = scenario_name,
          scaffold = scaffold_name,
          n_groups = ncol(L_work),
          singleton_n = sum(work_members == ncol(L_work)),
          label_agreement = mean(true_members[comparable] == work_members[comparable])
        )
      )
    }

    message(sprintf("Finished replicate %d/%d scenario %s", rep_id, n_reps, scenario_name))
  }
}

distance_summary <- summarize_distance(distance_results)
biology_summary <- summarize_biology(biology_results)

readr::write_csv(distance_results, file.path(out_dir, "graph_scaffold_distances.csv"))
readr::write_csv(distance_summary, file.path(out_dir, "graph_scaffold_distance_summary.csv"))
readr::write_csv(biology_results, file.path(out_dir, "graph_scaffold_biology.csv"))
readr::write_csv(biology_summary, file.path(out_dir, "graph_scaffold_biology_summary.csv"))
readr::write_csv(scaffold_meta, file.path(out_dir, "graph_scaffold_meta.csv"))
readr::write_csv(
  tibble::tibble(
    seed = seed_global,
    n_reps = n_reps,
    n_site = n_site,
    n_per_site = n_per_site,
    p = p,
    block_sizes = paste(block_sizes, collapse = ","),
    n_pc_covbat = n_pc_covbat,
    max_iter_gg = max_iter_gg
  ),
  file.path(out_dir, "run_meta.csv")
)

plot_df <- distance_summary |>
  dplyr::filter(metric %in% c("Frobenius", "Spectral", "Var")) |>
  dplyr::filter(method %in% c("ComBat", "CovBat", "GG-ComBat", "Blockwise ComBat")) |>
  dplyr::mutate(
    display = dplyr::if_else(is.na(scaffold), method, paste(method, scaffold, sep = ": ")),
    display = factor(display, levels = unique(display))
  )

distance_plot <- ggplot2::ggplot(plot_df, ggplot2::aes(x = display, y = median, fill = method)) +
  ggplot2::geom_col(width = 0.75) +
  ggplot2::facet_grid(metric ~ scenario, scales = "free_y") +
  ggplot2::coord_flip() +
  ggplot2::labs(x = NULL, y = "Median pairwise distance") +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(legend.position = "bottom")

ggplot2::ggsave(file.path(fig_dir, "graph_scaffold_distance_bars.pdf"), distance_plot, width = 10, height = 8)
ggplot2::ggsave(file.path(fig_dir, "graph_scaffold_distance_bars.png"), distance_plot, width = 10, height = 8, dpi = 300)

message("Wrote graph-scaffold sensitivity results to ", out_dir)
message("Wrote graph-scaffold sensitivity figures to ", fig_dir)
