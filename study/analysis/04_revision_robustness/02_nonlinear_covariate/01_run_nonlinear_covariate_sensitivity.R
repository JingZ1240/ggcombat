# Nonlinear covariate sensitivity analysis for AOAS revision.
#
# Reviewer 1 Major Comment 2 asks whether GG-ComBat's ComBat-style
# standardization can preserve nonlinear biological covariate effects, especially
# when covariates such as age are imbalanced by site. This script isolates that
# issue: the graph L is correctly specified, while the covariate design matrix is
# either misspecified (linear age only) or expanded to match a nonlinear effect.

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
  sub("~\\+~", " ", sub("^--file=", "", args_file[1]), fixed = FALSE)
} else {
  tryCatch(sys.frame(1)$ofile, error = function(e) NA_character_)
}
script_dir <- if (!is.na(script_path) && file.exists(script_path)) dirname(normalizePath(script_path)) else getwd()
root <- find_project_root(script_dir)

out_dir <- file.path(root, "study/analysis/04_revision_robustness/02_nonlinear_covariate/results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

load_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop("Missing required packages: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(lapply(packages, library, character.only = TRUE))
}

load_packages(c("dplyr", "tidyr", "tibble", "readr", "MASS", "neuroCombat", "CovBat"))

source(file.path(root, "R/covariance_utils.R"))
source(file.path(root, "R/standardize.R"))
source(file.path(root, "R/fit.R"))
source(file.path(root, "R/harmonize.R"))
source(file.path(root, "R/ggcombat.R"))
source(file.path(root, "R/evaluation_metrics.R"))

set.seed(as.integer(Sys.getenv("GGCOMBAT_NL_SEED", "20260915")))
n_reps <- as.integer(Sys.getenv("GGCOMBAT_NL_N_REPS", "100"))
n_site <- 3
n_per_site <- as.integer(Sys.getenv("GGCOMBAT_NL_N_PER_SITE", "300"))
block_sizes <- c(50, 50, 50, 50)
p <- sum(block_sizes)
n_pc_covbat <- as.integer(Sys.getenv("GGCOMBAT_NL_COVBAT_N_PC", "10"))
max_iter_gg <- as.integer(Sys.getenv("GGCOMBAT_NL_GG_MAX_ITER", "20"))

make_L <- function(block_sizes) {
  L <- matrix(0, nrow = sum(block_sizes), ncol = length(block_sizes))
  start <- 1
  for (k in seq_along(block_sizes)) {
    end <- start + block_sizes[[k]] - 1
    L[start:end, k] <- 1
    start <- end + 1
  }
  L
}

L_true <- make_L(block_sizes)

site_factor_covariances <- function() {
  list(
    matrix(c(1.00, 0.20, 0.10,
             0.20, 0.90, 0.15,
             0.10, 0.15, 1.10), 3, 3, byrow = TRUE),
    matrix(c(1.20, 0.35, 0.00,
             0.35, 1.00, 0.25,
             0.00, 0.25, 0.90), 3, 3, byrow = TRUE),
    matrix(c(0.85, 0.10, 0.30,
             0.10, 1.15, 0.20,
             0.30, 0.20, 1.00), 3, 3, byrow = TRUE)
  )
}

simulate_dataset <- function(rep_id, scenario) {
  set.seed(900000 + rep_id)
  batch <- factor(rep(seq_len(n_site), each = n_per_site), labels = paste0("Site", seq_len(n_site)))
  site_age_mean <- c(-0.8, 0.0, 0.8)
  age_c <- unlist(lapply(site_age_mean, function(mu) {
    pmax(pmin(stats::rnorm(n_per_site, mean = mu, sd = 0.45), 1.4), -1.4)
  }))
  age2 <- age_c^2 - mean(age_c^2)

  beta_age <- rep(0, p)
  beta_age[1:80] <- stats::runif(80, 0.25, 0.55)
  beta_age2 <- rep(0, p)
  if (scenario == "quadratic") {
    beta_age2[1:60] <- stats::runif(60, 0.35, 0.75)
  }

  gamma_blocks <- list(
    c(0.70, -0.25, 0.45, 0.10),
    c(-0.15, 0.55, -0.35, -0.10),
    c(-0.55, -0.20, 0.15, 0.25)
  )
  delta_blocks <- list(
    c(0.75, 1.00, 0.85, 0.95),
    c(1.15, 0.90, 1.10, 0.85),
    c(0.90, 1.20, 0.95, 1.05)
  )
  sigmas <- site_factor_covariances()

  Y_list <- vector("list", n_site)
  for (site in seq_len(n_site)) {
    idx <- ((site - 1) * n_per_site + 1):(site * n_per_site)
    factors <- MASS::mvrnorm(n_per_site, mu = rep(0, 3), Sigma = sigmas[[site]])
    latent <- matrix(0, nrow = n_per_site, ncol = p)
    latent[, 1:50] <- factors[, 1]
    latent[, 51:100] <- factors[, 2]
    latent[, 101:150] <- factors[, 3]

    gamma <- rep(gamma_blocks[[site]], block_sizes) + stats::rnorm(p, sd = 0.05)
    resid_sd <- rep(delta_blocks[[site]], block_sizes) * stats::runif(p, 0.85, 1.15)
    noise <- matrix(stats::rnorm(n_per_site * p), n_per_site, p)
    noise <- sweep(noise, 2, resid_sd, "*")

    bio <- age_c[idx] %o% beta_age + age2[idx] %o% beta_age2
    Y_list[[site]] <- bio + matrix(gamma, n_per_site, p, byrow = TRUE) + latent + noise
  }

  Y <- do.call(rbind, Y_list)
  colnames(Y) <- paste0("feature_", seq_len(p))
  list(
    Y = Y,
    batch = batch,
    age_c = age_c,
    age2 = age2,
    beta_age = beta_age,
    beta_age2 = beta_age2,
    scenario = scenario
  )
}

split_by_batch <- function(Y, batch) {
  split_idx <- split(seq_len(nrow(Y)), batch)
  lapply(split_idx, function(idx) Y[idx, , drop = FALSE])
}

oracle_residualize <- function(Y, dat, scenario) {
  cov_df <- data.frame(age_c = dat$age_c, age2 = dat$age2)
  if (scenario == "linear") {
    X <- stats::model.matrix(~ age_c, data = cov_df)
  } else {
    X <- stats::model.matrix(~ age_c + age2, data = cov_df)
  }
  fit <- stats::lm.fit(X, Y)
  Y - X %*% fit$coefficients
}

distance_summary <- function(Y, batch, dat, scenario, method, covariate_model) {
  residual <- oracle_residualize(Y, dat, scenario)
  pair_df <- evaluate_site_distances(split_by_batch(residual, batch))
  tibble::tibble(
    method = method,
    covariate_model = covariate_model,
    mean_distance = mean(pair_df$mean_distance),
    variance_distance = mean(pair_df$variance_distance),
    frobenius_covariance_distance = mean(pair_df$frobenius_covariance_distance),
    spectral_covariance_distance = mean(pair_df$spectral_covariance_distance)
  )
}

biology_metrics <- function(Y, dat, scenario, method, covariate_model) {
  cov_df <- data.frame(age_c = dat$age_c, age2 = dat$age2)
  if (scenario == "linear") {
    X <- stats::model.matrix(~ age_c, data = cov_df)
  } else {
    X <- stats::model.matrix(~ age_c + age2, data = cov_df)
  }
  coef <- stats::lm.fit(X, Y)$coefficients
  beta_age_hat <- coef["age_c", ]
  beta_age_rmse <- sqrt(mean((beta_age_hat - dat$beta_age)^2))
  beta_age_cor <- suppressWarnings(stats::cor(beta_age_hat, dat$beta_age))

  if (scenario == "linear") {
    beta_age2_rmse <- NA_real_
    beta_age2_cor <- NA_real_
  } else {
    beta_age2_hat <- coef["age2", ]
    beta_age2_rmse <- sqrt(mean((beta_age2_hat - dat$beta_age2)^2))
    beta_age2_cor <- suppressWarnings(stats::cor(beta_age2_hat, dat$beta_age2))
  }

  tibble::tibble(
    method = method,
    covariate_model = covariate_model,
    beta_age_rmse = beta_age_rmse,
    beta_age_cor = beta_age_cor,
    beta_age2_rmse = beta_age2_rmse,
    beta_age2_cor = beta_age2_cor
  )
}

run_combat_like <- function(dat, covariate_model, engine = c("combat", "covbat")) {
  engine <- match.arg(engine)
  cov_df <- data.frame(age_c = dat$age_c, age2 = dat$age2)
  mod <- if (covariate_model == "linear") {
    stats::model.matrix(~ age_c, data = cov_df)
  } else {
    stats::model.matrix(~ age_c + age2, data = cov_df)
  }

  if (engine == "combat") {
    out <- neuroCombat::neuroCombat(
      dat = t(dat$Y),
      batch = dat$batch,
      mod = mod,
      verbose = FALSE
    )
    return(t(out$dat.combat))
  }

  out <- CovBat::covbat(
    dat = t(dat$Y),
    bat = dat$batch,
    mod = mod,
    n.pc = n_pc_covbat
  )
  t(out$dat.covbat)
}

run_gg <- function(dat, covariate_model) {
  covariates <- if (covariate_model == "linear") {
    data.frame(age_c = dat$age_c)
  } else {
    data.frame(age_c = dat$age_c, age2 = dat$age2)
  }
  fit <- ggcombat(
    Y = dat$Y,
    batch = dat$batch,
    covariates = covariates,
    L = L_true,
    target = "average",
    singleton_group = ncol(L_true),
    variance_calibration = "empirical_target",
    singleton_calibration = "standardize",
    tol = 1e-5,
    max_iter = max_iter_gg,
    verbose = FALSE
  )
  fit$Y_harmonized
}

scenarios <- c("linear", "quadratic")
methods_to_run <- expand.grid(
  method = c("ComBat", "CovBat", "GG-ComBat"),
  covariate_model = c("linear", "expanded"),
  stringsAsFactors = FALSE
)

distance_rows <- list()
biology_rows <- list()

for (scenario in scenarios) {
  for (rep_id in seq_len(n_reps)) {
    message("Scenario=", scenario, " replicate=", rep_id, "/", n_reps)
    dat <- simulate_dataset(rep_id = rep_id, scenario = scenario)

    distance_rows[[length(distance_rows) + 1]] <- distance_summary(
      dat$Y, dat$batch, dat, scenario, "Raw", "none"
    ) |>
      dplyr::mutate(scenario = scenario, replicate = rep_id)
    biology_rows[[length(biology_rows) + 1]] <- biology_metrics(
      dat$Y, dat, scenario, "Raw", "none"
    ) |>
      dplyr::mutate(scenario = scenario, replicate = rep_id)

    for (m in seq_len(nrow(methods_to_run))) {
      method <- methods_to_run$method[[m]]
      covariate_model <- methods_to_run$covariate_model[[m]]

      Y_h <- tryCatch({
        if (method == "ComBat") {
          run_combat_like(dat, covariate_model, engine = "combat")
        } else if (method == "CovBat") {
          run_combat_like(dat, covariate_model, engine = "covbat")
        } else {
          run_gg(dat, covariate_model)
        }
      }, error = function(e) {
        warning("Method failed: ", method, " ", covariate_model, " scenario=", scenario,
                " rep=", rep_id, " error=", conditionMessage(e))
        NULL
      })

      if (is.null(Y_h)) next

      distance_rows[[length(distance_rows) + 1]] <- distance_summary(
        Y_h, dat$batch, dat, scenario, method, covariate_model
      ) |>
        dplyr::mutate(scenario = scenario, replicate = rep_id)
      biology_rows[[length(biology_rows) + 1]] <- biology_metrics(
        Y_h, dat, scenario, method, covariate_model
      ) |>
        dplyr::mutate(scenario = scenario, replicate = rep_id)
    }
  }
}

distance_df <- dplyr::bind_rows(distance_rows)
biology_df <- dplyr::bind_rows(biology_rows)

distance_summary_df <- distance_df |>
  tidyr::pivot_longer(
    cols = c(mean_distance, variance_distance, frobenius_covariance_distance, spectral_covariance_distance),
    names_to = "metric",
    values_to = "value"
  ) |>
  dplyr::group_by(scenario, method, covariate_model, metric) |>
  dplyr::summarise(
    n_reps = dplyr::n(),
    median = stats::median(value, na.rm = TRUE),
    mean = mean(value, na.rm = TRUE),
    sd = stats::sd(value, na.rm = TRUE),
    q25 = stats::quantile(value, 0.25, na.rm = TRUE),
    q75 = stats::quantile(value, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

biology_summary_df <- biology_df |>
  tidyr::pivot_longer(
    cols = c(beta_age_rmse, beta_age_cor, beta_age2_rmse, beta_age2_cor),
    names_to = "metric",
    values_to = "value"
  ) |>
  dplyr::group_by(scenario, method, covariate_model, metric) |>
  dplyr::summarise(
    n_reps = sum(!is.na(value)),
    median = stats::median(value, na.rm = TRUE),
    mean = mean(value, na.rm = TRUE),
    sd = stats::sd(value, na.rm = TRUE),
    q25 = stats::quantile(value, 0.25, na.rm = TRUE),
    q75 = stats::quantile(value, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::filter(n_reps > 0)

readr::write_csv(distance_df, file.path(out_dir, "nonlinear_covariate_distances.csv"))
readr::write_csv(biology_df, file.path(out_dir, "nonlinear_covariate_biology.csv"))
readr::write_csv(distance_summary_df, file.path(out_dir, "nonlinear_covariate_distance_summary.csv"))
readr::write_csv(biology_summary_df, file.path(out_dir, "nonlinear_covariate_biology_summary.csv"))

run_meta <- tibble::tibble(
  n_reps = n_reps,
  n_site = n_site,
  n_per_site = n_per_site,
  p = p,
  block_sizes = paste(block_sizes, collapse = ","),
  covbat_n_pc = n_pc_covbat,
  gg_max_iter = max_iter_gg,
  output_dir = normalizePath(out_dir)
)
readr::write_csv(run_meta, file.path(out_dir, "run_meta.csv"))

print(run_meta)
print(distance_summary_df |> dplyr::arrange(scenario, metric, method, covariate_model))
print(biology_summary_df |> dplyr::arrange(scenario, metric, method, covariate_model))
