# Clean entry point for parameter recovery and varying-dimension simulations.
#
# This script is intentionally not a full rewrite yet. The authoritative
# submitted notebook is `original/harmonization_sim3.Rmd`; this file records the
# formal run structure that should replace it after validation.

root <- normalizePath(file.path(getwd(), "../../.."))
source(file.path(root, "study/analysis/00_common/simulation_common.R"))

load_study_packages(c(
  "dplyr",
  "tibble",
  "purrr",
  "readr",
  "MASS",
  "MCMCpack",
  "neuroCombat",
  "CovBat",
  "ggplot2"
))
source_legacy_method(root)

run_np_parameter_recovery <- function(n_grid = c(50, 100, 250, 500),
                                      p_grid = c(100, 200, 300, 500),
                                      n_reps = 20,
                                      seed = 2026) {
  stop(
    "Not yet validated as a replacement for original/harmonization_sim3.Rmd. ",
    "Use the original notebook for submitted results until this script is fully ported.",
    call. = FALSE
  )
}

message("Parameter-recovery workflow skeleton loaded. Submitted source remains original/harmonization_sim3.Rmd.")
