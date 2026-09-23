# Summarize and visualize the operating-range simulation.

load_or_stop <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Missing required package: ", package, call. = FALSE)
  }
  suppressPackageStartupMessages(library(package, character.only = TRUE))
}

invisible(lapply(c("dplyr", "tidyr", "readr", "ggplot2", "tibble"), load_or_stop))

result_dir <- file.path(getwd(), "results")
figure_dir <- file.path(getwd(), "figures")
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

summary_path <- file.path(result_dir, "operating_range_distance_summary.csv")
if (!file.exists(summary_path)) {
  stop("Missing ", summary_path, ". Run 01_run_operating_range.R first.", call. = FALSE)
}

distance_summary <- readr::read_csv(summary_path, show_col_types = FALSE)

wide <- distance_summary |>
  dplyr::select(
    condition_id, additive_scale, variance_scale, covariance_scale,
    block_coherence, Method, Metric, median
  ) |>
  tidyr::pivot_wider(names_from = Method, values_from = median)

ratio_summary <- wide |>
  dplyr::transmute(
    condition_id, additive_scale, variance_scale, covariance_scale,
    block_coherence, Metric,
    gg_over_raw = `GG-ComBat` / Raw,
    gg_over_combat = `GG-ComBat` / ComBat,
    gg_over_covbat = `GG-ComBat` / CovBat
  )

overall_ratio_table <- ratio_summary |>
  dplyr::group_by(Metric) |>
  dplyr::summarise(
    median_gg_over_raw = stats::median(gg_over_raw),
    q25_gg_over_raw = stats::quantile(gg_over_raw, 0.25),
    q75_gg_over_raw = stats::quantile(gg_over_raw, 0.75),
    median_gg_over_combat = stats::median(gg_over_combat),
    q25_gg_over_combat = stats::quantile(gg_over_combat, 0.25),
    q75_gg_over_combat = stats::quantile(gg_over_combat, 0.75),
    median_gg_over_covbat = stats::median(gg_over_covbat),
    q25_gg_over_covbat = stats::quantile(gg_over_covbat, 0.25),
    q75_gg_over_covbat = stats::quantile(gg_over_covbat, 0.75),
    .groups = "drop"
  )

best_method_counts <- distance_summary |>
  dplyr::group_by(condition_id, Metric) |>
  dplyr::slice_min(median, n = 1, with_ties = FALSE) |>
  dplyr::ungroup() |>
  dplyr::count(Metric, Method, name = "n_conditions")

selected_conditions <- distance_summary |>
  dplyr::filter(condition_id %in% c(
    "a1_v1_c1_r1",
    "a2_v2_c2_r1",
    "a1_v1_c1_r0.25"
  )) |>
  dplyr::select(condition_id, Method, Metric, median, q25, q75) |>
  dplyr::arrange(condition_id, Metric, median)

readr::write_csv(ratio_summary, file.path(result_dir, "operating_range_ratio_summary.csv"))
readr::write_csv(overall_ratio_table, file.path(result_dir, "operating_range_overall_ratio_table.csv"))
readr::write_csv(best_method_counts, file.path(result_dir, "operating_range_best_method_counts.csv"))
readr::write_csv(selected_conditions, file.path(result_dir, "operating_range_selected_conditions.csv"))

plot_df <- ratio_summary |>
  dplyr::filter(Metric %in% c("Frobenius", "Spectral", "Var")) |>
  tidyr::pivot_longer(
    cols = c(gg_over_combat, gg_over_covbat),
    names_to = "comparison",
    values_to = "ratio"
  ) |>
  dplyr::mutate(
    comparison = dplyr::recode(
      comparison,
      gg_over_combat = "GG-ComBat / ComBat",
      gg_over_covbat = "GG-ComBat / CovBat"
    ),
    ratio = ratio
  ) |>
  dplyr::group_by(Metric, comparison, covariance_scale, block_coherence) |>
  dplyr::summarise(ratio = stats::median(ratio), .groups = "drop") |>
  dplyr::mutate(
    Metric = factor(Metric, levels = c("Frobenius", "Spectral", "Var")),
    block_coherence = factor(block_coherence, levels = sort(unique(block_coherence))),
    covariance_scale = factor(covariance_scale, levels = sort(unique(covariance_scale)))
  )

heatmap_plot <- ggplot2::ggplot(
  plot_df,
  ggplot2::aes(x = covariance_scale, y = block_coherence, fill = ratio)
) +
  ggplot2::geom_tile(color = "white", linewidth = 0.25) +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", ratio)), size = 2.8) +
  ggplot2::facet_grid(Metric ~ comparison) +
  ggplot2::scale_fill_gradient2(
    low = "#2166AC",
    mid = "#F7F7F7",
    high = "#B2182B",
    midpoint = 1,
    limits = c(0, max(1, max(plot_df$ratio, na.rm = TRUE))),
    name = "Ratio"
  ) +
  ggplot2::labs(
    x = "Covariance-effect scale",
    y = "Block coherence"
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    strip.background = ggplot2::element_rect(fill = "grey90", color = NA),
    legend.position = "right"
  )

ggplot2::ggsave(
  file.path(figure_dir, "operating_range_ratio_heatmap.pdf"),
  heatmap_plot,
  width = 7.5,
  height = 6.5
)
ggplot2::ggsave(
  file.path(figure_dir, "operating_range_ratio_heatmap.png"),
  heatmap_plot,
  width = 7.5,
  height = 6.5,
  dpi = 300
)

message("Wrote operating-range summaries to ", result_dir)
message("Wrote draft figures to ", figure_dir)
