# Generate manuscript-style figures from reproduced or submitted main simulation results.
#
# By default, this script reads submitted results and writes to
# `figures/reproduced/` rather than overwriting submitted figures.

root <- normalizePath(file.path(getwd(), "../../.."))
source(file.path(root, "study/analysis/00_common/simulation_common.R"))
load_study_packages(c("dplyr", "tidyr", "ggplot2", "readr"))

input_set <- Sys.getenv("GGCOMBAT_RESULT_SET", "submitted")
result_dir <- file.path(root, "study/analysis/01_main_simulation/results", input_set)
figure_dir <- file.path(root, "study/analysis/01_main_simulation/figures/reproduced")
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

method_colors <- c(
  "Raw" = "#7A7A7A",
  "ComBat" = "#3B6FB6",
  "CovBat" = "#4C9A8A",
  "GG-ComBat" = "#C05A4D"
)

distances <- readr::read_csv(file.path(result_dir, "distances.csv"), show_col_types = FALSE)
if (!"Method" %in% names(distances)) {
  # Legacy submitted CSV accidentally dropped Method; do not infer silently.
  stop(
    "The selected distances.csv does not contain Method. Use reproduced outputs ",
    "or repair submitted result provenance before regenerating figures.",
    call. = FALSE
  )
}

dist_summary <- distances |>
  dplyr::mutate(
    Method = dplyr::recode(Method, GGCombat = "GG-ComBat"),
    Method = factor(Method, levels = names(method_colors)),
    Metric = factor(Metric, levels = c("Mean", "Var", "Frobenius", "Spectral"))
  ) |>
  dplyr::group_by(Rep, Method, Metric) |>
  dplyr::summarise(Value = mean(Value), .groups = "drop")

p_dist <- ggplot2::ggplot(dist_summary, ggplot2::aes(x = Method, y = Value, fill = Method)) +
  ggplot2::geom_boxplot(outlier.size = 1, linewidth = 0.3, alpha = 0.85) +
  ggplot2::facet_wrap(~ Metric, scales = "free_y", nrow = 1) +
  ggplot2::scale_fill_manual(values = method_colors) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(
    text = ggplot2::element_text(family = "Times"),
    strip.background = ggplot2::element_rect(fill = "white", color = NA),
    strip.text = ggplot2::element_text(size = 12, face = "bold"),
    axis.title.x = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_text(size = 12),
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = ggplot2::element_text(size = 10),
    legend.position = "none",
    panel.grid.major = ggplot2::element_line(linewidth = 0.3, color = "gray90"),
    panel.grid.minor = ggplot2::element_blank(),
    panel.spacing = grid::unit(0.9, "lines"),
    plot.title = ggplot2::element_blank()
  ) +
  ggplot2::labs(y = "Distance")

ggplot2::ggsave(
  file.path(figure_dir, "fig_sim_distance.pdf"),
  p_dist,
  width = 7.0,
  height = 2.3,
  units = "in",
  device = "pdf"
)

message("Wrote reproduced distance figure to: ", figure_dir)
