#!/usr/bin/env Rscript
#
# Reads the per-sample CSV written by the MCP tool stratified_comparison()
# for a SINGLE scope (one field, N values, one extra_field/extra_value
# constraint or none) -- columns: value, rail_id, coverage -- and plots a
# boxplot with jittered individual points per value, same style as
# plot_tumor_vs_normal.R (boxplot+jitter, not violin, n annotated per
# group; log10 x-axis, since per-sample coverage typically spans orders of
# magnitude even when group medians look close -- checked, not assumed:
# confirmed on the KLK3/PRAD test case, where per-sample range was ~5200x
# despite group medians only spanning ~1.8x).
#
# Unlike plot_stratified_comparison.R (which compares two SCOPES --
# pan-cancer vs within-one-cancer-type -- for the same field/values, and
# needs two CSVs), this script is for a single stratified_comparison call:
# one CSV, N values of one field, in whatever scope that call used. It
# does NOT make a cancer-type-confound point; it just shows the values as
# called.
#
# Category order on the y-axis is the ORDER VALUES FIRST APPEAR IN THE CSV
# (i.e. the order they were passed to stratified_comparison), not sorted
# by magnitude -- appropriate when the values have a natural sequence
# (e.g. T-stage T2a < T2b < T2c < T3a < T3b < T4) that's more informative
# to preserve than a magnitude-sorted view would be.
#
# Usage:
#   Rscript plot_stratified_comparison_single_scope.R <csv_path> [title]
#
# csv_path is the per_sample_csv_path returned by stratified_comparison.
# title is used as the plot title; if omitted, a generic one is used.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript plot_stratified_comparison_single_scope.R <csv_path> [title]")
}
csv_path <- args[1]
plot_title <- if (length(args) >= 2) args[2] else "Per-sample coverage by value"

per_sample <- read_csv(csv_path, show_col_types = FALSE)

value_order <- unique(per_sample$value)  # first-appearance order, i.e. caller's own order

group_summary <- per_sample %>%
  group_by(value) %>%
  summarise(n = n(), coverage_median = median(coverage), coverage_avg = mean(coverage), .groups = "drop") %>%
  mutate(value = factor(value, levels = value_order)) %>%
  arrange(value)

cat("\n", plot_title, " (per-sample coverage, NOT library-size normalized -- see docs/plot-notes.md)\n\n", sep = "")
print(as.data.frame(group_summary), row.names = FALSE)
cat("\n")

plot_df <- per_sample %>%
  mutate(value = factor(value, levels = rev(value_order)))

label_df <- group_summary %>%
  mutate(value = factor(value, levels = rev(value_order)), label = paste0("n=", comma(n)))
group_max <- plot_df %>% group_by(value) %>% summarise(max_coverage = max(coverage), .groups = "drop")
label_df <- label_df %>% left_join(group_max, by = "value")

p <- ggplot(plot_df, aes(x = value, y = coverage)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, fill = "#2c7fb8", alpha = 0.25, color = "#2c7fb8") +
  geom_jitter(width = 0.15, height = 0, alpha = 0.25, size = 0.9, color = "#2c7fb8") +
  geom_text(data = label_df, aes(y = max_coverage, label = label), size = 3.5, hjust = -0.2) +
  coord_flip(clip = "off") +
  scale_y_log10(labels = comma, expand = expansion(mult = c(0.05, 0.3))) +
  labs(title = plot_title, x = NULL, y = "per-sample coverage (log10 scale)") +
  theme_minimal(base_size = 13) +
  theme(plot.margin = margin(5.5, 60, 5.5, 5.5))

out_dir <- file.path(dirname(csv_path), "..", "plots")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
safe_title <- gsub("[^A-Za-z0-9]+", "_", plot_title)
out_path <- file.path(out_dir, paste0(safe_title, ".png"))
ggsave(out_path, p, width = 8, height = 1.2 + 0.8 * length(value_order), dpi = 150)

cat("Saved plot to:", normalizePath(out_path), "\n")
