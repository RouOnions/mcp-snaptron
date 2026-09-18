#!/usr/bin/env Rscript
#
# Reads two per-sample CSVs written by the MCP tool stratified_comparison()
# -- columns: value, rail_id, coverage -- one pan-cancer (pooled) run and
# one constrained-to-one-cancer-type run of the SAME gene/field/values, and
# plots coverage_median as a grouped bar chart: value (e.g. Stage I / Stage
# IV) on the y-axis, scope (pan-cancer vs within one cancer type) as the
# grouping, side by side. Purpose: make the cancer-type confound visible --
# see docs/plot-notes.md for why the pan-cancer number alone is misleading.
#
# Usage:
#   Rscript plot_stratified_comparison.R <csv1> <label1> <csv2> <label2> [gene_name]
#
# csv1/csv2 are per_sample_csv_path values returned by stratified_comparison
# for the two scopes being compared; label1/label2 are how those scopes
# should be labeled on the chart (e.g. "Pan-cancer" and "Within KIRC").

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript plot_stratified_comparison.R <csv1> <label1> <csv2> <label2> [gene_name]")
}
csv1 <- args[1]
label1 <- args[2]
csv2 <- args[3]
label2 <- args[4]
gene_name <- if (length(args) >= 5) args[5] else sub("^stratified_comparison_([A-Za-z0-9]+)_.*$", "\\1", tools::file_path_sans_ext(basename(csv1)))

read_scope <- function(path, label) {
  read_csv(path, show_col_types = FALSE) %>%
    mutate(scope = label)
}

combined <- bind_rows(read_scope(csv1, label1), read_scope(csv2, label2))

summary_tbl <- combined %>%
  group_by(scope, value) %>%
  summarise(
    n = n(),
    coverage_median = median(coverage),
    coverage_avg = mean(coverage),
    .groups = "drop"
  )

cat("\n", gene_name, " by stage, pan-cancer vs within one cancer type (see docs/plot-notes.md)\n\n", sep = "")
print(as.data.frame(summary_tbl), row.names = FALSE)
cat("\n")

summary_tbl <- summary_tbl %>%
  mutate(
    scope = factor(scope, levels = c(label1, label2)),
    value = factor(value, levels = sort(unique(value)))
  )

p <- ggplot(summary_tbl, aes(x = value, y = coverage_median, fill = scope)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(
    aes(label = comma(coverage_median)),
    position = position_dodge(width = 0.7),
    vjust = -0.4, size = 3.2
  ) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.15))) +
  scale_fill_manual(values = c("#7570b3", "#1b9e77")) +
  labs(
    title = paste0(gene_name, " coverage_median by stage"),
    x = NULL,
    y = "coverage_median",
    fill = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top")

out_dir <- file.path(dirname(csv1), "..", "plots")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_path <- file.path(out_dir, paste0(gene_name, "_stratified_comparison.png"))
ggsave(out_path, p, width = 7, height = 5, dpi = 150)

cat("Saved plot to:", normalizePath(out_path), "\n")
