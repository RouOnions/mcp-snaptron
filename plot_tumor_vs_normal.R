#!/usr/bin/env Rscript
#
# Reads the per-sample CSV written by the MCP tool gene_tumor_vs_normal(gene)
# -- columns: group, rail_id, coverage -- and plots per-sample coverage per
# group (Primary Tumor / Solid Tissue Normal) as a boxplot with jittered
# individual points overlaid, ordered by coverage_median descending, with
# each group's n annotated, plus prints a plain-text summary table.
#
# Boxplot+jitter (not a violin), same reasoning as plot_tissue_expression.R:
# see docs/plot-notes.md.
#
# Usage:
#   Rscript plot_tumor_vs_normal.R <csv_path> [gene_name]
#
# csv_path is the per_sample_csv_path returned by the tool. gene_name is
# only used for the plot title/filename; if omitted it's guessed from the
# CSV filename. See docs/plot-notes.md for the normalization + cohort-size
# caveat that applies to this plot (deliberately not printed on the chart
# itself).

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript plot_tumor_vs_normal.R <csv_path> [gene_name]")
}
csv_path <- args[1]
gene_name <- if (length(args) >= 2) args[2] else sub("^gene_tumor_vs_normal_([A-Za-z0-9]+)_.*$", "\\1", tools::file_path_sans_ext(basename(csv_path)))

per_sample <- read_csv(csv_path, show_col_types = FALSE)

group_summary <- per_sample %>%
  group_by(group) %>%
  summarise(
    n = n(),
    coverage_median = median(coverage),
    coverage_avg = mean(coverage),
    .groups = "drop"
  ) %>%
  arrange(desc(coverage_median))

cat("\n", gene_name, " tumor vs normal (per-sample coverage, NOT library-size normalized -- see docs/plot-notes.md)\n\n", sep = "")
print(as.data.frame(group_summary), row.names = FALSE)
cat("\n")

plot_df <- per_sample %>%
  mutate(group = factor(group, levels = rev(group_summary$group)))

group_max <- plot_df %>%
  group_by(group) %>%
  summarise(max_coverage = max(coverage), .groups = "drop")

label_df <- group_summary %>%
  mutate(group = factor(group, levels = rev(group_summary$group))) %>%
  left_join(group_max, by = "group") %>%
  mutate(label = paste0("n=", comma(n)))

p <- ggplot(plot_df, aes(x = group, y = coverage)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, fill = "#d95f02", alpha = 0.25, color = "#d95f02") +
  geom_jitter(width = 0.15, height = 0, alpha = 0.15, size = 0.7, color = "#d95f02") +
  geom_text(data = label_df, aes(y = max_coverage, label = label), size = 3.5, hjust = -0.2) +
  coord_flip(clip = "off") +
  scale_y_log10(labels = comma, expand = expansion(mult = c(0.05, 0.3))) +
  labs(
    title = paste0(gene_name, " tumor vs normal (TCGA)"),
    x = NULL,
    y = "per-sample coverage (log10 scale)"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.margin = margin(5.5, 40, 5.5, 5.5))

out_dir <- file.path(dirname(csv_path), "..", "plots")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_path <- file.path(out_dir, paste0(gene_name, "_tumor_vs_normal.png"))
ggsave(out_path, p, width = 8, height = 4.5, dpi = 150)

cat("Saved plot to:", normalizePath(out_path), "\n")
