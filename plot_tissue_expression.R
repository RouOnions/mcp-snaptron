#!/usr/bin/env Rscript
#
# Reads the per-sample CSV written by the MCP tool
# gene_expression_across_tissues(gene) -- columns: tissue, rail_id, coverage
# -- and plots per-sample coverage per tissue as a boxplot with jittered
# individual points overlaid, ordered by coverage_median descending, plus
# prints a plain-text summary table.
#
# Boxplot+jitter (not a violin) is deliberate: several tissues here have
# many samples tied at the same discrete low value (see docs/plot-notes.md) and
# a few have very small n (e.g. Kidney n=11) -- a violin's kernel-density
# smoothing would misrepresent a pile of exact ties as a fake continuous
# shape, and would misrepresent a small-n tissue as a smooth density when
# there's almost no data behind it. Boxplot+jitter shows both honestly.
#
# Usage:
#   Rscript plot_tissue_expression.R <csv_path> [gene_name]
#
# csv_path is the per_sample_csv_path returned by the tool. gene_name is
# only used for the plot title/filename; if omitted it's guessed from the
# CSV filename. See docs/plot-notes.md for the normalization caveat that applies
# to this plot (deliberately not printed on the chart itself).

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript plot_tissue_expression.R <csv_path> [gene_name]")
}
csv_path <- args[1]
gene_name <- if (length(args) >= 2) args[2] else sub("^gene_expression_across_tissues_([A-Za-z0-9]+)_.*$", "\\1", tools::file_path_sans_ext(basename(csv_path)))

per_sample <- read_csv(csv_path, show_col_types = FALSE)

tissue_summary <- per_sample %>%
  group_by(tissue) %>%
  summarise(
    n = n(),
    coverage_median = median(coverage),
    coverage_avg = mean(coverage),
    .groups = "drop"
  ) %>%
  arrange(desc(coverage_median))

cat("\n", gene_name, " expression across GTEx tissues (per-sample coverage, NOT library-size normalized -- see docs/plot-notes.md)\n\n", sep = "")
print(as.data.frame(tissue_summary), row.names = FALSE)
cat("\n")

plot_df <- per_sample %>%
  mutate(tissue = factor(tissue, levels = rev(tissue_summary$tissue)))

p <- ggplot(plot_df, aes(x = tissue, y = coverage)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, fill = "#2c7fb8", alpha = 0.25, color = "#2c7fb8") +
  geom_jitter(width = 0.15, height = 0, alpha = 0.35, size = 0.9, color = "#2c7fb8") +
  coord_flip() +
  scale_y_log10(labels = comma) +
  labs(
    title = paste0(gene_name, " expression across GTEx tissues"),
    x = NULL,
    y = "per-sample coverage (log10 scale)"
  ) +
  theme_minimal(base_size = 13)

out_dir <- file.path(dirname(csv_path), "..", "plots")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_path <- file.path(out_dir, paste0(gene_name, "_tissue_expression.png"))
ggsave(out_path, p, width = 8, height = 5.5, dpi = 150)

cat("Saved plot to:", normalizePath(out_path), "\n")
