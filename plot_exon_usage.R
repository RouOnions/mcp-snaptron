#!/usr/bin/env Rscript
#
# Reads the per-sample CSV written by the MCP tool exon_usage(gene) --
# columns: exon, group, rail_id, coverage -- aggregates to per-exon
# tumor/normal coverage_median, computes the tumor/normal RATIO per exon
# (not two raw coverage lines -- the ratio is what makes cross-exon
# comparison legible when absolute coverage varies over orders of
# magnitude along the gene), and plots that ratio along genomic position
# (x-axis = exon start, ascending), LINEAR y-axis (not binned -- binning
# would merge nearby exons and could obscure a gradient that depends on
# per-exon resolution, e.g. MKI67's).
#
# Usage:
#   Rscript plot_exon_usage.R <csv_path> [gene_name]
#
# csv_path is the per_sample_csv_path returned by the tool. gene_name is
# only used for the plot title/filename; if omitted it's guessed from the
# CSV filename. See docs/plot-notes.md for the caveat that applies to this plot
# (deliberately not printed on the chart itself).

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript plot_exon_usage.R <csv_path> [gene_name]")
}
csv_path <- args[1]
gene_name <- if (length(args) >= 2) args[2] else sub("^exon_usage_([A-Za-z0-9]+)_.*$", "\\1", tools::file_path_sans_ext(basename(csv_path)))

per_sample <- read_csv(csv_path, show_col_types = FALSE)

exon_summary <- per_sample %>%
  group_by(exon, group) %>%
  summarise(n = n(), coverage_median = median(coverage), .groups = "drop") %>%
  pivot_wider(id_cols = exon, names_from = group, values_from = c(n, coverage_median)) %>%
  rename(
    tumor_n = `n_Primary Tumor`,
    normal_n = `n_Solid Tissue Normal`,
    tumor_median = `coverage_median_Primary Tumor`,
    normal_median = `coverage_median_Solid Tissue Normal`
  ) %>%
  mutate(
    start = as.numeric(sub("^chr[0-9XYM]+:([0-9]+)-.*$", "\\1", exon)),
    ratio = tumor_median / normal_median
  ) %>%
  arrange(start)

cat("\n", gene_name, " per-exon tumor/normal coverage_median ratio, by genomic position (see docs/plot-notes.md)\n\n", sep = "")
print(as.data.frame(exon_summary %>% select(exon, start, tumor_median, normal_median, ratio)), row.names = FALSE)
cat("\n")
cat("Largest ratio: ", exon_summary$exon[which.max(exon_summary$ratio)], " (", round(max(exon_summary$ratio), 2), "x)\n", sep = "")
cat("Smallest ratio:", exon_summary$exon[which.min(exon_summary$ratio)], "(", round(min(exon_summary$ratio), 2), "x)\n")
cat("\n")

p <- ggplot(exon_summary, aes(x = start, y = ratio)) +
  geom_hline(yintercept = 1, color = "grey50", linewidth = 0.3, linetype = "dashed") +
  geom_line(color = "#377eb8", linewidth = 0.6) +
  geom_point(color = "#377eb8", size = 2.2) +
  scale_x_continuous(labels = comma) +
  scale_y_continuous(labels = comma) +
  labs(
    title = paste0(gene_name, " per-exon tumor/normal coverage ratio"),
    x = "genomic position (exon start)",
    y = "tumor / normal coverage_median"
  ) +
  theme_minimal(base_size = 13)

out_dir <- file.path(dirname(csv_path), "..", "plots")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_path <- file.path(out_dir, paste0(gene_name, "_exon_usage.png"))
ggsave(out_path, p, width = 9, height = 5, dpi = 150)

cat("Saved plot to:", normalizePath(out_path), "\n")
