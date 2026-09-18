#!/usr/bin/env Rscript
#
# Reads a per-sample CSV (columns: exon, <group_col>, rail_id, coverage)
# built by fetching /tcgav2/exons or /gtexv2/exons once per value of an
# arbitrary field (e.g. clinical stage within one cancer type), and plots
# a heatmap: exon (genomic position order) x group value (ranked by peak
# coverage), color = log10(coverage_median). Same design as
# html_gene_dossier.py's section6_figure (GTEx exon-usage heatmap) --
# a heatmap, not a multi-line chart, since N exons x M groups is usually
# too many overlapping lines to read.
#
# This is a general-purpose version of that dossier section: it doesn't
# assume the grouping column is GTEx tissue -- it works for anything
# fetched into the same long-format shape (e.g. TCGA clinical stage,
# histology, etc., within one cancer type).
#
# Usage:
#   Rscript plot_exon_usage_by_field.R <csv_path> <group_column_name> [title]
#
# csv_path: the per-sample CSV. group_column_name: the name of the column
# to use as the heatmap's row grouping (e.g. "stage"). title: plot title;
# a generic one is used if omitted.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript plot_exon_usage_by_field.R <csv_path> <group_column_name> [title]")
}
csv_path <- args[1]
group_col <- args[2]
plot_title <- if (length(args) >= 3) args[3] else "Per-exon coverage by group"

per_sample <- read_csv(csv_path, show_col_types = FALSE)
per_sample$group <- per_sample[[group_col]]

summary_df <- per_sample %>%
  group_by(exon, group) %>%
  summarise(coverage_median = median(coverage), n = n(), .groups = "drop")

exon_order <- summary_df %>%
  distinct(exon, .keep_all = TRUE) %>%
  arrange(as.numeric(sub("^chr[0-9XYM]+:([0-9]+)-.*$", "\\1", exon))) %>%
  pull(exon)

# First-appearance order in the CSV (i.e. the order the caller fetched
# values in), not peak-ranked -- appropriate when groups have a natural
# progression (e.g. Stage I < II < III < IV), same choice
# plot_stratified_comparison_single_scope.R makes and for the same reason.
group_order <- unique(per_sample$group)

cat("\n", plot_title, " (per-exon coverage_median, NOT library-size normalized -- see docs/plot-notes.md)\n\n", sep = "")
n_by_group <- per_sample %>% distinct(group, rail_id) %>% count(group, name = "n_samples")
print(as.data.frame(n_by_group), row.names = FALSE)
cat("\n")

summary_df$exon <- factor(summary_df$exon, levels = exon_order)
summary_df$group <- factor(summary_df$group, levels = rev(group_order))

p <- ggplot(summary_df, aes(x = exon, y = group, fill = coverage_median)) +
  geom_tile() +
  scale_fill_gradient(low = "#f7fbff", high = "#08519c", trans = "log10", labels = scales::comma, name = "coverage_median\n(log10)") +
  labs(title = plot_title, x = "exon (genomic position order)", y = NULL) +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6), plot.margin = margin(5.5, 5.5, 5.5, 5.5))

out_dir <- file.path(dirname(csv_path), "..", "plots")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
safe_title <- gsub("[^A-Za-z0-9]+", "_", plot_title)
out_path <- file.path(out_dir, paste0(safe_title, ".png"))
ggsave(out_path, p, width = 12, height = 1.5 + 0.5 * length(group_order), dpi = 150)

cat("Saved plot to:", normalizePath(out_path), "\n")
