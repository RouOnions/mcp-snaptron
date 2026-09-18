#!/usr/bin/env Rscript
#
# Reads the per-junction CSV written by the MCP tool junction_usage(gene)
# -- long format, columns: junction, chromosome, start, end, strand,
# novelty_class, group, samples_count, cohort_total, pct_of_cohort,
# coverage_avg, coverage_median -- pivots it to one row per junction with
# tumor/normal pct_of_cohort and coverage_median side by side, optionally
# drops junctions below a minimum coverage_median in EITHER group (to
# exclude the ~1-read-per-sample "widely observed but trivially supported"
# junctions flagged in docs/plot-notes.md), ranks the remainder by
# |tumor_pct_of_cohort - normal_pct_of_cohort|, and plots the top N as a
# horizontal diverging bar chart: one bar per junction, bar length is the
# signed difference (positive = more common in tumor, negative = more
# common in normal), bar fill is novelty_class.
#
# Usage:
#   Rscript plot_junction_usage.R <csv_path> [gene_name] [top_n] [min_coverage_median]
#
# csv_path is the per_junction_csv_path returned by the tool. gene_name is
# only used for the plot title/filename; if omitted it's guessed from the
# CSV filename. top_n defaults to 10. min_coverage_median defaults to 0
# (no filtering -- reproduces the original unfiltered ranking). See
# docs/plot-notes.md for the normalization caveat that applies to this plot
# (deliberately not printed on the chart itself; the coverage_median
# threshold used IS shown in the title, since that's a stated parameter of
# this run, not an interpretive caveat).

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript plot_junction_usage.R <csv_path> [gene_name] [top_n] [min_coverage_median]")
}
csv_path <- args[1]
gene_name <- if (length(args) >= 2) args[2] else sub("^junction_usage_([A-Za-z0-9]+)_.*$", "\\1", tools::file_path_sans_ext(basename(csv_path)))
top_n <- if (length(args) >= 3) as.integer(args[3]) else 10L
min_coverage_median <- if (length(args) >= 4) as.numeric(args[4]) else 0

long_df <- read_csv(csv_path, show_col_types = FALSE)

wide_df <- long_df %>%
  select(junction, novelty_class, group, pct_of_cohort, coverage_median) %>%
  pivot_wider(
    id_cols = c(junction, novelty_class),
    names_from = group,
    values_from = c(pct_of_cohort, coverage_median),
    values_fill = 0
  ) %>%
  rename(
    tumor_pct = `pct_of_cohort_Primary Tumor`,
    normal_pct = `pct_of_cohort_Solid Tissue Normal`,
    tumor_coverage_median = `coverage_median_Primary Tumor`,
    normal_coverage_median = `coverage_median_Solid Tissue Normal`
  ) %>%
  mutate(pct_diff = tumor_pct - normal_pct)

n_before_filter <- nrow(wide_df)
wide_df <- wide_df %>%
  filter(tumor_coverage_median >= min_coverage_median, normal_coverage_median >= min_coverage_median) %>%
  arrange(desc(abs(pct_diff)))
n_after_filter <- nrow(wide_df)

top_junctions <- head(wide_df, top_n)

cat(
  "\n", gene_name, " top ", nrow(top_junctions), " differential junctions",
  " (coverage_median >= ", min_coverage_median, " required in BOTH groups -- ",
  n_after_filter, " of ", n_before_filter, " junctions passed; see docs/plot-notes.md)\n\n",
  sep = ""
)
print(
  as.data.frame(top_junctions %>% select(junction, novelty_class, tumor_pct, normal_pct, pct_diff, tumor_coverage_median, normal_coverage_median)),
  row.names = FALSE
)
cat("\n")

plot_df <- top_junctions %>%
  mutate(junction = factor(junction, levels = rev(junction)))

title_suffix <- if (min_coverage_median > 0) paste0(" (coverage_median ≥ ", min_coverage_median, ")") else ""

p <- ggplot(plot_df, aes(x = junction, y = pct_diff, fill = novelty_class)) +
  geom_col() +
  geom_hline(yintercept = 0, color = "grey40", linewidth = 0.3) +
  coord_flip() +
  scale_fill_brewer(palette = "Dark2") +
  labs(
    title = paste0(gene_name, " top ", nrow(plot_df), " differential junctions", title_suffix),
    x = NULL,
    y = "tumor_pct_of_cohort - normal_pct_of_cohort",
    fill = "novelty class"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")

out_dir <- file.path(dirname(csv_path), "..", "plots")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
suffix <- if (min_coverage_median > 0) paste0("_min", min_coverage_median) else ""
out_path <- file.path(out_dir, paste0(gene_name, "_junction_usage", suffix, ".png"))
ggsave(out_path, p, width = 9, height = 5.5, dpi = 150)

cat("Saved plot to:", normalizePath(out_path), "\n")
