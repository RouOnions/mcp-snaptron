"""
Standalone interactive HTML report generator for gene_expression_across_tissues.

Reads the per-sample CSV written by that MCP tool and produces ONE
self-contained HTML file (Plotly.js embedded INLINE, not loaded from a CDN,
so the file opens and renders with no internet connection and no external
dependencies at render time) with a boxplot + jittered per-sample points
chart, same visual design as plot_tissue_expression.R (boxplot+jitter, not
violin -- see docs/plot-notes.md for why; ordered by coverage_median descending;
log10 x-axis given the huge dynamic range).

This is an ADDITIONAL output format alongside the existing R/ggplot2 PNG
pipeline (plot_tissue_expression.R) -- it does not touch or replace it.
Unlike the PNG charts, this report is a standalone deliverable meant to be
opened and read on its own (e.g. emailed, shared), so it carries the
normalization caveat as visible text in the page body -- still not printed
ON the chart itself, same convention as the R plots, just present in the
surrounding report rather than only in docs/plot-notes.md.

Usage:
    python html_report_tissue_expression.py <csv_path> [gene_name]

csv_path is the per_sample_csv_path returned by
gene_expression_across_tissues. gene_name is only used for the report
title/filename; if omitted it's guessed from the CSV filename.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import pandas as pd
import plotly.graph_objects as go

NORMALIZATION_CAVEAT_HTML = (
    "Values are per-sample coverage, NOT library-size / depth normalized -- "
    "no AUC scaling has been applied. A gap between tissues of roughly "
    "&gt;10x is safe to call biological; a close gap may partly reflect "
    "sequencing depth differences rather than true expression differences."
)


def guess_gene_name(csv_path: Path) -> str:
    match = re.match(r"^gene_expression_across_tissues_([A-Za-z0-9]+)_", csv_path.stem)
    return match.group(1) if match else csv_path.stem


def build_report(csv_path: str, gene_name: str | None = None) -> Path:
    csv_path_obj = Path(csv_path)
    gene_name = gene_name or guess_gene_name(csv_path_obj)

    per_sample = pd.read_csv(csv_path_obj)

    tissue_summary = (
        per_sample.groupby("tissue")["coverage"]
        .agg(n="count", coverage_median="median", coverage_avg="mean")
        .reset_index()
        .sort_values("coverage_median", ascending=False)
    )
    tissue_order = tissue_summary["tissue"].tolist()

    fig = go.Figure()
    fig.add_trace(
        go.Box(
            x=per_sample["coverage"],
            y=per_sample["tissue"],
            orientation="h",
            boxpoints="all",
            jitter=0.5,
            pointpos=0,
            marker=dict(color="#2c7fb8", size=4, opacity=0.35),
            line=dict(color="#2c7fb8"),
            fillcolor="rgba(44,127,184,0.25)",
            hovertemplate="%{y}<br>coverage=%{x}<extra></extra>",
            showlegend=False,
        )
    )
    fig.update_layout(
        title=f"{gene_name} expression across GTEx tissues",
        xaxis_title="per-sample coverage (log10 scale)",
        xaxis_type="log",
        yaxis=dict(
            categoryorder="array",
            # Plotly draws the first category at the bottom of a vertical axis;
            # reverse so the highest-median tissue (first in tissue_order) ends
            # up at the top, matching plot_tissue_expression.R's layout.
            categoryarray=list(reversed(tissue_order)),
        ),
        template="plotly_white",
        height=100 + 40 * len(tissue_order),
        width=900,
        margin=dict(l=140),
    )

    summary_rows_html = "\n".join(
        f"<tr><td>{row.tissue}</td><td>{row.n}</td>"
        f"<td>{row.coverage_median:,.1f}</td><td>{row.coverage_avg:,.1f}</td></tr>"
        for row in tissue_summary.itertuples()
    )

    chart_html = fig.to_html(include_plotlyjs="inline", full_html=False, div_id="tissue-chart")

    page_html = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>{gene_name} expression across GTEx tissues</title>
<style>
  body {{ font-family: -apple-system, Segoe UI, Helvetica, Arial, sans-serif; margin: 0; padding: 24px; background: #fafafa; color: #222; }}
  h1 {{ font-size: 1.4rem; margin-bottom: 4px; }}
  .caveat {{ font-size: 0.85rem; color: #555; max-width: 800px; margin: 8px 0 24px; line-height: 1.4; }}
  table {{ border-collapse: collapse; margin-top: 16px; }}
  th, td {{ padding: 4px 12px; text-align: right; border-bottom: 1px solid #ddd; }}
  th:first-child, td:first-child {{ text-align: left; }}
  th {{ font-weight: 600; }}
</style>
</head>
<body>
<h1>{gene_name} expression across GTEx tissues</h1>
<p class="caveat">{NORMALIZATION_CAVEAT_HTML}</p>
{chart_html}
<table>
<thead><tr><th>tissue</th><th>n</th><th>coverage_median</th><th>coverage_avg</th></tr></thead>
<tbody>
{summary_rows_html}
</tbody>
</table>
</body>
</html>
"""

    out_dir = csv_path_obj.parent.parent / "reports"
    out_dir.mkdir(exist_ok=True)
    out_path = out_dir / f"{gene_name}_tissue_expression.html"
    out_path.write_text(page_html, encoding="utf-8")
    return out_path


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python html_report_tissue_expression.py <csv_path> [gene_name]")
        sys.exit(1)
    csv_arg = sys.argv[1]
    gene_arg = sys.argv[2] if len(sys.argv) >= 3 else None
    output_path = build_report(csv_arg, gene_arg)
    print(f"Saved report to: {output_path.resolve()}")
