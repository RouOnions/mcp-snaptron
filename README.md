# Snaptron MCP server

Built during MCP workshop at Institut Curie 2026 using claude code.
An MCP server wrapping the [Snaptron](https://snaptron.cs.jhu.edu) REST API
for cancer/tissue RNA-seq questions -- splice junctions,
gene/exon usage, tumor vs normal, tissue specificity across TCGA
(`tcgav2`) and GTEx (`gtexv2`). 

## Citation

Built on [Snaptron](https://snaptron.cs.jhu.edu):
Wilks C, Gaddipati P, Nellore A, Langmead B. "Snaptron: querying splicing
patterns across tens of thousands of RNA-seq samples." *Bioinformatics*
34(1):114–116 (2018). https://doi.org/10.1093/bioinformatics/btx547


## What's here

- **`snaptron_server.py`** -- the MCP server itself. 6 tools:
  `gene_expression_across_tissues`, `gene_tumor_vs_normal`,
  `junction_usage`, `exon_usage`, `stratified_comparison`,
  `generate_gene_dossier` (plus a `list_available_analyses` fallback).
  Each tool writes a tidy per-sample CSV to `outputs/` and, where a matching
  R script exists, auto-generates a PNG chart in `plots/` via `Rscript`.
- **`plot_*.R`** -- six standalone R/ggplot2 scripts, each reading a
  tool's output CSV and rendering one PNG. Callable directly:
  `Rscript plot_tissue_expression.R <csv_path> [gene_name]`.
  `plot_stratified_comparison_single_scope.R` (values as called, one
  scope) is wired into `stratified_comparison`'s auto-plot;
  `plot_stratified_comparison.R` (two scopes side by side, for the
  cancer-type-confound point) is not wired into any tool -- it takes two
  separate `stratified_comparison` calls' CSVs directly.
- **`html_gene_dossier.py`** -- generates one self-contained, interactive
  HTML "dossier" per gene (GTEx tissue ranking, TCGA pan-cancer and
  per-cancer-type tumor/normal, an auto-generated narrative, exon- and
  junction-level detail) into `reports/`. Long-running (~3-10 minutes
  depending on how much is already cached) -- see the module docstring.
- **`html_report_tissue_expression.py`** -- a lighter single-chart
  interactive HTML report (tissue expression only), same idea, faster.
- **`docs/plot-notes.md`** -- caveats and gotchas for each chart, kept off
  the charts themselves so they stay readable; read the relevant section
  before presenting a plot from this repo.

## Setup

### Python

```
pip install -r requirements.txt
```

Needs `mcp>=2.0` specifically -- v1's `FastMCP` was renamed to `MCPServer`
in v2 and this server uses the v2 API only (see `CLAUDE.md`).

### R

The `plot_*.R` scripts (and the auto-plot feature built into
`snaptron_server.py`'s tools) need R with:

```r
install.packages(c("readr", "dplyr", "tidyr", "ggplot2", "scales"))
```

`Rscript` must be reachable -- either on `PATH`, or in one of the standard
per-OS install locations `snaptron_server.py` falls back to (Windows:
`Program Files\R\R-*\bin`; Mac: Homebrew's `/opt/homebrew/bin` or
`/usr/local/bin`, or the official installer's
`/Library/Frameworks/R.framework/Resources/bin`; Linux: `/usr/bin` or
`/usr/local/bin`). If R isn't found anywhere, each tool still returns its
full JSON/CSV output -- only the PNG is skipped, with a note explaining how
to install R.

`plotly` (Python) is what powers the two `html_*.py` report generators --
there is no R-side `plotly` usage in this repo; the R scripts use
`ggplot2`.

## Running the server

```
python snaptron_server.py
```

Runs over stdio, the standard MCP transport -- point any MCP-compatible
client at it.

## Usage examples

Calling a tool (via any MCP client, or directly for local testing):

```python
import snaptron_server as s

# Where is this gene normally expressed? (all 31 confirmed GTEx tissues)
result = s.gene_expression_across_tissues_impl("LRRC10")
# -> {"tissue_summary": [{"tissue": "Heart", "samples_count": 933,
#      "coverage_median": 82099.0, ...}, ...], "per_sample_csv_path": "...",
#      "plot_png_path": "plots/LRRC10_tissue_expression.png", ...}

# Is this gene expressed differently in TCGA tumor vs normal?
result = s.gene_tumor_vs_normal_impl("MKI67")
# -> {"group_summary": [
#      {"group": "Primary Tumor", "samples_count": 9961, "coverage_median": 304248.0, ...},
#      {"group": "Solid Tissue Normal", "samples_count": 740, "coverage_median": 28420.0, ...}],
#      "plot_png_path": "plots/MKI67_tumor_vs_normal.png", ...}

# Does a clinical field (e.g. tumor stage) matter, within ONE cancer type
# (avoids the pan-cancer pooling confound -- see docs/plot-notes.md)?
result = s.stratified_comparison_impl(
    "MKI67", "cgc_case_pathologic_stage", ["Stage I", "Stage IV"],
    extra_field="gdc_cases.project.project_id", extra_value="KIRC",
)
```

Full multi-section HTML report for one gene (~3-10 minutes -- see
"Known limitations" below), either as the `generate_gene_dossier` MCP
tool, or directly from the shell:

```
python html_gene_dossier.py MKI67
```

Writes `reports/MKI67_dossier.html` -- GTEx tissue ranking, TCGA pan-cancer
and per-cancer-type tumor/normal, an auto-generated narrative, exon- and
junction-level detail, all in one self-contained file (opens with no
internet connection, no server needed). Rerunning on the same gene reuses
any matching CSVs already in `outputs/` instead of re-querying Snaptron
(see the module docstring for exactly what's cached and how).

## Known limitations

- **No library-size / depth normalization anywhere.** Values are per-sample
  coverage as Snaptron reports it, not AUC-scaled. A gap of roughly >10x
  between groups is safe to call biological; anything closer (below ~1.3x
  or above ~0.77x) should be treated as a close call, not a confirmed
  effect -- every tool's output carries this caveat explicitly.
- **`exon_usage` is coverage comparison, not differential exon usage.**
  It reports per-exon coverage in each group; it does not model relative
  inclusion rate the way DEXSeq-style differential exon usage would.
- **Not every clinical field applies to every cancer type.**
  `stratified_comparison`'s `field_nomenclature_caveat` flags this on
  every call, but concretely: TCGA-PRAD (prostate) does not populate
  `cgc_case_pathologic_stage` at all (every stage value, including
  substages, returns 0 samples) -- it's staged via `cgc_case_pathologic_t`
  (T-category) instead. Always check a field's `samples_count` for the
  specific cancer type you're using before trusting a result; see
  `CLAUDE.md`'s nomenclature note for the full reasoning.
- **`generate_gene_dossier` takes ~3-10 minutes**, dominated by a
  33-TCGA-project-code sweep in its per-cancer-type section. Warn a user
  before calling it rather than launching it on an offhand mention of a
  gene. Progress is cached incrementally (see `html_gene_dossier.py`'s
  docstring), so a rerun on the same gene is much faster.
