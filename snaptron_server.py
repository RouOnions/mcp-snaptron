"""
Snaptron MCP server.

Wraps the Snaptron REST API (https://snaptron.cs.jhu.edu) so an LLM agent can
ask cancer/tissue RNA-seq questions (gene expression across GTEx tissues,
tumor vs normal in TCGA, splice junction novelty) without having to know the
API's quirks itself. Those quirks are documented in CLAUDE.md in this repo
and are re-explained inline below wherever the code has to work around one.

Division of labour (see BUILD_PLAN.md):
    - This server only fetches, validates, parses, and tidies data. It
      writes big per-sample tables to CSV and returns a small summary +
      the CSV path.
    - Plotting (ggplot2) and human-readable table rendering happen on the
      calling agent's side, in R, reading the CSV this server writes.

Reading this file if you know R but are new to Python:
    - `def foo(x: str) -> dict:` is like an R function with a type hint on
      the argument and the return value; Python does not enforce these at
      runtime, they are just documentation (and the MCP SDK uses them to
      build the tool's schema for the calling LLM).
    - A dict `{"a": 1}` is like a named list `list(a = 1)` in R.
    - `pandas.DataFrame` is the closest thing Python has to an R data.frame.
"""

from __future__ import annotations

import io
import platform
import re
import shutil
import statistics
import subprocess
import urllib.parse
from datetime import datetime
from pathlib import Path

import pandas as pd
import requests

# mcp>=2.0 renamed FastMCP -> MCPServer (see the ModuleNotFoundError message
# you get if you try the old `from mcp.server.fastmcp import FastMCP` import
# with this SDK version). Checked live against the installed version
# (mcp 2.2.0) rather than assumed -- see the milestone-1 report.
from mcp.server.mcpserver import MCPServer

mcp = MCPServer("snaptron")

BASE_URL = "https://snaptron.cs.jhu.edu"

# Where per-sample CSVs get written. Kept inside the repo so the agent (and
# the user) can find them easily; each call writes a new timestamped file.
OUTPUT_DIR = Path(__file__).parent / "outputs"

# Rscript resolution, cross-platform. Checked PATH first (works everywhere R
# is installed and configured correctly); on this dev machine shutil.which
# returned None even with R installed, since the Windows R installer doesn't
# reliably add Rscript to PATH -- so PATH alone can't be trusted on any OS,
# hence the per-OS fallback locations below. Resolved once at import time;
# None if not found anywhere checked. This deliberately does NOT raise at
# import time even if R is completely missing: R only powers the OPTIONAL
# auto-plot feature (see run_plot_script), and crashing the whole MCP server
# -- all 6 tools, including the ones that don't touch R at all -- over a
# missing plotting dependency would be worse than degrading gracefully. The
# "fail with a clear error, not a silent None" requirement is satisfied at
# the point where the failure actually matters: every run_plot_script call
# returns an actionable message (install steps included) in its notes, not
# just a bare None.
def _find_rscript() -> str | None:
    which = shutil.which("Rscript") or shutil.which("Rscript.exe")
    if which:
        return which

    system = platform.system()
    candidates: list[Path] = []
    if system == "Windows":
        candidates.extend(sorted(Path(r"C:\Program Files\R").glob("R-*/bin/Rscript.exe"), reverse=True))
    elif system == "Darwin":
        candidates.extend(
            [
                Path("/opt/homebrew/bin/Rscript"),  # Homebrew, Apple Silicon
                Path("/usr/local/bin/Rscript"),  # Homebrew, Intel
                Path("/Library/Frameworks/R.framework/Resources/bin/Rscript"),  # official CRAN .pkg installer
            ]
        )
    else:  # Linux and other Unix-likes
        candidates.extend(
            [
                Path("/usr/bin/Rscript"),
                Path("/usr/local/bin/Rscript"),
            ]
        )

    for candidate in candidates:
        if candidate.exists():
            return str(candidate)
    return None


_RSCRIPT_NOT_FOUND_MESSAGE = (
    "Rscript not found -- checked PATH and standard R install locations for this OS "
    "(Windows: Program Files\\R\\R-*\\bin; Mac: /opt/homebrew/bin, /usr/local/bin, "
    "/Library/Frameworks/R.framework/Resources/bin; Linux: /usr/bin, /usr/local/bin). "
    "Install R from https://www.r-project.org/ and ensure Rscript is on PATH (or "
    "installed in one of the locations above) to enable automatic plot generation. "
    "This tool's JSON/CSV output is unaffected -- only the PNG plot is skipped."
)

_RSCRIPT_PATH = _find_rscript()
_PLOT_SCRIPT_DIR = Path(__file__).parent


def run_plot_script(script_name: str, args: list[str]) -> tuple[str | None, str]:
    """Run one of the plot_*.R scripts (each already prints "Saved plot to:
    <path>" on success -- see plot_*.R) and return (png_path_or_None, note).
    Never raises: a plotting failure must not break a tool's core JSON/CSV
    output, so every failure mode here returns None + an explanatory note
    instead."""
    if _RSCRIPT_PATH is None:
        return None, _RSCRIPT_NOT_FOUND_MESSAGE
    script_path = _PLOT_SCRIPT_DIR / script_name
    try:
        result = subprocess.run(
            [_RSCRIPT_PATH, str(script_path), *args],
            capture_output=True,
            text=True,
            timeout=60,
        )
    except subprocess.TimeoutExpired:
        return None, f"{script_name}: timed out after 60s -- plot skipped."
    if result.returncode != 0:
        return None, f"{script_name}: failed (exit {result.returncode}) -- {result.stderr.strip()[:300]}"
    match = re.search(r"Saved plot to: (.+)", result.stdout)
    if not match:
        return None, f"{script_name}: ran but did not report a saved plot path -- {result.stdout.strip()[:300]}"
    return match.group(1).strip(), f"Plot auto-generated via {script_name}."

# The full, live-confirmed panel of GTEx broad-tissue values (field SMTS) in
# gtexv2 -- NOT a guessed/curated subset. Derived by downloading the bulk
# per-sample metadata file (https://snaptron.cs.jhu.edu/data/gtexv2/samples.tsv
# -- a clean, header-matches-columns TSV, unlike the broken /samples API
# endpoint) and taking the distinct values of its SMTS column (19214 rows,
# 133 with no SMTS assigned -- excluded, not a tissue). All 31 values were
# then round-tripped against a live sfilter=SMTS:<value> query on GAPDH and
# checked against the bulk file's own per-value row count as a correctness
# check, not just an existence check -- see GTEX_SMTS_WORD_MATCH_CAVEAT for
# the one mismatch that check found.
GTEX_TISSUES = [
    "Adipose Tissue",
    "Adrenal Gland",
    "Bladder",
    "Blood",
    "Blood Vessel",
    "Bone Marrow",
    "Brain",
    "Breast",
    "Cervix Uteri",
    "Colon",
    "Esophagus",
    "Fallopian Tube",
    "Heart",
    "Kidney",
    "Liver",
    "Lung",
    "Muscle",
    "Nerve",
    "Ovary",
    "Pancreas",
    "Pituitary",
    "Prostate",
    "Salivary Gland",
    "Skin",
    "Small Intestine",
    "Spleen",
    "Stomach",
    "Testis",
    "Thyroid",
    "Uterus",
    "Vagina",
]

# CONFIRMED LIVE (not documented anywhere in Snaptron's own docs): sfilter on
# a text field matches whole WORDS/tokens within the field's value, not the
# full string. "Blood" is one word of "Blood Vessel", so
# sfilter=SMTS:Blood matches BOTH "Blood" (1048 samples in the bulk file)
# AND "Blood Vessel" (1398 samples) -- live samples_count for the "Blood"
# filter on GAPDH is 2446, exactly 1048+1398. Confirmed by testing all 31
# SMTS values above against their true bulk-file counts: this is the ONLY
# collision among them (checked live query count == bulk file count for
# every value; "Blo" alone matches nothing, ruling out plain substring
# matching -- it really is word-boundary, not prefix/substring). If a future
# tissue panel or other sfilter field introduces a value that is a whole
# word inside another value's text, re-run this same live
# query-count-vs-bulk-count check before trusting it -- don't assume this is
# the only instance across every field, only that it's the only one found
# in SMTS for GTEx broad tissue.
GTEX_SMTS_WORD_MATCH_CAVEAT = (
    "sfilter=SMTS:<value> matches whole words inside the SMTS field, not the "
    "exact string. 'Blood' therefore also catches 'Blood Vessel' samples at "
    "the raw API level. This tool corrects for it (see "
    "_correct_blood_smts_contamination) using ground-truth per-sample "
    "metadata from the bulk samples.tsv file, so the 'Blood' result below is "
    "clean -- but this caveat is surfaced because sfilter's word-match "
    "behavior is not documented by Snaptron and could bite a future query "
    "against SMTS or any other text field."
)

# Bulk per-sample metadata file: a clean alternative to the broken /samples
# API endpoint (see CLAUDE.md). Header count matches every row's field
# count -- unlike /samples?ids=..., which drops fields per-row wherever a
# sample lacks a value. Used ONLY to build a rail_id -> SMTS ground-truth
# lookup, to correct the one confirmed sfilter word-match collision above.
# Lazily downloaded once per process (24MB) and cached in memory.
GTEX_SAMPLES_METADATA_URL = "https://snaptron.cs.jhu.edu/data/gtexv2/samples.tsv"
_gtex_rail_id_to_smts_cache: dict[str, str] | None = None


def _get_gtex_rail_id_to_smts() -> dict[str, str]:
    """Return (and cache) a rail_id -> SMTS lookup built from the bulk gtexv2
    sample metadata file. Used to clean the "Blood" sfilter result of
    "Blood Vessel" contamination (see GTEX_SMTS_WORD_MATCH_CAVEAT)."""
    global _gtex_rail_id_to_smts_cache
    if _gtex_rail_id_to_smts_cache is None:
        response = requests.get(GTEX_SAMPLES_METADATA_URL, timeout=120)
        df = pd.read_csv(
            io.StringIO(response.text), sep="\t", dtype=str, usecols=["rail_id", "SMTS"]
        )
        _gtex_rail_id_to_smts_cache = dict(zip(df["rail_id"], df["SMTS"]))
    return _gtex_rail_id_to_smts_cache


def _correct_blood_smts_contamination(
    samples: list[tuple[str, int]]
) -> tuple[list[tuple[str, int]], str]:
    """Filter a "Blood"-sfiltered samples list down to rail_ids whose true
    SMTS (from the bulk metadata ground truth) is exactly "Blood", dropping
    the "Blood Vessel" samples the raw sfilter word-match let through, and
    recompute count/avg/median from the cleaned subset. Returns the cleaned
    samples list and an explanatory note."""
    rail_id_to_smts = _get_gtex_rail_id_to_smts()
    before = len(samples)
    cleaned = [(rail_id, coverage) for rail_id, coverage in samples if rail_id_to_smts.get(rail_id) == "Blood"]
    dropped = before - len(cleaned)
    note = (
        f"Blood: raw sfilter=SMTS:Blood matched {before} samples due to the "
        f"word-match collision with 'Blood Vessel' (see "
        f"GTEX_SMTS_WORD_MATCH_CAVEAT) -- dropped {dropped} true 'Blood Vessel' "
        f"samples using ground-truth metadata, leaving {len(cleaned)} true "
        f"'Blood' samples for the summary below."
    )
    return cleaned, note

NORMALIZATION_CAVEAT = (
    "Values are per-sample coverage (coverage_avg / coverage_median as "
    "reported by Snaptron), NOT library-size / depth normalized -- no AUC "
    "scaling has been applied (that would require per-sample library sizes, "
    "which live in Snaptron's /samples endpoint, and that endpoint's output "
    "is not safely parseable -- see CLAUDE.md). A large gap between groups "
    "(roughly >10x) is safe to call biological; a close gap may partly "
    "reflect sequencing depth differences rather than true expression "
    "differences."
)


# ---------------------------------------------------------------------------
# Shared fetch / validate / parse helpers
# ---------------------------------------------------------------------------


def _build_query(params: list[tuple[str, str]]) -> str:
    """Build a query string that matches Snaptron's confirmed-working format.

    Snaptron's own examples use a literal ':' inside filter values (e.g.
    `sfilter=cgc_sample_sample_type:Primary%20Tumor`) with only the space
    percent-encoded. `requests`' own `params=` dict encoding is not used
    here because it is not guaranteed to leave ':' unescaped, and this API
    has shown itself to be fragile around anything not tested exactly as
    documented -- so we replicate the confirmed-working encoding by hand.
    """
    parts = []
    for key, value in params:
        encoded_value = urllib.parse.quote(str(value), safe=":")
        parts.append(f"{key}={encoded_value}")
    return "&".join(parts)


def fetch_snaptron(
    compilation: str, endpoint: str, region: str, filters: list[tuple[str, str]] | None = None
) -> str:
    """GET one Snaptron endpoint and return the validated payload text.

    Two things this MUST do (see CLAUDE.md -- these are the load-bearing
    gotchas):

    1. Scan the body for a Python traceback. Snaptron returns HTTP 200 even
       for a malformed query and puts the traceback in the body instead of
       using a real error status code, so `response.status_code` cannot be
       trusted to detect failure.

    2. Only strip a leading plain-text status line (e.g. "200 OK\\n") if one
       is actually present. CLAUDE.md states every response has this
       prefix; live testing during this build showed that is only true for
       some endpoints (confirmed on /snaptron/registry) and NOT for the
       data endpoints this server actually calls (/genes, /snaptron,
       /exons all returned the TSV header as line 1, with no status-line
       prefix). Blindly stripping "line 1" on those endpoints would have
       silently eaten the real TSV header row. So: detect a status line by
       pattern (three digits + a space, e.g. "200 OK") rather than
       assuming it is always there.
    """
    query = _build_query([("regions", region)] + (filters or []))
    url = f"{BASE_URL}/{compilation}/{endpoint}?{query}"
    response = requests.get(url, timeout=60)
    body = response.text

    if "Traceback (most recent call last)" in body:
        raise RuntimeError(
            f"Snaptron returned an error for this query (HTTP status was "
            f"{response.status_code}, but the body contains a Python "
            f"traceback, which Snaptron uses instead of a real error "
            f"status). URL: {url}\n\nBody excerpt:\n{body[:1500]}"
        )

    first_line, _, rest = body.partition("\n")
    if re.match(r"^\d{3}\s", first_line):
        body = rest

    return body


def parse_tsv(text: str) -> pd.DataFrame:
    """Parse a Snaptron TSV payload into a DataFrame. Empty body -> empty df."""
    if not text.strip():
        return pd.DataFrame()
    return pd.read_csv(io.StringIO(text), sep="\t", dtype=str)


def split_packed_gene_column(df: pd.DataFrame) -> pd.DataFrame:
    """Split the packed `gene_id:gene_name:gene_type:bp_length` column.

    /genes and /exons responses pack four fields into one column instead of
    giving each its own. `regions=<GENE>` also returns every gene
    overlapping that gene's coordinate span, not just the named gene, so
    the caller still needs to filter rows down to gene_name == the gene
    they asked for.
    """
    packed_col = "gene_id:gene_name:gene_type:bp_length"
    parts = df[packed_col].str.split(":", n=3, expand=True)
    df = df.copy()
    df["gene_id"] = parts[0]
    df["gene_name"] = parts[1]
    df["gene_type"] = parts[2]
    df["bp_length"] = parts[3]
    return df


def parse_samples_column(samples_str: str) -> list[tuple[str, int]]:
    """Parse the packed `,railid:coverage,railid:coverage,...` samples column.

    This is how per-sample coverage is obtained without touching the
    broken /bases or /samples endpoints. The leading comma produces one
    empty token when split, which is silently skipped.
    """
    if not isinstance(samples_str, str) or not samples_str.strip():
        return []
    out: list[tuple[str, int]] = []
    for token in samples_str.split(","):
        token = token.strip()
        if not token:
            continue
        rail_id, _, coverage = token.partition(":")
        if not coverage:
            continue
        out.append((rail_id, int(coverage)))
    return out


def write_tidy_csv(rows: list[dict], filename_prefix: str) -> str | None:
    """Write per-sample rows to a timestamped CSV in OUTPUT_DIR. Returns the
    absolute path, or None if there were no rows to write."""
    if not rows:
        return None
    OUTPUT_DIR.mkdir(exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    csv_path = OUTPUT_DIR / f"{filename_prefix}_{timestamp}.csv"
    pd.DataFrame(rows).to_csv(csv_path, index=False)
    return str(csv_path.resolve())


# ---------------------------------------------------------------------------
# Tool 1: gene_expression_across_tissues
# ---------------------------------------------------------------------------


def gene_expression_across_tissues_impl(gene: str) -> dict:
    """Implementation, kept separate from the @mcp.tool wrapper below so it
    can be called and tested directly from a plain Python session (no MCP
    client needed) as well as through the MCP server.

    `gene` accepts either a gene symbol (e.g. "MKI67") or a raw region (e.g.
    "chr10:128096659-128099255") -- same as exon_usage (see
    _looks_like_region). With a region, client-side gene-symbol filtering is
    skipped (there is no symbol to match), and the returned "gene" field
    reports the actual gene symbol(s) found rather than the region string."""
    gene = gene.strip()
    is_region_query = _looks_like_region(gene)
    tissue_summary: list[dict] = []
    per_sample_rows: list[dict] = []
    notes: list[str] = []
    genes_found: set[str] = set()

    for tissue in GTEX_TISSUES:
        try:
            body = fetch_snaptron("gtexv2", "genes", gene, filters=[("sfilter", f"SMTS:{tissue}")])
        except RuntimeError as exc:
            notes.append(f"{tissue}: query failed -- {exc}")
            continue

        df = parse_tsv(body)
        if df.empty:
            notes.append(f"{tissue}: no rows returned")
            continue

        df = split_packed_gene_column(df)
        if is_region_query:
            matches = df
        else:
            matches = df[df["gene_name"].str.upper() == gene.upper()]
        if matches.empty:
            notes.append(
                f"{tissue}: gene '{gene}' not found among {len(df)} features "
                f"overlapping its coordinate span"
            )
            continue
        if len(matches) > 1:
            notes.append(
                f"{tissue}: {len(matches)} rows matched '{gene}' "
                f"(likely multiple gene model versions, or a region overlapping "
                f"more than one gene) -- using the first"
            )

        row = matches.iloc[0]
        genes_found.add(row["gene_name"])
        samples = parse_samples_column(row["samples"])

        if tissue == "Blood":
            samples, correction_note = _correct_blood_smts_contamination(samples)
            notes.append(correction_note)

        if not samples:
            notes.append(f"{tissue}: no samples left after parsing, skipping")
            continue

        coverages = [coverage for _, coverage in samples]
        samples_count = len(samples)
        coverage_avg = statistics.mean(coverages)
        coverage_median = statistics.median(coverages)

        tissue_summary.append(
            {
                "tissue": tissue,
                "samples_count": samples_count,
                "coverage_avg": coverage_avg,
                "coverage_median": coverage_median,
            }
        )

        for rail_id, coverage in samples:
            per_sample_rows.append({"tissue": tissue, "rail_id": rail_id, "coverage": coverage})

    tissue_summary.sort(key=lambda r: r["coverage_median"], reverse=True)

    resolved_gene = ", ".join(sorted(genes_found)) if is_region_query else gene
    if is_region_query and not genes_found:
        notes.append(f"{gene}: no annotated gene found overlapping this region")
    csv_filename_gene = re.sub(r"[^A-Za-z0-9._-]", "_", gene)  # region strings contain ':' -- not valid in Windows filenames
    plot_gene_label = resolved_gene if resolved_gene else csv_filename_gene

    csv_path = write_tidy_csv(per_sample_rows, f"gene_expression_across_tissues_{csv_filename_gene}")

    plot_png_path = None
    if csv_path:
        plot_png_path, plot_note = run_plot_script("plot_tissue_expression.R", [csv_path, plot_gene_label])
        notes.append(plot_note)

    return {
        "gene": resolved_gene or None,
        "query": gene,
        "compilation": "gtexv2",
        "tissue_summary": tissue_summary,
        "per_sample_csv_path": csv_path,
        "per_sample_csv_columns": ["tissue", "rail_id", "coverage"],
        "plot_png_path": plot_png_path,
        "notes": notes,
        "normalization_caveat": NORMALIZATION_CAVEAT,
        "smts_word_match_caveat": GTEX_SMTS_WORD_MATCH_CAVEAT,
    }


@mcp.tool()
def gene_expression_across_tissues(gene: str) -> dict:
    """Rank a gene's expression across the full panel of GTEx normal tissues.

    `gene` accepts a gene symbol (e.g. "MKI67") OR a raw region (e.g.
    "chr10:128096659-128099255", format chr:start-end) -- confirmed live
    both work. With a region, the result's "gene" field reports the actual
    gene symbol(s) found overlapping it (not the region string), and
    "query" echoes back exactly what was passed in.

    Queries GTEx (compilation gtexv2) for the given gene symbol once per
    tissue in the full, live-confirmed SMTS (broad tissue) panel -- all 31
    distinct values found in gtexv2's bulk sample metadata, not a curated
    subset (see GTEX_TISSUES). That's 31 sequential HTTP requests, so this
    call is slower than a single query -- budget roughly 1.5-4s per tissue
    (~70-100s total; see CLAUDE.md's speed note).

    Returns a per-tissue summary (samples_count, coverage_avg,
    coverage_median) ranked by coverage_median, plus the absolute path to a
    CSV of every individual sample's coverage for further analysis/plotting
    in R. Use this to answer "where is gene X normally expressed?" or to
    check tissue-specificity of a gene.

    coverage_avg/coverage_median here are computed locally from each
    tissue's per-sample coverage values (not read directly off the API
    response) so that the "Blood" result can be corrected for a confirmed
    API bug: sfilter=SMTS:Blood word-matches "Blood Vessel" too and silently
    pools both tissues' samples together. See smts_word_match_caveat in the
    result for detail. All other 30 tissues were validated to have no such
    collision.
    """
    return gene_expression_across_tissues_impl(gene)


# ---------------------------------------------------------------------------
# Tool 2: gene_tumor_vs_normal
# ---------------------------------------------------------------------------

TCGA_SAMPLE_TYPE_FIELD = "cgc_sample_sample_type"
# Exact GDC vocabulary strings -- case-sensitive, confirmed in CLAUDE.md.
TCGA_GROUPS = ["Primary Tumor", "Solid Tissue Normal"]

# tcgav2 whole-cohort sample counts per sample-type group. There is no
# Snaptron endpoint that returns this directly (see CLAUDE.md). Confirmed by
# querying gene_tumor_vs_normal on four broadly/ubiquitously-expressed genes
# (GAPDH, ACTB, MKI67, TP53): every one of them returned samples_count 9961
# for Primary Tumor and 740 for Solid Tissue Normal. A gene-specific count
# would vary across genes of very different expression breadth; getting the
# identical number four times only makes sense if this IS the group's total
# cohort size. Used by junction_usage to normalize a junction's usage as a
# fraction of its group's cohort, instead of a raw (cohort-size-confounded)
# sample count.
TCGA_COHORT_TOTALS = {
    "Primary Tumor": 9961,
    "Solid Tissue Normal": 740,
}


def _fetch_tcga_gene_row(gene: str, sample_type_value: str, is_region_query: bool = False) -> tuple[dict | None, list[str]]:
    """Fetch /tcgav2/genes filtered to one TCGA sample-type group and return
    the row matching the requested gene symbol (or None + a note explaining
    why not). `is_region_query`: same as exon_usage's _looks_like_region --
    when `gene` is a raw chr:start-end region rather than a symbol, skip the
    gene-name filter entirely (there is no symbol to match)."""
    notes: list[str] = []
    try:
        body = fetch_snaptron(
            "tcgav2", "genes", gene, filters=[("sfilter", f"{TCGA_SAMPLE_TYPE_FIELD}:{sample_type_value}")]
        )
    except RuntimeError as exc:
        notes.append(f"{sample_type_value}: query failed -- {exc}")
        return None, notes

    df = parse_tsv(body)
    if df.empty:
        notes.append(f"{sample_type_value}: no rows returned")
        return None, notes

    df = split_packed_gene_column(df)
    matches = df if is_region_query else df[df["gene_name"].str.upper() == gene.upper()]
    if matches.empty:
        notes.append(
            f"{sample_type_value}: gene '{gene}' not found among {len(df)} features "
            f"overlapping its coordinate span"
        )
        return None, notes
    if len(matches) > 1:
        notes.append(
            f"{sample_type_value}: {len(matches)} rows matched '{gene}' "
            f"(likely multiple gene model versions, or a region overlapping more "
            f"than one gene) -- using the first"
        )

    return matches.iloc[0].to_dict(), notes


def gene_tumor_vs_normal_impl(gene: str) -> dict:
    """Implementation, kept separate from the @mcp.tool wrapper for direct
    testability (see gene_expression_across_tissues_impl for why).

    `gene` accepts either a gene symbol (e.g. "MKI67") or a raw region (e.g.
    "chr10:128096659-128099255") -- same as exon_usage (see
    _looks_like_region)."""
    gene = gene.strip()
    is_region_query = _looks_like_region(gene)
    group_summary: list[dict] = []
    per_sample_rows: list[dict] = []
    notes: list[str] = []
    genes_found: set[str] = set()

    for sample_type_value in TCGA_GROUPS:
        row, row_notes = _fetch_tcga_gene_row(gene, sample_type_value, is_region_query)
        notes.extend(row_notes)
        if row is None:
            continue
        genes_found.add(row["gene_name"])

        try:
            samples_count = int(row["samples_count"])
            coverage_avg = float(row["coverage_avg"])
            coverage_median = float(row["coverage_median"])
        except (ValueError, TypeError):
            notes.append(f"{sample_type_value}: could not parse coverage fields, skipping")
            continue

        group_summary.append(
            {
                "group": sample_type_value,
                "samples_count": samples_count,
                "coverage_avg": coverage_avg,
                "coverage_median": coverage_median,
            }
        )

        for rail_id, coverage in parse_samples_column(row["samples"]):
            per_sample_rows.append({"group": sample_type_value, "rail_id": rail_id, "coverage": coverage})

    cohort_note = None
    by_group = {g["group"]: g["samples_count"] for g in group_summary}
    if "Primary Tumor" in by_group and "Solid Tissue Normal" in by_group and by_group["Solid Tissue Normal"] > 0:
        ratio = by_group["Primary Tumor"] / by_group["Solid Tissue Normal"]
        cohort_note = (
            f"Cohort sizes are imbalanced: {by_group['Primary Tumor']} tumor vs "
            f"{by_group['Solid Tissue Normal']} normal samples ({ratio:.1f}x more tumor). "
            f"This is exactly why coverage_avg/coverage_median (per-sample) are reported "
            f"instead of coverage_sum, which would be dominated by cohort size rather than "
            f"biology."
        )

    resolved_gene = ", ".join(sorted(genes_found)) if is_region_query else gene
    if is_region_query and not genes_found:
        notes.append(f"{gene}: no annotated gene found overlapping this region")
    csv_filename_gene = re.sub(r"[^A-Za-z0-9._-]", "_", gene)  # region strings contain ':' -- not valid in Windows filenames
    plot_gene_label = resolved_gene if resolved_gene else csv_filename_gene

    csv_path = write_tidy_csv(per_sample_rows, f"gene_tumor_vs_normal_{csv_filename_gene}")

    plot_png_path = None
    if csv_path:
        plot_png_path, plot_note = run_plot_script("plot_tumor_vs_normal.R", [csv_path, plot_gene_label])
        notes.append(plot_note)

    return {
        "gene": resolved_gene or None,
        "query": gene,
        "compilation": "tcgav2",
        "group_summary": group_summary,
        "cohort_size_note": cohort_note,
        "per_sample_csv_path": csv_path,
        "per_sample_csv_columns": ["group", "rail_id", "coverage"],
        "plot_png_path": plot_png_path,
        "notes": notes,
        "normalization_caveat": NORMALIZATION_CAVEAT,
    }


@mcp.tool()
def gene_tumor_vs_normal(gene: str) -> dict:
    """Compare a gene's expression between TCGA primary tumor and matched
    solid tissue normal samples (pan-cancer, compilation tcgav2).

    `gene` accepts a gene symbol (e.g. "MKI67") OR a raw region (e.g.
    "chr10:128096659-128099255", format chr:start-end) -- confirmed live
    both work. With a region, the result's "gene" field reports the actual
    gene symbol(s) found overlapping it (not the region string), and
    "query" echoes back exactly what was passed in.

    Queries /tcgav2/genes once per group (Primary Tumor, Solid Tissue
    Normal) via the cgc_sample_sample_type field. Returns a two-group
    summary (samples_count, coverage_avg, coverage_median), a note on the
    cohort-size imbalance between the groups, plus the absolute path to a
    CSV of every individual sample's coverage for further analysis/plotting
    in R. Use this to answer "is gene X expressed differently in tumor vs
    normal tissue?". Note: coverage_median can floor at an identical low
    value for both groups if the gene is not meaningfully expressed in
    either -- report coverage_avg alongside it.
    """
    return gene_tumor_vs_normal_impl(gene)


# ---------------------------------------------------------------------------
# Tool 3: junction_usage
# ---------------------------------------------------------------------------


def _classify_junction_novelty(left_annotated: str, right_annotated: str, annotated: str) -> str:
    """Classify a junction's novelty from left_annotated/right_annotated/annotated.

    left_annotated/right_annotated are "0" (that splice site is novel) or a
    comma-separated list of annotation-source codes (that site is known).
    This is richer than a flat annotated==0 check (see CLAUDE.md):
        - both sites novel            -> "both-ends-novel"
        - exactly one site novel      -> "one-end-novel"
        - both sites known, but this
          exact pairing is not        -> "novel-combination"
        - both known, pairing known   -> "known"
    """
    left_novel = left_annotated == "0"
    right_novel = right_annotated == "0"
    if left_novel and right_novel:
        return "both-ends-novel"
    if left_novel or right_novel:
        return "one-end-novel"
    if annotated == "0":
        return "novel-combination"
    return "known"


def _fetch_tcga_junctions(gene: str, sample_type_value: str, novel_only: bool) -> tuple[pd.DataFrame, list[str]]:
    """Fetch /tcgav2/snaptron (junctions) filtered to one TCGA sample-type
    group, optionally restricted to novel junctions server-side."""
    notes: list[str] = []
    filters = [("sfilter", f"{TCGA_SAMPLE_TYPE_FIELD}:{sample_type_value}")]
    if novel_only:
        filters.append(("rfilter", "annotated:0"))

    try:
        body = fetch_snaptron("tcgav2", "snaptron", gene, filters=filters)
    except RuntimeError as exc:
        notes.append(f"{sample_type_value}: query failed -- {exc}")
        return pd.DataFrame(), notes

    df = parse_tsv(body)
    if df.empty:
        notes.append(f"{sample_type_value}: no junctions returned")
    return df, notes


def _fetch_gtex_junctions(gene: str, tissue: str) -> tuple[pd.DataFrame, list[str]]:
    """Fetch /gtexv2/snaptron (junctions) filtered to one GTEx SMTS tissue."""
    notes: list[str] = []
    try:
        body = fetch_snaptron("gtexv2", "snaptron", gene, filters=[("sfilter", f"SMTS:{tissue}")])
    except RuntimeError as exc:
        notes.append(f"{tissue}: query failed -- {exc}")
        return pd.DataFrame(), notes

    df = parse_tsv(body)
    if df.empty:
        notes.append(f"{tissue}: no junctions returned")
    return df, notes


def _junction_usage_tissue(gene: str, notes: list[str]) -> dict:
    per_tissue_df: dict[str, pd.DataFrame] = {}
    for tissue in GTEX_TISSUES:
        df, tissue_notes = _fetch_gtex_junctions(gene, tissue)
        notes.extend(tissue_notes)
        per_tissue_df[tissue] = df

    long_rows: list[dict] = []
    junction_data: dict[str, dict] = {}

    for tissue, df in per_tissue_df.items():
        if df.empty:
            continue
        for _, row in df.iterrows():
            try:
                start = int(row["start"])
                end = int(row["end"])
            except (ValueError, TypeError):
                continue

            samples = parse_samples_column(row["samples"])
            if tissue == "Blood":
                samples, correction_note = _correct_blood_smts_contamination(samples)
                notes.append(f"{row['chromosome']}:{start}-{end}: {correction_note}")
            if not samples:
                continue

            coverages = [coverage for _, coverage in samples]
            samples_count = len(samples)
            coverage_avg = statistics.mean(coverages)
            coverage_median = statistics.median(coverages)

            junction_key = f"{row['chromosome']}:{start}-{end}"
            novelty_class = _classify_junction_novelty(
                row.get("left_annotated", ""), row.get("right_annotated", ""), row.get("annotated", "")
            )

            long_rows.append(
                {
                    "junction": junction_key,
                    "chromosome": row["chromosome"],
                    "start": start,
                    "end": end,
                    "strand": row.get("strand"),
                    "novelty_class": novelty_class,
                    "group": tissue,
                    "samples_count": samples_count,
                    "coverage_avg": coverage_avg,
                    "coverage_median": coverage_median,
                }
            )

            entry = junction_data.setdefault(
                junction_key,
                {
                    "junction": junction_key,
                    "chromosome": row["chromosome"],
                    "start": start,
                    "end": end,
                    "strand": row.get("strand"),
                    "novelty_class": novelty_class,
                    "tissue_summary": [],
                },
            )
            entry["tissue_summary"].append(
                {
                    "tissue": tissue,
                    "samples_count": samples_count,
                    "coverage_avg": coverage_avg,
                    "coverage_median": coverage_median,
                }
            )

    for entry in junction_data.values():
        entry["tissue_summary"].sort(key=lambda t: t["coverage_median"], reverse=True)

    junction_rows = sorted(junction_data.values(), key=lambda e: e["start"])
    ranked_by_peak_tissue = sorted(
        junction_data.values(),
        key=lambda e: e["tissue_summary"][0]["coverage_median"] if e["tissue_summary"] else 0,
        reverse=True,
    )

    csv_path = write_tidy_csv(long_rows, f"junction_usage_tissue_{gene}")
    notes.append(
        "plot_png_path is None: plot_junction_usage.R only supports tumor_vs_normal "
        "CSVs (it pivots on 'Primary Tumor'/'Solid Tissue Normal' group values) -- no "
        "plot script exists yet for tissue-mode junction usage."
    )

    return {
        "gene": gene,
        "comparison": "tissue",
        "compilation": "gtexv2",
        "total_junctions_found": len(junction_rows),
        "top_junctions_by_peak_tissue_coverage": ranked_by_peak_tissue[:25],
        "per_junction_csv_path": csv_path,
        "plot_png_path": None,
        "per_junction_csv_columns": [
            "junction",
            "chromosome",
            "start",
            "end",
            "strand",
            "novelty_class",
            "group",
            "samples_count",
            "coverage_avg",
            "coverage_median",
        ],
        "notes": notes,
        "normalization_caveat": (
            "Each junction's tissue_summary is ranked by coverage_median descending across all "
            "31 confirmed GTEx tissues (same panel and Blood/Blood Vessel sfilter word-match "
            "correction as gene_expression_across_tissues -- see GTEX_SMTS_WORD_MATCH_CAVEAT). "
            "Unlike tumor_vs_normal, no pct_of_cohort correction is applied here -- there is no "
            "established per-tissue cohort-total analogous to TCGA_COHORT_TOTALS for junctions, "
            "so samples_count is a raw count and tissue cohort sizes differ substantially (e.g. "
            "Brain n~1239 vs Kidney n~11) -- treat a tissue's samples_count as informative context, "
            "not a normalized fraction. coverage_median values are NOT library-size normalized -- "
            + NORMALIZATION_CAVEAT
        ),
    }


def junction_usage_impl(gene: str, comparison: str = "tumor_vs_normal", novel_only: bool = False) -> dict:
    """Implementation, kept separate from the @mcp.tool wrapper for direct
    testability (see gene_expression_across_tissues_impl for why).

    `comparison`: "tumor_vs_normal" (default, tcgav2) is unchanged -- ranks
    junctions by pct_of_cohort_diff between TCGA tumor and normal. "tissue"
    (gtexv2) mirrors exon_usage's tissue generalization -- each junction
    gets a nested tissue_summary list across all 31 confirmed GTEx tissues,
    ranked by coverage_median, including the same Blood/Blood Vessel
    sfilter correction tool 1 applies (that collision is a property of
    sfilter itself, not of /genes specifically -- see
    GTEX_SMTS_WORD_MATCH_CAVEAT). novel_only is only meaningful for
    tumor_vs_normal (tissue mode does not restrict server-side by
    annotated:0 -- novelty_class is still reported per junction either way)."""
    if comparison not in ("tumor_vs_normal", "tissue"):
        raise ValueError(f"Unsupported comparison '{comparison}' -- must be 'tumor_vs_normal' or 'tissue'.")

    gene = gene.strip()

    if comparison == "tissue":
        notes: list[str] = [
            "Junctions are every /snaptron row overlapping the gene's coordinate "
            "span (regions=<gene> resolves to a span, same as /genes and /exons "
            "-- see CLAUDE.md). Unlike /genes and /exons, the junction schema has "
            "no gene-name column to filter on client-side, so for genes with "
            "tightly neighboring genes a small number of junctions from that "
            "neighbor may be included."
        ]
        return _junction_usage_tissue(gene, notes)

    notes = [
        "Junctions are every /snaptron row overlapping the gene's coordinate "
        "span (regions=<gene> resolves to a span, same as /genes and /exons "
        "-- see CLAUDE.md). Unlike /genes and /exons, the junction schema has "
        "no gene-name column to filter on client-side, so for genes with "
        "tightly neighboring genes a small number of junctions from that "
        "neighbor may be included."
    ]

    per_group_df: dict[str, pd.DataFrame] = {}
    for sample_type_value in TCGA_GROUPS:
        df, group_notes = _fetch_tcga_junctions(gene, sample_type_value, novel_only)
        notes.extend(group_notes)
        per_group_df[sample_type_value] = df

    # Long-format tidy rows for the CSV: one row per junction per group it
    # appeared in.
    long_rows: list[dict] = []
    # Joined-across-groups rows, keyed by "chrom:start-end", for the ranked
    # differential-usage summary.
    joined: dict[str, dict] = {}

    group_prefix = {"Primary Tumor": "tumor", "Solid Tissue Normal": "normal"}

    for sample_type_value, df in per_group_df.items():
        if df.empty:
            continue
        cohort_total = TCGA_COHORT_TOTALS[sample_type_value]
        prefix = group_prefix[sample_type_value]

        for _, row in df.iterrows():
            try:
                samples_count = int(row["samples_count"])
                coverage_avg = float(row["coverage_avg"])
                coverage_median = float(row["coverage_median"])
                start = int(row["start"])
                end = int(row["end"])
            except (ValueError, TypeError):
                continue

            pct_of_cohort = samples_count / cohort_total
            junction_key = f"{row['chromosome']}:{start}-{end}"
            novelty_class = _classify_junction_novelty(
                row.get("left_annotated", ""), row.get("right_annotated", ""), row.get("annotated", "")
            )

            long_rows.append(
                {
                    "junction": junction_key,
                    "chromosome": row["chromosome"],
                    "start": start,
                    "end": end,
                    "strand": row.get("strand"),
                    "novelty_class": novelty_class,
                    "group": sample_type_value,
                    "samples_count": samples_count,
                    "cohort_total": cohort_total,
                    "pct_of_cohort": round(pct_of_cohort, 4),
                    "coverage_avg": coverage_avg,
                    "coverage_median": coverage_median,
                }
            )

            entry = joined.setdefault(
                junction_key,
                {
                    "junction": junction_key,
                    "chromosome": row["chromosome"],
                    "start": start,
                    "end": end,
                    "strand": row.get("strand"),
                    "novelty_class": novelty_class,
                    "tumor_samples_count": 0,
                    "tumor_pct_of_cohort": 0.0,
                    "tumor_coverage_median": None,
                    "normal_samples_count": 0,
                    "normal_pct_of_cohort": 0.0,
                    "normal_coverage_median": None,
                },
            )
            entry[f"{prefix}_samples_count"] = samples_count
            entry[f"{prefix}_pct_of_cohort"] = round(pct_of_cohort, 4)
            entry[f"{prefix}_coverage_median"] = coverage_median

    junction_rows = list(joined.values())
    for entry in junction_rows:
        entry["pct_of_cohort_diff"] = round(entry["tumor_pct_of_cohort"] - entry["normal_pct_of_cohort"], 4)
    junction_rows.sort(key=lambda r: abs(r["pct_of_cohort_diff"]), reverse=True)

    novelty_class_counts: dict[str, int] = {}
    for entry in junction_rows:
        novelty_class_counts[entry["novelty_class"]] = novelty_class_counts.get(entry["novelty_class"], 0) + 1

    csv_path = write_tidy_csv(long_rows, f"junction_usage_{gene}")

    plot_png_path = None
    if csv_path:
        plot_png_path, plot_note = run_plot_script("plot_junction_usage.R", [csv_path, gene])
        notes.append(plot_note)

    return {
        "gene": gene,
        "comparison": comparison,
        "novel_only": novel_only,
        "compilation": "tcgav2",
        "cohort_totals": TCGA_COHORT_TOTALS,
        "total_junctions_found": len(junction_rows),
        "novelty_class_counts": novelty_class_counts,
        "top_differential_junctions": junction_rows[:25],
        "per_junction_csv_path": csv_path,
        "plot_png_path": plot_png_path,
        "per_junction_csv_columns": [
            "junction",
            "chromosome",
            "start",
            "end",
            "strand",
            "novelty_class",
            "group",
            "samples_count",
            "cohort_total",
            "pct_of_cohort",
            "coverage_avg",
            "coverage_median",
        ],
        "notes": notes,
        "normalization_caveat": (
            "pct_of_cohort = a junction's samples_count divided by its group's "
            "whole-cohort total (9961 tumor / 740 normal -- see "
            "TCGA_COHORT_TOTALS), i.e. the fraction of that group's samples "
            "using this junction at all. This corrects for the 13.5x cohort "
            "size imbalance between groups, unlike a raw samples_count "
            "comparison. coverage_median values themselves are still NOT "
            "library-size normalized -- " + NORMALIZATION_CAVEAT + " "
            "Ranking by pct_of_cohort_diff alone can surface junctions that "
            "are widely observed but only ~1 read deep per sample "
            "(coverage_median 1.0) -- in a highly-expressed gene sequenced "
            "across thousands of samples, that pattern is often sequencing "
            "noise rather than a real splicing event. Check coverage_median "
            "alongside pct_of_cohort_diff before treating a top-ranked "
            "junction as biologically meaningful."
        ),
    }


@mcp.tool()
def junction_usage(gene: str, comparison: str = "tumor_vs_normal", novel_only: bool = False) -> dict:
    """Compare splice junction usage for a gene, either between TCGA tumor
    and normal (comparison="tumor_vs_normal", default, compilation tcgav2)
    or across the full 31-tissue GTEx panel (comparison="tissue",
    compilation gtexv2, same tissue list and Blood/Blood Vessel sfilter
    correction as gene_expression_across_tissues), including junction
    novelty classification either way.

    Queries /{compilation}/snaptron once per group (2 groups for
    tumor_vs_normal, 31 tissues for tissue) for the gene's coordinate span,
    keys junctions by "chrom:start-end". Each junction is classified as
    "both-ends-novel", "one-end-novel", "novel-combination" (both splice
    sites individually known, but not previously seen spliced together), or
    "known", using left_annotated/right_annotated/annotated.

    tumor_vs_normal: usage is reported as pct_of_cohort (samples_count /
    that group's whole tcgav2 cohort total: 9961 tumor, 740 normal) rather
    than a raw sample count, since raw counts would otherwise be dominated
    by the ~13.5x larger tumor cohort. Junctions are ranked by the absolute
    difference in pct_of_cohort between groups. Set novel_only=True to
    restrict server-side to junctions Snaptron has not annotated at all
    (rfilter=annotated:0) -- this still includes all three novel subtypes
    above except "known". Returns the top 25 most differentially-used
    junctions inline (top_differential_junctions).

    tissue: each junction gets a nested tissue_summary list (all 31
    tissues, ranked by coverage_median) instead of a tumor/normal diff --
    there is no established per-tissue cohort-total to normalize against,
    so samples_count here is a raw count, not pct_of_cohort. novel_only is
    not applied server-side in this mode. Returns the top 25 junctions by
    peak tissue coverage inline (top_junctions_by_peak_tissue_coverage).

    Both modes also return the absolute path to a CSV with every junction
    found (one row per junction per group it appeared in) for further
    analysis/plotting in R.
    """
    return junction_usage_impl(gene, comparison=comparison, novel_only=novel_only)


# ---------------------------------------------------------------------------
# Tool 4: stratified_comparison
# ---------------------------------------------------------------------------

# Used as the reference gene for deriving a cohort's total sample count under
# an arbitrary field:value filter (see stratified_comparison_impl). GAPDH is
# broadly expressed enough that its samples_count under any real field:value
# filter should equal that value's whole-cohort size, the same logic already
# used (and confirmed on GAPDH/ACTB/MKI67/TP53) to derive TCGA_COHORT_TOTALS.
HOUSEKEEPING_GENE_DEFAULT = "GAPDH"


# Confirmed live (registry-scan-then-verify, same pattern as the other
# fields): discriminates cleanly by TCGA project/cancer-type, using the bare
# project code with no "TCGA-" prefix (e.g. "BRCA", not "TCGA-BRCA" --
# "TCGA-BRCA" returned 0 rows when tried). GAPDH samples_count: BRCA=1256,
# LUAD=601, LUSC=555. Combining it with another sfilter (e.g. stage) via a
# second repeated sfilter param confirmed ANDs correctly (Stage I alone:
# 1403, BRCA alone: 1256, Stage I AND BRCA together: 105).
# NOTE: a `study` field also exists in the registry and looks like a
# shorter alternative, but a live sfilter query on it triggered a
# ChunkedEncodingError (connection dropped) -- same failure signature as
# the already-documented broken endpoints. Left untested further; use
# gdc_cases.project.project_id instead.
TCGA_PROJECT_FIELD = "gdc_cases.project.project_id"


def _fetch_tcga_gene_row_multi(gene: str, filters: list[tuple[str, str]], label: str) -> tuple[dict | None, list[str]]:
    """Fetch /tcgav2/genes filtered by one or more field:value sfilter
    constraints (ANDed together via repeated sfilter params -- confirmed
    live, see TCGA_PROJECT_FIELD comment above) and return the row matching
    the requested gene symbol (or None + a note explaining why not). `label`
    is only used to make notes human-readable."""
    notes: list[str] = []
    sfilters = [("sfilter", f"{field}:{value}") for field, value in filters]
    try:
        body = fetch_snaptron("tcgav2", "genes", gene, filters=sfilters)
    except RuntimeError as exc:
        notes.append(f"{label}: query failed -- {exc}")
        return None, notes

    df = parse_tsv(body)
    if df.empty:
        notes.append(f"{label}: no rows returned")
        return None, notes

    df = split_packed_gene_column(df)
    matches = df[df["gene_name"].str.upper() == gene.upper()]
    if matches.empty:
        notes.append(
            f"{label}: gene '{gene}' not found among {len(df)} features "
            f"overlapping its coordinate span"
        )
        return None, notes
    if len(matches) > 1:
        notes.append(
            f"{label}: {len(matches)} rows matched gene symbol '{gene}' "
            f"(likely multiple gene model versions) -- using the first"
        )

    return matches.iloc[0].to_dict(), notes


def stratified_comparison_impl(
    gene: str,
    field: str,
    values: list[str],
    housekeeping_gene: str = HOUSEKEEPING_GENE_DEFAULT,
    extra_field: str | None = None,
    extra_value: str | None = None,
) -> dict:
    """Implementation, kept separate from the @mcp.tool wrapper for direct
    testability (see gene_expression_across_tissues_impl for why)."""
    gene = gene.strip()
    field = field.strip()
    has_extra_constraint = bool(extra_field and extra_value)
    notes: list[str] = []
    cohort_total_cache: dict[str, int | None] = {}
    group_summary: list[dict] = []
    per_sample_rows: list[dict] = []

    def build_filters(value: str) -> list[tuple[str, str]]:
        filters = [(field, value)]
        if has_extra_constraint:
            filters.append((extra_field, extra_value))
        return filters

    def get_cohort_total(value: str) -> int | None:
        # Cached per value within this call only -- field/values (and any
        # extra_field/extra_value) are caller-specified per BUILD_PLAN's
        # instruction not to hardcode totals the way TCGA_COHORT_TOTALS does
        # for tumor/normal. Applying the SAME extra constraint to the
        # housekeeping-gene query as to the target gene is what makes
        # pct_of_cohort mean "fraction of THIS constrained cohort" rather
        # than silently falling back to the pan-cancer total.
        if value in cohort_total_cache:
            return cohort_total_cache[value]
        hk_row, hk_notes = _fetch_tcga_gene_row_multi(
            housekeeping_gene, build_filters(value), f"cohort total for {value!r}"
        )
        notes.extend(hk_notes)
        total = None
        if hk_row is not None:
            try:
                total = int(hk_row["samples_count"])
            except (ValueError, TypeError):
                notes.append(f"{value}: could not parse {housekeeping_gene} samples_count for cohort total")
        cohort_total_cache[value] = total
        return total

    for value in values:
        cohort_total = get_cohort_total(value)

        row, row_notes = _fetch_tcga_gene_row_multi(gene, build_filters(value), value)
        notes.extend(row_notes)
        if row is None:
            continue

        try:
            samples_count = int(row["samples_count"])
            coverage_avg = float(row["coverage_avg"])
            coverage_median = float(row["coverage_median"])
        except (ValueError, TypeError):
            notes.append(f"{value}: could not parse coverage fields, skipping")
            continue

        pct_of_cohort = samples_count / cohort_total if cohort_total else None

        group_summary.append(
            {
                "value": value,
                "samples_count": samples_count,
                "cohort_total": cohort_total,
                "pct_of_cohort": round(pct_of_cohort, 4) if pct_of_cohort is not None else None,
                "coverage_avg": coverage_avg,
                "coverage_median": coverage_median,
            }
        )

        for rail_id, coverage in parse_samples_column(row["samples"]):
            per_sample_rows.append({"value": value, "rail_id": rail_id, "coverage": coverage})

    csv_prefix = f"stratified_comparison_{gene}_{field}"
    if has_extra_constraint:
        csv_prefix += f"_within_{extra_value}"
    csv_prefix = csv_prefix.replace(".", "_").replace(" ", "_")
    csv_path = write_tidy_csv(per_sample_rows, csv_prefix)

    plot_png_path = None
    if csv_path:
        plot_title = f"{gene} by {field}" + (f" within {extra_value}" if has_extra_constraint else "")
        plot_png_path, plot_note = run_plot_script(
            "plot_stratified_comparison_single_scope.R", [csv_path, plot_title]
        )
        notes.append(plot_note)
    notes.append(
        "plot_stratified_comparison_single_scope.R (used for plot_png_path above) shows the "
        "values AS CALLED, in one scope -- it does not make the pan-cancer-vs-within-one-"
        "cancer-type confound point plot_stratified_comparison.R does; that script still "
        "needs two separate stratified_comparison calls (one unconstrained, one with "
        "extra_field/extra_value) passed to it directly, it is not wired into this tool."
    )

    pan_cancer_pooling_caveat = (
        "PAN-CANCER POOLING: no extra_field/extra_value constraint was supplied, so each "
        f"value of '{field}' pools every TCGA cancer type together (e.g. one value of "
        "this field means that value across ALL cancer types at once, not within a "
        "single cancer type). Since a gene's expression can vary by cancer type "
        f"independent of '{field}', and different cancer types can have different "
        f"distributions across '{field}' values, a difference seen here may be driven by "
        "which cancer types happen to be more common at each value rather than by "
        f"'{field}' itself. To rule this out, pass extra_field='{TCGA_PROJECT_FIELD}' and "
        "extra_value=<a specific project code, e.g. 'BRCA'> to compare within one cancer "
        "type."
        if not has_extra_constraint
        else (
            f"Constrained to {extra_field}={extra_value}, so this comparison is WITHIN one "
            "TCGA cancer type/project and is not confounded by cancer-type composition "
            f"the way an unconstrained '{field}' comparison would be."
        )
    )

    return {
        "gene": gene,
        "field": field,
        "values": values,
        "extra_constraint": {"field": extra_field, "value": extra_value} if has_extra_constraint else None,
        "compilation": "tcgav2",
        "housekeeping_gene_used_for_cohort_totals": housekeeping_gene,
        "group_summary": group_summary,
        "per_sample_csv_path": csv_path,
        "per_sample_csv_columns": ["value", "rail_id", "coverage"],
        "plot_png_path": plot_png_path,
        "notes": notes,
        "cancer_type_confound_caveat": pan_cancer_pooling_caveat,
        "field_nomenclature_caveat": (
            f"Not all cancer types share the same clinical-field nomenclature or even use "
            f"the same field at all -- '{field}' values that work for one cancer type can "
            f"return n=0 for every value in another. Confirmed live: TCGA-PRAD (prostate) "
            f"does not populate cgc_case_pathologic_stage at all (every 'Stage I'..'Stage IV' "
            f"value, including substages, returned 0 samples) -- prostate cancer in TCGA is "
            f"staged via cgc_case_pathologic_t (T-category) instead, not an overall AJCC "
            f"stage grouping. Always check each value's samples_count for the SPECIFIC "
            f"cancer type before trusting a result, per CLAUDE.md's nomenclature note -- "
            f"don't assume a field/value pair that works for one cancer type carries over "
            f"to another."
        ),
        "normalization_caveat": (
            f"cohort_total for each value is derived live in this call by querying "
            f"{housekeeping_gene} filtered to that value (and to extra_field/extra_value, "
            "if given) and reading its samples_count (the same method used to determine "
            "the hardcoded tumor/normal totals used elsewhere in this server) -- it is "
            "NOT hardcoded, since field and values here are caller-specified and may not "
            "be the tumor/normal split. pct_of_cohort = samples_count / cohort_total. "
            "coverage_avg/coverage_median are NOT library-size normalized -- " + NORMALIZATION_CAVEAT
        ),
    }


@mcp.tool()
def stratified_comparison(
    gene: str,
    field: str,
    values: list[str],
    housekeeping_gene: str = HOUSEKEEPING_GENE_DEFAULT,
    extra_field: str | None = None,
    extra_value: str | None = None,
) -> dict:
    """Compare a gene's TCGA expression across an arbitrary set of values of
    a clinical/pathology field (e.g. tumor stage, histology), generalizing
    gene_tumor_vs_normal beyond the fixed tumor/normal split.

    Queries /tcgav2/genes once per value in `values`, filtered via
    sfilter=<field>:<value>. `field` must be an exact tcgav2 registry field
    name (see /snaptron/registry). Do not guess a field name from how it
    sounds -- verify it actually discriminates with a real query first; some
    registry fields are narrow per-cancer-type technical fields rather than
    general clinical variables (see CLAUDE.md's SMTSPAX note for why this
    matters). cgc_case_pathologic_stage (values like "Stage I", "Stage IA",
    "Stage IV") is confirmed working -- but NOT for every cancer type: e.g.
    TCGA-PRAD (prostate) does not populate this field at all (every stage
    value returns 0 samples), since prostate cancer is conventionally
    staged via cgc_case_pathologic_t (T-category: T2a, T2b, T2c, T3a, T3b,
    T4, also confirmed working) rather than an overall AJCC stage grouping.
    Not all cancer types share the same clinical-field nomenclature --
    check each value's samples_count for the SPECIFIC cancer type you're
    using before trusting a result; don't assume a field/value pair that
    works for one cancer type carries over to another (see this result's
    own field_nomenclature_caveat, and CLAUDE.md).

    IMPORTANT -- cancer-type confound: without extra_field/extra_value, each
    value of `field` pools every TCGA cancer type together (e.g. "Stage I"
    means Stage I across ALL cancer types at once). A gene's expression can
    vary by cancer type independent of `field`, and cancer types differ in
    their distribution across `field` values, so a pan-cancer difference can
    be driven by cancer-type composition rather than `field` itself. Pass
    extra_field="gdc_cases.project.project_id" (confirmed working; use the
    bare project code, e.g. extra_value="BRCA", NOT "TCGA-BRCA") to restrict
    the whole comparison to one cancer type and remove this confound. This
    extra constraint is applied to BOTH the target gene query and the
    housekeeping-gene cohort-total query, so pct_of_cohort still means
    "fraction of this constrained cohort," not the pan-cancer total.

    Unlike gene_tumor_vs_normal, cohort totals per value are NOT hardcoded:
    each value's whole-cohort sample count is derived live, in this same
    call, by querying `housekeeping_gene` (default GAPDH) filtered to that
    value (and extra_field/extra_value, if given) and reading its
    samples_count -- the same method originally used to determine the
    tumor/normal totals. This is cached per value within one call but never
    across calls, since field/values can be anything the caller passes.

    Returns a per-value summary (samples_count, cohort_total, pct_of_cohort,
    coverage_avg, coverage_median), the absolute path to a per-sample CSV,
    and (via plot_stratified_comparison_single_scope.R) an auto-generated
    boxplot+jitter PNG of the values as called, ordered by the sequence
    `values` was passed in (not sorted by magnitude -- appropriate when the
    values have a natural progression, e.g. T-stage). The result also
    carries a cancer_type_confound_caveat field describing whether this
    call is pan-cancer pooled or constrained to one cancer type, and a
    field_nomenclature_caveat noting that field/value support varies by
    cancer type.
    """
    return stratified_comparison_impl(
        gene, field, values, housekeeping_gene=housekeeping_gene, extra_field=extra_field, extra_value=extra_value
    )


# ---------------------------------------------------------------------------
# Tool 5: exon_usage
# ---------------------------------------------------------------------------

# /exons returns the SAME 18-field row shape as /genes (confirmed live on
# TP53 and MKI67), but its header is malformed: the live response declares
# only 1 field ("DataSource:Type") while every data row has 18 tab-separated
# fields (confirmed independent of this code via raw curl). Header-based
# parsing (pandas.read_csv, what parse_tsv/split_packed_gene_column use for
# /genes and /snaptron) silently misreads this -- pandas treats the extra
# data columns as an implicit index, collapsing everything into 1-2 garbage
# columns. This is why an earlier recon pass concluded "/exons sfilter
# doesn't work": it was checking row count (which never changes under a
# filter on /genes either -- only the per-row aggregate stats do) using a
# parser that was silently broken on this endpoint regardless. Re-tested
# properly: sfilter DOES filter /exons correctly (samples_count tracked the
# known 9961/740 tumor/normal cohort totals exactly, and matched /genes'
# samples_count exactly for the same gene under gdc_cases.project.project_id
# filtering). Fix: parse /exons positionally (split on tab, discard the
# broken header line, use fixed field indices) instead of via pandas.
EXONS_EXPECTED_FIELD_COUNT = 18
# 0 DataSource:Type, 1 snaptron_id, 2 chromosome, 3 start, 4 end, 5 length,
# 6 strand, 7-9 blank/NA, 10 blank, 11 packed id, 12 samples,
# 13 samples_count, 14 coverage_sum, 15 coverage_avg, 16 coverage_median,
# 17 compilation_id -- same layout as /genes, confirmed by position.


def _fetch_tcga_exons_rows(gene: str, sample_type_value: str) -> tuple[list[list[str]], list[str]]:
    """Fetch /tcgav2/exons filtered to one TCGA sample-type group and return
    positionally-parsed rows (see EXONS_EXPECTED_FIELD_COUNT comment above
    for why positional, not pandas/header-based)."""
    notes: list[str] = []
    try:
        body = fetch_snaptron(
            "tcgav2", "exons", gene, filters=[("sfilter", f"{TCGA_SAMPLE_TYPE_FIELD}:{sample_type_value}")]
        )
    except RuntimeError as exc:
        notes.append(f"{sample_type_value}: query failed -- {exc}")
        return [], notes

    lines = [line for line in body.split("\n") if line.strip()]
    if not lines:
        notes.append(f"{sample_type_value}: no rows returned")
        return [], notes

    data_lines = lines[1:]  # line 0 is the known-broken /exons header -- always discard
    rows = [line.split("\t") for line in data_lines]
    bad_rows = [r for r in rows if len(r) != EXONS_EXPECTED_FIELD_COUNT]
    if bad_rows:
        notes.append(
            f"{sample_type_value}: {len(bad_rows)} of {len(rows)} /exons rows did not have "
            f"the expected {EXONS_EXPECTED_FIELD_COUNT} fields -- dropped rather than guessed at"
        )
    rows = [r for r in rows if len(r) == EXONS_EXPECTED_FIELD_COUNT]
    return rows, notes


def _split_packed_exon_id(packed: str) -> tuple[str, str, str, str]:
    """Split /exons' packed identifier column. Unlike /genes' 4-field
    gene_id:gene_name:gene_type:bp_length, /exons repeats gene_id: 5 fields,
    gene_id:gene_id:gene_name:gene_type:bp_length (confirmed live on TP53 and
    MKI67). Returns (gene_id, gene_name, gene_type, bp_length), gene_id
    deduplicated."""
    parts = packed.split(":", 4)
    gene_id, _gene_id_repeated, gene_name, gene_type, bp_length = parts
    return gene_id, gene_name, gene_type, bp_length


REGION_PATTERN = re.compile(r"^chr[0-9A-Za-z]+:\d+-\d+$")


def _looks_like_region(value: str) -> bool:
    """True if `value` is a raw chr:start-end region rather than a gene
    symbol. Snaptron's regions= param accepts both forms; this server needs
    to know which one it got, because client-side gene-name filtering (to
    drop neighboring genes' features from a symbol query -- see CLAUDE.md)
    is meaningless and actively wrong to apply to a literal region query."""
    return bool(REGION_PATTERN.match(value))


def _fetch_gtex_exons_rows(gene: str, tissue: str) -> tuple[list[list[str]], list[str]]:
    """Fetch /gtexv2/exons filtered to one GTEx SMTS tissue and return
    positionally-parsed rows -- same broken-header workaround as
    _fetch_tcga_exons_rows (see EXONS_EXPECTED_FIELD_COUNT comment)."""
    notes: list[str] = []
    try:
        body = fetch_snaptron("gtexv2", "exons", gene, filters=[("sfilter", f"SMTS:{tissue}")])
    except RuntimeError as exc:
        notes.append(f"{tissue}: query failed -- {exc}")
        return [], notes

    lines = [line for line in body.split("\n") if line.strip()]
    if not lines:
        notes.append(f"{tissue}: no rows returned")
        return [], notes

    data_lines = lines[1:]  # line 0 is the known-broken /exons header -- always discard
    rows = [line.split("\t") for line in data_lines]
    bad_rows = [r for r in rows if len(r) != EXONS_EXPECTED_FIELD_COUNT]
    if bad_rows:
        notes.append(
            f"{tissue}: {len(bad_rows)} of {len(rows)} /exons rows did not have "
            f"the expected {EXONS_EXPECTED_FIELD_COUNT} fields -- dropped rather than guessed at"
        )
    rows = [r for r in rows if len(r) == EXONS_EXPECTED_FIELD_COUNT]
    return rows, notes


def _region_note(is_region_query: bool) -> str:
    base = (
        "Exons are every /exons row overlapping the gene's coordinate span "
        "whose packed gene_name matches exactly (regions=<gene> resolves to "
        "a span and can include neighboring genes' exons -- see CLAUDE.md); "
        "those are filtered out client-side the same way /genes handles it."
    )
    if is_region_query:
        return base + (
            " (A raw chr:start-end region was passed instead of a gene symbol -- "
            "gene-name filtering is skipped, every exon in the region is kept.)"
        )
    return base


def _exon_usage_tumor_vs_normal(gene: str, is_region_query: bool, notes: list[str]) -> dict:
    per_group_rows: dict[str, list[list[str]]] = {}
    for sample_type_value in TCGA_GROUPS:
        rows, group_notes = _fetch_tcga_exons_rows(gene, sample_type_value)
        notes.extend(group_notes)
        per_group_rows[sample_type_value] = rows

    group_prefix = {"Primary Tumor": "tumor", "Solid Tissue Normal": "normal"}
    exon_data: dict[str, dict] = {}
    per_sample_rows: list[dict] = []
    genes_found: set[str] = set()

    for sample_type_value, rows in per_group_rows.items():
        prefix = group_prefix[sample_type_value]
        for r in rows:
            chromosome, start, end, strand = r[2], r[3], r[4], r[6]
            packed_id, samples_str = r[11], r[12]
            gene_id, gene_name, gene_type, _bp_length = _split_packed_exon_id(packed_id)
            if not is_region_query and gene_name.upper() != gene.upper():
                continue  # neighboring gene's exon overlapping the same span -- skip
            genes_found.add(gene_name)

            try:
                samples_count = int(r[13])
                coverage_avg = float(r[15])
                coverage_median = float(r[16])
                start_i, end_i = int(start), int(end)
            except (ValueError, TypeError):
                notes.append(f"{sample_type_value}: could not parse a row for exon {chromosome}:{start}-{end}, skipping")
                continue

            exon_key = f"{chromosome}:{start_i}-{end_i}"
            entry = exon_data.setdefault(
                exon_key,
                {
                    "exon": exon_key,
                    "chromosome": chromosome,
                    "start": start_i,
                    "end": end_i,
                    "strand": strand,
                    "gene_id": gene_id,
                    "gene_name": gene_name,
                    "gene_type": gene_type,
                    "tumor_samples_count": 0,
                    "tumor_coverage_avg": None,
                    "tumor_coverage_median": None,
                    "normal_samples_count": 0,
                    "normal_coverage_avg": None,
                    "normal_coverage_median": None,
                },
            )
            entry[f"{prefix}_samples_count"] = samples_count
            entry[f"{prefix}_coverage_avg"] = coverage_avg
            entry[f"{prefix}_coverage_median"] = coverage_median

            for rail_id, coverage in parse_samples_column(samples_str):
                per_sample_rows.append({"exon": exon_key, "group": sample_type_value, "rail_id": rail_id, "coverage": coverage})

    exon_summary = sorted(exon_data.values(), key=lambda e: e["start"])

    cohort_note = None
    if all(per_group_rows.get(sample_type_value) for sample_type_value in TCGA_GROUPS):
        cohort_note = (
            f"Cohort sizes are imbalanced: ~{TCGA_COHORT_TOTALS['Primary Tumor']} tumor vs "
            f"~{TCGA_COHORT_TOTALS['Solid Tissue Normal']} normal samples in tcgav2 overall "
            f"(a given exon's own samples_count is usually close to but not always identical "
            f"to these whole-cohort totals, since not every sample has coverage at every "
            f"exon). This is why coverage_avg/coverage_median (per-sample) are reported "
            f"instead of coverage_sum."
        )

    resolved_gene = ", ".join(sorted(genes_found)) if is_region_query else gene
    if is_region_query and not genes_found:
        notes.append(f"{gene}: no annotated gene found overlapping this region")

    csv_filename_gene = re.sub(r"[^A-Za-z0-9._-]", "_", gene)  # region strings contain ':' -- not valid in Windows filenames
    csv_path = write_tidy_csv(per_sample_rows, f"exon_usage_{csv_filename_gene}")

    plot_png_path = None
    if csv_path:
        plot_png_path, plot_note = run_plot_script("plot_exon_usage.R", [csv_path, csv_filename_gene])
        notes.append(plot_note)

    return {
        "gene": resolved_gene or None,
        "query": gene,
        "comparison": "tumor_vs_normal",
        "compilation": "tcgav2",
        "exon_summary": exon_summary,
        "cohort_size_note": cohort_note,
        "per_sample_csv_path": csv_path,
        "per_sample_csv_columns": ["exon", "group", "rail_id", "coverage"],
        "plot_png_path": plot_png_path,
        "notes": notes,
        "normalization_caveat": NORMALIZATION_CAVEAT,
    }


def _exon_usage_tissue(gene: str, is_region_query: bool, notes: list[str]) -> dict:
    per_tissue_rows: dict[str, list[list[str]]] = {}
    for tissue in GTEX_TISSUES:
        rows, tissue_notes = _fetch_gtex_exons_rows(gene, tissue)
        notes.extend(tissue_notes)
        per_tissue_rows[tissue] = rows

    exon_data: dict[str, dict] = {}
    per_sample_rows: list[dict] = []
    genes_found: set[str] = set()

    for tissue, rows in per_tissue_rows.items():
        for r in rows:
            chromosome, start, end, strand = r[2], r[3], r[4], r[6]
            packed_id, samples_str = r[11], r[12]
            gene_id, gene_name, gene_type, _bp_length = _split_packed_exon_id(packed_id)
            if not is_region_query and gene_name.upper() != gene.upper():
                continue  # neighboring gene's exon overlapping the same span -- skip
            genes_found.add(gene_name)

            try:
                start_i, end_i = int(start), int(end)
            except (ValueError, TypeError):
                notes.append(f"{tissue}: could not parse a row for exon {chromosome}:{start}-{end}, skipping")
                continue

            samples = parse_samples_column(samples_str)
            if tissue == "Blood":
                samples, correction_note = _correct_blood_smts_contamination(samples)
                notes.append(f"{chromosome}:{start_i}-{end_i}: {correction_note}")
            if not samples:
                continue

            coverages = [coverage for _, coverage in samples]
            samples_count = len(samples)
            coverage_avg = statistics.mean(coverages)
            coverage_median = statistics.median(coverages)

            exon_key = f"{chromosome}:{start_i}-{end_i}"
            entry = exon_data.setdefault(
                exon_key,
                {
                    "exon": exon_key,
                    "chromosome": chromosome,
                    "start": start_i,
                    "end": end_i,
                    "strand": strand,
                    "gene_id": gene_id,
                    "gene_name": gene_name,
                    "gene_type": gene_type,
                    "tissue_summary": [],
                },
            )
            entry["tissue_summary"].append(
                {
                    "tissue": tissue,
                    "samples_count": samples_count,
                    "coverage_avg": coverage_avg,
                    "coverage_median": coverage_median,
                }
            )

            for rail_id, coverage in samples:
                per_sample_rows.append({"exon": exon_key, "group": tissue, "rail_id": rail_id, "coverage": coverage})

    for entry in exon_data.values():
        entry["tissue_summary"].sort(key=lambda t: t["coverage_median"], reverse=True)

    exon_summary = sorted(exon_data.values(), key=lambda e: e["start"])

    resolved_gene = ", ".join(sorted(genes_found)) if is_region_query else gene
    if is_region_query and not genes_found:
        notes.append(f"{gene}: no annotated gene found overlapping this region")

    csv_filename_gene = re.sub(r"[^A-Za-z0-9._-]", "_", gene)
    csv_path = write_tidy_csv(per_sample_rows, f"exon_usage_tissue_{csv_filename_gene}")
    notes.append(
        "plot_png_path is None: plot_exon_usage.R only supports tumor_vs_normal "
        "CSVs (it pivots on 'Primary Tumor'/'Solid Tissue Normal' group values) -- no "
        "plot script exists yet for tissue-mode exon usage."
    )

    return {
        "gene": resolved_gene or None,
        "query": gene,
        "comparison": "tissue",
        "compilation": "gtexv2",
        "exon_summary": exon_summary,
        "per_sample_csv_path": csv_path,
        "per_sample_csv_columns": ["exon", "group", "rail_id", "coverage"],
        "plot_png_path": None,
        "notes": notes,
        "normalization_caveat": NORMALIZATION_CAVEAT,
    }


def exon_usage_impl(gene: str, comparison: str = "tumor_vs_normal") -> dict:
    """Implementation, kept separate from the @mcp.tool wrapper for direct
    testability (see gene_expression_across_tissues_impl for why).

    `gene` accepts either a gene symbol (e.g. "MKI67") or a raw region (e.g.
    "chr10:128096659-128099255") -- confirmed live that Snaptron's regions=
    resolves both. When a raw region is given, client-side gene-symbol
    filtering is skipped (there is no symbol to match), and the returned
    "gene" field reports the actual gene symbol(s) found among the exons in
    that region, not the region string itself.

    `comparison`: "tumor_vs_normal" (default, tcgav2) mirrors
    gene_tumor_vs_normal -- one entry per exon with tumor_*/normal_* fields.
    "tissue" (gtexv2) mirrors gene_expression_across_tissues -- one entry
    per exon with a nested tissue_summary list (all 31 confirmed GTEx
    tissues, ranked by coverage_median descending), including the same
    Blood/Blood Vessel sfilter word-match correction tool 1 applies (see
    GTEX_SMTS_WORD_MATCH_CAVEAT) -- that collision is a property of sfilter
    itself, not of /genes specifically, so it applies here too."""
    if comparison not in ("tumor_vs_normal", "tissue"):
        raise ValueError(f"Unsupported comparison '{comparison}' -- must be 'tumor_vs_normal' or 'tissue'.")

    gene = gene.strip()
    is_region_query = _looks_like_region(gene)
    notes: list[str] = [
        _region_note(is_region_query),
        "This is exon-level COVERAGE comparison (per-exon samples_count / "
        "coverage_avg / coverage_median), NOT DEXSeq-grade differential exon "
        "usage (relative inclusion rate accounting for overall gene "
        "expression change) -- see BUILD_PLAN.md.",
    ]

    if comparison == "tumor_vs_normal":
        return _exon_usage_tumor_vs_normal(gene, is_region_query, notes)
    return _exon_usage_tissue(gene, is_region_query, notes)


@mcp.tool()
def exon_usage(gene: str, comparison: str = "tumor_vs_normal") -> dict:
    """Compare a gene's per-exon coverage, either between TCGA tumor and
    normal (comparison="tumor_vs_normal", default, compilation tcgav2) or
    across the full 31-tissue GTEx panel (comparison="tissue", compilation
    gtexv2, same tissue list and Blood-collision correction as
    gene_expression_across_tissues).

    `gene` accepts a gene symbol (e.g. "MKI67") OR a raw region (e.g.
    "chr10:128096659-128099255", format chr:start-end) -- confirmed live
    both work. With a region, the result's "gene" field reports the actual
    gene symbol(s) found overlapping it (not the region string), and
    "query" echoes back exactly what was passed in.

    Queries /{compilation}/exons once per group (2 groups for
    tumor_vs_normal, 31 tissues for tissue), parsed POSITIONALLY rather
    than via header (the /exons endpoint returns a malformed header -- one
    declared field vs 18 actual data fields per row -- confirmed live; see
    the comment above EXONS_EXPECTED_FIELD_COUNT). Returns one entry per
    exon (ordered by genomic position): tumor_vs_normal gives each exon
    flat tumor_*/normal_* fields; tissue gives each exon a nested
    tissue_summary list ranked by coverage_median. Also returns the
    absolute path to a per-sample CSV (columns: exon, group, rail_id,
    coverage) for further analysis/plotting in R.

    This is exon-level COVERAGE comparison, not DEXSeq-grade differential
    exon usage (which would need relative-inclusion-rate modeling
    accounting for overall gene expression change) -- label any output
    accordingly.
    """
    return exon_usage_impl(gene, comparison=comparison)


# ---------------------------------------------------------------------------
# Tool 6: generate_gene_dossier
# ---------------------------------------------------------------------------

# Thin wrapper around html_gene_dossier.py's build_dossier() -- that script
# is NOT rewritten or modified here, just exposed through call_tool the same
# way the other 5 tools are. html_gene_dossier.py does `import
# snaptron_server as s` itself; when this file is run as __main__ (the
# normal way an MCP server starts), sys.modules only has this module under
# the key "__main__", not "snaptron_server", so that import would otherwise
# load a SECOND, separate copy of this whole file (a second MCPServer
# instance, harmless but wasteful). Registering this module under its own
# name first makes html_gene_dossier's import resolve to the already-running
# instance instead.
import sys as _sys  # noqa: E402

_sys.modules.setdefault("snaptron_server", _sys.modules[__name__])


@mcp.tool()
def generate_gene_dossier(gene: str) -> dict:
    """Generate a full multi-section HTML gene dossier: GTEx tissue ranking,
    TCGA pan-cancer tumor/normal, per-cancer-type tumor/normal across all 33
    TCGA project codes, an auto-generated narrative summary, and exon-level
    detail for the cancer types that turn out well-powered. Wraps
    html_gene_dossier.py's build_dossier() unchanged.

    *** TAKES ROUGHLY 10 MINUTES *** -- confirmed live on HLX: 596.9s end to
    end, dominated by ~130 sequential API queries in the per-cancer-type
    sweep (33 project codes). THE CALLING AGENT MUST WARN THE USER AND GET
    THEIR GO-AHEAD BEFORE CALLING THIS TOOL -- do not launch it on an
    offhand mention of a gene the way a quick lookup tool would be called;
    treat it like any other slow, resource-consuming action that deserves
    confirmation first. Progress is saved incrementally to a scratch file,
    so an interrupted run resumes rather than restarting from scratch on a
    retry -- but the first attempt still has to pay the full ~10 minutes
    somewhere.

    Returns the absolute path to the generated self-contained HTML file
    (Plotly embedded inline, opens with no internet connection) and the
    measured generation time in seconds.
    """
    import time as _time

    import html_gene_dossier  # imported here, not at module level, since it in turn imports this module

    t0 = _time.time()
    out_path = html_gene_dossier.build_dossier(gene)
    elapsed = _time.time() - t0

    return {
        "gene": gene.strip(),
        "dossier_html_path": str(out_path.resolve()),
        "generation_time_seconds": round(elapsed, 1),
    }


# ---------------------------------------------------------------------------
# Fallback / status tool
# ---------------------------------------------------------------------------


@mcp.tool()
def list_available_analyses() -> dict:
    """List the analysis tools this server exposes and their build status."""
    return {
        "tools": [
            {
                "name": "gene_expression_across_tissues",
                "status": "available",
                "description": "Rank a gene's expression across a fixed panel of GTEx tissues.",
            },
            {
                "name": "gene_tumor_vs_normal",
                "status": "available",
                "description": "Compare a gene's expression in TCGA tumor vs normal samples.",
            },
            {
                "name": "junction_usage",
                "status": "available",
                "description": "Compare splice junction usage (incl. novel junctions) between tumor and normal.",
            },
            {
                "name": "stratified_comparison",
                "status": "available",
                "description": "Compare a gene's TCGA expression across arbitrary values of a clinical field (e.g. tumor stage).",
            },
            {
                "name": "exon_usage",
                "status": "available",
                "description": "Compare a gene's per-exon coverage in TCGA tumor vs normal samples (coverage comparison, not DEXSeq-grade differential usage).",
            },
            {
                "name": "generate_gene_dossier",
                "status": "available",
                "description": "Full multi-section HTML gene dossier (tissue, pan-cancer, per-cancer-type, narrative, exon detail). TAKES ~10 MINUTES -- warn the user before calling.",
            },
        ]
    }


if __name__ == "__main__":
    mcp.run()
