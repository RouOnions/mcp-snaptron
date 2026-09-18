"""
Gene dossier generator: one gene name in, one self-contained interactive
HTML report out. Mirrors the scope of the ad-hoc HLX report from an earlier
session (see docs/plot-notes.md's "HLX per-cancer-type" section) as a proper,
reusable script instead of a one-off exploration.

Sections, in order:
  1. GTEx tissue ranking (gene_expression_across_tissues, all 31 tissues)
  2. TCGA pan-cancer tumor vs normal (gene_tumor_vs_normal)
  3. Per-cancer-type tumor vs normal (stratified_comparison looped over all
     33 TCGA project codes), flagging which have a usable matched-normal n
     (>=10) and which don't
  4. An auto-generated plain-language narrative from the numbers above
  5. Exon-level detail (reusing exon_usage's building blocks, extended with
     a project_id constraint the shared tool doesn't have) for ONLY the
     cancer types flagged well-powered in section 3
  6. GTEx exon usage across all 31 tissues (exon_usage(comparison="tissue"))
     -- a heatmap (exon x tissue, color = coverage_median), not a multi-line
     chart, since N exons x 31 tissues would be too many overlapping lines
  7. Junction usage, both TCGA tumor_vs_normal (pan-cancer) and GTEx tissue
     mode (junction_usage for both), filtered to "well-powered" junctions
     (coverage_median >= MIN_NORMAL_N, same threshold and rationale as
     section 3/5's cancer-type selection -- see docs/plot-notes.md's junction
     min_coverage_median finding) so the charts show real signal, not the
     ~1-read noise floor that dominates an unfiltered junction ranking

This script calls the MCP tool implementation functions directly (imports
snaptron_server.py as a plain module) -- the same already-tested functions
used throughout this project, not a reimplementation. It does not modify
snaptron_server.py or any plot_*.R script.

Usage:
    python html_gene_dossier.py <gene_name>

Long-running: section 3 alone is ~130 sequential API queries (33 project
codes x 2 sample-type values x 2 queries per value -- one for the
housekeeping-gene cohort total, one for the actual gene). Progress is
saved incrementally to a scratch JSON file after every project code, so an
interrupt or crash does not lose completed queries -- rerunning the script
resumes from where it left off. The scratch directory is a plain
subdirectory of this script's own location (Path(__file__).parent /
"scratch"), NOT /tmp -- that assumption broke on this Windows setup in an
earlier session; a relative path next to the script works on any OS,
mirroring the same pattern snaptron_server.py already uses for OUTPUT_DIR.

Caching across reruns: sections 3 and 5 already skip re-fetching on a
rerun for the same gene (their scratch JSON persists after a completed run
and is not cleared, so a full rerun's loops see every code already
"completed" and do nothing). Sections 1, 2, 6, and 7 had NO such caching
-- checked live, confirmed they always re-hit the API -- so
_load_cached_csv() was added: each of those sections first looks in
outputs/ for a CSV from ANY earlier run (this gene + this tool +
comparison mode), and if one exists, re-aggregates the summary numbers
from it via pandas instead of re-querying Snaptron. Not day-scoped -- any
prior matching CSV counts, not just "from earlier today" -- since the
point is avoiding a redundant identical query, not enforcing freshness.
"""

from __future__ import annotations

import json
import re
import sys
import time
from pathlib import Path

import math

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

sys.path.insert(0, str(Path(__file__).parent))
import snaptron_server as s  # noqa: E402

SCRATCH_DIR = Path(__file__).parent / "scratch"
SCRATCH_DIR.mkdir(exist_ok=True)

REPORTS_DIR = Path(__file__).parent / "reports"
REPORTS_DIR.mkdir(exist_ok=True)

# Standard 33 TCGA pan-cancer-atlas project codes.
TCGA_PROJECT_CODES = [
    "ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM",
    "HNSC", "KICH", "KIRC", "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC",
    "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD",
    "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM",
]

# Same threshold used in the earlier ad-hoc HLX per-cancer-type analysis
# (docs/plot-notes.md: "16 of 33 types kept, requiring >=10 matched normal samples").
MIN_NORMAL_N = 10

CHART_COLOR = "#2c7fb8"
UP_COLOR = "#d95f02"    # tumor higher than normal
DOWN_COLOR = "#1f78b4"  # normal higher than tumor


def _load_cached_csv(prefix: str, gene: str) -> Path | None:
    """Most recent outputs/ CSV matching this EXACT tool+comparison+gene
    combination from ANY earlier run, or None.

    `prefix` must be the literal filename prefix write_tidy_csv used for
    this tool+comparison (e.g. "junction_usage" for tumor_vs_normal vs
    "junction_usage_tissue" for tissue mode -- comparison mode is baked in
    as part of the literal prefix string, not a separate field).

    Matched with an ANCHORED regex, not Path.glob's "*_{gene}_*.csv" --
    glob's wildcard cannot tell where one tool+comparison's own infix ends
    and `gene` begins, so prefix="junction_usage" with gene="tissue" would
    otherwise glob-match "junction_usage_tissue_MKI67_<ts>.csv" (a
    DIFFERENT gene's tissue-mode file), silently serving tissue data to a
    tumor_vs_normal caller. The anchored pattern requires `gene` to appear
    immediately after `prefix_`, and the remainder to be exactly
    write_tidy_csv's own timestamp format (YYYYMMDD_HHMMSS) -- no room for
    an unexpected infix like "_tissue_" to sneak in between prefix and gene
    either way, confirmed by test (see the dossier caching notes)."""
    pattern = re.compile(rf"^{re.escape(prefix)}_{re.escape(gene)}_\d{{8}}_\d{{6}}\.csv$")
    candidates = sorted(p for p in s.OUTPUT_DIR.glob("*.csv") if pattern.match(p.name))
    return candidates[-1] if candidates else None


# ---------------------------------------------------------------------------
# Section 1: GTEx tissue ranking
# ---------------------------------------------------------------------------


def section1_tissue(gene: str) -> dict:
    cached = _load_cached_csv("gene_expression_across_tissues", gene)
    if cached:
        print(f"[1/7] GTEx tissue expression -- reusing cached CSV {cached.name}", file=sys.stderr)
        per_sample = pd.read_csv(cached)
        tissue_summary = (
            per_sample.groupby("tissue")["coverage"]
            .agg(samples_count="count", coverage_median="median", coverage_avg="mean")
            .reset_index()
            .to_dict("records")
        )
        return {"gene": gene, "tissue_summary": tissue_summary, "per_sample_csv_path": str(cached)}

    print(f"[1/7] GTEx tissue expression ({len(s.GTEX_TISSUES)} tissues)...", file=sys.stderr)
    return s.gene_expression_across_tissues_impl(gene)


def section1_figure(result: dict, gene: str) -> go.Figure:
    tissue_summary = sorted(result["tissue_summary"], key=lambda t: t["coverage_median"], reverse=True)
    tissue_order = [t["tissue"] for t in tissue_summary]
    per_sample = pd.read_csv(result["per_sample_csv_path"])

    fig = go.Figure()
    fig.add_trace(
        go.Box(
            x=per_sample["coverage"],
            y=per_sample["tissue"],
            orientation="h",
            boxpoints="all",
            jitter=0.5,
            pointpos=0,
            marker=dict(color=CHART_COLOR, size=3, opacity=0.3),
            line=dict(color=CHART_COLOR),
            fillcolor="rgba(44,127,184,0.25)",
            hovertemplate="%{y}<br>coverage=%{x}<extra></extra>",
            showlegend=False,
        )
    )
    fig.update_layout(
        title=f"{gene} expression across GTEx tissues",
        xaxis_title="per-sample coverage (log10 scale)",
        xaxis_type="log",
        yaxis=dict(categoryorder="array", categoryarray=list(reversed(tissue_order))),
        template="plotly_white",
        height=100 + 26 * len(tissue_order),
        margin=dict(l=140),
    )
    return fig


# ---------------------------------------------------------------------------
# Section 2: TCGA pan-cancer tumor vs normal
# ---------------------------------------------------------------------------


def section2_pan_cancer(gene: str) -> dict:
    cached = _load_cached_csv("gene_tumor_vs_normal", gene)
    if cached:
        print(f"[2/7] TCGA pan-cancer tumor vs normal -- reusing cached CSV {cached.name}", file=sys.stderr)
        per_sample = pd.read_csv(cached)
        group_summary = (
            per_sample.groupby("group")["coverage"]
            .agg(samples_count="count", coverage_median="median", coverage_avg="mean")
            .reset_index()
            .to_dict("records")
        )
        return {"gene": gene, "group_summary": group_summary, "per_sample_csv_path": str(cached)}

    print("[2/7] TCGA pan-cancer tumor vs normal...", file=sys.stderr)
    return s.gene_tumor_vs_normal_impl(gene)


def section2_figure(result: dict, gene: str) -> go.Figure:
    per_sample = pd.read_csv(result["per_sample_csv_path"])
    group_order = [g["group"] for g in sorted(result["group_summary"], key=lambda g: g["coverage_median"], reverse=True)]

    fig = go.Figure()
    fig.add_trace(
        go.Box(
            x=per_sample["coverage"],
            y=per_sample["group"],
            orientation="h",
            boxpoints="all",
            jitter=0.5,
            pointpos=0,
            marker=dict(color=UP_COLOR, size=3, opacity=0.2),
            line=dict(color=UP_COLOR),
            fillcolor="rgba(217,95,2,0.2)",
            hovertemplate="%{y}<br>coverage=%{x}<extra></extra>",
            showlegend=False,
        )
    )
    fig.update_layout(
        title=f"{gene} tumor vs normal (TCGA pan-cancer)",
        xaxis_title="per-sample coverage (log10 scale)",
        xaxis_type="log",
        yaxis=dict(categoryorder="array", categoryarray=list(reversed(group_order))),
        template="plotly_white",
        height=260,
        margin=dict(l=160),
    )
    return fig


# ---------------------------------------------------------------------------
# Section 3: per-cancer-type tumor vs normal (the expensive loop)
# ---------------------------------------------------------------------------


def section3_per_cancer_type(gene: str) -> list[dict]:
    scratch_path = SCRATCH_DIR / f"dossier_{gene}_per_cancer_type.json"
    results: list[dict] = []
    completed = set()

    if scratch_path.exists():
        results = json.loads(scratch_path.read_text())
        completed = {r["project"] for r in results}
        print(f"[3/7] Resuming per-cancer-type loop: {len(completed)}/{len(TCGA_PROJECT_CODES)} already done", file=sys.stderr)
    else:
        print(f"[3/7] Per-cancer-type tumor vs normal ({len(TCGA_PROJECT_CODES)} project codes)...", file=sys.stderr)

    for i, code in enumerate(TCGA_PROJECT_CODES, 1):
        if code in completed:
            continue
        t0 = time.time()
        try:
            r = s.stratified_comparison_impl(
                gene,
                "cgc_sample_sample_type",
                ["Primary Tumor", "Solid Tissue Normal"],
                extra_field=s.TCGA_PROJECT_FIELD,
                extra_value=code,
            )
            by_value = {g["value"]: g for g in r["group_summary"]}
            tumor = by_value.get("Primary Tumor")
            normal = by_value.get("Solid Tissue Normal")
            entry = {
                "project": code,
                "tumor_n": tumor["samples_count"] if tumor else 0,
                "tumor_median": tumor["coverage_median"] if tumor else None,
                "normal_n": normal["samples_count"] if normal else 0,
                "normal_median": normal["coverage_median"] if normal else None,
                "error": None,
            }
        except Exception as exc:  # noqa: BLE001 -- a single bad project code must not abort a ~130-query run
            entry = {
                "project": code,
                "tumor_n": 0,
                "tumor_median": None,
                "normal_n": 0,
                "normal_median": None,
                "error": str(exc)[:300],
            }
        entry["usable"] = bool(entry["normal_n"] >= MIN_NORMAL_N and entry["tumor_n"] >= MIN_NORMAL_N)
        if entry["usable"] and entry["normal_median"] and entry["normal_median"] > 0 and entry["tumor_median"] is not None:
            entry["ratio"] = entry["tumor_median"] / entry["normal_median"]
        else:
            entry["ratio"] = None
        results.append(entry)
        completed.add(code)

        # Incremental save -- crash/interrupt safety, per project instructions.
        scratch_path.write_text(json.dumps(results, indent=2))
        elapsed = time.time() - t0
        print(
            f"  [{i}/{len(TCGA_PROJECT_CODES)}] {code}: tumor n={entry['tumor_n']} normal n={entry['normal_n']} "
            f"usable={entry['usable']} ratio={entry['ratio']} ({elapsed:.1f}s)",
            file=sys.stderr,
        )

    return results


def section3_figure(results: list[dict], gene: str) -> go.Figure:
    usable = [r for r in results if r["usable"] and r["ratio"] is not None]
    usable.sort(key=lambda r: r["ratio"])

    colors = [UP_COLOR if r["ratio"] >= 1 else DOWN_COLOR for r in usable]
    fig = go.Figure()
    fig.add_trace(
        go.Bar(
            x=[r["ratio"] for r in usable],
            y=[r["project"] for r in usable],
            orientation="h",
            marker=dict(color=colors),
            text=[f"{r['ratio']:.2f}x (n={r['tumor_n']}/{r['normal_n']})" for r in usable],
            textposition="outside",
            hovertemplate="%{y}: %{x:.2f}x<extra></extra>",
        )
    )
    fig.add_vline(x=1, line_dash="dash", line_color="grey")
    fig.update_layout(
        title=f"{gene} tumor/normal ratio by TCGA cancer type (well-powered only, n≥10 both groups)",
        xaxis_title="tumor / normal coverage_median",
        template="plotly_white",
        height=100 + 26 * len(usable),
        margin=dict(l=80, r=140),
        showlegend=False,
    )
    return fig


# ---------------------------------------------------------------------------
# Section 4: auto-generated narrative
# ---------------------------------------------------------------------------


def section4_narrative(gene: str, pan_cancer: dict, per_type: list[dict]) -> str:
    usable = [r for r in per_type if r["usable"] and r["ratio"] is not None]
    not_usable = [r for r in per_type if not r["usable"]]

    if not usable:
        return (
            f"<p>No TCGA cancer type had both a tumor and matched-normal cohort of at least "
            f"{MIN_NORMAL_N} samples for {gene}, so no per-cancer-type comparison could be made. "
            f"{len(not_usable)} of {len(per_type)} project codes were checked and excluded for "
            f"insufficient n.</p>"
        )

    highest = max(usable, key=lambda r: r["ratio"])
    lowest = min(usable, key=lambda r: r["ratio"])

    pan_by_group = {g["group"]: g for g in pan_cancer["group_summary"]}
    pan_tumor = pan_by_group.get("Primary Tumor")
    pan_normal = pan_by_group.get("Solid Tissue Normal")
    pan_ratio = None
    if pan_tumor and pan_normal and pan_normal["coverage_median"]:
        pan_ratio = pan_tumor["coverage_median"] / pan_normal["coverage_median"]

    up_types = [r for r in usable if r["ratio"] >= 1.5]
    down_types = [r for r in usable if r["ratio"] <= 0.67]
    masking_note = ""
    if up_types and down_types:
        masking_note = (
            f" Pan-cancer pooling can obscure this: {len(up_types)} well-powered cancer type(s) "
            f"show tumor at least 1.5x higher than normal ({', '.join(r['project'] for r in up_types)}) "
            f"while {len(down_types)} show tumor at least 1.5x LOWER than normal "
            f"({', '.join(r['project'] for r in down_types)}) -- opposing effects in different cancer "
            f"types that a single pan-cancer number would average out or flatten, not resolve."
        )
    elif pan_ratio is not None and 0.8 <= pan_ratio <= 1.25 and (highest["ratio"] >= 2 or lowest["ratio"] <= 0.5):
        masking_note = (
            f" The pan-cancer comparison alone looks close to flat (~{pan_ratio:.2f}x), which would "
            f"understate the per-cancer-type range actually observed ({lowest['ratio']:.2f}x-{highest['ratio']:.2f}x) "
            f"-- check individual cancer types before concluding {gene} is unchanged in tumor vs normal."
        )

    return (
        f"<p><b>{gene}</b>: of {len(per_type)} TCGA cancer types checked, {len(usable)} had a usable "
        f"matched-normal cohort (n&ge;{MIN_NORMAL_N} in both tumor and normal); "
        f"{len(not_usable)} were excluded for insufficient n ({', '.join(r['project'] for r in not_usable)}). "
        f"Among the {len(usable)} usable types, the largest tumor/normal ratio is "
        f"<b>{highest['project']}</b> at {highest['ratio']:.2f}x (tumor median {highest['tumor_median']:,.1f}, "
        f"n={highest['tumor_n']}, vs normal median {highest['normal_median']:,.1f}, n={highest['normal_n']}), "
        f"and the smallest is <b>{lowest['project']}</b> at {lowest['ratio']:.2f}x (tumor median "
        f"{lowest['tumor_median']:,.1f}, n={lowest['tumor_n']}, vs normal median {lowest['normal_median']:,.1f}, "
        f"n={lowest['normal_n']}).{masking_note} These are coverage_median ratios, per-sample but NOT "
        f"library-size normalized -- treat ratios below ~1.3x or above ~0.77x as a close call rather than "
        f"a confirmed effect, per this project's existing normalization caveat.</p>"
    )


# ---------------------------------------------------------------------------
# Section 5: exon-level detail for well-powered cancer types only
# ---------------------------------------------------------------------------


def _fetch_exons_for_project(gene: str, sample_type_value: str, project_code: str) -> list[list[str]]:
    """Reuses exon_usage's building blocks (fetch_snaptron, positional
    parsing, EXONS_EXPECTED_FIELD_COUNT) but adds a project_id constraint
    that exon_usage_impl itself does not support -- same approach the
    earlier HLX ad-hoc analysis used (see docs/plot-notes.md), kept local to
    this script rather than added to the shared tool's public contract."""
    filters = [
        ("sfilter", f"{s.TCGA_SAMPLE_TYPE_FIELD}:{sample_type_value}"),
        ("sfilter", f"{s.TCGA_PROJECT_FIELD}:{project_code}"),
    ]
    body = s.fetch_snaptron("tcgav2", "exons", gene, filters=filters)
    lines = [line for line in body.split("\n") if line.strip()]
    if len(lines) < 2:
        return []
    rows = [line.split("\t") for line in lines[1:]]
    return [r for r in rows if len(r) == s.EXONS_EXPECTED_FIELD_COUNT]


def section5_exon_detail(gene: str, per_type: list[dict]) -> dict[str, list[dict]]:
    usable_codes = [r["project"] for r in per_type if r["usable"]]
    print(f"[5/7] Exon-level detail for {len(usable_codes)} well-powered cancer type(s): {usable_codes}", file=sys.stderr)

    scratch_path = SCRATCH_DIR / f"dossier_{gene}_exon_by_type.json"
    all_results: dict[str, list[dict]] = {}
    if scratch_path.exists():
        all_results = json.loads(scratch_path.read_text())

    for code in usable_codes:
        if code in all_results:
            continue
        t0 = time.time()
        try:
            tumor_rows = _fetch_exons_for_project(gene, "Primary Tumor", code)
            normal_rows = _fetch_exons_for_project(gene, "Solid Tissue Normal", code)
        except Exception as exc:  # noqa: BLE001 -- one bad project code must not abort the loop
            print(f"  [{code}] query failed -- {str(exc)[:200]}", file=sys.stderr)
            all_results[code] = []
            scratch_path.write_text(json.dumps(all_results, indent=2))
            continue

        exon_data: dict[str, dict] = {}
        for rows, prefix in [(tumor_rows, "tumor"), (normal_rows, "normal")]:
            for r in rows:
                gene_id, gene_name, _gene_type, _bp = s._split_packed_exon_id(r[11])
                if gene_name.upper() != gene.upper():
                    continue
                try:
                    start_i, end_i = int(r[3]), int(r[4])
                    coverage_median = float(r[16])
                except (ValueError, TypeError):
                    continue
                key = f"{r[2]}:{start_i}-{end_i}"
                entry = exon_data.setdefault(key, {"exon": key, "start": start_i})
                entry[f"{prefix}_coverage_median"] = coverage_median

        exon_list = sorted(exon_data.values(), key=lambda e: e["start"])
        for e in exon_list:
            tm, nm = e.get("tumor_coverage_median"), e.get("normal_coverage_median")
            e["ratio"] = (tm / nm) if (tm is not None and nm not in (None, 0)) else None

        all_results[code] = exon_list
        scratch_path.write_text(json.dumps(all_results, indent=2))
        print(f"  [{code}] {len(exon_list)} exons ({time.time() - t0:.1f}s)", file=sys.stderr)

    return all_results


def section5_figure(exon_by_type: dict[str, list[dict]], gene: str) -> go.Figure | None:
    """Small multiples -- one panel per cancer type, same x-axis (genomic
    position) and same shared y-axis range across all panels -- rather than
    overlapping colored lines on one plot. A distinct cancer type's
    exon-level pattern should be visible at a glance without having to
    discriminate line colors across many overlapping series."""
    codes = [c for c, exons in exon_by_type.items() if any(e["ratio"] is not None for e in exons)]
    if not codes:
        return None

    per_code_exons = {code: [e for e in exon_by_type[code] if e["ratio"] is not None] for code in codes}
    all_ratios = [e["ratio"] for exons in per_code_exons.values() for e in exons]
    y_min, y_max = min(all_ratios), max(all_ratios)
    y_pad = (y_max - y_min) * 0.08 or 0.5
    y_range = [max(0, y_min - y_pad), y_max + y_pad]  # a ratio can't be negative

    n = len(codes)
    cols = math.ceil(math.sqrt(n))
    rows = math.ceil(n / cols)

    fig = make_subplots(
        rows=rows,
        cols=cols,
        subplot_titles=codes,
        shared_xaxes=True,
        shared_yaxes=True,
        horizontal_spacing=0.04,
        vertical_spacing=0.08 if rows > 1 else 0.15,
    )

    for i, code in enumerate(codes):
        row, col = divmod(i, cols)
        row, col = row + 1, col + 1
        exons = per_code_exons[code]
        fig.add_trace(
            go.Scatter(
                x=[e["start"] for e in exons],
                y=[e["ratio"] for e in exons],
                mode="lines+markers",
                line=dict(color=UP_COLOR),
                marker=dict(size=4),
                showlegend=False,
                hovertemplate=f"{code}" + "<br>pos=%{x}<br>ratio=%{y:.2f}x<extra></extra>",
            ),
            row=row,
            col=col,
        )
        fig.add_hline(y=1, line_dash="dash", line_color="grey", line_width=1, row=row, col=col)

    fig.update_yaxes(range=y_range)
    fig.update_xaxes(showticklabels=False)  # genomic coordinates are long and repeat identically per panel -- declutter
    fig.update_layout(
        title=f"{gene} per-exon tumor/normal ratio, by cancer type (well-powered types only, shared y-axis)",
        template="plotly_white",
        height=220 * rows,
        margin=dict(t=80),
    )
    fig.update_annotations(font_size=12)
    return fig


# ---------------------------------------------------------------------------
# Section 6: GTEx exon usage across all 31 tissues (heatmap)
# ---------------------------------------------------------------------------


def section6_gtex_exon(gene: str) -> dict:
    cached = _load_cached_csv("exon_usage_tissue", gene)
    if cached:
        print(f"[6/7] GTEx exon usage -- reusing cached CSV {cached.name}", file=sys.stderr)
        return {"gene": gene, "per_sample_csv_path": str(cached)}

    print(f"[6/7] GTEx exon usage ({len(s.GTEX_TISSUES)} tissues)...", file=sys.stderr)
    return s.exon_usage_impl(gene, comparison="tissue")


def section6_figure(result: dict, gene: str) -> go.Figure | None:
    per_sample = pd.read_csv(result["per_sample_csv_path"])
    if per_sample.empty:
        return None

    summary = per_sample.groupby(["exon", "group"])["coverage"].median().reset_index()

    starts = summary["exon"].str.extract(r":(\d+)-")[0].astype(int)
    exon_order = summary.assign(start=starts).drop_duplicates("exon").sort_values("start")["exon"].tolist()
    tissue_peak = summary.groupby("group")["coverage"].max().sort_values(ascending=False)
    tissue_order = tissue_peak.index.tolist()

    pivot = summary.pivot(index="group", columns="exon", values="coverage").reindex(index=tissue_order, columns=exon_order)
    z = np.log10(pivot.values.astype(float))  # log10: same huge-dynamic-range reason as the other coverage charts here

    fig = go.Figure(
        data=go.Heatmap(
            z=z,
            x=exon_order,
            y=tissue_order,
            colorscale="Blues",
            colorbar=dict(title="log10(coverage_median)"),
            customdata=pivot.values,
            hovertemplate="%{y} / %{x}<br>coverage_median=%{customdata}<extra></extra>",
        )
    )
    fig.update_layout(
        title=f"{gene} per-exon coverage across GTEx tissues",
        xaxis_title="exon (genomic position order)",
        yaxis_title="tissue (ranked by peak coverage)",
        template="plotly_white",
        height=120 + 20 * len(tissue_order),
        margin=dict(b=140),
    )
    fig.update_xaxes(tickangle=45)
    return fig


# ---------------------------------------------------------------------------
# Section 7: Junction usage -- TCGA tumor/normal + GTEx tissue, well-powered only
# ---------------------------------------------------------------------------


def section7_junction_tumor_normal(gene: str) -> dict:
    cached = _load_cached_csv("junction_usage", gene)
    if cached:
        print(f"[7a/7] TCGA junction usage (tumor vs normal) -- reusing cached CSV {cached.name}", file=sys.stderr)
        return {"gene": gene, "per_junction_csv_path": str(cached)}

    print("[7a/7] TCGA junction usage (tumor vs normal)...", file=sys.stderr)
    return s.junction_usage_impl(gene, comparison="tumor_vs_normal")


def section7_tumor_normal_figure(result: dict, gene: str) -> go.Figure | None:
    long_df = pd.read_csv(result["per_junction_csv_path"])
    if long_df.empty:
        return None

    wide = long_df.pivot_table(
        index=["junction", "novelty_class"], columns="group", values=["pct_of_cohort", "coverage_median"]
    )
    wide.columns = [f"{a}_{b}" for a, b in wide.columns]
    wide = wide.reset_index()
    for col in [
        "pct_of_cohort_Primary Tumor",
        "pct_of_cohort_Solid Tissue Normal",
        "coverage_median_Primary Tumor",
        "coverage_median_Solid Tissue Normal",
    ]:
        if col not in wide.columns:
            wide[col] = 0.0
    wide = wide.fillna(0.0)

    # "well-powered" junction -- same MIN_NORMAL_N threshold/rationale as sections
    # 3/5's cancer-type selection, applied here to coverage_median instead of n
    # (see docs/plot-notes.md's junction min_coverage_median finding: an unfiltered
    # ranking surfaces ~1-read noise, not real signal).
    powered = wide[
        (wide["coverage_median_Primary Tumor"] >= MIN_NORMAL_N)
        & (wide["coverage_median_Solid Tissue Normal"] >= MIN_NORMAL_N)
    ].copy()
    if powered.empty:
        return None

    powered["diff"] = powered["pct_of_cohort_Primary Tumor"] - powered["pct_of_cohort_Solid Tissue Normal"]
    powered = powered.reindex(powered["diff"].abs().sort_values(ascending=False).index).head(25)
    powered = powered.sort_values("diff")

    colors = [UP_COLOR if d >= 0 else DOWN_COLOR for d in powered["diff"]]
    fig = go.Figure()
    fig.add_trace(
        go.Bar(
            x=powered["diff"],
            y=powered["junction"],
            orientation="h",
            marker=dict(color=colors),
            text=powered["novelty_class"],
            textposition="outside",
            hovertemplate="%{y}: %{x:.4f}<extra></extra>",
        )
    )
    fig.add_vline(x=0, line_dash="dash", line_color="grey")
    fig.update_layout(
        title=f"{gene} tumor/normal junction usage (well-powered only, coverage_median≥{MIN_NORMAL_N} both groups)",
        xaxis_title="tumor_pct_of_cohort - normal_pct_of_cohort",
        template="plotly_white",
        height=100 + 26 * len(powered),
        margin=dict(l=220, r=140),
        showlegend=False,
    )
    return fig


def section7_junction_tissue(gene: str) -> dict:
    cached = _load_cached_csv("junction_usage_tissue", gene)
    if cached:
        print(f"[7b/7] GTEx junction usage (tissue) -- reusing cached CSV {cached.name}", file=sys.stderr)
        return {"gene": gene, "per_junction_csv_path": str(cached)}

    print(f"[7b/7] GTEx junction usage (tissue, {len(s.GTEX_TISSUES)} tissues)...", file=sys.stderr)
    return s.junction_usage_impl(gene, comparison="tissue")


def section7_tissue_figure(result: dict, gene: str) -> go.Figure | None:
    long_df = pd.read_csv(result["per_junction_csv_path"])
    if long_df.empty:
        return None

    # well-powered junction: at least one tissue with coverage_median >= MIN_NORMAL_N
    peak_by_junction = long_df.groupby("junction")["coverage_median"].max()
    powered_junctions = peak_by_junction[peak_by_junction >= MIN_NORMAL_N].sort_values(ascending=False).head(25).index
    if len(powered_junctions) == 0:
        return None

    filtered = long_df[long_df["junction"].isin(powered_junctions)]
    pivot = filtered.pivot_table(index="group", columns="junction", values="coverage_median", aggfunc="median")
    junction_order = list(powered_junctions)
    tissue_peak = pivot.max(axis=1).sort_values(ascending=False)
    pivot = pivot.reindex(index=tissue_peak.index, columns=junction_order)

    z = np.log10(pivot.values.astype(float))
    fig = go.Figure(
        data=go.Heatmap(
            z=z,
            x=junction_order,
            y=tissue_peak.index.tolist(),
            colorscale="Oranges",
            colorbar=dict(title="log10(coverage_median)"),
            customdata=pivot.values,
            hovertemplate="%{y} / %{x}<br>coverage_median=%{customdata}<extra></extra>",
        )
    )
    fig.update_layout(
        title=f"{gene} junction usage across GTEx tissues (well-powered only, coverage_median≥{MIN_NORMAL_N} in ≥1 tissue)",
        xaxis_title="junction (ranked by peak coverage)",
        yaxis_title="tissue (ranked by peak coverage)",
        template="plotly_white",
        height=120 + 20 * len(tissue_peak),
        margin=dict(b=160),
    )
    fig.update_xaxes(tickangle=45)
    return fig


# ---------------------------------------------------------------------------
# Assembly
# ---------------------------------------------------------------------------


def build_dossier(gene: str) -> Path:
    t_start = time.time()
    gene = gene.strip()

    r1 = section1_tissue(gene)
    r2 = section2_pan_cancer(gene)
    r3 = section3_per_cancer_type(gene)
    narrative_html = section4_narrative(gene, r2, r3)
    r5 = section5_exon_detail(gene, r3)
    r6 = section6_gtex_exon(gene)
    r7a = section7_junction_tumor_normal(gene)
    r7b = section7_junction_tissue(gene)

    fig1 = section1_figure(r1, gene)
    fig2 = section2_figure(r2, gene)
    fig3 = section3_figure(r3, gene)
    fig5 = section5_figure(r5, gene)
    fig6 = section6_figure(r6, gene)
    fig7a = section7_tumor_normal_figure(r7a, gene)
    fig7b = section7_tissue_figure(r7b, gene)

    chart1_html = fig1.to_html(include_plotlyjs="inline", full_html=False, div_id="chart1")
    chart2_html = fig2.to_html(include_plotlyjs=False, full_html=False, div_id="chart2")
    chart3_html = fig3.to_html(include_plotlyjs=False, full_html=False, div_id="chart3")
    chart5_html = fig5.to_html(include_plotlyjs=False, full_html=False, div_id="chart5") if fig5 else "<p><i>No well-powered cancer types -- no exon-level section.</i></p>"
    chart6_html = fig6.to_html(include_plotlyjs=False, full_html=False, div_id="chart6") if fig6 else "<p><i>No exon data returned -- no GTEx exon heatmap.</i></p>"
    chart7a_html = (
        fig7a.to_html(include_plotlyjs=False, full_html=False, div_id="chart7a")
        if fig7a
        else f"<p><i>No junction cleared the coverage_median&ge;{MIN_NORMAL_N} floor in both TCGA groups -- no tumor/normal junction chart.</i></p>"
    )
    chart7b_html = (
        fig7b.to_html(include_plotlyjs=False, full_html=False, div_id="chart7b")
        if fig7b
        else f"<p><i>No junction cleared the coverage_median&ge;{MIN_NORMAL_N} floor in any GTEx tissue -- no tissue junction heatmap.</i></p>"
    )

    not_usable = [r for r in r3 if not r["usable"]]
    not_usable_html = ", ".join(
        f"{r['project']} (tumor n={r['tumor_n']}, normal n={r['normal_n']})" for r in not_usable
    ) or "none"

    caveat_html = (
        "Values are per-sample coverage, NOT library-size / depth normalized -- no AUC scaling has "
        "been applied. A gap of roughly &gt;10x is safe to call biological; a closer gap may partly "
        "reflect sequencing depth rather than true expression differences. Section 3/5 ratios use "
        f"coverage_median and require n&ge;{MIN_NORMAL_N} in both tumor and normal to be shown at all. "
        f"Sections 6/7 heatmaps and the section 7 bar chart are similarly restricted to "
        f"coverage_median&ge;{MIN_NORMAL_N} -- an unfiltered ranking mostly surfaces ~1-read noise, "
        f"not real signal (see docs/plot-notes.md)."
    )

    elapsed = time.time() - t_start

    page_html = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>{gene} gene dossier</title>
<style>
  body {{ font-family: -apple-system, Segoe UI, Helvetica, Arial, sans-serif; margin: 0; padding: 24px; background: #fafafa; color: #222; max-width: 1100px; }}
  h1 {{ font-size: 1.6rem; margin-bottom: 4px; }}
  h2 {{ font-size: 1.15rem; margin-top: 40px; border-bottom: 1px solid #ddd; padding-bottom: 4px; }}
  .caveat {{ font-size: 0.85rem; color: #555; max-width: 900px; margin: 8px 0 24px; line-height: 1.4; }}
  .narrative {{ font-size: 0.95rem; line-height: 1.5; background: #fff; border: 1px solid #ddd; border-radius: 6px; padding: 16px; max-width: 900px; }}
  .excluded {{ font-size: 0.85rem; color: #777; margin-top: 8px; }}
  .meta {{ font-size: 0.8rem; color: #888; }}
</style>
</head>
<body>
<h1>{gene} gene dossier</h1>
<p class="caveat">{caveat_html}</p>

<h2>1. GTEx tissue expression</h2>
{chart1_html}

<h2>2. TCGA pan-cancer tumor vs normal</h2>
{chart2_html}

<h2>3. Per-cancer-type tumor vs normal</h2>
{chart3_html}
<p class="excluded">Excluded for insufficient matched-normal n (&lt;{MIN_NORMAL_N}): {not_usable_html}</p>

<h2>4. Summary</h2>
<div class="narrative">{narrative_html}</div>

<h2>5. Exon-level detail (well-powered cancer types only)</h2>
{chart5_html}

<h2>6. GTEx exon usage across tissues</h2>
{chart6_html}

<h2>7a. TCGA junction usage (tumor vs normal)</h2>
{chart7a_html}

<h2>7b. GTEx junction usage across tissues</h2>
{chart7b_html}

<p class="meta">Generated in {elapsed:.1f}s.</p>
</body>
</html>
"""

    out_path = REPORTS_DIR / f"{gene}_dossier.html"
    out_path.write_text(page_html, encoding="utf-8")
    print(f"\nSaved dossier to: {out_path.resolve()}", file=sys.stderr)
    print(f"Total time: {elapsed:.1f}s", file=sys.stderr)
    return out_path


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python html_gene_dossier.py <gene_name>")
        sys.exit(1)
    build_dossier(sys.argv[1])
