# Plot notes

Caveats that used to be printed on the chart itself (title/subtitle) live
here instead, one entry per plot type, so figures stay uncluttered but the
caveat is still on record. Read the relevant entry before presenting a
plot from this repo.

## Tissue expression (`plot_tissue_expression.R`, tool: `gene_expression_across_tissues`)

Per-sample coverage (and coverage_median/coverage_avg) are NOT
library-size / depth normalized -- no AUC scaling has been applied. A gap
between tissues of roughly >10x is safe to call biological; a close gap
may partly reflect sequencing depth differences rather than true
expression differences.

Several low-expressing tissues show a dense cluster of samples tied at the
same low discrete coverage value (e.g. 152 in the LRRC10 test case), with a
real tail of higher individual values above it -- confirmed by checking
the actual per-sample distribution, not assumed. This is not a hard floor
across all samples; it's the single most common value, consistent with it
representing roughly one read's worth of coverage over the region, with
multiples of it (2 reads, 3 reads, ...) as the next most common values.
The boxplot+jitter shows this honestly (a dense low cluster plus a
scattered tail); a violin plot was deliberately avoided here because its
kernel-density smoothing would (a) turn a pile of exact ties into a fake
continuous shape, and (b) misrepresent tissues with very small n (e.g.
Kidney, n=11 in the LRRC10 test case) as a smooth density when there is
almost no data behind it.

## Stratified comparison (`plot_stratified_comparison.R`, tool: `stratified_comparison`)

Bars are coverage_median (per-sample coverage, NOT library-size normalized
-- same caveat as tissue expression above), grouped by scope: pan-cancer
(pooled across all TCGA cancer types) vs constrained to one cancer type
via extra_field/extra_value. The pan-cancer bars are shown deliberately
alongside the within-cancer-type bars, not because pan-cancer pooling is a
valid comparison on its own -- it isn't (see stratified_comparison's own
cancer_type_confound_caveat) -- but specifically to make the size of that
confound visible: MKI67 Stage I vs Stage IV pan-cancer looked like a ~4.5x
difference (70,473 vs 320,286), but within TCGA-KIRC alone (n=300/n=104)
it's ~1.3x (81,430 vs 109,672.5). Read the pan-cancer bars as "what you'd
wrongly conclude without controlling for cancer type," not as a result in
their own right.

## Tumor vs normal (`plot_tumor_vs_normal.R`, tool: `gene_tumor_vs_normal`)

Same per-sample coverage / not-library-size-normalized caveat as tissue
expression above. The two groups' cohort sizes are very unequal (9,961
Primary Tumor vs 740 Solid Tissue Normal in tcgav2 -- a ~13.5x imbalance),
which is why the chart annotates n for each group directly rather than
leaving the reader to guess it from point density; do not read relative
jitter-cloud density as a proxy for relative group size without checking
the n labels, since the lower per-point alpha on the tumor group (needed
so 9,961 overlapping points stay legible) also reduces how solid that
cloud looks. Boxplot+jitter, not violin, same reasoning as tissue
expression (exact-tie piles and n honesty).

## Junction usage (`plot_junction_usage.R`, tool: `junction_usage`)

Bar length is tumor_pct_of_cohort - normal_pct_of_cohort (fraction of each
group's whole cohort using that junction at all, not a raw sample count --
see junction_usage's own normalization_caveat for why). Bars are colored by
novelty_class so novel/known splicing can be read at a glance.

IMPORTANT, checked directly rather than assumed: for the TP53 top-10 test
case, every one of the top 5 ranked junctions has coverage_median = 1.0 in
BOTH tumor and normal. These are junctions detected in many samples
(hundreds to thousands, given TP53's very high expression and the cohort
sizes involved) but supported by only about one read per sample where
detected. Ranking by pct_of_cohort_diff alone surfaces "widely observed"
junctions, not necessarily "strongly used" ones -- a chart of the top
differential junctions by this metric is a map of splice-site usage
breadth, not evidence of a robust alternative-splicing event, unless
coverage_median is also checked and is not sitting at the ~1-read floor.
Check coverage_median (in the CSV, not shown on this chart) before
presenting any specific bar here as a meaningful splicing difference.

`plot_junction_usage.R` takes an optional 4th argument, min_coverage_median
(default 0 = no filtering), which drops any junction below that
coverage_median in EITHER group before ranking. Confirmed on TP53
(coverage_median >= 10 in both groups): only 11 of 2,526 junctions pass at
all. The ranking flips character completely -- all 11 are "known"
junctions (no novel ones survive at this depth in this test case), each
used in 97-100% of samples in BOTH groups, with small NEGATIVE differences
(normal slightly higher than tumor, roughly 0.4-2.0 percentage points).
This reads as real signal in the sense that it clears the noise floor
(coverage_median 22-71, not ~1), but the differences themselves are small
-- consistent with a modest, not dramatic, shift toward slightly less
consistent use of TP53's canonical junctions in tumor vs normal, not
evidence of tumor-specific novel splicing (that signal, if present, is
apparently below this coverage threshold in this gene/cohort and was only
visible -- unreliably -- in the unfiltered, noise-dominated ranking
above).

Confirmed on MET (coverage_median >= 10 in both groups): only 19 of 1,724
junctions pass. All 19 are "known" (no novel/novel-combination junctions
survive at this depth), all 19 NEGATIVE (normal higher than tumor), with a
gap of roughly 6-15 percentage points -- larger than TP53's 0.4-2.0pp at
the same threshold, and coverage_median 12-47 in the survivors, comfortably
above the ~1-read floor. This is a real signal by the same test used for
TP53 (clears the floor, unanimous direction).

BUT: the >=10 filter that makes this ranking trustworthy is also why it
says nothing about MET exon 14 skipping, MET's actual well-known
cancer-relevant splicing event. Exon 14 skipping is typically a MINORITY
isoform within any single tumor sample (most reads at that locus still
come from the canonical, exon-14-included transcript, even in a sample
where skipping is present and clinically significant) -- exactly the kind
of low-coverage-per-sample, rare-in-the-aggregate signal this filter is
designed to exclude. So this filtered ranking and its plot show "canonical
MET splicing is slightly less consistent, genome-wide-style, in tumor than
normal" -- a real but modest effect across MET's ordinary junctions -- NOT
evidence about exon 14 skipping specifically, which would need to be
looked up as a targeted, individual junction (not surfaced by a
coverage-filtered top-N ranking) rather than inferred from this chart.
State this plainly if `plots/MET_junction_usage_min10.png` comes up in the
demo -- it is easy to mistake for "we found the exon 14 skipping signal"
because MET is the gene most people associate with that specific event,
and this plot is not that.

## Exon usage (`plot_exon_usage.R`, tool: `exon_usage`)

Y-axis is the tumor/normal coverage_median RATIO per exon (aggregated from
the per-sample CSV), not two raw coverage lines -- the ratio is what stays
legible when absolute per-exon coverage spans orders of magnitude along the
gene (as it does for both test genes below). LINEAR y-axis (not log10, and
not binned -- binning would merge nearby exons and obscure a gradient that
depends on per-exon resolution). Checked before committing to linear: the
ratio range for both test genes (1.04x-12.86x, ~12x total spread at most)
is modest enough that the largest point (MKI67's 12.86x) does not dominate
or flatten the rest of the chart -- the declining MKI67 gradient stays
clearly visible on a linear scale, arguably more intuitively than log10
since vertical distance now matches the actual ratio difference. Revisit
this choice if a future gene's ratio range is much wider (e.g. one exon at
100x+ while others sit near 1x) -- that could reintroduce the domination
problem linear scale is currently avoiding. X-axis is exon start position
(genomic coordinate, not exon rank/order), with a dashed reference line at
ratio=1 (no tumor/normal difference).

Same per-sample-coverage / not-library-size-normalized caveat as the other
coverage plots -- see tissue expression above.

Two test cases, deliberately contrasting:
- **MKI67**: ratio is NOT flat -- it forms a real gradient along the gene,
  from 12.86x (chr10:128096659-128099255, the exon nearest one end) down to
  1.79x (chr10:128126193-128126235, a 42bp exon near the other end),
  decreasing roughly monotonically in between. This is a concrete example
  of what exon-level resolution adds beyond gene_tumor_vs_normal's single
  ~10.7x whole-gene number: the tumor/normal effect is not uniform across
  MKI67's exons, it is strongest at one end and weakest at the other.
- **TP53**: ratio is flat by comparison, 1.04x-1.44x across all 39 exons,
  consistent with the modest ~1.3x whole-gene effect already found. No
  exon stands out the way MKI67's do.
Don't over-read the MKI67 gradient shape as a specific biological
mechanism (e.g. don't claim it's a particular alternative-splicing or
3'-bias artifact) without independent evidence -- the plot shows WHERE the
tumor/normal effect concentrates along the gene, not WHY.

## HLX per-cancer-type tumor/normal + exon usage (ad hoc, no dedicated tool/plot script yet)

Produced by hand-looping `gene_tumor_vs_normal`-style queries per TCGA
project code and a custom `/exons` fetch (same positional-parsing approach
as `exon_usage`, extended with a `gdc_cases.project.project_id` sfilter) --
not run through a dedicated R plot script, so recorded here instead of
under a `plot_*.R` heading. Underlying data saved to `outputs/`:
`stratified_comparison_HLX_gdc_cases_project_project_id_tumor_vs_normal_20260918_105913.csv`
(33 TCGA project codes, tumor+normal coverage_avg/median/n) and
`exon_usage_HLX_KIRC_20260918_105913.csv` /
`exon_usage_HLX_LUSC_20260918_105913.csv` /
`exon_usage_HLX_COAD_20260918_105913.csv` (per-region coverage, tumor vs
normal, for the three flagged cancer types).

**Finding: pan-cancer pooling hid a real, strong effect instead of
inflating a fake one -- the mirror image of the MKI67 Stage I/IV lesson
above.** `gene_tumor_vs_normal` on HLX pan-cancer looks flat (median
33,665 tumor vs 34,772 normal, n=9,961/740, ratio ~0.97x). Looping the
same tumor-vs-normal comparison per TCGA project code (16 of 33 types
kept, requiring >=10 matched normal samples) shows HLX is actually
strongly dysregulated, in OPPOSITE directions in different cancer types:
~5x HIGHER in tumor in KIRC (116,978 vs 23,536, n=543/72) but ~3-3.6x
LOWER in tumor in LUSC (33,314 vs 120,002, n=504/51) and COAD (16,295 vs
49,866, n=503/41). All three are well-powered, not thin-n artifacts. The
opposing signs cancel out under pan-cancer pooling, producing a null
result that is itself an artifact of pooling -- same underlying lesson as
the MKI67/stage confound (pooling across cancer types confounds any
per-gene comparison) but manifesting as a false negative here instead of
the earlier false positive-looking inflation. Always check per-cancer-type
before trusting (or dismissing) a pan-cancer `gene_tumor_vs_normal` result.

Exon-level follow-up on the three flagged cancer types: HLX is annotated
as a single coding exon in humans, so Snaptron's `/exons` endpoint splits
it into 6 sub-regions that are UTR/transcript-boundary segments, not
intron-exon splice boundaries. In all three cancer types, most regions
scale roughly proportionally with the whole-gene tumor/normal change (not
independent differential exon usage). One exception in all three: a 30bp
micro-region (chr1:220881069-220881099) whose share of total gene
coverage moves opposite to the rest of the gene -- relatively depleted in
KIRC tumor, relatively enriched in LUSC/COAD tumor -- but its absolute
coverage is only in the tens of reads (vs hundreds of thousands for the
whole gene), so this reads as a low-confidence lead, not a confirmed
splicing event without deeper, targeted follow-up.
