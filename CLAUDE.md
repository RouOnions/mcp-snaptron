# Snaptron API -- confirmed findings from pre-hackathon recon

Context: MCP server wrapping the Snaptron REST API
(https://snaptron.cs.jhu.edu) for cancer/tissue RNA-seq questions
(splice junctions, novel junctions, gene/exon usage, tumor vs normal,
tissue specificity). Everything below was empirically tested against
the live server on tcgav2 and gtexv2, not taken from docs alone.

## Non-negotiable response-handling rules

1. **The plain-text status line is ENDPOINT-SPECIFIC, not universal**
   (corrected after live testing). Only `/snaptron/registry` prefixes the
   body with `200 OK\n` before the payload. `/genes`, `/snaptron`, `/exons`
   and even traceback responses do NOT -- they start directly with the TSV
   header (or with `Traceback`). Detect a status line with a pattern (first
   line matches `^\d{3}\s+\S`, e.g. `200 OK`) and strip it ONLY if present.
   NEVER blindly strip line 1 -- on `/genes` that silently deletes the real
   `DataSource:Type...` header row and corrupts every parse.
2. **A malformed query returns HTTP 200 with a Python traceback in the
   body.** The HTTP status code cannot be trusted to detect failure.
   Every call must scan the body text for `Traceback (most recent call
   last)` before treating the response as data. This is the single most
   important gotcha -- do not skip it. (Confirmed live.)
3. Compilation is a **URL path segment**, not a query parameter:
   `https://snaptron.cs.jhu.edu/{compilation}/{endpoint}?...`
   e.g. `/tcgav2/snaptron?regions=TP53`, NOT `/snaptron?compilation=tcgav2`.
   (Confirmed live.)

## Tech stack note (MCP Python SDK)

The official `mcp` package is now **v2.x** (`pip install mcp` pulls 2.x). In
v2, `FastMCP` was renamed to `MCPServer` and moved:
`from mcp.server.mcpserver import MCPServer` (NOT `mcp.server.fastmcp`, which
was removed in v2). The `@mcp.tool()` decorator is unchanged. Note there is
also a SEPARATE standalone `fastmcp` package (v4.x, `from fastmcp import
FastMCP`) -- different thing; this project uses the official `mcp` one.

## Confirmed working endpoints (tcgav2)

| Endpoint | Purpose | Status |
|---|---|---|
| `/snaptron/registry` | lists compilations + their metadata fields | works |
| `/{comp}/snaptron?regions=X` | junctions | works, clean schema |
| `/{comp}/genes?regions=X` | gene-level quantification | works, clean schema -- CONFIRMED live end-to-end (tissue-expression tool built on it) |
| `/{comp}/exons?regions=X` | exon-level | works BUT sfilter did NOT reduce row count when tested -- do not assume it supports tumor/normal filtering like genes/junctions does. Re-test before relying on it. |
| `/{comp}/bases?regions=chr:start-end` | per-base coverage | **AVOID -- see below** |
| `/{comp}/samples?ids=...` | sample metadata lookup | **AVOID -- see below** |

## Confirmed filter syntax

- Sample-metadata filter: `sfilter=<field>:<value>` (colon, not `=`).
  e.g. `sfilter=cgc_sample_sample_type:Primary%20Tumor` -- always
  URL-encode the value (`%20` for spaces).
- Row/feature filter: `rfilter=<field>:<value>`, e.g. `rfilter=annotated:0`.
- Numeric operators go BEFORE the colon: `rfilter=samples_count>:5`
  (means >=5, confirmed inclusive). `samples_count>=5` (operator+equals
  combined) crashes with a traceback.
- Multiple filters: repeat the param (`&rfilter=a:b&rfilter=c:d`),
  combines as logical AND.
- The `filters=` parameter (as opposed to `rfilter=`/`sfilter=`) is
  broken -- crashes with a traceback. Don't use it.

## Confirmed grouping fields

- **TCGA tumor vs normal**: field `cgc_sample_sample_type`, values
  `Primary Tumor` / `Solid Tissue Normal` (exact GDC vocabulary,
  case-sensitive). Discrimination verified: different values produce
  different row counts AND different sample IDs.
- **GTEx tissue**: field `SMTS` (broad tissue, e.g. "Heart", "Liver") --
  use this one. `SMTSD` is detailed tissue (e.g. "Heart - Left
  Ventricle") if finer granularity is needed. Do NOT use `SMTSPAX` or
  `SMTSISCH` -- these look like tissue fields in a naive keyword search
  but are QC/technical covariates (fixation time, ischemic time), not
  tissue names. This exact mistake caused a server-side connection drop
  when tested.
  - **Full tissue list (31 values), live-confirmed, not guessed**: the
    registry lists field names/types but NOT enumerated values for text
    fields. Get the real list from the bulk metadata file
    `https://snaptron.cs.jhu.edu/data/gtexv2/samples.tsv` (~24MB, 19,214
    rows) -- read its `SMTS` column, take `value_counts()`. Then
    round-trip EVERY value with a live `sfilter=SMTS:<value>` query
    (e.g. on GAPDH) and compare the returned `samples_count` to the bulk
    file's row count for that value -- this is the check that catches
    collisions (see next point). The 31: Adipose Tissue, Adrenal Gland,
    Bladder, Blood, Blood Vessel, Bone Marrow, Brain, Breast, Cervix
    Uteri, Colon, Esophagus, Fallopian Tube, Heart, Kidney, Liver, Lung,
    Muscle, Nerve, Ovary, Pancreas, Pituitary, Prostate, Salivary Gland,
    Skin, Small Intestine, Spleen, Stomach, Testis, Thyroid, Uterus,
    Vagina.
  - **`SMTS:Blood` collision -- CONFIRMED, one-time fix applied**:
    `sfilter=SMTS:Blood` word-matches `Blood Vessel` too, silently
    pooling both tissues (live count 2,446 = true Blood 1,048 + Blood
    Vessel 1,398). Confirmed via `sfilter=SMTS:Blo` returning zero rows
    (partial-word doesn't match, but `Blood` as a whole word inside
    `Blood Vessel` does). Fix: cross-reference the raw `SMTS:Blood`
    result against the bulk metadata file locally to strip out the
    Blood Vessel contamination, then recompute avg/median from the
    clean subset (verified this local recompute matches the server's
    own numbers exactly for uncontaminated tissues -- no new drift
    introduced). All other 30 of 31 values checked clean via the same
    live-vs-bulk-file round-trip, including the smallest cohorts
    (Fallopian Tube n=9, Cervix Uteri n=19, Bladder n=21) where a
    collision would be most visible -- Blood is the only one. This check
    only covers SMTS itself; other fields (TCGA fields, SMTSD) have NOT
    been checked for the same word-match behavior -- re-verify if used.

## Known-broken -- do not build tools on these

- **`/bases` (per-base coverage over arbitrary regions)**: broken header
  (only names ~2 of ~11,373 per-sample columns), sample identity is not
  recoverable from the response, payload explodes with region size
  (10kb region = 232MB, ~16s), and `sids=` parameter to narrow it either
  crashes or drops the connection (`ChunkedEncodingError`) on every test
  run. Four independent failure modes. Treat as non-viable for an MCP
  tool. Worth stating as a documented negative finding, not silently
  ignoring.
- **`/samples?ids=...`**: `sids=` param returns empty (wrong param name).
  `ids=` works but the response has a column-count mismatch between
  header (~939 fields) and data rows (~932) that is NOT a simple
  positional offset -- several fields are dropped per-row wherever a
  sample lacks a value, while the header lists the full superset. Do
  not consume this endpoint's output by column position. The `fields=`
  parameter meant to request a narrow, safe subset also crashes.
- **`sfilter=study:...` on tcgav2 `/genes`**: triggers a connection drop
  (`ChunkedEncodingError: Response ended prematurely`) -- same failure
  signature as `/bases`/`/samples`. `study` looked like a plausible
  shorthand for cancer-type/project (5-char field name), but is broken.
  Use `gdc_cases.project.project_id` instead (confirmed working -- see
  below), with the BARE project code as the value (`BRCA`, not
  `TCGA-BRCA` -- the prefixed form returns 0 rows).

## Confirmed clinical/stratification fields (tcgav2)

- **Cancer type / project**: `gdc_cases.project.project_id`, bare codes
  (`BRCA`, `LUAD`, `KIRC`, etc, NOT `TCGA-` prefixed). Confirmed working;
  essential for constraining any pan-cancer field (like stage) to one
  cancer type -- see caveat below.
- **Pathologic stage**: `cgc_case_pathologic_stage` (values like
  `Stage I`, `Stage IV`) and `gdc_cases.diagnoses.tumor_stage` (lowercase
  values, near-duplicate GDC harmonization, similar but NOT identical
  counts -- pick one, don't assume interchangeable).
- **IMPORTANT -- pan-cancer pooling confound**: querying stage (or any
  clinical field) WITHOUT also constraining `project_id` pools all ~33
  TCGA cancer types together. This confounds cancer-type composition
  with the field of interest -- e.g. MKI67 Stage I vs Stage IV showed
  ~4.5x pan-cancer, which shrank to ~1.3x once constrained to one cancer
  type (KIRC), because cancer types differ hugely in baseline expression
  AND in their stage distribution. ALWAYS constrain clinical-field
  comparisons to one `project_id` for a valid comparison. When doing so,
  the housekeeping-gene cohort-total query must ALSO be constrained to
  the same project, or the "fraction of cohort" math silently reverts to
  meaning "fraction of pan-cancer total."
- **Stage nomenclature/coding varies by cancer type**: not all cancer
  types are staged I-IV (some use grade, some other systems); some code
  only substages (Stage IIA, IIIB...) with few or no samples matching
  the bare "Stage I"/"Stage IV" string. ALWAYS check n for both compared
  values WITHIN the target cancer type before trusting a result -- do not
  assume a field/value pair that works pan-cancer, or for one cancer
  type, carries over to another. Example: BRCA Stage IV bare-string n=22
  (too thin to trust) despite BRCA otherwise being a well-populated
  project -- likely because many BRCA Stage IV cases are coded as
  substages instead.
  - **Confirmed more extreme case: PRAD (prostate) doesn't use overall
    stage AT ALL.** Checked live: `cgc_case_pathologic_stage` returns 0
    samples for EVERY value (Stage I through Stage IV, including every
    substage) when constrained to TCGA-PRAD -- not thin, genuinely zero.
    Prostate cancer is conventionally staged by pathologic T-category
    instead; `cgc_case_pathologic_t` is confirmed working and populated
    for PRAD (T2a n=14, T2b n=11, T2c n=198, T3a n=175, T3b n=140, T4
    n=13). Lesson generalizes beyond "check n" to "check the field
    returns ANY data at all for this cancer type" -- some cancer types
    don't use the field you'd naively reach for, not just thinly.

## Gotchas that look like bugs but aren't

- `regions=<GENE_SYMBOL>` resolves to that gene's coordinate span and
  returns **every** gene/junction/exon overlapping that span, not just
  the named gene (e.g. `regions=TP53` on `/genes` returns 3 rows: TP53
  + 2 neighboring genes). Must filter client-side on the gene symbol
  inside the packed id column.
- `/genes` and `/exons` pack 4 values into one column:
  `gene_id:gene_name:gene_type:bp_length` -- split on `:` client-side.
- `/exons` header has three literal `NA` column names (inherited
  placeholder slots from the junction schema) -- always empty, ignore.
- `left_annotated`/`right_annotated` on `/snaptron` (junctions) are
  NOT booleans -- each is either `"0"` or a comma-separated list of
  annotation-source codes (e.g. `gC19` = Gencode v19). This lets you
  distinguish one-end-novel vs both-ends-novel vs "known sites used in
  a novel combination," which is richer than a flat annotated==0 check.

## Normalization -- read before making any tumor/normal or tissue claim

- `coverage_sum` is a raw sum -- dominated by cohort size (more samples
  in a group -> bigger sum, regardless of biology). Do not compare raw
  sums across groups of different size (e.g. tcgav2 has ~9-13x more
  Primary Tumor samples than Solid Tissue Normal for most genes).
- `coverage_avg`/`coverage_median` divide by sample count, removing the
  cohort-size confound. These are NOT depth/library-size normalized
  (no AUC scaling applied) -- a large gap (>10x) is safely biological, a
  close call between groups may partly reflect sequencing depth
  differences. State this caveat in any output rather than overclaiming.
- No endpoint found that returns total cohort size per sample-type
  value directly. A whole-gene query's own samples_count is a usable
  proxy for the region of interest.
- MEDIAN vs MEAN behaviour (observed live on LRRC10 across GTEx tissues):
  non-expressing tissues floor at an identical `coverage_median` ~152,
  almost certainly ~one read's worth of coverage over the region. So
  `coverage_median` is a clean ON/OFF discriminator (expressed vs silent)
  but loses all resolution among low tissues (they all read ~152);
  `coverage_avg` keeps gradation but is outlier-sensitive (e.g. a few
  high samples inflate it). Report both. Expect a lowly-expressed gene
  to floor in BOTH groups of a tumor/normal comparison.

## Confirmed-viable MVP tool scope (recon-tested end to end)

1. **Gene-level tumor vs normal** (`/genes` + `sfilter` on
   `cgc_sample_sample_type`) -- most robust, coverage_avg/median already
   per-sample-normalized.
2. **Junction differential + novelty** (`/snaptron` + `sfilter` for
   tumor/normal, `rfilter=annotated:0` for novel junctions, and
   left/right_annotated for one-end vs both-ends novel).
3. **Tissue specificity across GTEx** (`/genes` on `gtexv2` + `sfilter`
   on `SMTS`, looping the full 31-tissue panel -- see SMTS section above
   -- ranking by coverage_avg/median). Validated on two contrasting
   cases: LRRC10 (cardiac-restricted) -- Heart 82,099 median vs
   everything else floored at ~150-152 (one-read noise floor), a clean
   isolated signal. ADAR (ubiquitously expressed) -- NO floor anywhere,
   every tissue in the hundred-thousands-to-millions range, top tissue
   (Blood, 1,906,809 median, post Blood-Vessel-collision fix) only
   ~10-20% above the next few (Spleen, Pituitary, Uterus, Lung) -- a
   modest top-of-a-continuum result, not a tissue-restriction claim.
   Confirms the tool handles both a sharply-restricted gene and a
   broadly-expressed one correctly. Full 31-tissue panel costs ~90s per
   call (~2.9s/tissue, linear, no batching) vs ~25-30s for the earlier
   10-tissue version -- budget for it.

Exon-level tumor/normal and any `/bases`-based region coverage are NOT
confirmed -- treat as future work / open questions, not MVP scope,
unless re-tested live at the event with more time.

4. **Stratified comparison, generic** (`/genes` + `sfilter` on any
   confirmed clinical field, e.g. `cgc_case_pathologic_stage`, optionally
   constrained by `gdc_cases.project.project_id` to one cancer type --
   see "Confirmed clinical/stratification fields" above). Cohort totals
   for arbitrary field/values are NOT hardcoded like tumor/normal --
   derive live per value via a housekeeping-gene query, same method used
   to discover the original tumor/normal totals. Validated on MKI67
   Stage I vs Stage IV: pan-cancer pooled showed ~4.5x, shrank to ~1.3x
   within KIRC alone (well-powered, n=300/104) -- the pan-cancer number
   was substantially inflated by cancer-type composition, not a clean
   stage effect. BRCA was tried first but rejected for the headline
   number (Stage IV n=22, too thin to trust) -- always check n within the
   target cancer type before trusting a stratified result, per the
   nomenclature caveat above.

## Speed (single sequential requests, not load-tested)

Typical query: ~1.5-4s regardless of region size (no fast path for
small regions; filtering happens server-side after full pull, not as
pushdown). `/bases` scales badly with region length (avoid, see above).
An agent looping many regions serially will accumulate real wall-clock
time -- budget for it, don't assume sub-second responses.

## Tool design principle for the MCP server itself

Don't expose a generic passthrough (raw compilation/field/filter
params) and expect the calling agent to pick the right catalogue --
that pushes two days of debugging onto every future call, and it's how
the SMTSPAX mistake above happened. Instead, bake the compilation and
field choice into distinct, named, well-described tools
(`compare_tumor_vs_normal`, `compare_across_tissues`, etc.) so the tool
name itself answers "which catalogue" before the agent has to guess.
Keep one clearly-labeled fallback/exploratory tool for anything not
covered, rather than making the curated tools also do double duty.
