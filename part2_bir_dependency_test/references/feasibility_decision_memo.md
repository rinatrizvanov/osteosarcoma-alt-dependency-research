# Feasibility Decision Memo — TCAB1/ALT Osteosarcoma Project

**Author:** Rinat Rizvanov
**Date:** 2026-05-17
**Status:** Decision required (advisor input)

## Finding

The feasibility audit (notebooks 04–05) and RNA-availability diagnostics
(notebooks 06–07) establish that **neither aim is estimable from the data
currently in `data/raw/`**:

- **Aim 1** (TCAB1-disrupted vs intact, differential expression): 2 disrupted
  samples have file-backed RNA, 0 intact do. NOT ESTIMABLE.
- **Aim 2** (ALT classifier): 0 ALT-labeled samples have file-backed RNA. NOT
  ESTIMABLE.

## Root cause

Every integrated cohort is **label-only in its public release**: de Nonneville
(ALT/methylation table), Lawlor (SV/IHC), Chen 2014 PCGP (SV/FISH),
Taglialatela (cell-line ALT list). The only expression matrices on disk
(DepMap, 1,719 cell lines; TARGET-OS, 88 tumors) carry no architecture or ALT
labels. Labels and expression live in disjoint datasets.

## This is a designed catch, not a pipeline failure

The Stage 0 diagnostics (06, 07) are read-only checks built specifically to
detect this before any analysis was written. They prevented an artifactual
result. The pipeline behaved correctly.

## Recovery routes (not yet exhausted)

| Route | Source | Effort | Dependency |
|---|---|---|---|
| A. de Nonneville expression | ArrayExpress E-MTAB-11743 (public) | Medium | None (but verify assay type — RNA vs methylation) |
| B. Taglialatela/ALTitude → DepMap RNA | Cell-line-name → ModelID join, data already on disk | Low | None |
| C. Lawlor PDX RNA | Internal lab data | Low if it exists | Internal access (PI/postdoc) |

Notebook 07 only matched exact sample IDs; it did not test name-based DepMap
mapping (Route B) or external public deposits (Route A). Notebook 08 will test
these before any decision is final.

## Recommendation

1. Pursue Route B immediately (no external dependency; notebook 08).
2. Download and inspect E-MTAB-11743 to confirm Route A's assay type.
3. Request Lawlor internal RNA (Route C) — the highest-value route if available.
4. If all routes fail: pivot to a methods/feasibility paper documenting the
   reusable multi-source labeling pipeline and the public-data gap. This remains
   a publishable, portfolio-worthy contribution.

## Decision needed from advisor

- Can the lab share RNA-seq for the Lawlor 2026 GCCRI PDX models?
- Is reframing Aim 1 around an expression-proxy contrast acceptable if
  architecture-confirmed disruption is not recoverable?

  ## Update (notebook 08): Route B partially rescues Aim 2

Notebook 08 found 37 osteosarcoma DepMap cell lines with file-backed RNA,
of which ~31 have ALTitude curated ALT calls not yet merged into the manifest.
ALTitude DepMap calls are defined as `alt_silver` in label_spec_v1.0 and are
the spec-sanctioned fallback when gold ALT counts are insufficient (they are).

- **Aim 2:** likely VIABLE as a single-platform DepMap signature-score analysis
  using ALTitude `alt_silver` labels. Class balance pending notebook 09.
- **Aim 1:** still NOT ESTIMABLE via Route B. Requires Route A/C or a reframe
  to functional/expression-proxy TCAB1 status.

Revised recommendation: pursue Aim 2 via ALTitude-DepMap as the primary
runnable analysis; keep Aim 1 dependent on Route C (Lawlor RNA) or reframed.

## Route A resolved (manual ArrayExpress check, 2026-05-18)

E-MTAB-11743 is an Affymetrix gene-expression MICROARRAY (study type
"transcription profiling by array", array design A-GEOD-19803), NOT RNA-seq
and NOT methylation. 12 samples deposited (manifest has 22 de Nonneville
ALT-labeled samples — only ≤12 are recoverable, mapping via SDRF needed).

Implication: microarray expression cannot be pooled with DepMap RNA-seq into
one matrix (incompatible measurement technologies). Route A is therefore NOT
a primary-cohort expansion. Its legitimate uses:
  (a) independent cross-platform validation of a DepMap-derived ALT signature
      (primary tumor, orthogonal platform — strengthens generalizability claim);
  (b) small standalone within-platform replication sub-analysis.

Route A is downgraded from "potential primary cohort" to "external validation
cohort, pending notebook 09 producing a signature to validate."

## Update (notebook 09 corrected, 2026-05-19)

Notebook 09's first verdict (NOT ESTIMABLE) was a column-binding bug — the
code read `Curated.Marker` (ATRX-mutation feature) instead of `ALT_pred_class`
(ALTitude's ALT classification). Diagnostic confirmed the correct field, and
sanity checks against literature ground truth (U2OS=ALT+, SJSA1=ALT-, etc.)
pass.

Corrected Aim 2 inventory (31 RNA-backed OS DepMap cell lines):

- Block 1 (gold-only, experimental): 12 ALT+ vs 7 ALT-, both arms with
  published C-circle/APB/TRAP evidence. Within-cohort viable.
- Block 2 (gold+silver, sensitivity): 18 ALT+ vs 7 ALT-, adding 6 ALTitude
  prediction-only ALT+ lines. Silver inflates positive class only.
- Non-OS flagged and excluded: TM31 (glioblastoma). SKLU1 mapped but not in
  DepMap RNA matrix.

VERDICT: Aim 2 WITHIN-COHORT VIABLE on Block 1 (primary) with Block 2 as
sensitivity. Single-platform analysis (DepMap log2(TPM+1)), no cross-cohort
confounding. Methodology per label_spec_v1.0 §2: regularized signature score
with leave-one-out CV, not a calibrated binary classifier.

Aim 1 status unchanged: NOT ESTIMABLE without Route C (Lawlor PDX RNA) or
reframe to a TCAB1-expression-quantile contrast within DepMap RNA.

## Update (manifest direction-aware, 2026-05-20)

Manifest now carries `alt_call_direction` (positive/negative/unknown), populated
for all informative ALT rows. Silver lines split: 6 positive / 5 negative.
Audit re-run (notebooks 04/05) with corrected manifest:
  [paste Pass B and Pass A verdict summaries]

Aim 2 final modeling design (notebook 11):
- Primary strict contrast: DepMap_26Q1, RNA-backed, gold labels only.
  12 ALT+ vs 8 ALT- (20 samples, single platform, all peer-reviewed).
- Pragmatic sensitivity: same cohort + RNA filter, gold + silver labels.
  18 ALT+ vs 13 ALT- (31 samples).
- Arm split keyed on `alt_call_direction`, NOT `alt_status`.
- de Nonneville microarray (E-MTAB-11743): EXTERNAL VALIDATION ONLY,
  not pooled into training (different measurement technology).

Aim 1 status unchanged: NOT ESTIMABLE without Route C (Lawlor PDX RNA).

## Audit re-run with direction-aware counts (2026-05-30, kernel-restart-clean)

After the manifest direction fix (notebook 10) and the audit code fix to key
Aim 2 arm counts on `alt_call_direction` instead of `alt_status`:

| Definition | Aim2 + | Aim2 - | Aim2 verdict          | Aim1 primary | Aim1 TP53 |
|------------|--------|--------|------------------------|--------------|-----------|
| strict     | 60     | 19     | GO (penalized only)    | MARGINAL     | GO        |
| pragmatic  | 66     | 24     | GO (unpenalized)       | MARGINAL     | GO        |
| broad      | 66     | 24     | GO (unpenalized)       | MARGINAL     | GO        |

Peer-reviewed sensitivity (excluding Lawlor 2026 preprint):
| Definition | Aim2 + | Aim2 - | Aim1 disrupted | Aim1 verdict |
|------------|--------|--------|----------------|--------------|
| strict     | 39     | 16     | 4              | NO-GO        |
| pragmatic  | 45     | 21     | 4              | NO-GO        |
| broad      | 45     | 21     | 5              | NO-GO        |

Aim 2 modeling design (notebook 11):
- Primary: strict, gold-only, DepMap-only single-platform = 12 ALT+ / 8 ALT-.
- Pragmatic sensitivity: + silver, DepMap-only = 18 ALT+ / 13 ALT-.
- Arm split keyed on `alt_call_direction`, NOT `alt_status`.
- Pass B `cohort × alt_call_direction` tables: NO dominance warnings under any
  definition; DepMap is best-balanced cohort (8 neg / 12 pos under strict).
- Aim 2 cross-cohort 60/19 (strict) and 66/24 (pragmatic) are documented as
  the audit's potential universe; the analytical universe is DepMap-only to
  preserve single-platform cleanliness.
- de Nonneville microarray (E-MTAB-11743): EXTERNAL VALIDATION ONLY, not
  pooled into training.

Aim 1 status: MARGINAL primary, GO TP53-restricted secondary, NOT ESTIMABLE
under Stage 0 RNA-backing diagnostic. Pending Lawlor PDX RNA via PI request
(Route C) or deliberate quantile-based reframe.

Aim 2 modeling proceeds in notebook 11.