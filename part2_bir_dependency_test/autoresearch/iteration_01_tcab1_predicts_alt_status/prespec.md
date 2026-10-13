# Iteration 01 pre-specification — FROZEN

**Hypothesis H1.** TCAB1/WRAP53 loss predicts ALT **status**, not dependency within ALT.
**Status:** frozen before any exposure–outcome join was performed.
**Supersedes:** `prespec.DRAFT.md` (committed `07220fd`, before any scan report existed).
**Disclosures:** `DISCLOSURES.md` (committed before this file).
**Scans informing it:** `scan_exposure.md`, `scan_metadata.md`, `scan_literature.md`.

---

## 1. Hypothesis and directional prediction

**H1.** 17p13.1/TCAB1 disruption is an ALT-**initiating** lesion: it predicts
whether an osteosarcoma line is ALT+ at all, and carries little or no
information about which recombination machinery an already-ALT+ line depends on.

**Prediction.** ALT-positive osteosarcoma lines have **lower** WRAP53
log2(TPM+1) than ALT-negative lines. With δ = Cliff's delta for
(ALT+ vs ALT−) on WRAP53 expression, the prediction is **δ < 0**, so
`predicted_sign = -1`.

**Mechanistic basis.** Keegan et al. (2026) *BMC Cancer*,
https://doi.org/10.1186/s12885-026-16250-5: SVs in TP53 intron 1 inactivate
TP53 and TCAB1 together via the head-to-head 1,419-nt overlap at 17p13.1,
present in "approximately 40% of ALT positive osteosarcoma tumors", framed as
"an **early** and previously unrecognized event contributing to the
**activation** of the ALT pathway." In-text figure: SVs in TCAB1 in ~35% of
ALT+ samples "and none of the ALT negative samples."
`Literature/INDEX.md` Prespec Conflict **#8**.

**Not Stage 1 renamed.** Stage 1 asked, within ALT+ lines, whether WRAP53
expression modulates a BIR-minus-HR Chronos composite. H1 asks, across ALT+ and
ALT− lines, whether WRAP53 expression discriminates ALT status. Different
cohort, contrast and outcome. **No gene-effect value is read in this iteration;
`src.dependency_stats` is never imported and a gate check asserts so.**

## 2. What the literature scan changed, before freezing

`scan_literature.md` establishes that the supported half of H1 is already
published and that the half under test here is not supported but actively
undercut. Recorded now, so the result cannot be re-framed afterwards:

1. **Keegan's claim is tumour-level and SV/protein-level, not cell-line
   mRNA-level.** Two cell lines appear in the paper, both as controls (SJSA1,
   HEK293). The only ALT+/ALT− contrast is 21 TARGET-OS tumours (10+/11−) with
   ALT status from a computational predictor, not a tissue assay.
2. **Keegan's own RNA-seq contradicts the mRNA proxy.** For TCAB1-SV samples
   disrupting the α/γ promoters but sparing β, "TCAB1 RNA expression levels were
   **similar to normal controls**" (Supplemental Fig. 3). The lesion need not
   lower the transcript.
3. **The two canonical ALT+ osteosarcoma lines are WRAP53-HIGH.** Mahmoudi et al.
   (2011) https://doi.org/10.1038/cddis.2010.90 Fig. 1a lists U2OS and Saos-2
   among WRAP53-**overexpressing** lines. This makes the counter-direction a
   live, literature-grounded alternative rather than a nuisance.
4. **Low WRAP53 does not force loss of telomerase.** Wang et al. (2022)
   https://doi.org/10.1002/1878-0261.13324 report TCAB1-deleted lines with
   significantly reduced expression but "**no correlation between TCAB1 protein
   levels and telomerase activity**" (HL60: negligible TCAB1, highest activity).
   The causal chain TCAB1-low → telomerase-null → ALT has a broken first link.
5. **Three independent feature-selection exercises already discarded WRAP53.**
   Mentions of WRAP53/TCAB1: Armendáriz-Castillo 28-gene ALT signature **0**,
   Bennett ALTitude **0**, Taglialatela **0**. The 28-gene signature instead
   includes **CCT5** — the chaperonin setting TCAB1 *protein* (Jiang 2025). The
   field captures this axis at the protein layer. This repo's own Aim 2 agrees:
   WRAP53 coefficient exactly 0.000000 in all three cohorts, bootstrap
   frequency 0.17/0.12, LOOCV fold frequency 0.0000
   (`references/notebook_13_report.md`).
6. **TCAB1 loss is a passenger of TP53 SVs and is insufficient alone.** Keegan:
   "loss of TCAB1 may be a consequence of TP53 SVs"; TCAB1 SVs occurred only in
   TP53-SV tumours; TP53 SVs in 50% of ALT+ **and 36% of ALT−**; and "there were
   two ALT negative tumors with SVs in TCAB1." In a uniformly TP53-disrupted
   osteosarcoma panel the upstream exposure is close to non-discriminating by
   construction.

**Consequence accepted in advance:** the prior probability of H1 in its
mRNA form is low, and this iteration is run as a *diagnostic* with a
pre-registered asymmetry (§8), not as a confirmatory test.

## 3. Known label–exposure pairs — 4 of 31

Per `DISCLOSURES.md`:

| # | Lines | How known |
|---|---|---|
| D1 | ACH-002834 (OS384, gold+), ACH-003185 (NRH-OS4, silver+) — both WRAP53 = 0 | Stage 1 prespec §4d, recorded before any Stage 1 outcome was read |
| D2 | ACH-000364 (U2OS, +), ACH-000748 (SJSA1, −) | ground-truth anchors in `src/` docstrings, joinable to the exposure table |

Guards: **R3** drops the D1 pair; **R3b** drops all four. If the primary and
either guard disagree in direction, **the guard is reported as the primary
interpretation.** The main agent has not read the per-line exposure table at the
time of this freeze (`DISCLOSURES.md` D3).

## 4. Cohort — filter trace fixed before the join

From the frozen `references/master_manifest.csv`; ALT labels never re-derived.

| Step | N |
|---|---|
| `cohort == "DepMap_26Q1"` | 35 |
| `alt_call_direction ∈ {positive, negative}` (4 `alt_unknown` dropped) | 31 |
| present in the expression matrix | **31** (the RNA filter costs zero) |
| → arms | **18 ALT+ / 13 ALT−** |

Reproduces Stage 1's 16/9/1 exactly when the CRISPR filter is added, which
validates the trace. CRISPR screening is **not** required here, gaining 6 lines
(+24%) and improving balance from 1.78:1 to 1.38:1: ACH-002178 (NY, gold+),
ACH-002689 (OS525, gold+), ACH-000613 (HOS, high-conf−), ACH-001813 (C242,
silver−), ACH-001817 (C393, silver−), ACH-002779 (NCC-ESOS1-C1, silver−).

**Primary cohort is patient-collapsed: 18 ALT+ / 12 ALT− (N = 30).**
ACH-000613 (HOS) and ACH-001001 (143B) share `PatientID = PT-epgMPc` — 143B is
a HOS derivative — and both fall in the ALT− arm, so 2 of 13 negatives are one
patient. **Rule fixed in advance: retain ACH-001001 (143B), the Stage 1 cohort
member; drop ACH-000613 (HOS), the newly gained line.** The uncollapsed 18v13 is
carried as a pre-specified sensitivity (R7).

Anchors asserted before analysis, aborting on failure: ACH-000364 (U2OS) =
positive, ACH-000748 (SJSA1) = negative.

## 5. Exposure

WRAP53 log2(TPM+1), direct quantification, inheriting Stage 1 notebook 02's
`direct` strandedness branch and its documented limitation that Broad never
states library strandedness. **No new exposure variable is derived**, so
anti-bias rule (b) imposes no additional notebook.

From `scan_exposure.md` (exposure-side only, no ALT join): over the 34
RNA-backed osteosarcoma lines the distribution is **not bimodal** — one
continuous mode spanning 2.041–5.888 plus a discrete **2-line spike at exactly
0.0**, with no low-but-nonzero values. The largest interior gap (2.041 units) is
3.3× the next largest and 1.58 SD. ρ(WRAP53, WRAP53 CN) = +0.554 (p = 8.30e-4);
ρ(WRAP53, TP53) = +0.244 (p = 0.164) — both independently reproducing Stage 1's
reported values.

**Design consequence.** The only naturally separable low group has n = 2, below
`stratify.MIN_CELL_SIZE = 4`; every percentile cut (0.25/0.33/0.50 →
3.078/3.336/3.603) falls inside the dense body and is arbitrary. **The primary
is therefore threshold-free** — it groups by the ALT label and treats WRAP53 as
continuous. No cutpoint is chosen anywhere in the primary. Where a cut is used
descriptively, the inequality is pinned to strict `<` (the tertile cut 3.335580
is an observed value, so `<` gives 11/23 and `<=` gives 12/22).

## 6. Primary analysis

- **Measured variable:** WRAP53 log2(TPM+1), continuous, per line.
- **Grouping:** ALT-positive vs ALT-negative (the frozen label).
- **Test:** two-sided Mann-Whitney U. `method="exact"` if the pooled sample has
  no ties, else `method="asymptotic"`, recording `mwu_method`. **A tie is
  expected** — the two WRAP53 = 0 lines are tied with each other — so the
  asymptotic fallback is the likely path. Inherited from
  `dependency_stats.compare_groups`.
- **Effect size:** Cliff's δ with a 2000-resample bootstrap 95% CI, plus the
  Hodges-Lehmann median shift in log2(TPM+1) units.
- **Discrimination:** AUROC reported as **(1 − δ)/2** — a monotone transform of
  δ, not a second test. "Predicts status" is a discrimination claim, and the
  classifier under test is "**lower** WRAP53 ⇒ ALT-positive", so δ = −1 (every
  ALT+ line below every ALT− line) maps to AUROC = 1. See the §14 erratum.
- **Multiplicity:** none within the iteration (one primary test). The binding
  correction is the campaign-level BH across all iterations' primaries; **the
  success criterion is experiment-level q < 0.1**, not this iteration's raw p.

## 7. Pre-specified secondary and robustness analyses

| # | Check | n | Expectation if H1 true |
|---|---|---|---|
| R1 | Gold/high-confidence labels only (drop `alt_silver`) | 12v8 | same direction |
| R2 | Restrict to the 25 Stage 1 CRISPR-screened lines | 16v9 | same direction; shows the result is not an artefact of cohort expansion |
| R3 | Exclude the D1 pair (ACH-002834, ACH-003185) | 16v12 | same direction — known-information guard |
| R3b | Exclude all four known pairs (D1 + D2) | 15v11 | same direction — known-information guard |
| R4 | WRAP53 gene-level copy number as a strand-agnostic exposure | ≤30 | same direction |
| R5 | TP53 log2(TPM+1) as the exposure. **Descriptive, non-gating** | 18v12 | a 17p13.1 co-deletion architecture predicts a same-direction TP53 effect; a WRAP53-specific effect predicts a weaker one |
| R6 | Leave-one-line-out on the primary | 30 refits | no sign flip |
| R7 | Uncollapsed cohort (restore ACH-000613/HOS) | 18v13 | same direction |
| R8 | **Leave-one-depositor-out** (see §9) | varies | no sign flip |
| R9 | Bottom-tertile enrichment for ALT+, Fisher exact. **Descriptive only** | — | the literature predicts a *mixture*, so any signal sits in the lower tail |

R5 is deliberately non-gating: the two genes overlap head-to-head, so they are
not independent, and gating on TP53 would penalise H1 for a correlation the
locus architecture guarantees.

R9 is descriptive because the natural cut cannot reach significance: with both
WRAP53 = 0 lines in the ALT+ arm, the most extreme attainable 2×2 table gives
Fisher two-sided p = 0.497 (one-sided 0.329 = C(18,2)/C(31,2)). **No
configuration of the natural cut can reach p < 0.05.**

## 8. Criteria — and a pre-registered asymmetry

Evaluated mechanically by `src.verdict.evaluate_iteration` with
`predicted_sign = -1`. The six outcomes partition the (direction, p, |δ|) space
exhaustively, so no result can fall into a gap as Stage 1's did.

Support conditions — all must pass for `supported` rather than `partial_support`:

| Key | Condition |
|---|---|
| `primary` | two-sided p < 0.05 **and** δ ≤ −0.330 |
| `r1_gold_only` | R1 same direction (δ < 0) |
| `r3_drop_known` | R3 same direction (δ < 0) |
| `r3b_drop_all_known` | R3b same direction (δ < 0) |
| `r6_loo` | no sign flip on any leave-one-line-out refit |
| `r8_loDo` | no sign flip on any leave-one-depositor-out refit |
| `discrimination` | AUROC ≥ 0.70, i.e. δ ≤ −0.400 (see §14 erratum) |

**The asymmetry, fixed in advance.** Because power differs sharply by direction
and by what each direction would mean:

- A **non-significant result in the predicted direction** is *uninformative*.
  Per §9 the study has 27% power at Keegan's own 35% prevalence even under
  perfect transcript separation, and ~10% at the realistic attenuated
  magnitude. Such a result must be reported as uninformative about H1, and
  **must not** be offered as evidence against it.
- A **counter-directional result** (ALT+ lines WRAP53-*higher*) at |δ| ≥ 0.435
  *is* detectable and *is* interpretable: it would corroborate the
  measurement-layer objection (§2 items 2–5) and Mahmoudi's observation that
  U2OS and Saos-2 are WRAP53-overexpressing. It would be evidence that the mRNA
  proxy does not carry the lesion, **not** evidence that TCAB1 is unrelated to
  ALT — Keegan's tumour-level SV claim is untouched by this design either way.

A `supported` verdict is **not** a campaign positive until it clears
experiment-wide BH at q < 0.1.

## 9. Power, stated in advance

Computed by exhaustive enumeration of the exact Mann-Whitney null
(`src.verdict.min_detectable_delta`) and by simulation under the literature's
mixture alternative. Minimum detectable |δ| (two-sided exact p < 0.05):

| Cohort | n | MDE \|δ\| |
|---|---|---|
| **Primary, patient-collapsed** | 18v12 | **0.435** |
| R7 uncollapsed | 18v13 | 0.427 |
| R2 Stage 1 lines | 16v9 | 0.486 |
| R1 gold-only | 12v8 | 0.542 |

**The mixture alternative.** The literature predicts a mixture, not a location
shift: a fraction φ of ALT+ lines carry the lesion. If k of n ALT+ lines are
transcript-low and rank below every ALT− line while the rest are exchangeable,
then **δ = −k/n exactly** (verified by simulation). Power at 18v12:

| φ = k/18 | δ | power |
|---|---|---|
| 0.167 (k=3) | −0.167 | 8.7% |
| 0.222 (k=4) | −0.222 | 12.7% |
| 0.333 (k=6) ≈ **Keegan 35%** | −0.333 | **27.4%** |
| 0.389 (k=7) ≈ **Keegan 40%** | −0.389 | **39.1%** |
| 0.444 (k=8) | −0.444 | 53.5% |
| 0.500 (k=9) | −0.500 | 71.0% |

**So:** at the literature's own central prevalence, and assuming *perfect*
transcript separation — the ceiling, which Keegan's own RNA-seq (§2 item 2)
says is not met — this design has **27% power**. At the realistic attenuated
magnitude (|δ| < 0.2) it has **under 13%**. This is recorded before the result
is seen, and it is the reason §8 treats the two directions asymmetrically.

*Correction recorded:* an earlier working note in this session asserted the
design was "futile" on the grounds that the literature's point ceiling (0.350)
sits below the detection floor (0.435). That framing was wrong — δ is a random
variable, so a k = 6 mixture still clears the floor 27% of the time. The design
is badly underpowered, not futile.

## 10. Minimum-cell guard

If either arm has n < 4 after §4, H1 is declared not estimable, logged to
`NULL_SPACE.md` as `not-estimable` with the numbers, and the iteration closes.
Inherited from `stratify.MIN_CELL_SIZE = 4`. Not expected to fire at 18v12.

## 11. Confounders that cannot be controlled — declared, not adjusted

From `scan_metadata.md`. These are stated now because none can be fixed later:

1. **Depositor clusters with arm.** UCSF/Sweet-Cordero 5+:1−; Oslo NRH 3+:1−;
   Montefiore/Gorlick 0+:3−. `SampleCollectionSite = soft_tissue`,
   `GrowthPattern = Mixed` and the media IDs are **source aliases for the same
   lines**, not independent axes. 3 of the 4 gained negatives come from two
   investigator blocks, so expanding the cohort increases institutional
   imbalance of the negative arm. **R8 (leave-one-depositor-out) is the guard,
   and it is a support condition**, not an afterthought.
2. **Batch confounding is unfalsifiable with the data on hand.** `Model.csv` has
   no date, era or onboarding column, and `data/raw/depmap/` has no
   `OmicsProfiles.csv` — so platform, library prep, read length, pipeline
   version and profiling date are unknowable in this repo. For an
   RNA-**expression** contrast that is precisely the confounder one would want
   to exclude. **A significant WRAP53 difference in either direction cannot be
   distinguished from a batch artefact**, and the memo must say so whatever the
   result.
3. **Prior treatment, stage and grade cannot be adjusted for or even
   described** — 14 `Model.csv` columns are 100% missing, including `Stage`,
   `PatientTumorGrade` and `PatientTreatmentType`; `PatientTreatmentStatus` is
   85.7% missing.
4. **All 31 labels trace to one upstream source** (ALTitude 2024: 19 curated +
   11 prediction-only + 1 low-signal C-circle cross-check). Tier is a perfect
   function of provenance, so `alt_status` and `source_confidence` are
   redundant and **no independent label replication exists**. R1 therefore tests
   label tier, not label independence.
5. **`label_spec_v1.0.md` §2 is a placeholder** — its body reads
   `[Copy from v0.9 Section 2 — no changes]` and v0.9 is absent from the
   repository. Tier definitions are reconstructed from
   `references/label_rules.py`. The authoritative written definition of the ALT
   tiers does not exist in this repo.
6. **Forbidden variables.** `sv_architecture_auto` must not be used: it is
   artefactual, 33 of 35 rows reading "Both-single-SV" from a single CN segment
   spanning both WRAP53 and TP53 (`label_spec_v1.0.md` §6.1).
   `telomerase_status` must not be used as a control: it near-separates against
   the arm (8/13 negatives telomerase-positive, 0/18 positives), making it a
   label proxy rather than a covariate.

## 12. Rule (c) — mandatory subtype secondary does not bind

Anti-bias rule (c) binds only when a **composite** endpoint is used. This
iteration's endpoint is a single continuous variable, so it does not bind.

Recorded for later iterations: **no defensible type I vs type II ALT marker
exists in DepMap 26Q1.** The paper that profiles the two subtypes directly
(Lee et al. 2026 *Mol Cells*, https://doi.org/10.1016/j.mocell.2026.100405)
could not be retrieved — `Literature/INDEX.md` lists it as the highest-value
gap, supplementary files only.

## 13. Two unretrieved papers, noted not chased

`scan_literature.md` flags two papers that would sharpen this iteration and are
not in the corpus. Per the loop's missing-data rule they are **noted, not
pursued**:

- https://doi.org/10.1016/j.cancergen.2026.06.001 — "Beyond ATRX: the
  heterogeneous basis of ALT **in osteosarcoma**" (2026). The corpus has no
  OS-specific ALT-heterogeneity source.
- https://doi.org/10.1080/15384101.2025.2609646 — "TCAB1: a key promoter of
  tumour growth through telomere maintenance and senescence evasion" (2026).
  Directly on the exposure.

## 14. Changelog

| Version | Date | Notes |
|---|---|---|
| draft | 2026-10-12 | `prespec.DRAFT.md`, committed `07220fd` before any scan report existed |
| **v1.0.1 erratum** | 2026-10-12 | **Pre-analysis correction of a transcription error, committed before any exposure–outcome join.** v1.0 §6 wrote the discrimination metric as AUROC = (δ + 1)/2. With δ defined in this repo as P(ALT+ > ALT−) − P(ALT+ < ALT−) (`dependency_stats.cliffs_delta`) and the prediction being *lower* WRAP53 in ALT+, that formula is inverted: it maps the ideal case δ = −1 to AUROC = 0. It also made the `discrimination` support condition (AUROC ≥ 0.70) require δ ≥ +0.400, which directly contradicts the `primary` condition (δ ≤ −0.330) and rendered the `supported` verdict **unreachable**. The correct transform is AUROC = (1 − δ)/2, under which AUROC ≥ 0.70 ⟺ δ ≤ −0.400 — consistent with `primary`. Nothing else changed: not the test, the direction, the cohort, the criteria thresholds or the power analysis. The alternative — running as frozen and reporting a mechanically unsatisfiable condition as a "failure" — would have been worse than a documented fix. |
| **v1.0 frozen** | 2026-10-12 | Scan-informed. Primary moved to the patient-collapsed 18v12; threshold-free primary justified by the exposure's shape; R3b, R7, R8, R9 added; §8 asymmetry, §9 power table and §11 uncontrollable confounders recorded. Frozen before any exposure–outcome join. |
