# Label Spec v1.0 — TCAB1/ALT Osteosarcoma Project

**Status:** Frozen 2026-05-15 after feasibility audit (notebooks 04 + 05).
Changes after this version require a new version number and changelog entry.

## Audit verdict (binding)

The feasibility audit (Pass A and Pass B) recommends **pragmatic as the primary
definition set** for both aims, with the following actual N per group:

| Aim | Definition | Disrupted | Intact | Verdict |
|---|---|---|---|---|
| Aim 1 primary | Pragmatic | 8 | 32 | MARGINAL (below strict ≥10; clears ≥8 floor) |
| Aim 1 TP53-restricted | Pragmatic | — | 12 total | GO (clears ≥8 floor) |
| Aim 2 (signature) | Pragmatic | 48 pos | 11 neg | Signature scoring only (negatives <20) |

The strict definition was not selected because its disrupted count (8) is identical
to pragmatic but neither clears the ≥10 strict primary threshold. The broad definition
adds only 1 disrupted sample (9 total) and introduces 16 uncertain DepMap cell lines
that contribute to neither arm of the primary contrast.

## 1. TCAB1 Status (unchanged from v0.9)

### 1a. Genomic status (`tcab1_genomic_status`)
[Copy exactly from v0.9 Section 1a — no changes]

### 1b. Functional status (`tcab1_functional_status`)
[Copy exactly from v0.9 Section 1b — no changes]

### 1c. Primary analysis labels (updated for v1.0)

For the Aim 1 primary contrast, derive:
- `tcab1_overall_status = candidate_disrupted`: `disrupted_high_conf` OR architecture
  call `Both-single-SV` / `Both-independent-SVs` / `TCAB1-only` (any WRAP53 hit),
  AND not contradicted by `functional_intact_supported`.
- `tcab1_overall_status = candidate_intact`: architecture call `Neither` / `TP53-only`
  (WRAP53 not hit), OR DepMap cell line with `functional_intact_supported` and no
  genomic disruption evidence, AND not contradicted by `functional_loss_supported`.

All other combinations: `tcab1_overall_status = unknown`, excluded from primary contrast.

## 2. ALT Status (unchanged from v0.9)

[Copy from v0.9 Section 2 — no changes]

**Aim 2 implementation in v1.0:** With N=48 positive / N=11 negative, train a
regularized signature score (elastic net or LASSO logistic regression with strong
regularization) rather than a binary classifier with confidence intervals. Validation
by stratified leave-one-cohort-out cross-validation. Report AUROC and AUPRC with
confidence intervals from CV folds. Do NOT report calibrated probabilities at small N.

## 3. Required Auxiliary Columns (unchanged from v0.9)

[Copy from v0.9 Section 3 — no changes]

## 4. Analysis Design (v1.0 — frozen)

### Primary contrast (Aim 1)
`candidate_disrupted` (N=8) vs `candidate_intact` (N=32), pragmatic definition.
Method: DESeq2 with TP53 status and cohort as covariates.
Multiple testing: BH-corrected, FDR < 0.05.
Effect size threshold: |log2FC| > 1.0 for headline gene list; report full ranked
list as supplementary.

### Required sensitivity analyses (run alongside primary, report both)
- **Peer-reviewed labels only:** drops disrupted from 8 → 4. Report the gene
  intersection between peer-reviewed-only and full-cohort results.
- **TP53-restricted (secondary GO contrast):** primary contrast within TP53-altered
  samples (N=12 total). This controls for the strongest confounder in osteosarcoma.
  If the TCAB1 effect disappears here, it is a TP53-SV effect — report as a
  meaningful negative result.
- **Cohort-stratified:** report effect estimates by cohort to assess batch effects.
- **Broad definition exploratory:** apply primary contrast to broad universe
  (N=9 disrupted + N=32 intact), report alongside primary.

### Aim 2 (signature score)
N=48 alt_gold_positive + N=11 alt_gold_negative across 4 cohorts.
Method: elastic net logistic regression with cohort-stratified leave-one-out CV.
Report performance with 95% CI from CV folds.
Apply trained signature score to DepMap OS cell lines as exploratory output
(no labels available there).

### Aim 3 (single-cell projection)
Deferred. Plan after Aim 1 + Aim 2 signatures are characterized.

## 5. Go/No-Go Thresholds (record of v0.9 thresholds; not adjusted post-hoc)

[Copy from v0.9 Section 6 verbatim — these were set before the audit and remain
as the criteria. The audit verdict above documents that the project clears the
floor (≥8 per group) but not the preferred minimum (≥10 per group).]

## 6. Data Limitations (binding for paper Methods section)

These limitations must be reported in any deliverable produced from this project:

1. **DepMap 26Q1 lacks structural variant calls.** Architecture calls for DepMap
   cell lines are derived from copy-number segments and treated as `unknown`
   for primary architecture analysis. Three CN-derived `disrupted_high_conf`
   calls (Saos-2, NRH-OS1, NRH-OS4) are documented with caveats in
   `data/manifests/master_manifest.csv` notes column.

2. **CCLE 2019 SvABA data has minimal overlap with modern DepMap OS cohort.**
   Integration attempted (notebook 02); only 1 of 35 DepMap OS lines overlaps,
   yielding zero locus-window SVs.

3. **Half of Aim 1 disrupted signal comes from preprint sources** (Lawlor 2026
   PDX + Taglialatela 2026). Peer-reviewed sensitivity analysis is mandatory
   in all Aim 1 reporting.

4. **Aim 2 negative class is small (N=11).** Signature scoring with cross-validation
   is appropriate; binary classification with calibrated probabilities is not.

## 7. Changelog

| Version | Date | Notes |
|---|---|---|
| v0.9 | 2026-04-29 | Initial spec for feasibility audit |
| v1.0 | 2026-05-15 | Frozen after feasibility audit. Pragmatic definition selected as primary. Aim 1 marginal (N=8 vs 32), Aim 2 signature-scoring (N=48 vs 11), TP53-restricted secondary GO (N=12). |