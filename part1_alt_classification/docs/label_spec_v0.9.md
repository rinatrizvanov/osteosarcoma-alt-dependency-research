# Label Spec v0.9 — TCAB1/ALT Osteosarcoma Project

**Status:** Provisional. All thresholds and category boundaries are subject to revision after the
feasibility audit. Frozen as v1.0 only after sample-size and confounding tables are reviewed.

**Versioning rule:** Any change after v1.0 freeze requires a new version number and a written
rationale in the changelog. Thresholds may not be adjusted after differential expression results
are seen.

---

## 1. TCAB1 Status

### 1a. Genomic status (`tcab1_genomic_status`)

| Category | Criteria |
|---|---|
| `disrupted_high_conf` | Breakpoint inside WRAP53 gene body (coding exons or promoter/5′ UTR), OR focal deletion removing ≥1 coding exon, OR rearrangement whose reconstructed derivative chromosome predicts a truncated or fusion WRAP53 transcript. Biallelic or hemizygous-on-deleted-background preferred when allele information is available. |
| `candidate_regulatory_disruption` | SV within 5 kb of WRAP53 gene body without direct coding-sequence damage. Reserved for sensitivity analyses only. |
| `intact` | No focal deletion affecting WRAP53 coding sequence, no breakpoint involving the WRAP53 locus, no evidence of transcript truncation or fusion. Broad chromosome-level aneuploidy without focal involvement does not exclude this category. |
| `unknown` | Distal SVs of unclear effect, low WGS coverage at the locus, or assay not performed. |

### 1b. Functional status (`tcab1_functional_status`)

| Category | Criteria |
|---|---|
| `functional_loss_supported` | Low TCAB1 RNA (bottom quartile within cohort) AND at least one corroborating evidence: low protein, promoter hypermethylation, or genomic disruption. |
| `functional_loss_possible` | Low TCAB1 RNA without corroborating evidence. |
| `functional_intact_supported` | RNA in upper three quartiles AND no evidence suggestive of loss. |
| `unknown` | RNA not available, or conflicting evidence (e.g., genomic disruption with normal RNA). |

### 1c. Primary analysis labels

For the Aim 1 primary contrast, derive two columns from the genomic and functional status:

- `tcab1_primary_disrupted` = `disrupted_high_conf` AND no contradictory functional evidence
  (i.e., not `functional_intact_supported`).
- `tcab1_primary_intact` = genomic `intact` AND no contradictory functional evidence
  (i.e., not `functional_loss_supported`).

All other combinations excluded from the primary contrast.

Also store a continuous `tcab1_confidence` score (0–1) combining genomic and functional evidence
weights, for sensitivity analyses at varying confidence cutoffs.

---

## 2. ALT Status

| Category | Criteria |
|---|---|
| `alt_gold_positive` | Two independent lines of evidence: C-circle positive + APB positive, OR C-circle positive + telomere FISH ultrabright foci, OR any single ALT assay positive + TRAP-negative / TERT-low. |
| `alt_negative_high_conf` | Strong evidence against ALT: TRAP-positive with TERT expression, OR ALTitude-negative with corroborating telomere features, OR multiple ALT assays negative. |
| `alt_silver` | ALTitude calls on DepMap cell lines (validated computational labels). Treat as near-gold for Aim 2 training when gold counts are insufficient. |
| `alt_proxy` | Single computational call without orthogonal validation (e.g., own WGS-feature calls on TARGET-OS). For exploratory analyses only. |
| `alt_unknown` | Partial evidence, conflicting assays, or no direct assay performed. |

**Aim 2 training preference order:**
1. `alt_gold_positive` vs `alt_negative_high_conf` only (preferred if N ≥ 20 per class).
2. Add `alt_silver` to expand training set if gold-only is too small.
3. `alt_proxy` reserved for external validation, never for training.

---

## 3. Required Auxiliary Columns

Track separately and do not collapse into status labels:

- `tp53_status`: `altered` / `not_altered` / `unknown`
- `telomerase_status`: `positive` / `negative` / `unknown` (based on TRAP, TERT expression,
  ALTitude features)
- `sv_architecture_17p13_1`: `TP53-only` / `TCAB1-only` / `Both-single-SV` /
  `Both-independent-SVs` / `Neither` / `not_applicable`
- `cohort`: source dataset
- `purity_quartile`: when tumor purity is available
- `notes`: free-text reasoning for manually curated calls
- `alt_source`: citation for ALT calls drawn from literature

The SV architecture column directly tests the TCAB1 preprint's central model that a single
structural event co-disables TP53 and TCAB1; it cannot be reconstructed from gene-level status
alone.

---

## 4. Analysis Design

### Primary contrast
`tcab1_primary_disrupted` vs `tcab1_primary_intact`, all eligible RNA-seq samples.

### Secondary contrasts
- `tcab1_primary_disrupted` vs `tcab1_primary_intact`, telomerase-negative samples only.
- `functional_loss_supported` vs `functional_intact_supported`.
- Primary contrast within `alt_gold_positive` samples only.

### Required sensitivity analyses
- Primary contrast restricted to TP53-altered samples only. Holds TP53 status constant. If the
  TCAB1 effect disappears here, the apparent "TCAB1 effect" is a TP53-SV effect — this is a
  meaningful negative result, not a project failure.
- Primary contrast at varying `tcab1_confidence` thresholds (e.g., > 0.7, > 0.8, > 0.9).
- Primary contrast including `candidate_regulatory_disruption` cases as disrupted.
- Cohort-stratified analysis to assess batch confounding.

---

## 5. Parallel Definition Sets for the Feasibility Audit

Pre-specify three definition sets and apply all three before choosing the primary:

| Set | TCAB1 disrupted | ALT positive |
|---|---|---|
| Strict | `disrupted_high_conf` only, no contradictory functional evidence | `alt_gold_positive` only |
| Pragmatic | `disrupted_high_conf` + `candidate_regulatory_disruption` with supporting low RNA | `alt_gold_positive` + `alt_silver` |
| Broad | Any SV in WRAP53 locus OR bottom-quartile RNA | All non-`alt_unknown` positive calls including proxy |

**Selection rule:** the primary set is the strictest definition that clears the go/no-go thresholds
(Section 6). The remaining two become sensitivity analyses, reported alongside the primary.

This rule is fixed before counts are seen and prevents post-hoc threshold manipulation.

---

## 6. Go/No-Go Thresholds

*Set before the audit, applied after.*

### Aim 1 primary contrast
- Minimum: **10 per group** after all inclusion filters.
- Preferred: 15–20 per group.
- Below 8 per group: do not freeze that definition as primary.

### Aim 1 TP53-restricted sensitivity
- Minimum: **8 per group**.
- Below 8: keep as exploratory/descriptive.

### Aim 1 telomerase-negative secondary
- Below 8 per group: report as exploratory only.
- 8–12 per group: secondary analysis.
- 15+ per group: full secondary analysis with formal statistics.

### Aim 2 classifier
- Minimum: **20 per class** for unpenalized models.
- 15 per class acceptable with heavy regularization (LASSO, nested CV) and feature count ≤ 10.
- Below 15 per class: present as signature scoring with cross-validation, not "classifier."

### Cohort balance (Pass B)
- Within the chosen primary definition, no single cohort may contribute >70% of either class.
- If a cohort dominates, either stratify by cohort in all analyses or restrict to balanced cohorts.

---

## 7. Audit Workflow (Order of Operations)

1. **Build manifest** — auto-fill metadata, RNA/WGS/WES availability, TP53 status, TCAB1 RNA
   levels, telomerase indicators.
2. **Manual curation** — SV architecture, TCAB1 genomic status, gold ALT labels from literature
   supplements, conflict cases.
3. **Pass B first: confounding census** — run cohort × status contingency tables before tuning
   thresholds. If cohort dominates fatally, redesign before continuing.
4. **Apply strict / pragmatic / broad in parallel** — generate counts under all three definition
   sets.
5. **Pass A: availability census** — counts for each aim under each definition.
6. **Decision** — apply Section 6 thresholds and Section 5 selection rule.
7. **Freeze v1.0** — write final spec, tag repo, begin DE pipeline.

*Pass B precedes Pass A because fatal confounding cannot be fixed by tuning labels.*

---

## 8. Inclusion / Exclusion Rules

**Include:**
- Osteosarcoma samples (primary, metastatic, recurrent — track separately).
- DepMap osteosarcoma cell lines for Aim 2 (with cohort flag for batch handling).

**Exclude from primary analyses:**
- Samples with `tcab1_genomic_status = unknown` AND `tcab1_functional_status = unknown`.
- Samples with `alt_unknown` for any analysis requiring ALT status.
- Samples with RNA-seq quality failures (define per cohort during QC).

**Track but flag:**
- Pediatric vs adult osteosarcoma (potential effect modifier).
- Pre- vs post-chemotherapy biopsies (treatment can alter expression).

---

## 9. What This Spec Does Not Yet Decide

*Deliberately deferred until after the audit:*

- Final RNA normalization strategy across cohorts (depends on which cohorts pool successfully).
- Final feature set for the Aim 2 classifier (depends on N and feature stability).
- Whether to include cell lines in Aim 1 DE or restrict to tumors.
- Single-cell projection method for Aim 3 (depends on signature characteristics from Aim 1).

---

## 10. Changelog

| Version | Notes |
|---|---|
| v0.9 | Initial spec for feasibility audit. |
| v1.0 | To be created after audit. Will document which definition set was selected as primary, the actual N per group, and any adjustments to thresholds or auxiliary columns. |
