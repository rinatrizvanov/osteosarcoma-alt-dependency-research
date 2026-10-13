# Notebook 06 — Aim 1 RNA Inventory

_Generated: 2026-05-17T19:45:58Z_


Read-only diagnostic. Determines whether the Aim 1 primary contrast (pragmatic `tcab1_overall_status` — `candidate_disrupted` vs `candidate_intact`) is statistically estimable from the RNA available on disk, or whether TCAB1 status is confounded with cohort / platform / sample type.


## Inputs

- Manifest: `data/manifests/master_manifest.csv` (115 rows)
- DepMap RNA: `data/raw/depmap/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv` (1719 unique ModelIDs)
- TARGET-OS TPM: `data/raw/target_os/TARGET-OS.star_tpm.tsv.gz` (88 sample columns; not joined to manifest as a cohort)

## Primary contrast composition

Total contrast rows: **40** (disrupted=8, intact=32).

Arm × cohort:

```
arm                disrupted  intact
cohort                              
DepMap_26Q1                2       0
Lawlor_2026_GCCRI          4       0
PCGP_2014                  2      32
```

## Per-sample RNA-backing inventory

| sample_id | cohort | arm | sample_type | manifest_has_rna | file_backed_rna | expression_source | data_type |
|---|---|---|---|---|---|---|---|
| ACH-003171 | DepMap_26Q1 | disrupted | cell_line | True | True | OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv | log2(TPM+1) protein-coding |
| ACH-003185 | DepMap_26Q1 | disrupted | cell_line | True | True | OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv | log2(TPM+1) protein-coding |
| OS33 | Lawlor_2026_GCCRI | disrupted | pdx | False | False | Lawlor 2026 ships SV/IHC labels only; no expression matrix in data/raw/. | unknown |
| OS37 | Lawlor_2026_GCCRI | disrupted | pdx | False | False | Lawlor 2026 ships SV/IHC labels only; no expression matrix in data/raw/. | unknown |
| OS39R | Lawlor_2026_GCCRI | disrupted | pdx | False | False | Lawlor 2026 ships SV/IHC labels only; no expression matrix in data/raw/. | unknown |
| 0560-000 | Lawlor_2026_GCCRI | disrupted | pdx | False | False | Lawlor 2026 ships SV/IHC labels only; no expression matrix in data/raw/. | unknown |
| SJOS001103_D1 | PCGP_2014 | intact | tumor | False | False | Chen 2014 PCGP cohort ships TP53 SV / ATRX / telomere-FISH tables only; no expression matrix in data/raw/. | unknown |
| SJOS001105_D1 | PCGP_2014 | intact | tumor | False | False | Chen 2014 PCGP cohort ships TP53 SV / ATRX / telomere-FISH tables only; no expression matrix in data/raw/. | unknown |
| SJOS001106_D1 | PCGP_2014 | intact | tumor | False | False | Chen 2014 PCGP cohort ships TP53 SV / ATRX / telomere-FISH tables only; no expression matrix in data/raw/. | unknown |
| SJOS001107_M1 | PCGP_2014 | intact | tumor | False | False | Chen 2014 PCGP cohort ships TP53 SV / ATRX / telomere-FISH tables only; no expression matrix in data/raw/. | unknown |
| SJOS001107_M2 | PCGP_2014 | intact | tumor | False | False | Chen 2014 PCGP cohort ships TP53 SV / ATRX / telomere-FISH tables only; no expression matrix in data/raw/. | unknown |
| SJOS001112_M2 | PCGP_2014 | intact | tumor | False | False | Chen 2014 PCGP cohort ships TP53 SV / ATRX / telomere-FISH tables only; no expression matrix in data/raw/. | unknown |
| SJOS001117_D1 | PCGP_2014 | intact | tumor | False | False | Chen 2014 PCGP cohort ships TP53 SV / ATRX / telomere-FISH tables only; no expression matrix in data/raw/. | unknown |
| SJOS001118_D1 | PCGP_2014 | intact | tumor | False | False | Chen 2014 PCGP cohort ships TP53 SV / ATRX / telomere-FISH tables only; no expression matrix in data/raw/. | unknown |
| SJOS001119_D1 | PCGP_2014 | intact | tumor | False | False | Chen 2014 PCGP cohort ships TP53 SV / ATRX / telomere-FISH tables only; no expression matrix in data/raw/. | unknown |
| SJOS001120_D1 | PCGP_2014 | intact | tumor | False | False | Chen 2014 PCGP cohort ships TP53 SV / ATRX / telomere-FISH tables only; no expression matrix in data/raw/. | unknown |
| SJOS001123_D1 | PCGP_2014 | intact | tumor | False | False | Chen 2014 PCGP cohort ships TP53 SV / ATRX / telomere-FISH tables only; no expression matrix in data/raw/. | unknown |
| SJOS001125_D1 | PCGP_2014 | intact | tumor | False | False | Chen 2014 PCGP cohort ships TP53 SV / ATRX / telomere-FISH tables only; no expression matrix in data/raw/. | unknown |
| SJOS001126_D1 | PCGP_2014 | intact | tumor | False | False | Chen 2014 PCGP cohort ships TP53 SV / ATRX / telomere-FISH tables only; no expression matrix in data/raw/. | unknown |
| SJOS001127_D1 | PCGP_2014 | intact | tumor | False | False | Chen 2014 PCGP cohort ships TP53 SV / ATRX / telomere-FISH tables only; no expression matrix in data/raw/. | unknown |
| SJOS001_M | PCGP_2014 | intact | tumor | False | False | Chen 2014 PCGP cohort ships TP53 SV / ATRX / telomere-FISH tables only; no expression matrix in data/raw/. | unknown |
| SJOS002_D | PCGP_2014 | disrupted | tumor | False | False | Chen 2014 PCGP cohort ships TP53 SV / ATRX / telomere-FISH tables only; no expression matrix in data/raw/. | unknown |
| SJOS003_D | PCGP_2014 | intact | tumor | False | False | Chen 2014 PCGP cohort ships TP53 SV / ATRX / telomere-FISH tables only; no expression matrix in data/raw/. | unknown |
| SJOS004_D | PCGP_2014 | intact | tumor | False | False | Chen 2014 PCGP cohort ships TP53 SV / ATRX / telomere-FISH tables only; no expression matrix in data/raw/. | unknown |
| SJOS005_D | PCGP_2014 | intact | tumor | False | False | Chen 2014 PCGP cohort ships TP53 SV / ATRX / telomere-FISH tables only; no expression matrix in data/raw/. | unknown |
| SJOS006_D | PCGP_2014 | intact | tumor | False | False | Chen 2014 PCGP cohort ships TP53 SV / ATRX / telomere-FISH tables only; no expression matrix in data/raw/. | unknown |
| SJOS007_D | PCGP_2014 | intact | tumor | False | False | Chen 2014 PCGP cohort ships TP53 SV / ATRX / telomere-FISH tables only; no expression matrix in data/raw/. | unknown |
| SJOS008_D | PCGP_2014 | intact | tumor | False | False | Chen 2014 PCGP cohort ships TP53 SV / ATRX / telomere-FISH tables only; no expression matrix in data/raw/. | unknown |
| SJOS009_D | PCGP_2014 | intact | tumor | False | False | Chen 2014 PCGP cohort ships TP53 SV / ATRX / telomere-FISH tables only; no expression matrix in data/raw/. | unknown |
| SJOS010_D | PCGP_2014 | intact | tumor | False | False | Chen 2014 PCGP cohort ships TP53 SV / ATRX / telomere-FISH tables only; no expression matrix in data/raw/. | unknown |
| SJOS010_M | PCGP_2014 | intact | tumor | False | False | Chen 2014 PCGP cohort ships TP53 SV / ATRX / telomere-FISH tables only; no expression matrix in data/raw/. | unknown |
| SJOS011_D | PCGP_2014 | intact | tumor | False | False | Chen 2014 PCGP cohort ships TP53 SV / ATRX / telomere-FISH tables only; no expression matrix in data/raw/. | unknown |
| SJOS012_D | PCGP_2014 | intact | tumor | False | False | Chen 2014 PCGP cohort ships TP53 SV / ATRX / telomere-FISH tables only; no expression matrix in data/raw/. | unknown |
| SJOS013_D | PCGP_2014 | intact | tumor | False | False | Chen 2014 PCGP cohort ships TP53 SV / ATRX / telomere-FISH tables only; no expression matrix in data/raw/. | unknown |
| SJOS014_D | PCGP_2014 | disrupted | tumor | False | False | Chen 2014 PCGP cohort ships TP53 SV / ATRX / telomere-FISH tables only; no expression matrix in data/raw/. | unknown |
| SJOS015_D | PCGP_2014 | intact | tumor | False | False | Chen 2014 PCGP cohort ships TP53 SV / ATRX / telomere-FISH tables only; no expression matrix in data/raw/. | unknown |
| SJOS016_D | PCGP_2014 | intact | tumor | False | False | Chen 2014 PCGP cohort ships TP53 SV / ATRX / telomere-FISH tables only; no expression matrix in data/raw/. | unknown |
| SJOS017_D | PCGP_2014 | intact | tumor | False | False | Chen 2014 PCGP cohort ships TP53 SV / ATRX / telomere-FISH tables only; no expression matrix in data/raw/. | unknown |
| SJOS018_D | PCGP_2014 | intact | tumor | False | False | Chen 2014 PCGP cohort ships TP53 SV / ATRX / telomere-FISH tables only; no expression matrix in data/raw/. | unknown |
| SJOS019_D | PCGP_2014 | intact | tumor | False | False | Chen 2014 PCGP cohort ships TP53 SV / ATRX / telomere-FISH tables only; no expression matrix in data/raw/. | unknown |

## Cross-tabulations (file-backed RNA only)

### `cohort` × arm

| cohort | disrupted |
|---|---|
| DepMap_26Q1 | 2 |

### `sample_type` × arm

| sample_type | disrupted |
|---|---|
| cell_line | 2 |

### `data_type` × arm

| data_type | disrupted |
|---|---|
| log2(TPM+1) protein-coding | 2 |

## Confounding assessment (file-backed RNA only)

- File-backed totals: **disrupted=2**, **intact=0**, total=2.
- arm `disrupted` (2 sample(s)):
    - by cohort: DepMap_26Q1: 2 (100%)
    - by sample_type: cell_line: 2 (100%)
    - by data_type: log2(TPM+1) protein-coding: 2 (100%)
- arm `intact` (0 sample(s)):
    - (no file-backed samples)

- Cohorts with both arms file-backed: `[]`
- Cohorts with both arms file-backed under the same data_type: `[]`
- `has_within_cohort_co_occurrence`: **False**
- `has_within_cohort_co_occurrence_same_dtype`: **False**

## AIM 1 DESIGN VERDICT

**Category: NOT ESTIMABLE AS DEFINED**

Supporting counts:

- File-backed RNA disrupted samples: 2
- File-backed RNA intact samples: 0
- Cohorts contributing to disrupted arm (file-backed): ['DepMap_26Q1']
- Cohorts contributing to intact arm (file-backed): []
- Cohorts carrying both arms (same data_type): []

Justification:

> File-backed RNA is available for 2 disrupted sample(s) (from cohort(s) DepMap_26Q1) and 0 intact sample(s) (from cohort(s) (none)). TCAB1 status is perfectly or near-perfectly confounded with cohort/platform/sample_type, so the Aim 1 primary contrast as currently defined is not statistically estimable. The contrast must be redefined (e.g. require RNA-backed cohorts on both sides) or restricted (e.g. drop arms without file-backed RNA) before any DE analysis can proceed.
