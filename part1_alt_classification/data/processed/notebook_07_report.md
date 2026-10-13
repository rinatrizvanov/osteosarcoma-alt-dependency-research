# Notebook 07 — Aim 2 RNA Inventory

_Generated: 2026-05-17T20:43:05Z_


Read-only diagnostic. Determines whether the Aim 2 ALT classifier contrast (`alt_status` ∈ `{alt_gold_positive, alt_negative_high_conf}`) is statistically estimable from the RNA available on disk, or whether ALT status is confounded with cohort / platform / sample type.


## Inputs

- Manifest: `data/manifests/master_manifest.csv` (115 rows)
- DepMap RNA: `data/raw/depmap/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv` (1719 unique ModelIDs)
- TARGET-OS TPM: `data/raw/target_os/TARGET-OS.star_tpm.tsv.gz` (88 sample columns; not joined to manifest as a cohort)

## Primary contrast composition

Total contrast rows: **59** (positive=48, negative=11).

Arm × cohort:

```
arm                 alt_gold_positive  alt_negative_high_conf
cohort                                                       
Lawlor_2026_GCCRI                  16                       1
PCGP_2014                          11                       2
Taglialatela_2026                   5                       2
de_Nonneville_2022                 16                       6
```

## Per-sample RNA-backing inventory

| sample_id | cohort | alt_status | sample_type | manifest_has_rna | file_backed_rna | expression_source | data_type |
|---|---|---|---|---|---|---|---|
| deNonneville_2523.0 | de_Nonneville_2022 | alt_gold_positive | tumor | False | False | de Nonneville 2022 ships a tumor ALT-call table (EV3) only; no expression matrix in data/raw/. | unknown |
| deNonneville_4435.0 | de_Nonneville_2022 | alt_gold_positive | tumor | False | False | de Nonneville 2022 ships a tumor ALT-call table (EV3) only; no expression matrix in data/raw/. | unknown |
| deNonneville_4680.0 | de_Nonneville_2022 | alt_gold_positive | tumor | False | False | de Nonneville 2022 ships a tumor ALT-call table (EV3) only; no expression matrix in data/raw/. | unknown |
| deNonneville_7695.0 | de_Nonneville_2022 | alt_gold_positive | tumor | False | False | de Nonneville 2022 ships a tumor ALT-call table (EV3) only; no expression matrix in data/raw/. | unknown |
| deNonneville_9496.0 | de_Nonneville_2022 | alt_gold_positive | tumor | False | False | de Nonneville 2022 ships a tumor ALT-call table (EV3) only; no expression matrix in data/raw/. | unknown |
| deNonneville_8781.0 | de_Nonneville_2022 | alt_gold_positive | tumor | False | False | de Nonneville 2022 ships a tumor ALT-call table (EV3) only; no expression matrix in data/raw/. | unknown |
| deNonneville_8299.0 | de_Nonneville_2022 | alt_gold_positive | tumor | False | False | de Nonneville 2022 ships a tumor ALT-call table (EV3) only; no expression matrix in data/raw/. | unknown |
| deNonneville_4462.0 | de_Nonneville_2022 | alt_gold_positive | tumor | False | False | de Nonneville 2022 ships a tumor ALT-call table (EV3) only; no expression matrix in data/raw/. | unknown |
| deNonneville_1941.0 | de_Nonneville_2022 | alt_gold_positive | tumor | False | False | de Nonneville 2022 ships a tumor ALT-call table (EV3) only; no expression matrix in data/raw/. | unknown |
| deNonneville_2254.0 | de_Nonneville_2022 | alt_gold_positive | tumor | False | False | de Nonneville 2022 ships a tumor ALT-call table (EV3) only; no expression matrix in data/raw/. | unknown |
| deNonneville_3103.0 | de_Nonneville_2022 | alt_gold_positive | tumor | False | False | de Nonneville 2022 ships a tumor ALT-call table (EV3) only; no expression matrix in data/raw/. | unknown |
| deNonneville_4001.0 | de_Nonneville_2022 | alt_gold_positive | tumor | False | False | de Nonneville 2022 ships a tumor ALT-call table (EV3) only; no expression matrix in data/raw/. | unknown |
| deNonneville_4078.0 | de_Nonneville_2022 | alt_gold_positive | tumor | False | False | de Nonneville 2022 ships a tumor ALT-call table (EV3) only; no expression matrix in data/raw/. | unknown |
| deNonneville_4591.0 | de_Nonneville_2022 | alt_gold_positive | tumor | False | False | de Nonneville 2022 ships a tumor ALT-call table (EV3) only; no expression matrix in data/raw/. | unknown |
| deNonneville_4665.0 | de_Nonneville_2022 | alt_gold_positive | tumor | False | False | de Nonneville 2022 ships a tumor ALT-call table (EV3) only; no expression matrix in data/raw/. | unknown |
| deNonneville_5136.0 | de_Nonneville_2022 | alt_gold_positive | tumor | False | False | de Nonneville 2022 ships a tumor ALT-call table (EV3) only; no expression matrix in data/raw/. | unknown |
| deNonneville_7855.0 | de_Nonneville_2022 | alt_negative_high_conf | tumor | False | False | de Nonneville 2022 ships a tumor ALT-call table (EV3) only; no expression matrix in data/raw/. | unknown |
| deNonneville_5090.0 | de_Nonneville_2022 | alt_negative_high_conf | tumor | False | False | de Nonneville 2022 ships a tumor ALT-call table (EV3) only; no expression matrix in data/raw/. | unknown |
| deNonneville_5664.0 | de_Nonneville_2022 | alt_negative_high_conf | tumor | False | False | de Nonneville 2022 ships a tumor ALT-call table (EV3) only; no expression matrix in data/raw/. | unknown |
| deNonneville_6500.0 | de_Nonneville_2022 | alt_negative_high_conf | tumor | False | False | de Nonneville 2022 ships a tumor ALT-call table (EV3) only; no expression matrix in data/raw/. | unknown |
| deNonneville_6537.0 | de_Nonneville_2022 | alt_negative_high_conf | tumor | False | False | de Nonneville 2022 ships a tumor ALT-call table (EV3) only; no expression matrix in data/raw/. | unknown |
| deNonneville_6699.0 | de_Nonneville_2022 | alt_negative_high_conf | tumor | False | False | de Nonneville 2022 ships a tumor ALT-call table (EV3) only; no expression matrix in data/raw/. | unknown |
| SJSA1 | Lawlor_2026_GCCRI | alt_negative_high_conf | pdx | False | False | Lawlor 2026 ships SV/IHC labels only; no expression matrix in data/raw/. | unknown |
| OS33 | Lawlor_2026_GCCRI | alt_gold_positive | pdx | False | False | Lawlor 2026 ships SV/IHC labels only; no expression matrix in data/raw/. | unknown |
| OS37 | Lawlor_2026_GCCRI | alt_gold_positive | pdx | False | False | Lawlor 2026 ships SV/IHC labels only; no expression matrix in data/raw/. | unknown |
| OS39R | Lawlor_2026_GCCRI | alt_gold_positive | pdx | False | False | Lawlor 2026 ships SV/IHC labels only; no expression matrix in data/raw/. | unknown |
| OS41 | Lawlor_2026_GCCRI | alt_gold_positive | pdx | False | False | Lawlor 2026 ships SV/IHC labels only; no expression matrix in data/raw/. | unknown |
| OS42 | Lawlor_2026_GCCRI | alt_gold_positive | pdx | False | False | Lawlor 2026 ships SV/IHC labels only; no expression matrix in data/raw/. | unknown |
| OS46 | Lawlor_2026_GCCRI | alt_gold_positive | pdx | False | False | Lawlor 2026 ships SV/IHC labels only; no expression matrix in data/raw/. | unknown |
| OS48 | Lawlor_2026_GCCRI | alt_gold_positive | pdx | False | False | Lawlor 2026 ships SV/IHC labels only; no expression matrix in data/raw/. | unknown |
| OS49 | Lawlor_2026_GCCRI | alt_gold_positive | pdx | False | False | Lawlor 2026 ships SV/IHC labels only; no expression matrix in data/raw/. | unknown |
| OS50 | Lawlor_2026_GCCRI | alt_gold_positive | pdx | False | False | Lawlor 2026 ships SV/IHC labels only; no expression matrix in data/raw/. | unknown |
| OS52 | Lawlor_2026_GCCRI | alt_gold_positive | pdx | False | False | Lawlor 2026 ships SV/IHC labels only; no expression matrix in data/raw/. | unknown |
| OS55 | Lawlor_2026_GCCRI | alt_gold_positive | pdx | False | False | Lawlor 2026 ships SV/IHC labels only; no expression matrix in data/raw/. | unknown |
| OS58 | Lawlor_2026_GCCRI | alt_gold_positive | pdx | False | False | Lawlor 2026 ships SV/IHC labels only; no expression matrix in data/raw/. | unknown |
| 0522-000 | Lawlor_2026_GCCRI | alt_gold_positive | pdx | False | False | Lawlor 2026 ships SV/IHC labels only; no expression matrix in data/raw/. | unknown |
| 0560-000 | Lawlor_2026_GCCRI | alt_gold_positive | pdx | False | False | Lawlor 2026 ships SV/IHC labels only; no expression matrix in data/raw/. | unknown |
| 1922 | Lawlor_2026_GCCRI | alt_gold_positive | pdx | False | False | Lawlor 2026 ships SV/IHC labels only; no expression matrix in data/raw/. | unknown |
| 1981 | Lawlor_2026_GCCRI | alt_gold_positive | pdx | False | False | Lawlor 2026 ships SV/IHC labels only; no expression matrix in data/raw/. | unknown |
| SAOS2 | Taglialatela_2026 | alt_gold_positive | cell_line | False | False | Taglialatela 2026 ships a cell-line list with ALT calls only; no expression matrix in data/raw/. | unknown |
| G292 | Taglialatela_2026 | alt_gold_positive | cell_line | False | False | Taglialatela 2026 ships a cell-line list with ALT calls only; no expression matrix in data/raw/. | unknown |
| OS384 | Taglialatela_2026 | alt_gold_positive | cell_line | False | False | Taglialatela 2026 ships a cell-line list with ALT calls only; no expression matrix in data/raw/. | unknown |
| HOS | Taglialatela_2026 | alt_negative_high_conf | cell_line | False | False | Taglialatela 2026 ships a cell-line list with ALT calls only; no expression matrix in data/raw/. | unknown |
| SJSA1 | Taglialatela_2026 | alt_negative_high_conf | cell_line | False | False | Taglialatela 2026 ships a cell-line list with ALT calls only; no expression matrix in data/raw/. | unknown |
| SKLU1 | Taglialatela_2026 | alt_gold_positive | cell_line | False | False | Taglialatela 2026 ships a cell-line list with ALT calls only; no expression matrix in data/raw/. | unknown |
| TM31 | Taglialatela_2026 | alt_gold_positive | cell_line | False | False | Taglialatela 2026 ships a cell-line list with ALT calls only; no expression matrix in data/raw/. | unknown |
| SJOS001_M | PCGP_2014 | alt_gold_positive | tumor | False | False | Chen 2014 PCGP cohort ships TP53 SV / ATRX / telomere-FISH tables only; no expression matrix in data/raw/. | unknown |
| SJOS002_D | PCGP_2014 | alt_gold_positive | tumor | False | False | Chen 2014 PCGP cohort ships TP53 SV / ATRX / telomere-FISH tables only; no expression matrix in data/raw/. | unknown |
| SJOS004_D | PCGP_2014 | alt_negative_high_conf | tumor | False | False | Chen 2014 PCGP cohort ships TP53 SV / ATRX / telomere-FISH tables only; no expression matrix in data/raw/. | unknown |
| SJOS005_D | PCGP_2014 | alt_gold_positive | tumor | False | False | Chen 2014 PCGP cohort ships TP53 SV / ATRX / telomere-FISH tables only; no expression matrix in data/raw/. | unknown |
| SJOS006_D | PCGP_2014 | alt_gold_positive | tumor | False | False | Chen 2014 PCGP cohort ships TP53 SV / ATRX / telomere-FISH tables only; no expression matrix in data/raw/. | unknown |
| SJOS008_D | PCGP_2014 | alt_gold_positive | tumor | False | False | Chen 2014 PCGP cohort ships TP53 SV / ATRX / telomere-FISH tables only; no expression matrix in data/raw/. | unknown |
| SJOS009_D | PCGP_2014 | alt_gold_positive | tumor | False | False | Chen 2014 PCGP cohort ships TP53 SV / ATRX / telomere-FISH tables only; no expression matrix in data/raw/. | unknown |
| SJOS012_D | PCGP_2014 | alt_negative_high_conf | tumor | False | False | Chen 2014 PCGP cohort ships TP53 SV / ATRX / telomere-FISH tables only; no expression matrix in data/raw/. | unknown |
| SJOS014_D | PCGP_2014 | alt_gold_positive | tumor | False | False | Chen 2014 PCGP cohort ships TP53 SV / ATRX / telomere-FISH tables only; no expression matrix in data/raw/. | unknown |
| SJOS015_D | PCGP_2014 | alt_gold_positive | tumor | False | False | Chen 2014 PCGP cohort ships TP53 SV / ATRX / telomere-FISH tables only; no expression matrix in data/raw/. | unknown |
| SJOS016_D | PCGP_2014 | alt_gold_positive | tumor | False | False | Chen 2014 PCGP cohort ships TP53 SV / ATRX / telomere-FISH tables only; no expression matrix in data/raw/. | unknown |
| SJOS017_D | PCGP_2014 | alt_gold_positive | tumor | False | False | Chen 2014 PCGP cohort ships TP53 SV / ATRX / telomere-FISH tables only; no expression matrix in data/raw/. | unknown |
| SJOS019_D | PCGP_2014 | alt_gold_positive | tumor | False | False | Chen 2014 PCGP cohort ships TP53 SV / ATRX / telomere-FISH tables only; no expression matrix in data/raw/. | unknown |

## Cross-tabulations (file-backed RNA only)

### `cohort` × alt_status

(no rows)

### `sample_type` × alt_status

(no rows)

### `data_type` × alt_status

(no rows)

## Confounding assessment (file-backed RNA only)

- File-backed totals: **positive=0**, **negative=0**, total=0.
- arm `alt_gold_positive` (0 sample(s)):
    - (no file-backed samples)
- arm `alt_negative_high_conf` (0 sample(s)):
    - (no file-backed samples)

- Cohorts with both arms file-backed: `[]`
- Cohorts with both arms file-backed under the same data_type: `[]`
- `has_within_cohort_co_occurrence`: **False**
- `has_within_cohort_co_occurrence_same_dtype`: **False**

## AIM 2 DESIGN VERDICT

**Category: NOT ESTIMABLE AS DEFINED**

Supporting counts:

- File-backed RNA ALT-positive samples: 0
- File-backed RNA ALT-negative samples: 0
- Cohorts contributing to ALT-positive arm (file-backed): []
- Cohorts contributing to ALT-negative arm (file-backed): []
- Cohorts carrying both arms (same data_type): []

Justification:

> File-backed RNA is available for 0 ALT-positive sample(s) (from cohort(s) (none)) and 0 ALT-negative sample(s) (from cohort(s) (none)). ALT status is perfectly or near-perfectly confounded with cohort / platform / sample_type, so the Aim 2 contrast as currently defined is not statistically estimable. The contrast must be redefined (e.g. require RNA-backed cohorts on both sides) or restricted (e.g. drop arms without file-backed RNA) before any ALT classifier can be trained.
