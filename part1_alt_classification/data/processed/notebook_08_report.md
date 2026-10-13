# Notebook 08 — Recovery Route Inventory

_Generated: 2026-05-18T16:43:12Z_


Read-only diagnostic testing three RNA-recovery routes not evaluated in notebooks 06/07 (exact `sample_id` matching only).


## Route B — DepMap cell-line-name join

- Candidates: **38** (Taglialatela + ALTitude DepMap panel overlap).

- Mapped to ModelID: **38**

- In DepMap RNA matrix: **37**

- ALT-positive + RNA-backed: **4**

- ALT-negative + RNA-backed: **2**


### Route B mapping table

| cell_line_name | mapped_model_id | in_depmap_rna | alt_status | source | cohort |

|---|---|---|---|---|---|

| G292CLONEA141B1 | ACH-000082 | True | alt_unknown | pending_literature_merge | DepMap_26Q1 |
| MG63 | ACH-000359 | True | alt_unknown | pending_literature_merge | DepMap_26Q1 |
| U2OS | ACH-000364 | True | alt_unknown | pending_literature_merge | DepMap_26Q1 |
| SAOS2 | ACH-000410 | True | alt_unknown | pending_literature_merge | DepMap_26Q1 |
| HOS | ACH-000613 | True | alt_unknown | pending_literature_merge | DepMap_26Q1 |
| SJSA1 | ACH-000748 | True | alt_unknown | pending_literature_merge | DepMap_26Q1 |
| 143B | ACH-001001 | True | alt_unknown | pending_literature_merge | DepMap_26Q1 |
| SYO1 | ACH-001275 | True | alt_unknown | pending_literature_merge | DepMap_26Q1 |
| HUO9 | ACH-001526 | True | alt_unknown | pending_literature_merge | DepMap_26Q1 |
| CAL72 | ACH-001715 | True | alt_unknown | pending_literature_merge | DepMap_26Q1 |
| C242 | ACH-001813 | True | alt_unknown | pending_literature_merge | DepMap_26Q1 |
| OS252 | ACH-001814 | True | alt_unknown | pending_literature_merge | DepMap_26Q1 |
| C393 | ACH-001817 | True | alt_unknown | pending_literature_merge | DepMap_26Q1 |
| NOS1 | ACH-002067 | True | alt_unknown | pending_literature_merge | DepMap_26Q1 |
| HSOS1 | ACH-002069 | True | alt_unknown | pending_literature_merge | DepMap_26Q1 |
| HUO3N1 | ACH-002142 | True | alt_unknown | pending_literature_merge | DepMap_26Q1 |
| NY | ACH-002178 | True | alt_unknown | pending_literature_merge | DepMap_26Q1 |
| CCLFPEDS0019T | ACH-002433 | True | alt_unknown | pending_literature_merge | DepMap_26Q1 |
| OS186 | ACH-002470 | True | alt_unknown | pending_literature_merge | DepMap_26Q1 |
| OS052 | ACH-002471 | True | alt_unknown | pending_literature_merge | DepMap_26Q1 |
| OS457 | ACH-002688 | True | alt_unknown | pending_literature_merge | DepMap_26Q1 |
| OS525 | ACH-002689 | True | alt_unknown | pending_literature_merge | DepMap_26Q1 |
| OS526 | ACH-002690 | True | alt_unknown | pending_literature_merge | DepMap_26Q1 |
| NCCESOS1C1 | ACH-002779 | True | alt_unknown | pending_literature_merge | DepMap_26Q1 |
| OS384 | ACH-002834 | True | alt_unknown | pending_literature_merge | DepMap_26Q1 |
| NRHOS1 | ACH-003171 | True | alt_unknown | pending_literature_merge | DepMap_26Q1 |
| NRHGCT1 | ACH-003176 | True | alt_unknown | pending_literature_merge | DepMap_26Q1 |
| NRHOS2 | ACH-003178 | True | alt_unknown | pending_literature_merge | DepMap_26Q1 |
| NRHOS4 | ACH-003185 | True | alt_unknown | pending_literature_merge | DepMap_26Q1 |
| C337 | ACH-003573 | True | alt_unknown | pending_literature_merge | DepMap_26Q1 |
| 531MII | ACH-003587 | True | alt_unknown | pending_literature_merge | DepMap_26Q1 |
| SAOS2 | ACH-000410 | True | alt_gold_positive | Taglialatela2026 | Taglialatela_2026 |
| G292 | ACH-000082 | True | alt_gold_positive | Taglialatela2026 | Taglialatela_2026 |
| OS384 | ACH-002834 | True | alt_gold_positive | Taglialatela2026 | Taglialatela_2026 |
| HOS | ACH-000613 | True | alt_negative_high_conf | Taglialatela2026 | Taglialatela_2026 |
| SJSA1 | ACH-000748 | True | alt_negative_high_conf | Lawlor et al. bioRxiv 2026.01.02.691864|Taglialatela2026 | Taglialatela_2026 |
| SKLU1 | ACH-000309 | False | alt_gold_positive | Taglialatela2026 | Taglialatela_2026 |
| TM31 | ACH-000592 | True | alt_gold_positive | Taglialatela2026 | Taglialatela_2026 |

## Route A — de Nonneville (E-MTAB-11743)

- Manifest samples with ALT labels: **22** (positive=16, negative=6).

- Expression deposited at ArrayExpress **E-MTAB-11743** (not downloaded in this notebook).

- **Human verification required:** confirm assay type is RNA-seq, not methylation array.


```
TODO — Route A (de Nonneville / E-MTAB-11743)
[ ] Download E-MTAB-11743 from ArrayExpress (do NOT assume assay type).
[ ] Open sample metadata: confirm modality is RNA-seq (not 450K methylation alone).
[ ] Map ArrayExpress sample names to these manifest sample_ids:
      ALT-positive (16): deNonneville_1941.0, deNonneville_2254.0, deNonneville_2523.0, deNonneville_3103.0, deNonneville_4001.0, deNonneville_4078.0, deNonneville_4435.0, deNonneville_4462.0, deNonneville_4591.0, deNonneville_4665.0, deNonneville_4680.0, deNonneville_5136.0, deNonneville_7695.0, deNonneville_8299.0, deNonneville_8781.0, deNonneville_9496.0
      ALT-negative (6): deNonneville_5090.0, deNonneville_5664.0, deNonneville_6500.0, deNonneville_6537.0, deNonneville_6699.0, deNonneville_7855.0
[ ] If RNA-seq confirmed: add expression file to data/raw/ and register in EXPRESSION_SOURCES.
[ ] Re-run notebooks 07 and 08 after integration.
```

## Route C — Lawlor internal PDX RNA

- Lawlor PDX rows in manifest: **17**

- With ALT and/or architecture labels: **17** (ALT-positive=16, ALT-negative=1, architecture-labeled=4).

- RNA is **not public**; requires internal lab access.


```
TODO — Route C (Lawlor internal PDX RNA — not public)
[ ] Request matched RNA-seq (raw or processed counts/TPM) from the Lawlor/GCCRI lab.
[ ] Samples needing matched RNA (17 PDX models):
      All Lawlor sample_ids: 0522-000, 0560-000, 1922, 1981, OS33, OS37, OS39R, OS41, OS42, OS46, OS48, OS49, OS50, OS52, OS55, OS58, SJSA1
[ ] Priority for Aim 2 (ALT labels):
      ALT-positive (16): 0522-000, 0560-000, 1922, 1981, OS33, OS37, OS39R, OS41, OS42, OS46, OS48, OS49, OS50, OS52, OS55, OS58
      ALT-negative (1): SJSA1
[ ] Priority for Aim 1 (architecture / TCAB1):
      Architecture-labeled: 0560-000, OS33, OS37, OS39R
[ ] On receipt: place files under data/raw/ (new subfolder) and register in EXPRESSION_SOURCES.
[ ] Re-run notebooks 06–08 after integration.
```

## Combined — Route-B-only best case

Assumes Route B succeeds; Routes A and C remain pending.


### Aim 1 (TCAB1 pragmatic contrast)

- File-backed RNA samples: **2** (disrupted=2, intact=0).

| cohort | disrupted |
|---|---|
| DepMap_26Q1 | 2 |

- Verdict: **NOT ESTIMABLE AS DEFINED**

> File-backed RNA is available for 2 disrupted sample(s) (from cohort(s) DepMap_26Q1) and 0 intact sample(s) (from cohort(s) (none)). TCAB1 status is perfectly or near-perfectly confounded with cohort/platform/sample_type, so the Aim 1 primary contrast as currently defined is not statistically estimable. The contrast must be redefined (e.g. require RNA-backed cohorts on both sides) or restricted (e.g. drop arms without file-backed RNA) before any DE analysis can proceed.


### Aim 2 (ALT classifier contrast)

- File-backed RNA samples: **7** (positive=4, negative=3).

| cohort | alt_gold_positive | alt_negative_high_conf |
|---|---|---|
| Lawlor_2026_GCCRI | 0 | 1 |
| Taglialatela_2026 | 4 | 2 |

- Cohorts with both arms (same data_type): `[('Taglialatela_2026', 'log2(TPM+1) protein-coding')]`

- Verdict: **NOT ESTIMABLE AS DEFINED**

> File-backed RNA is available for 4 ALT-positive sample(s) (from cohort(s) Taglialatela_2026) and 3 ALT-negative sample(s) (from cohort(s) Lawlor_2026_GCCRI, Taglialatela_2026). ALT status is perfectly or near-perfectly confounded with cohort / platform / sample_type, so the Aim 2 contrast as currently defined is not statistically estimable. The contrast must be redefined (e.g. require RNA-backed cohorts on both sides) or restricted (e.g. drop arms without file-backed RNA) before any ALT classifier can be trained.


## RECOVERY VERDICT

**Route B (achievable now, zero external dependency):** cell-line join recovers **37** DepMap RNA profiles (4 ALT-positive, 2 ALT-negative among candidates with gold ALT calls). Within the Aim 2 primary contrast, Route B alone yields **7** file-backed samples (4 positive, 3 negative) — verdict **NOT ESTIMABLE AS DEFINED** (within-cohort minimum = 3 per arm). Aim 1 remains **NOT ESTIMABLE AS DEFINED** (2 disrupted, 0 intact file-backed).


**Routes A and C (pending):** de Nonneville (`E-MTAB-11743`, 22 labeled tumors) and Lawlor (17 labeled PDX models) would add the largest labeled cohorts if expression is obtained and sample IDs mapped. Neither was counted in the Route-B-only verdict above.


**Recommendation:**

1. Proceed with Route B integration for Taglialatela cell lines (6/7 map to DepMap RNA; SKLU1 maps to ACH-000309 but is absent from the expression matrix).

2. Inspect E-MTAB-11743 assay metadata before downloading (Route A).

3. Request Lawlor PDX RNA from the lab (Route C) — highest leverage for both aims.

