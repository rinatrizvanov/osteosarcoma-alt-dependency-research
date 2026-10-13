# Exposure scan — WRAP53 (TCAB1) in DepMap 26Q1 osteosarcoma

Iteration 01, exposure side only. Run **blind to ALT status**: no ALT label was
loaded, joined, or inspected at any point in producing this file. WRAP53 is the
exposure variable here; cutpoints below are characterised *before* any
relationship to the outcome has been seen.

- Date of scan: 2026-10-12
- Release: DepMap 26Q1
- Interpreter: project venv (`.venv/bin/python`); pandas 2.3.3, numpy 2.5.4, scipy 1.18.1
- Loaders reused verbatim (no reimplementation): `src.cohort.osteosarcoma_models`,
  `src.cohort.rna_model_ids`, `src.stratify.osteosarcoma_rna_models`,
  `src.stratify.load_expression`, `src.stratify.load_wrap53_cn`
- Sources read: `data/raw/depmap/Model.csv`,
  `data/raw/depmap/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv`,
  `data/raw/depmap/OmicsCNGeneWGS.csv`

---

## 1. Membership steps (exact N)

| Step | N | Lost at this step |
|---|---|---|
| `Model.csv` rows matching `'steosarcoma'` (OncotreePrimaryDisease or OncotreeSubtype) | **35** | — |
| ∩ ModelIDs present in the expression matrix (`rna_model_ids()`) | **34** | 1 — `ACH-001094` (HTK-, Osteosarcoma) |
| of those 34, WRAP53 gene-level CN available (non-NaN in `OmicsCNGeneWGS.csv`) | **33** | 1 — `ACH-000613` (HOS) absent from the WGS CN file |

The 34-line RNA-backed set is the reference set used for every statistic and
cutpoint below. CN coverage is 33/34 (97.1%); the one CN gap is a missing row in
`OmicsCNGeneWGS.csv`, not a NaN value. No CRISPR screen membership was probed —
`CRISPRGeneEffect.csv` was never opened.

---

## 2. Distribution of WRAP53 log2(TPM+1), n = 34

| Statistic | Value |
|---|---|
| n | 34 |
| min | 0.000000 |
| 25th percentile (Q1) | 3.077646 |
| median (Q2) | 3.602946 |
| 75th percentile (Q3) | 4.558937 |
| max | 5.887913 |
| mean | 3.602633 |
| SD (ddof = 1) | 1.295000 |
| IQR | 1.481291 |
| **count at exactly 0.0** | **2** |

Supporting shape statistics: skewness −1.060, excess kurtosis +1.461,
Shapiro–Wilk W = 0.9195 (p = 0.0157). Both the negative skew and the normality
rejection are driven by the two exact zeros; among the 32 non-zero lines the
distribution is min 2.041, max 5.888, mean 3.828, SD 0.946.

Counts below low thresholds: `< 0.1` → 2 lines; `< 0.5` → 2 lines. That is, there
are no intermediate low-but-nonzero values between 0.000 and 2.041 — see §6.

---

## 3. Per-line table, sorted ascending by WRAP53 expression

No ALT column is present, by design.

| # | ModelID | StrippedCellLineName | CellLineName | WRAP53 log2(TPM+1) | WRAP53 CN | TP53 log2(TPM+1) |
|---:|---|---|---|---:|---:|---:|
| 1 | ACH-002834 | OS384 | OS384 | 0.000000 | 0.0000895 | 0.000000 |
| 2 | ACH-003185 | NRHOS4 | NRH-OS4 | 0.000000 | 0.0000075 | 0.105981 |
| 3 | ACH-002142 | HUO3N1 | HuO-3N1 | 2.041187 | 0.591413 | 1.290525 |
| 4 | ACH-001818 | C396 | C396 | 2.116691 | 0.927663 | 4.603757 |
| 5 | ACH-001712 | HS860T | Hs 860.T | 2.476716 | 0.959032 | 4.637402 |
| 6 | ACH-001817 | C393 | C393 | 2.495810 | 0.948905 | 5.079182 |
| 7 | ACH-003171 | NRHOS1 | NRH-OS1 | 2.748135 | 0.703361 | 2.531332 |
| 8 | ACH-001814 | OS252 | OS252 | 2.940866 | 0.660346 | 0.897178 |
| 9 | ACH-001813 | C242 | C242 | 3.066105 | 0.920135 | 3.960224 |
| 10 | ACH-003176 | NRHGCT1 | NRH-GCT1 | 3.112270 | 1.006233 | 6.317411 |
| 11 | ACH-002470 | OS186 | OS186 | 3.247747 | 0.942557 | 2.390417 |
| 12 | ACH-002433 | CCLFPEDS0019T | CCLF_PEDS_0019_T | 3.335580 | 1.311482 | 6.541247 |
| 13 | ACH-003178 | NRHOS2 | NRH-OS2 | 3.426647 | 0.661215 | 1.056812 |
| 14 | ACH-001715 | CAL72 | CAL-72 | 3.458638 | 0.686669 | 0.290298 |
| 15 | ACH-002778 | NCCOS1X2C1 | NCC-OS1-X2-C1 | 3.462052 | 1.690069 | 1.508246 |
| 16 | ACH-003573 | C337 | C337 | 3.483829 | 0.784511 | 0.458391 |
| 17 | ACH-002688 | OS457 | OS457 | 3.488528 | 0.711127 | 4.089312 |
| 18 | ACH-000748 | SJSA1 | SJSA-1 | 3.717364 | 0.868129 | 6.082545 |
| 19 | ACH-002779 | NCCESOS1C1 | NCC-ESOS1-C1 | 3.851769 | 0.968657 | 6.464547 |
| 20 | ACH-000364 | U2OS | U-2 OS | 4.100406 | 0.700637 | 4.993584 |
| 21 | ACH-000359 | MG63 | MG-63 | 4.145083 | 1.161826 | 2.343769 |
| 22 | ACH-001275 | SYO1 | SYO-1 | 4.275134 | 0.981565 | 5.733338 |
| 23 | ACH-000410 | SAOS2 | Saos-2 | 4.286105 | 0.842855 | 5.007913 |
| 24 | ACH-001526 | HUO9 | HuO9 | 4.300614 | 1.090083 | 4.533165 |
| 25 | ACH-003587 | 531MII | 531MII | 4.415177 | 1.393285 | 0.537147 |
| 26 | ACH-002067 | NOS1 | NOS-1 | 4.606857 | 1.494031 | 5.930165 |
| 27 | ACH-002690 | OS526 | OS526 | 4.618644 | 1.026855 | 1.228166 |
| 28 | ACH-002471 | OS052 | OS052 | 4.627875 | 0.906692 | 5.426916 |
| 29 | ACH-001001 | 143B | 143B | 4.785828 | 1.237176 | 6.133198 |
| 30 | ACH-000613 | HOS | HOS | 4.800227 | *n/a* | 6.800059 |
| 31 | ACH-002689 | OS525 | OS525 | 4.877960 | 1.739669 | 1.554225 |
| 32 | ACH-002178 | NY | NY | 5.022417 | 0.040120 | 0.102668 |
| 33 | ACH-002069 | HSOS1 | HS-Os-1 | 5.269336 | 2.357259 | 5.934514 |
| 34 | ACH-000082 | G292CLONEA141B1 | G-292, clone A141B1 | 5.887913 | 2.339662 | 2.474945 |

Notes on the table, factual only:
- Ties at 0.000000 (rows 1–2) are broken by ModelID.
- `ACH-000613` (HOS) has no row in `OmicsCNGeneWGS.csv`; CN is reported as *n/a*
  rather than imputed.
- `ACH-002178` (NY) carries near-zero WRAP53 CN (0.0401) with high expression
  (5.022) — a CN/expression discordance; `ACH-002069` and `ACH-000082` carry the
  highest CN (~2.34–2.36).
- TP53 expression is included as a pre-specified covariate/neighbour readout only.
  It spans 0.000 (`ACH-002834`) to 6.800 (`ACH-000613`).

---

## 4. Candidate cutpoint rules and resulting group sizes (34-line reference set)

| Rule | Numeric cutpoint | n low | n high | Low arm = rows (table §3) |
|---|---:|---:|---:|---|
| (a) exactly-zero vs non-zero | `WRAP53 == 0.0` (cut = 0.000000) | **2** | **32** | 1–2 |
| (b) bottom tertile, 33.3rd percentile (`< cut`) | **3.335580** | **11** | **23** | 1–11 |
| (c) median split (`< cut`) | **3.602946** | **17** | **17** | 1–17 |
| (d) bottom quartile, 25th percentile (`< cut`) | **3.077646** | **9** | **25** | 1–9 |

Conventions, stated so the counts are reproducible: quantiles are
`numpy.quantile` with the default linear interpolation; the low arm is
`value < cut` (strict), which matches the existing `src.stratify` convention.
Four mechanical notes on where each cut lands relative to the observed values:

- (b) The tertile cut 3.335580 is *itself* an observed value (`ACH-002433`), so a
  strict `<` places that line in the high arm; `<=` would give 12/22.
- (d) The quartile cut 3.077646 is *not* an observed value: it is the
  interpolated point between rows 9 and 10 (3.066105, 3.112270). Strict `<` and
  `<=` therefore both give 9 low / 25 high, and the split is unambiguous.
- (c) The median 3.602946 is likewise not an observed value — it is the midpoint
  of rows 17 and 18 (3.488528, 3.717364) — so `<` and `<=` both give 17/17.
- Only rule (b) is convention-sensitive. Any downstream code should fix the
  inequality explicitly rather than rely on it.

Group sizes for rules (b)–(d) are driven entirely by percentile position, not by
any gap in the data (see §6). Rule (a) is the only rule whose boundary coincides
with a discontinuity in the observed values.

---

## 5. Spearman correlations across the reference set

| Pair | n | rho | p |
|---|---:|---:|---:|
| WRAP53 log2(TPM+1) vs WRAP53 gene-level CN | 33 | **+0.5537** | **8.30 × 10⁻⁴** |
| WRAP53 log2(TPM+1) vs TP53 log2(TPM+1) | 34 | **+0.2443** | **0.1637** |

Computed with `scipy.stats.spearmanr`, two-sided, pairwise-complete (the CN
correlation drops `ACH-000613`, which has no CN row). The CN association is
moderate and positive: WRAP53 expression tracks WRAP53 copy number, so the
exposure is partly dosage-driven and CN is not an independent confirmation
channel. The WRAP53–TP53 expression correlation is weak and not significant at
the conventional 0.05 level on n = 34.

---

## 6. Shape of the distribution: zero-inflated spike, not a low-expression mode

Judged from the sorted values in §3, the gap structure is:

| Rank boundary | Lower value | Upper value | Gap |
|---|---:|---:|---:|
| 2 → 3 | 0.000000 | 2.041187 | **2.041187** |
| 33 → 34 | 5.269336 | 5.887913 | 0.618577 |
| 4 → 5 | 2.116691 | 2.476716 | 0.360025 |
| 6 → 7 | 2.495810 | 2.748135 | 0.252325 |
| 19 → 20 | 3.851769 | 4.100406 | 0.248637 |
| 32 → 33 | 5.022417 | 5.269336 | 0.246694 |

**Verdict stated plainly.** The distribution is **not** bimodal in the sense of
two populated modes. It is a single, roughly continuous mode over 2.04–5.89 plus a
**discrete two-line spike at exactly 0.0**. The only real discontinuity is the
2.04-unit gap between the two zeros and the third-lowest line — 3.3× larger than
the next-largest gap and 1.58 SD of the full distribution. Inside the non-zero
range the values are dense and gap-free: the largest internal gap is 0.619, and
that one is in the upper tail, not at any candidate low/high boundary. Shapiro–Wilk
rejects normality (p = 0.0157), but that rejection is produced by the two zeros,
not by a second cluster of low-but-nonzero lines.

Consequences for separability, stated without reference to any outcome:

- A "WRAP53-low" group defined as **exactly zero** is naturally and
  unambiguously separable. It is a discrete, data-defined category — plausibly a
  floor/undetected-transcript or homozygous-loss category, consistent with both
  zero-expression lines also carrying essentially zero WRAP53 CN
  (8.95 × 10⁻⁵ and 7.50 × 10⁻⁶). But it has only **n = 2**, below the existing
  `src.stratify.MIN_CELL_SIZE` of 4.
- Any **percentile-based** "WRAP53-low" group — tertile, median, or quartile — is
  **arbitrary**. Those cuts fall inside the dense continuous body of the
  distribution, where adjacent lines differ by ~0.01–0.25 log2 units, so a
  cutpoint moved by a fraction of a unit reassigns lines. Nothing in the data
  marks a boundary at 3.08, 3.34 or 3.60; the resulting group sizes are
  consequences of the chosen percentile and the strict-inequality convention, not
  of structure in the exposure.
- There is therefore a direct trade-off that cannot be resolved from the exposure
  alone: the only non-arbitrary cut yields n = 2, and every cut that yields a
  usable arm size is arbitrary. This should be settled by pre-specification
  before the outcome is unblinded, not after.

---

## 7. LIMITATIONS OF THIS SCAN

1. **n = 34.** All quantiles, the SD, and both correlations rest on 34
   observations. Percentile cutpoints have wide sampling uncertainty at this size;
   a single added or removed line can shift the tertile across an observed value
   and move a line between arms. No confidence intervals on the cutpoints are
   reported here, and none should be inferred.
2. **The two zeros are not mechanistically resolved.** log2(TPM+1) = 0.000000 is
   consistent with true absence of transcript, with a homozygous deletion, and
   with a quantification floor or low library depth. The near-zero CN for both
   lines is consistent with loss, but gene-level CN from `OmicsCNGeneWGS.csv` was
   not checked against the segment-level file, and no RNA depth or QC metric was
   inspected. The zero category is treated as an observed value, not a validated
   biological state.
3. **Single expression modality, single transcript level.** Only
   `OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv` was read, at gene level.
   No isoform, protein, promoter-methylation, or mutation data was consulted, so
   this scan cannot speak to whether low WRAP53 mRNA corresponds to low TCAB1
   protein or to a loss-of-function allele. A line with normal mRNA and a damaging
   variant would be invisible here.
4. **CN coverage is 33/34** and CN was not imputed, so the expression-vs-CN
   correlation excludes `ACH-000613` (HOS). CN is also only moderately correlated
   with expression (rho = 0.554), so it corroborates rather than validates the
   exposure.
5. **CN/expression discordances are unexplained.** `ACH-002178` (NY) has WRAP53
   CN 0.0401 with expression 5.022. This was noted, not investigated.
6. **No normalisation for confounders.** Values are raw log2(TPM+1); nothing was
   adjusted for growth pattern, media, passage, sex, primary-vs-metastasis,
   tumour purity, or the OS-subtype heterogeneity visible in `Model.csv`
   (`ACH-003176` is a giant-cell-tumour-derived line, `ACH-001275` and
   `ACH-002779` carry non-conventional subtype annotations). Whether any of these
   covary with WRAP53 was not tested.
7. **Only two covariates were loaded** (WRAP53 CN, TP53 expression). No
   multivariable or partial-correlation analysis was performed, and the single
   WRAP53–TP53 test is unadjusted for multiplicity.
8. **No formal modality test.** The unimodal-plus-spike judgement in §6 is a
   descriptive reading of sorted values and gaps. No Hartigan dip test, kernel
   density estimate, or mixture-model fit was run, and n = 34 would give any such
   test low power regardless.
9. **Blindness constrains what can be checked.** By design this scan cannot
   assess whether the reference set is balanced on the outcome, whether the
   exposure distribution differs between outcome groups, or whether any cutpoint
   is informative. Those are outcome-side questions and are out of scope here by
   construction.
10. **One release, no replication.** DepMap 26Q1 only. No independent
    osteosarcoma expression cohort was used to check that the zero spike or the
    quantile positions reproduce.

---

## 8. Blindness attestation — prohibited files NOT opened

I confirm that in producing this file I did **not** read, load, parse, join,
stream, `grep`, or otherwise inspect the contents of any of the following:

- `references/master_manifest.csv` — **not opened**
- `data/processed/cohort_alt_sets.csv` — **not opened**
- `data/processed/tcab1_strata.csv` and
  `data/processed/tcab1_strata_alt_negative.csv` — **not opened**
- any other file or column carrying `alt_status` or `alt_call_direction`
  (including `data/processed/primary_result.json`,
  `data/processed/composite_scores.csv`, `data/processed/per_gene_results.csv`,
  `data/processed/stage1_verdict.json`, `data/processed/robustness_table.csv`,
  and every `notebook_0*_report.md`) — **not opened**
- `data/raw/depmap/CRISPRGeneEffect.csv` — **not opened**; `src.dependency_stats`
  — **not imported**
- `docs/stage1_outcome_memo.md` — **not opened**
- `README.md` — **not opened**

Functions called were restricted to `osteosarcoma_models()`, `rna_model_ids()`,
`osteosarcoma_rna_models()`, `load_expression(..., symbols=["WRAP53","TP53"])`
and `load_wrap53_cn(...)`. None of these touch the manifest or any ALT label;
`src.cohort.build_cohort()`, which does join the manifest, was **not** called.
`src.cohort` and `src.stratify` were imported as modules — importing binds
`MANIFEST_PATH` as a path constant only and reads no file.

No existing file was modified. This document is the only file written into the
repository by this scan; the scan script itself was written to the session
scratchpad, outside the repository.
