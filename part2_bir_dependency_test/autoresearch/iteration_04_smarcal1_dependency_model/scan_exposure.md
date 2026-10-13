# Exposure scan — iteration 04 candidate exposures (DepMap 26Q1 osteosarcoma)

Exposure side only, run **blind to the outcome**. Iteration 04's outcome is SMARCAL1
Chronos gene effect at the level of individual cell lines; nothing in this file is
derived from a gene-effect value. `src.dependency_stats` was **not imported**,
`CRISPRGeneEffect.csv` was **not read** except as a membership probe
(`src.cohort.crispr_model_ids()`, which streams the first field of each row and reads
no gene-effect column). `autoresearch/iteration_04_smarcal1_dependency_model/DISCLOSURES.md`
was read first; nothing it declares as already-seen is used here, because this file
contains no outcome quantity at all.

- Date of scan: 2026-10-12
- Release: DepMap 26Q1 (every input sha256-verified against `src.utils.EXPECTED_SHA256`)
- Interpreter: `uv run python`, project venv, Python 3.12; pandas 2.3.3, numpy 2.5.4, scipy 1.18.1
- Loaders reused, not reimplemented: `src.stratify.load_expression`,
  `src.stratify.expression_gene_columns`, `src.stratify.osteosarcoma_rna_models`,
  `src.cohort.osteosarcoma_models`, `src.cohort.rna_model_ids`,
  `src.cohort.crispr_model_ids`, `src.cohort.build_cohort` + `assert_all`
- Sources read: `Model.csv`, `OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv`,
  `OmicsCNGeneWGS.csv`, `OmicsCNSegmentsWGS.csv`, `references/master_manifest.csv`
  (ALT labels, used for set membership only — see the b0 note below)
- Scan scripts were written to the session scratchpad, outside the repository. This
  document is the only file this scan writes into the repo.

**Anti-bias rule b0 (iteration 01).** No table pairing ModelIDs with ALT labels was
printed to a terminal at any point. The ALT+ subset below is reported as distribution
summaries and ID-withheld sorted vectors; the only ModelIDs named anywhere in this
file are the four ATRX copy-number lines that iteration 02 already adjudicated by name.

**Anti-bias rule b (iteration 01).** Every threshold, gap and shape judgement is made
on the 34 RNA-backed osteosarcoma lines, the widest non-outcome superset. The
16-line ALT+ figures in §4 are reported so the prespec knows the *realised* spread of
each exposure in the analysis set; they are not used to choose a cutpoint.

### Candidate exposures

20 genes. Nineteen from the SMARCAL1/ALT replication-stress mechanism (Taglialatela 2026;
Broderick 2023; Cox/Flynn 2016), plus SMARCAL1 expression itself as an exposure
distinct from its gene effect (Liu 2023; Brosnan-Cashman 2021 report SMARCAL1 loss in
ALT+ tumours). The twentieth, **`ZRANB3`**, was added after this scan began, at the
literature check's request: a SNF2-family fork remodeler functionally paralogous to
SMARCAL1 but named nowhere in the mechanism motivating the primary, to be registered as
a **predicted-null negative-control exposure**. It is characterised here on exactly the
same footing as the other 19, with no special treatment, and §6 reports its correlation
with the nominated primary (`BLM`) and first alternate (`PRIMPOL`) explicitly — a
negative control collinear with the primary is not a control.

Candidates: `PRIMPOL`, `BLM`, `RPA1`, `RPA2`, `RPA3`, `POLA1`, `TOP3A`, `RMI1`, `RMI2`, `PML`, `EXD2`, `MUS81`, `SLX4`, `FANCM`, `WRN`, `ATRX`, `DAXX`, `TERT`, `SMARCAL1`, `ZRANB3`.

All 20 resolve to exactly one column in both the expression matrix and
`OmicsCNGeneWGS.csv`; no symbol was ambiguous or absent.

---

## 0. Membership (exact N)

| Step | N |
|---|---|
| `Model.csv` rows matching `'steosarcoma'` | 35 |
| ∩ present in the expression matrix = **reference set** (`osteosarcoma_rna_models()`) | **34** |
| of those, ALT+ and CRISPR-screened (frozen manifest, `alt_call_direction == 'positive'`) | **16** |
| ALT− and CRISPR-screened | 9 |
| reference-set lines with a row in `OmicsCNGeneWGS.csv` | 33 (1 missing: `ACH-000613`) |
| reference-set lines with rows in `OmicsCNSegmentsWGS.csv` | 33 (same line missing) |
| ALT+16 lines with gene-level CN **and** segment coverage | **16 / 16** |

`build_cohort()` passed `assert_all()` — ground-truth anchors, set sizes (16/9/1) and
tier counts all as pre-specified. All 16 ALT+ lines are inside the 34-line reference set,
so the reference set is a true superset of the analysis set.

---

## 1. Distribution of log2(TPM+1) over the 34 RNA-backed lines

| Gene | n | n = 0 | n < 0.5 | n < 1 | mean | SD | median | Q1 | Q3 | IQR | min | max |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| **PRIMPOL** | 34 | 0 | 0 | 0 | 3.695 | 0.710 | 3.704 | 3.159 | 4.211 | 1.052 | 2.187 | 5.085 |
| **BLM** | 34 | 0 | 1 | 5 | 3.207 | 1.268 | 3.575 | 3.093 | 4.086 | 0.994 | 0.462 | 4.793 |
| **RPA1** | 34 | 0 | 0 | 0 | 6.142 | 0.707 | 6.168 | 5.628 | 6.688 | 1.060 | 4.731 | 7.503 |
| **RPA2** | 34 | 0 | 0 | 0 | 6.131 | 0.755 | 6.131 | 5.653 | 6.581 | 0.929 | 4.371 | 8.060 |
| **RPA3** | 34 | 0 | 0 | 0 | 5.865 | 0.988 | 6.079 | 5.335 | 6.544 | 1.209 | 3.642 | 7.493 |
| **POLA1** | 34 | 0 | 0 | 0 | 3.631 | 0.862 | 3.908 | 3.374 | 4.272 | 0.898 | 1.167 | 4.576 |
| **TOP3A** | 34 | 0 | 0 | 0 | 4.042 | 1.012 | 4.156 | 3.201 | 4.748 | 1.547 | 1.995 | 6.268 |
| **RMI1** | 34 | 0 | 0 | 0 | 4.246 | 0.918 | 4.398 | 3.853 | 4.722 | 0.869 | 1.615 | 6.173 |
| **RMI2** | 34 | 0 | 0 | 1 | 3.990 | 1.041 | 4.109 | 3.470 | 4.494 | 1.024 | 0.959 | 5.892 |
| **PML** | 34 | 0 | 0 | 0 | 4.323 | 0.748 | 4.477 | 3.838 | 4.821 | 0.983 | 2.618 | 5.762 |
| **EXD2** | 34 | 0 | 1 | 1 | 3.662 | 0.936 | 3.855 | 3.389 | 4.171 | 0.782 | 0.096 | 5.164 |
| **MUS81** | 34 | 0 | 0 | 0 | 4.679 | 0.868 | 4.555 | 4.000 | 5.165 | 1.165 | 3.456 | 6.933 |
| **SLX4** | 34 | 0 | 0 | 7 | 1.601 | 0.627 | 1.417 | 1.158 | 1.941 | 0.783 | 0.771 | 3.139 |
| **FANCM** | 34 | 0 | 0 | 1 | 2.491 | 0.847 | 2.523 | 1.977 | 3.139 | 1.162 | 0.989 | 3.848 |
| **WRN** | 34 | 0 | 0 | 0 | 3.355 | 0.598 | 3.335 | 2.972 | 3.616 | 0.644 | 2.076 | 4.643 |
| **ATRX** | 34 | 0 | 0 | 0 | 4.622 | 0.947 | 4.669 | 4.211 | 5.261 | 1.050 | 1.719 | 6.882 |
| **DAXX** | 34 | 0 | 0 | 0 | 5.559 | 0.508 | 5.396 | 5.181 | 6.012 | 0.831 | 4.604 | 6.495 |
| **TERT** | 34 | 15 | 27 | 30 | 0.389 | 0.846 | 0.013 | 0.000 | 0.162 | 0.162 | 0.000 | 3.397 |
| **SMARCAL1** | 34 | 0 | 1 | 1 | 3.580 | 0.871 | 3.809 | 3.337 | 4.063 | 0.726 | 0.014 | 4.825 |
| **ZRANB3** | 34 | 0 | 0 | 0 | 2.398 | 0.444 | 2.422 | 2.122 | 2.690 | 0.568 | 1.356 | 3.249 |

Quartiles are `numpy.quantile` with default linear interpolation; SD uses `ddof = 1`.

Only one candidate has any line at exactly 0.000000: **TERT**, with
15 of 34 at zero and 27 of 34 below 0.5. Two further genes have a single
sub-1.0 line — **SMARCAL1** (min 0.0141) and **EXD2** (min 0.0962) — and
**BLM** has 5 and **SLX4** 7 lines below 1.0. Everything else is bounded away from the
quantification floor across the whole reference set.

---

## 2. Shape: unimodal or bimodal, and the gap diagnostic

The diagnostic is the one iteration 01 used on WRAP53 (where a 2.041-unit gap at 3.3×
the next-largest gap and 1.58 SD identified a discrete two-line spike): sort the 34
values, take all 33 consecutive differences, and report the largest, its position
(`k below` = how many lines fall below it), the gap divided by the SD, and the ratio to
the next-largest gap. Because a largest gap at `k = 1` or `k = 33` is a single-point
outlier rather than a mode boundary, the **balanced gap** is also reported: the largest
gap with at least 4 lines on *both* sides, which is the only kind that could define a
usable low group under `src.stratify.MIN_CELL_SIZE = 4`.

| Gene | largest gap | gap / SD | × next gap | k below | balanced gap (≥4 / ≥4) | / SD | × next eligible | k below | low-group size at the balanced gap | skew | ex. kurt | Shapiro p |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| **PRIMPOL** | 0.494 | 0.69 | 1.72 | 1 | 0.287 | 0.40 | 1.35 | 23 | 23 | 0.07 | -0.64 | 7.6e-01 |
| **BLM** | 0.815 | 0.64 | 1.07 | 6 | 0.815 | 0.64 | 1.07 | 6 | 6 | -1.10 | 0.02 | 1.7e-04 |
| **RPA1** | 0.421 | 0.60 | 1.77 | 1 | 0.238 | 0.34 | 1.13 | 16 | 16 | 0.04 | -0.91 | 5.0e-01 |
| **RPA2** | 0.822 | 1.09 | 1.31 | 33 | 0.259 | 0.34 | 1.24 | 28 | 28 | 0.11 | 0.42 | 1.0e+00 |
| **RPA3** | 0.471 | 0.48 | 1.52 | 2 | 0.297 | 0.30 | 1.00 | 4 | 4 | -0.73 | -0.13 | 6.7e-02 |
| **POLA1** | 0.582 | 0.67 | 1.08 | 3 | 0.541 | 0.63 | 1.94 | 4 | 4 | -1.51 | 1.77 | 1.1e-04 |
| **TOP3A** | 0.622 | 0.61 | 1.56 | 32 | 0.338 | 0.33 | 1.92 | 13 | 13 | 0.05 | -0.44 | 7.3e-01 |
| **RMI1** | 1.075 | 1.17 | 2.11 | 1 | 0.509 | 0.55 | 1.22 | 4 | 4 | -0.60 | 1.28 | 2.3e-01 |
| **RMI2** | 1.253 | 1.20 | 3.35 | 1 | 0.374 | 0.36 | 1.11 | 4 | 4 | -0.63 | 1.10 | 4.0e-01 |
| **PML** | 0.473 | 0.63 | 1.20 | 2 | 0.220 | 0.29 | 1.11 | 11 | 11 | -0.53 | -0.13 | 3.5e-01 |
| **EXD2** | 1.846 | 1.97 | 3.80 | 1 | 0.359 | 0.38 | 1.64 | 4 | 4 | -1.83 | 5.55 | 4.9e-04 |
| **MUS81** | 0.440 | 0.51 | 1.03 | 31 | 0.362 | 0.42 | 1.16 | 29 | 29 | 0.91 | 0.49 | 3.4e-02 |
| **SLX4** | 0.349 | 0.56 | 1.47 | 33 | 0.238 | 0.38 | 1.35 | 28 | 28 | 0.89 | -0.08 | 1.0e-02 |
| **FANCM** | 0.233 | 0.28 | 1.16 | 27 | 0.233 | 0.28 | 1.31 | 27 | 27 | -0.03 | -1.02 | 2.7e-01 |
| **WRN** | 0.373 | 0.62 | 1.36 | 2 | 0.270 | 0.45 | 1.46 | 27 | 27 | 0.17 | 0.01 | 6.2e-01 |
| **ATRX** | 1.276 | 1.35 | 1.04 | 33 | 0.429 | 0.45 | 1.72 | 5 | 5 | -0.73 | 2.05 | 8.5e-02 |
| **DAXX** | 0.282 | 0.56 | 1.16 | 21 | 0.282 | 0.56 | 2.16 | 21 | 21 | 0.25 | -1.07 | 5.8e-02 |
| **TERT** | 1.381 | 1.63 | 1.82 | 32 | 0.757 | 0.90 | 2.29 | 30 | 30 | 2.66 | 6.68 | 2.9e-09 |
| **SMARCAL1** | 1.835 | 2.11 | 2.18 | 1 | 0.247 | 0.28 | 1.67 | 7 | 7 | -2.29 | 7.87 | 4.3e-05 |
| **ZRANB3** | 0.275 | 0.62 | 1.43 | 1 | 0.192 | 0.43 | 1.03 | 4 | 4 | -0.18 | -0.08 | 9.7e-01 |

### Sorted values, 34 lines, ModelIDs withheld (rule b0)

```
PRIMPOL   2.19 2.68 2.72 2.82 2.97 2.98 2.98 2.99 3.16 3.17 3.25 3.37 3.42 3.54 3.54 3.69 3.69 3.71 3.72 3.73 3.83 3.85 3.85 4.14 4.17 4.23 4.29 4.38 4.59 4.66 4.69 4.71 4.85 5.09
BLM       0.46 0.57 0.79 0.81 0.96 1.28 2.09 2.32 3.09 3.11 3.21 3.28 3.30 3.31 3.35 3.39 3.47 3.68 3.77 3.84 3.94 3.95 3.99 4.00 4.04 4.10 4.15 4.16 4.16 4.19 4.26 4.50 4.73 4.79
RPA1      4.73 5.15 5.22 5.24 5.25 5.35 5.49 5.59 5.62 5.65 5.66 5.70 5.73 5.76 5.92 5.92 6.16 6.17 6.31 6.37 6.44 6.50 6.58 6.60 6.61 6.71 6.73 6.73 6.81 7.02 7.13 7.16 7.29 7.50
RPA2      4.37 5.00 5.03 5.16 5.37 5.41 5.47 5.54 5.64 5.69 5.70 5.79 5.92 5.93 6.07 6.08 6.12 6.14 6.20 6.24 6.38 6.42 6.42 6.43 6.58 6.58 6.60 6.69 6.95 6.96 7.09 7.17 7.24 8.06
RPA3      3.64 3.72 4.19 4.19 4.49 4.79 5.04 5.23 5.28 5.50 5.54 5.60 5.75 5.78 5.82 6.00 6.00 6.15 6.28 6.31 6.34 6.36 6.36 6.40 6.53 6.55 6.64 6.77 6.83 6.86 6.88 6.89 7.21 7.49
POLA1     1.17 1.61 1.65 2.23 2.77 3.05 3.06 3.25 3.33 3.50 3.50 3.50 3.63 3.76 3.82 3.86 3.90 3.92 3.93 3.95 4.03 4.05 4.15 4.20 4.25 4.28 4.29 4.29 4.32 4.33 4.40 4.42 4.47 4.58
TOP3A     2.00 2.37 2.77 2.77 2.83 2.97 3.13 3.14 3.19 3.23 3.34 3.35 3.52 3.86 4.03 4.09 4.15 4.17 4.20 4.30 4.40 4.56 4.63 4.64 4.66 4.78 4.80 4.94 4.95 5.03 5.04 5.34 5.96 6.27
RMI1      1.61 2.69 2.69 2.86 3.37 3.41 3.62 3.77 3.82 3.94 3.97 3.99 4.03 4.28 4.32 4.33 4.39 4.41 4.43 4.44 4.49 4.54 4.61 4.66 4.71 4.73 4.75 4.78 4.79 5.21 5.24 5.51 5.82 6.17
RMI2      0.96 2.21 2.30 2.59 2.96 3.12 3.45 3.46 3.47 3.48 3.57 3.62 3.74 3.81 3.90 4.00 4.08 4.14 4.19 4.21 4.25 4.36 4.36 4.39 4.39 4.53 4.81 4.83 4.91 5.24 5.28 5.51 5.66 5.89
PML       2.62 2.71 3.18 3.36 3.39 3.49 3.59 3.64 3.81 3.93 3.94 4.15 4.25 4.30 4.32 4.43 4.47 4.48 4.61 4.63 4.64 4.68 4.69 4.73 4.75 4.84 4.87 4.90 4.92 5.12 5.17 5.24 5.37 5.76
EXD2      0.10 1.94 2.18 2.66 3.02 3.11 3.11 3.12 3.34 3.53 3.56 3.61 3.63 3.69 3.72 3.75 3.85 3.86 3.89 3.90 4.00 4.01 4.01 4.02 4.13 4.19 4.20 4.21 4.34 4.45 4.62 4.66 4.92 5.16
MUS81     3.46 3.54 3.59 3.75 3.83 3.87 3.88 3.89 4.00 4.01 4.25 4.26 4.31 4.33 4.47 4.52 4.55 4.56 4.57 4.66 4.67 4.76 4.85 4.92 5.12 5.18 5.23 5.23 5.31 5.67 5.99 6.43 6.51 6.93
SLX4      0.77 0.83 0.95 0.97 0.97 0.98 1.00 1.12 1.15 1.17 1.18 1.22 1.26 1.27 1.38 1.39 1.41 1.42 1.47 1.52 1.63 1.65 1.73 1.78 1.91 1.95 1.97 2.12 2.36 2.54 2.70 2.72 2.79 3.14
FANCM     0.99 1.07 1.18 1.38 1.42 1.48 1.64 1.81 1.97 2.00 2.02 2.03 2.08 2.16 2.25 2.32 2.49 2.55 2.59 2.68 2.71 2.86 2.88 3.00 3.01 3.18 3.28 3.51 3.52 3.55 3.66 3.74 3.84 3.85
WRN       2.08 2.21 2.58 2.65 2.83 2.90 2.92 2.93 2.96 3.01 3.08 3.09 3.13 3.15 3.26 3.27 3.33 3.34 3.37 3.37 3.39 3.46 3.50 3.50 3.60 3.62 3.72 3.99 4.00 4.17 4.29 4.36 4.37 4.64
ATRX      1.72 2.95 3.22 3.34 3.50 3.93 4.05 4.09 4.21 4.23 4.25 4.50 4.53 4.56 4.62 4.63 4.67 4.67 4.77 4.79 4.93 4.94 5.18 5.18 5.21 5.28 5.32 5.39 5.42 5.46 5.55 5.59 5.61 6.88
DAXX      4.60 4.85 4.89 4.98 5.02 5.15 5.16 5.17 5.17 5.21 5.26 5.27 5.30 5.32 5.33 5.33 5.36 5.44 5.45 5.46 5.55 5.83 5.85 5.97 6.01 6.01 6.06 6.13 6.17 6.27 6.30 6.32 6.33 6.50
TERT      0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.01 0.01 0.02 0.02 0.02 0.03 0.05 0.09 0.11 0.14 0.17 0.50 0.53 0.77 0.92 1.68 1.70 3.08 3.40
SMARCAL1  0.01 1.85 2.69 2.74 2.89 3.01 3.06 3.31 3.31 3.41 3.44 3.45 3.53 3.54 3.65 3.73 3.80 3.82 3.83 3.86 3.88 3.91 3.98 3.98 3.98 4.09 4.14 4.17 4.21 4.28 4.35 4.48 4.51 4.83
ZRANB3    1.36 1.63 1.71 1.75 1.95 1.98 2.00 2.11 2.11 2.15 2.16 2.22 2.26 2.30 2.39 2.40 2.41 2.44 2.48 2.49 2.50 2.50 2.55 2.65 2.67 2.70 2.72 2.77 2.80 2.82 3.01 3.07 3.21 3.25
```

### Read plainly

**One gene is genuinely bimodal: BLM.** Its low group is
0.46 0.57 0.79 0.81 0.96 1.28, then a 0.815 gap to 2.09 — and the next-largest gap
(0.763) sits two lines further up, so the region between 1.28 and 3.09 holds only two
lines. That is **a natural low group of 6 (or 8, taking the second gap) out of 34**, the
only balanced gap in the panel that falls at a sparse region rather than inside a dense
body. It is also the only gene whose largest *and* balanced gap coincide.

**Three genes are unimodal with a single extreme low outlier**, which is a different
thing and must not be read as a low group:

- **SMARCAL1** — one line at 0.0141, then 1.849: a 1.835-unit gap at
  2.11 SD and 2.18× the next gap, with **n = 1** below it. Above that line the
  distribution is dense and unimodal (balanced gap only 0.247 = 0.28 SD).
- **EXD2** — one line at 0.0962, then 1.942: gap 1.846 = 1.97 SD, 3.80× the next,
  n = 1 below. This single point produces EXD2's excess kurtosis of
  5.55 and its Shapiro rejection.
- **RMI2** and **RMI1** — one line each at 0.96 and 1.61, gaps of 1.253 and
  1.075 with n = 1 below; dense thereafter.

**Two genes are unimodal with a single extreme high outlier**: ATRX (top value 6.882,
gap 1.276 at k = 33) and RPA2 (8.060, gap 0.822 at k = 33).

**TERT is neither unimodal nor bimodal — it is a zero spike with a thin tail.**
15 lines at exactly 0, 27 below 0.5, then 6 lines scattered from 0.50 to 3.40. Its
largest gap (1.381) is at k = 32, i.e. between the two highest lines. Shapiro W =
0.532, p = 2.9e-09.

**The remaining 13 genes are unimodal with no interior structure.** Their balanced gaps
run 0.192–0.541 units, i.e. 0.28–0.63 SD, with ratios to the next eligible
gap of 1.00–2.16. Eleven of the 13 are below 1.5, i.e. the widest interior gap is
indistinguishable from the second widest, so nothing in the data marks a boundary and
any low group defined on them would be a percentile artefact of the kind iteration 01
warned against. The two exceptions are **POLA1** (ratio 1.94 at 0.63 SD, the 4-line low
tail discussed above and in §5 — X-linked, so it tracks hemizygosity) and **DAXX**
(ratio 2.16 but a gap of only 0.282 units, on the narrowest distribution in the panel;
see §3). Neither is a mode boundary of the BLM kind.

**ZRANB3**, the nominated negative control, belongs squarely in that last group and is
the most featureless distribution in the panel: no zeros, skew
-0.18, excess kurtosis -0.08, Shapiro p = 0.97, largest gap
0.275 (0.62 SD, 1.43× the next) sitting at k = 1, balanced gap
0.192 (0.43 SD). All 34 values are distinct and all 16 ALT+ values are
distinct — no ties anywhere — which matters for a rank-based control (see the
negative-control note at the end of the ranking).

**Does any gene have a natural low group of size ≥ 4?** Over the 34-line reference set:
**yes, exactly one — BLM (n = 6, or 8 at the second gap).** POLA1 comes closest among the
rest, with 4 lines below a 0.541-unit gap (0.63 SD), but see §5: POLA1 is X-linked and
that low tail tracks hemizygosity and sex, not a focal lesion. No other gene has a
≥ 4-line group behind any gap wider than 0.63 SD.

---

## 3. Dynamic range and usability as a continuous exposure

Range and CV over the 34-line reference set, and the same over the 16 ALT+ lines that
would actually carry the regression. `floor` counts lines at or near the quantification
floor. A continuous exposure needs (i) no zero-inflation, (ii) spread large relative to
measurement noise, (iii) no single point carrying the fit.

| Gene | range 34 | SD 34 | CV 34 | range 16 | SD 16 | CV 16 | IQR 16 | n = 0 (16) | n < 1 (16) | dominant single point? | flag |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|
| **PRIMPOL** | 2.90 | 0.710 | 0.192 | 2.17 | 0.761 | 0.209 | 1.321 | 0 | 0 | no | — |
| **BLM** | 4.33 | 1.268 | 0.395 | 3.83 | 0.895 | 0.249 | 0.814 | 0 | 1 | no | bimodal over 34; only 2 low lines in the 16 |
| **RPA1** | 2.77 | 0.707 | 0.115 | 2.01 | 0.651 | 0.106 | 0.980 | 0 | 0 | no | narrow (CV 0.106) |
| **RPA2** | 3.69 | 0.755 | 0.123 | 2.59 | 0.652 | 0.102 | 0.640 | 0 | 0 | yes — 1 high line, gap 0.973 (2.67× next) | narrow (CV 0.102) |
| **RPA3** | 3.85 | 0.988 | 0.168 | 3.30 | 0.722 | 0.114 | 0.573 | 0 | 0 | yes — 1 low line, gap 1.348 (2.25× next) | spread in the 16 is one low line |
| **POLA1** | 3.41 | 0.862 | 0.237 | 1.52 | 0.503 | 0.131 | 0.809 | 0 | 0 | no | X-linked; narrow in the 16 |
| **TOP3A** | 4.27 | 1.012 | 0.250 | 3.30 | 0.989 | 0.225 | 1.264 | 0 | 0 | no | — |
| **RMI1** | 4.56 | 0.918 | 0.216 | 1.99 | 0.572 | 0.124 | 0.538 | 0 | 0 | no | narrow |
| **RMI2** | 4.93 | 1.041 | 0.261 | 2.93 | 0.760 | 0.196 | 0.732 | 0 | 0 | no | — |
| **PML** | 3.14 | 0.748 | 0.173 | 3.05 | 0.835 | 0.196 | 1.164 | 0 | 0 | no | — |
| **EXD2** | 5.07 | 0.936 | 0.256 | 4.83 | 1.128 | 0.299 | 0.714 | 0 | 1 | yes — 1 low line, gap 2.567 (3.78× next) | **single-point outlier inside the 16** |
| **MUS81** | 3.48 | 0.868 | 0.185 | 3.05 | 0.837 | 0.179 | 1.258 | 0 | 0 | yes — 1 high line, gap 0.832 (1.52× next) | — |
| **SLX4** | 2.37 | 0.627 | 0.392 | 2.37 | 0.678 | 0.414 | 0.790 | 0 | 4 | no | near the floor (median 1.47) |
| **FANCM** | 2.86 | 0.847 | 0.340 | 2.43 | 0.792 | 0.273 | 1.095 | 0 | 0 | no | — |
| **WRN** | 2.57 | 0.598 | 0.178 | 2.57 | 0.662 | 0.200 | 0.526 | 0 | 0 | no | — |
| **ATRX** | 5.16 | 0.947 | 0.205 | 5.16 | 1.204 | 0.264 | 1.171 | 0 | 0 | no | outlier-driven at both ends; X-linked |
| **DAXX** | 1.89 | 0.508 | 0.091 | 1.34 | 0.468 | 0.084 | 0.796 | 0 | 0 | no | **narrowest range in the panel** |
| **TERT** | 3.40 | 0.846 | 2.172 | 0.17 | 0.054 | 1.792 | 0.022 | 8 | 16 | no | **zero-inflated / degenerate** |
| **SMARCAL1** | 4.81 | 0.871 | 0.243 | 2.13 | 0.607 | 0.167 | 0.873 | 0 | 0 | yes — 1 high line, gap 0.541 (1.81× next) | — |
| **ZRANB3** | 1.89 | 0.444 | 0.185 | 1.58 | 0.427 | 0.180 | 0.431 | 0 | 0 | no | **smallest SD in the panel**; nominated negative control |

**Flagged as unusable or marginal on range/zero-inflation grounds alone:**

- **TERT** cannot serve. Over the 16 ALT+ lines it is 8 exact zeros, an IQR of
  0.022, a range of 0.168 and a maximum of 0.168. The entire analysis set sits
  inside the bottom 5% of the gene's reference range. This is a constant with noise, not
  a variable. (Biologically unsurprising — the set is ALT+ — but the consequence is
  statistical, not biological.)
- **DAXX** has the smallest range and the smallest CV in the panel: 1.34 and 0.084
  over the 16 (versus 3.30 / 0.225 for TOP3A). It is not degenerate, but a slope
  would be estimated across a 1.3-unit window, so it is the weakest continuous candidate
  that is not outright unusable.
- **EXD2** has a single line at 0.096 with the next at 2.663 *inside the ALT+ set*.
  As a continuous exposure it is effectively a 1-vs-15 indicator; one line would carry
  the fit (excess kurtosis 7.89 in the 16).
- **SLX4** is the lowest-expressed non-degenerate candidate (median 1.474 in the 16,
  4 of 16 below 1.0, 7 of 34 below 1.0). Its CV is high (0.414) precisely
  because the mean is near the floor, where log2(TPM+1) is least reliable.

No candidate has near-zero variance in the strict sense; the smallest SD over the 16 is
TERT's 0.054 (degenerate) and then DAXX's 0.468.

---

## 4. Restricted to the 16 ALT+ CRISPR-screened lines

Exposure side only — this touches no outcome. The 16 were identified from the frozen
manifest (`alt_call_direction == 'positive'`) intersected with `crispr_model_ids()`,
which reads row labels and no gene-effect value. **ModelIDs are withheld throughout this
section** (rule b0); sorted vectors are given so the prespec can see the realised shape.

| Gene | n | n = 0 | mean | SD | median | Q1 | Q3 | IQR | min | max | range | largest gap | / SD | k below | balanced gap | k below | Shapiro p |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| **PRIMPOL** | 16 | 0 | 3.648 | 0.761 | 3.480 | 2.988 | 4.308 | 1.321 | 2.680 | 4.853 | 2.17 | 0.454 | 0.60 | 10 | 0.454 | 10 | 1.2e-01 |
| **BLM** | 16 | 0 | 3.591 | 0.895 | 3.805 | 3.341 | 4.155 | 0.814 | 0.962 | 4.793 | 3.83 | 1.361 | 1.52 | 1 | 0.220 | 6 | 5.9e-03 |
| **RPA1** | 16 | 0 | 6.155 | 0.651 | 6.049 | 5.654 | 6.635 | 0.980 | 5.152 | 7.163 | 2.01 | 0.291 | 0.45 | 13 | 0.265 | 10 | 2.5e-01 |
| **RPA2** | 16 | 0 | 6.391 | 0.652 | 6.293 | 6.035 | 6.674 | 0.640 | 5.470 | 8.060 | 2.59 | 0.973 | 1.49 | 15 | 0.365 | 12 | 3.3e-01 |
| **RPA3** | 16 | 0 | 6.310 | 0.722 | 6.346 | 6.116 | 6.690 | 0.573 | 4.194 | 7.493 | 3.30 | 1.348 | 1.87 | 1 | 0.190 | 12 | 2.4e-02 |
| **POLA1** | 16 | 0 | 3.855 | 0.503 | 3.935 | 3.455 | 4.264 | 0.809 | 3.053 | 4.576 | 1.52 | 0.394 | 0.78 | 6 | 0.394 | 6 | 1.7e-01 |
| **TOP3A** | 16 | 0 | 4.394 | 0.989 | 4.347 | 3.705 | 4.969 | 1.264 | 2.967 | 6.268 | 3.30 | 0.631 | 0.64 | 4 | 0.631 | 4 | 6.4e-01 |
| **RMI1** | 16 | 0 | 4.612 | 0.572 | 4.427 | 4.310 | 4.848 | 0.538 | 3.824 | 5.819 | 1.99 | 0.486 | 0.85 | 12 | 0.486 | 12 | 2.5e-01 |
| **RMI2** | 16 | 0 | 3.883 | 0.760 | 3.855 | 3.459 | 4.191 | 0.732 | 2.587 | 5.512 | 2.93 | 0.856 | 1.13 | 14 | 0.256 | 6 | 5.8e-01 |
| **PML** | 16 | 0 | 4.257 | 0.835 | 4.379 | 3.537 | 4.701 | 1.164 | 2.708 | 5.762 | 3.05 | 0.569 | 0.68 | 5 | 0.569 | 5 | 9.1e-01 |
| **EXD2** | 16 | 0 | 3.779 | 1.128 | 4.009 | 3.652 | 4.366 | 0.714 | 0.096 | 4.923 | 4.83 | 2.567 | 2.27 | 1 | 0.169 | 9 | 4.2e-04 |
| **MUS81** | 16 | 0 | 4.666 | 0.837 | 4.565 | 3.969 | 5.228 | 1.258 | 3.456 | 6.507 | 3.05 | 0.832 | 0.99 | 15 | 0.546 | 9 | 3.8e-01 |
| **SLX4** | 16 | 0 | 1.638 | 0.678 | 1.474 | 1.128 | 1.919 | 0.790 | 0.771 | 3.139 | 2.37 | 0.583 | 0.86 | 13 | 0.211 | 5 | 1.6e-01 |
| **FANCM** | 16 | 0 | 2.902 | 0.792 | 3.091 | 2.435 | 3.530 | 1.095 | 1.422 | 3.848 | 2.43 | 0.478 | 0.60 | 4 | 0.478 | 4 | 1.8e-01 |
| **WRN** | 16 | 0 | 3.301 | 0.662 | 3.297 | 2.999 | 3.525 | 0.526 | 2.076 | 4.643 | 2.57 | 0.693 | 1.05 | 2 | 0.130 | 7 | 6.5e-01 |
| **ATRX** | 16 | 0 | 4.564 | 1.204 | 4.616 | 4.168 | 5.340 | 1.171 | 1.719 | 6.882 | 5.16 | 1.276 | 1.06 | 15 | 0.314 | 7 | 6.8e-01 |
| **DAXX** | 16 | 0 | 5.575 | 0.468 | 5.444 | 5.246 | 6.042 | 0.796 | 4.980 | 6.318 | 1.34 | 0.468 | 1.00 | 11 | 0.468 | 11 | 2.8e-02 |
| **TERT** | 16 | 8 | 0.030 | 0.054 | 0.003 | 0.000 | 0.022 | 0.022 | 0.000 | 0.168 | 0.17 | 0.063 | 1.15 | 13 | 0.012 | 12 | 2.8e-05 |
| **SMARCAL1** | 16 | 0 | 3.626 | 0.607 | 3.594 | 3.237 | 4.110 | 0.873 | 2.691 | 4.825 | 2.13 | 0.541 | 0.89 | 15 | 0.299 | 4 | 8.4e-01 |
| **ZRANB3** | 16 | 0 | 2.369 | 0.427 | 2.405 | 2.112 | 2.542 | 0.431 | 1.631 | 3.214 | 1.58 | 0.271 | 0.63 | 2 | 0.206 | 12 | 9.4e-01 |

### Sorted values, 16 ALT+ lines, ModelIDs withheld

```
PRIMPOL   2.68 2.72 2.82 2.98 2.99 3.16 3.25 3.42 3.54 3.71 4.17 4.29 4.38 4.69 4.71 4.85
BLM       0.96 2.32 3.30 3.31 3.35 3.47 3.68 3.77 3.84 3.94 3.95 4.15 4.16 4.19 4.26 4.79
RPA1      5.15 5.35 5.62 5.65 5.66 5.70 5.73 5.92 6.17 6.31 6.58 6.60 6.73 7.02 7.13 7.16
RPA2      5.47 5.64 5.70 5.93 6.07 6.08 6.14 6.20 6.38 6.43 6.58 6.58 6.95 6.96 7.09 8.06
RPA3      4.19 5.54 6.00 6.00 6.15 6.28 6.31 6.34 6.36 6.53 6.55 6.64 6.83 6.86 6.89 7.49
POLA1     3.05 3.06 3.25 3.33 3.50 3.50 3.90 3.92 3.95 4.15 4.20 4.25 4.29 4.33 4.42 4.58
TOP3A     2.97 3.14 3.19 3.23 3.86 4.03 4.20 4.30 4.40 4.63 4.80 4.94 5.04 5.34 5.96 6.27
RMI1      3.82 3.94 3.99 4.28 4.32 4.33 4.39 4.41 4.44 4.66 4.71 4.73 5.21 5.24 5.51 5.82
RMI2      2.59 2.96 3.12 3.45 3.46 3.48 3.74 3.81 3.90 4.00 4.08 4.19 4.21 4.39 5.24 5.51
PML       2.71 3.18 3.36 3.39 3.59 4.15 4.25 4.32 4.43 4.47 4.48 4.64 4.87 5.12 5.37 5.76
EXD2      0.10 2.66 3.34 3.53 3.69 3.85 3.89 4.00 4.02 4.19 4.21 4.34 4.45 4.62 4.66 4.92
MUS81     3.46 3.75 3.87 3.88 4.00 4.01 4.33 4.56 4.57 5.12 5.18 5.23 5.23 5.31 5.67 6.51
SLX4      0.77 0.95 0.97 1.00 1.17 1.38 1.39 1.42 1.52 1.63 1.73 1.91 1.95 2.54 2.72 3.14
FANCM     1.42 1.64 2.00 2.08 2.55 2.59 2.68 3.00 3.18 3.28 3.51 3.52 3.55 3.74 3.84 3.85
WRN       2.08 2.21 2.90 2.96 3.01 3.08 3.13 3.26 3.33 3.37 3.46 3.50 3.60 4.00 4.29 4.64
ATRX      1.72 2.95 3.50 4.05 4.21 4.23 4.25 4.56 4.67 4.93 5.18 5.32 5.39 5.59 5.61 6.88
DAXX      4.98 5.02 5.17 5.21 5.26 5.32 5.33 5.44 5.45 5.46 5.55 6.01 6.13 6.27 6.30 6.32
TERT      0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.01 0.01 0.02 0.02 0.03 0.09 0.14 0.17
SMARCAL1  2.69 2.74 2.89 3.01 3.31 3.44 3.53 3.54 3.65 3.80 3.83 4.09 4.17 4.21 4.28 4.83
ZRANB3    1.63 1.71 1.98 2.11 2.11 2.15 2.30 2.40 2.41 2.44 2.48 2.49 2.70 2.77 3.01 3.21
```

### What changes between the 34 and the 16

This is the part the prespec needs, because three of the four interesting features of
the reference distributions **do not survive restriction to the analysis set**:

- **BLM's low group collapses.** 6 of 34 reference lines sit below 1.28; only
  1 of the 16 ALT+ lines does. In the 16, BLM is 0.962, 2.323, then 14 lines
  packed into 1.496 units. The one real bimodality in the panel is almost entirely
  outside the analysis set.
- **SMARCAL1's near-zero line is not in the ALT+ set.** min over 34 = 0.0141;
  min over 16 = 2.691. Within the 16, SMARCAL1 is the best-behaved variable in
  the panel: range 2.691–4.825, SD 0.607, skew 0.11, excess kurtosis
  -0.52, Shapiro p = 0.84, no gap above 0.299 with 4 lines either side.
- **POLA1's low tail is not in the ALT+ set either.** min over 34 = 1.167, min over
  16 = 3.053; the 16-line range shrinks to 1.52, the narrowest of any
  non-degenerate candidate except DAXX.
- **EXD2's outlier *is* in the ALT+ set** (0.096, next 2.663) — the one case where the
  reference-set pathology transfers intact.

Genes whose realised 16-line spread is widest, and therefore have the most signal to
regress against anything: **ATRX** (SD 1.204, range 5.16), **EXD2** (SD 1.128, range 4.83), **TOP3A** (SD 0.989, range 3.30), **BLM** (SD 0.895, range 3.83), **MUS81** (SD 0.837, range 3.05).

The only genes with an interior balanced gap in the 16 that is both ≥ 0.45 units and
≥ 2× the next eligible gap — i.e. the only ones that could support a *data-defined*
binary split with both arms ≥ 4 — are:

| Gene | balanced gap | / SD | × next eligible | arms (low / high) | boundary |
|---|---:|---:|---:|---:|---|
| **TOP3A** | 0.631 | 0.64 | 2.68 | 4 / 12 | 3.232 → 3.863 |
| **PML** | 0.569 | 0.68 | 2.52 | 5 / 11 | 3.586 → 4.155 |
| **RMI1** | 0.486 | 0.85 | 2.29 | 12 / 4 | 4.727 → 5.212 |
| **DAXX** | 0.468 | 1.00 | 4.09 | 11 / 5 | 5.545 → 6.014 |
| **PRIMPOL** | 0.454 | 0.60 | 2.66 | 10 / 6 | 3.714 → 4.167 |

These are offered as a description of the exposure's shape, not as a recommendation:
every one is a cut through a near-continuous body, and none is a discontinuity of the
WRAP53-zero kind. For comparison, BLM's genuine mode boundary in the 34-line reference
set sits at 0.815 units and 0.64 SD with only two lines either side of it inside the 16.

---

## 5. Copy number: apparent loss, and whether segments corroborate it

Binding rule from iteration 02: a gene-level value in `OmicsCNGeneWGS.csv` must be
validated against `OmicsCNSegmentsWGS.csv` before it is used, because three of four
apparent homozygous ATRX deletions were refuted that way. Iteration 03 §4.3 established
the cohort-wide prior: all 16 ALT+ lines carry 2–105 focal sub-2 Mb near-zero segments
(median 11 per line), so a gene-level zero in this cohort is neither rare nor credible.

### 5.1 Gene-level CN over the reference set (33 lines with CN) and the ALT+16

| Gene | min 34 | median 34 | max 34 | <0.1 | <0.3 | <0.6 | >1.5 | min 16 | <0.1 | <0.3 | <0.6 | >1.5 | expr vs own CN (rho, 33 lines) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| **PRIMPOL** | 0.3588 | 0.802 | 1.508 | 0 | 0 | 7 | 1 | 0.3588 | 0 | 0 | 7 | 1 | +0.592 (p = 2.8e-04) |
| **BLM** | 0.8531 | 1.124 | 2.312 | 0 | 0 | 0 | 8 | 0.8531 | 0 | 0 | 0 | 7 | +0.485 (p = 4.3e-03) |
| **RPA1** | 0.6177 | 0.999 | 1.841 | 0 | 0 | 0 | 3 | 0.6680 | 0 | 0 | 0 | 1 | +0.548 (p = 9.5e-04) |
| **RPA2** | 0.8945 | 1.085 | 1.707 | 0 | 0 | 0 | 2 | 0.8945 | 0 | 0 | 0 | 2 | +0.570 (p = 5.4e-04) |
| **RPA3** | 0.5462 | 1.165 | 1.738 | 0 | 0 | 1 | 5 | 0.5462 | 0 | 0 | 1 | 4 | +0.659 (p = 3.0e-05) |
| **POLA1** | 0.3159 | 0.611 | 1.183 | 0 | 0 | 16 | 0 | 0.3689 | 0 | 0 | 7 | 0 | +0.392 (p = 2.4e-02) |
| **TOP3A** | 0.7103 | 1.265 | 3.042 | 0 | 0 | 0 | 11 | 0.7103 | 0 | 0 | 0 | 8 | +0.722 (p = 2.2e-06) |
| **RMI1** | 0.6390 | 1.065 | 2.342 | 0 | 0 | 0 | 3 | 0.6390 | 0 | 0 | 0 | 2 | +0.652 (p = 3.9e-05) |
| **RMI2** | 0.6636 | 0.951 | 1.430 | 0 | 0 | 0 | 0 | 0.6636 | 0 | 0 | 0 | 0 | +0.408 (p = 1.9e-02) |
| **PML** | 0.4632 | 1.010 | 1.586 | 0 | 0 | 1 | 2 | 0.4632 | 0 | 0 | 1 | 1 | +0.158 (p = 3.8e-01) |
| **EXD2** | 0.6125 | 1.082 | 1.825 | 0 | 0 | 0 | 5 | 0.8393 | 0 | 0 | 0 | 3 | +0.573 (p = 4.9e-04) |
| **MUS81** | 0.8492 | 1.067 | 2.370 | 0 | 0 | 0 | 5 | 0.8492 | 0 | 0 | 0 | 2 | +0.514 (p = 2.2e-03) |
| **SLX4** | 0.4848 | 0.942 | 1.659 | 0 | 0 | 1 | 1 | 0.4848 | 0 | 0 | 1 | 0 | +0.276 (p = 1.2e-01) |
| **FANCM** | 0.7021 | 1.139 | 1.854 | 0 | 0 | 0 | 5 | 0.7482 | 0 | 0 | 0 | 4 | +0.534 (p = 1.4e-03) |
| **WRN** | 0.4278 | 0.991 | 1.462 | 0 | 0 | 2 | 0 | 0.4278 | 0 | 0 | 2 | 0 | +0.453 (p = 8.2e-03) |
| **ATRX** | 0.0002 | 0.583 | 1.073 | 4 | 4 | 18 | 0 | 0.0002 | 4 | 4 | 8 | 0 | -0.123 (p = 5.0e-01) |
| **DAXX** | 0.6873 | 1.117 | 1.761 | 0 | 0 | 0 | 5 | 0.9073 | 0 | 0 | 0 | 2 | +0.450 (p = 8.6e-03) |
| **TERT** | 0.8194 | 1.066 | 2.074 | 0 | 0 | 0 | 2 | 0.8194 | 0 | 0 | 0 | 0 | +0.229 (p = 2.0e-01) |
| **SMARCAL1** | 0.5259 | 0.966 | 1.296 | 0 | 0 | 2 | 0 | 0.5346 | 0 | 0 | 1 | 0 | +0.716 (p = 2.8e-06) |
| **ZRANB3** | 0.5199 | 0.963 | 1.354 | 0 | 0 | 1 | 0 | 0.5199 | 0 | 0 | 1 | 0 | +0.256 (p = 1.5e-01) |

**Only one candidate shows any apparent homozygous or focal loss anywhere in the
reference set: ATRX**, in 4 lines, all 4 inside the ALT+16. This reproduces iterations
02 and 03 exactly (4 lines below 0.1, 4 below 0.3, 18 of 33 below 0.6). Every other
candidate's minimum gene-level CN across all 33 covered lines is ≥ 0.316 (POLA1), and
11 of the 20 never fall below 0.6 at all. **There is no second deletion event in this
panel to validate or refute.**

The sub-0.6 counts for ATRX (18/33) and POLA1 (16/33) are **X-linkage, not loss**: both
genes are on chrX (ATRX Xq21.1, POLA1 Xp22.11) and `OmicsCNGeneWGS.csv` is already
ploidy-normalised, so a hemizygous male line reads ~0.5. Checked directly against
`Model.csv` `Sex` over the 34 (17 Male / 13 Female / 4 Unknown):

| Gene | modality | Male median | Female median | Mann–Whitney p |
|---|---|---:|---:|---:|
| ATRX | gene-level CN | 0.581 | 0.760 | 0.032 |
| ATRX | log2(TPM+1) | 4.771 | 5.181 | 0.277 |
| POLA1 | gene-level CN | 0.611 | 0.887 | 0.127 |
| POLA1 | log2(TPM+1) | 3.504 | 4.202 | **0.017** |

POLA1 **expression** differs by sex at p = 0.017 over the 34 lines. Any iteration using
POLA1 expression as its exposure inherits sex as a confounder. ATRX expression does not
(p = 0.277).

### 5.2 Segment corroboration

Iteration 03 §4.2 noted that this repository holds no gene annotation, so the segment
check is a manual per-locus exercise. GRCh38 intervals were therefore hand-entered for
all 20 candidates (ATRX's taken verbatim from iteration 02's memo,
`chrX:77,504,878–77,786,233`) and then **validated**, by correlating each gene's
gene-level CN against the length-weighted mean `RelCopyRatio` of the segments spanning
its hand-entered interval, across the 33 CN-covered lines. A wrong interval would
decorrelate; every one is tight:

| Gene | GRCh38 interval used | Spearman(gene CN, spanning-segment ratio) | mean # spanning segments | min raw spanning-segment ratio over 33 lines | lines with any spanning segment < 0.1 |
|---|---|---:|---:|---:|---:|
| **PRIMPOL** | chr4:184,588,000–184,631,000 | +0.980 | 1.00 | 0.37184 | 0 |
| **BLM** | chr15:90,717,000–90,816,000 | +0.987 | 1.00 | 0.80701 | 0 |
| **RPA1** | chr17:1,733,000–1,803,000 | +0.983 | 1.03 | 0.52813 | 0 |
| **RPA2** | chr1:27,890,000–27,916,000 | +0.961 | 1.00 | 0.84297 | 0 |
| **RPA3** | chr7:7,636,000–7,730,000 | +0.973 | 1.03 | 0.52597 | 0 |
| **POLA1** | chrX:24,693,000–24,961,000 | +0.979 | 1.00 | 0.28678 | 0 |
| **TOP3A** | chr17:18,254,000–18,318,000 | +0.983 | 1.00 | 0.66920 | 0 |
| **RMI1** | chr9:83,264,000–83,290,000 | +0.970 | 1.00 | 0.63682 | 0 |
| **RMI2** | chr16:11,343,000–11,371,000 | +0.960 | 1.03 | 0.63463 | 0 |
| **PML** | chr15:73,994,000–74,048,000 | +0.982 | 1.03 | 0.43377 | 0 |
| **EXD2** | chr14:69,159,000–69,207,000 | +0.978 | 1.00 | 0.60756 | 0 |
| **MUS81** | chr11:65,834,000–65,841,000 | +0.940 | 1.00 | 0.83986 | 0 |
| **SLX4** | chr16:3,581,000–3,609,000 | +0.972 | 1.00 | 0.49957 | 0 |
| **FANCM** | chr14:45,133,000–45,215,000 | +0.987 | 1.06 | 0.67708 | 0 |
| **WRN** | chr8:31,033,000–31,176,000 | +0.986 | 1.00 | 0.42396 | 0 |
| **ATRX** | chrX:77,504,878–77,786,233 | +0.854 | 1.18 | 0.00000 | 1 |
| **DAXX** | chr6:33,318,000–33,322,000 | +0.944 | 1.00 | 0.66032 | 0 |
| **TERT** | chr5:1,253,000–1,296,000 | +0.955 | 1.03 | 0.78576 | 0 |
| **SMARCAL1** | chr2:216,413,000–216,492,000 | +0.944 | 1.03 | 0.48838 | 0 |
| **ZRANB3** | chr2:135,170,000–135,403,000 | +0.980 | 1.12 | 0.58799 | 0 |

Every interval validates at rho ≥ +0.85 (19 of 20 at ≥ +0.94, ZRANB3 at +0.980), and the mean number of
spanning segments is ~1 for every gene, so these loci sit inside single broad segments
in essentially every line. ATRX is the lowest at +0.854 — exactly because it is the one
gene where gene-level and segment-level disagree, which is the finding below.

**The raw segment-level minima confirm the gene-level picture independently.** Across
all 33 lines and all 20 loci, the only spanning segment anywhere below 0.1 is a single
ATRX one. No candidate has a focal near-zero segment that the gene-level file *missed*,
so there is no hidden deletion either.

### 5.3 Adjudication of the four apparent ATRX losses

Independently re-derived here; it reproduces iteration 02's verdict line for line.
ModelIDs are named only because iteration 02's memo already named them.

| Line | gene-level ATRX CN | spanning segment(s) | verdict |
|---|---:|---|---|
| `ACH-000364` | 0.004586 | 66,388,000–96,077,000 (29689 kb) ratio 0.670960 | **refuted** — locus sits inside one broad segment at near-hemizygous ratio |
| `ACH-001715` | 0.000182 | 55,687,000–77,555,000 (21868 kb) ratio 0.445461; 77,555,000–77,736,000 (181 kb) ratio 0.000001; 77,736,000–82,726,000 (4990 kb) ratio 0.436619 | **corroborated** — genuine focal intragenic deletion |
| `ACH-002067` | 0.000552 | 63,757,000–155,983,000 (92226 kb) ratio 0.406464 | **refuted** — locus sits inside one broad segment at near-hemizygous ratio |
| `ACH-002688` | 0.000347 | 65,430,000–112,810,000 (47380 kb) ratio 0.423750 | **refuted** — locus sits inside one broad segment at near-hemizygous ratio |

`ACH-001715` (CAL72 in iteration 02's memo) carries a 181 kb segment at
`RelCopyRatio = 1 × 10⁻⁶` immediately inside the ATRX interval, flanked within the same
interval by segments at 0.445 and 0.437 — a real focal intragenic deletion. The other
three lines' ATRX loci fall inside single segments of 29.7, 92.2 and 47.4 Mb at ratios
0.671, 0.406 and 0.424: hemizygous-looking arm- or chromosome-scale states, not
deletions. **1 of 4 corroborated, 3 of 4 refuted** — identical to iteration 02.

**Consequence for this iteration.** There is no usable copy-number exposure in this
candidate panel. The only validated event has n = 1. Every other CN axis that would
clear `MIN_CELL_SIZE = 4` is a threshold through the middle of a near-diploid continuum,
which iteration 03 §4.1 already ruled out as an arbitrary cut with no biological event
behind it. Copy number remains useful here only as a *corroborating* channel for an
expression exposure, and even then weakly: expression-vs-own-CN rho runs from +0.16
(PML) to +0.72 (TOP3A), with SMARCAL1 at +0.716 — so for the dosage-sensitive genes, CN
is not an independent confirmation of expression but largely the same variable.

---

## 6. Collinearity among the candidate exposures

Spearman rho, expression vs expression only, over the 34-line reference set. No gene
effect is involved. `scipy.stats.spearmanr`, two-sided, n = 34 for every cell,
190 unordered pairs.

| | PRIMPOL | BLM | RPA1 | RPA2 | RPA3 | POLA1 | TOP3A | RMI1 | RMI2 | PML | EXD2 | MUS81 | SLX4 | FANCM | WRN | ATRX | DAXX | TERT | SMARCAL1 | ZRANB3 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| **PRIMPOL** | — | +0.27 | +0.40 | +0.44 | +0.32 | +0.04 | +0.26 | +0.39 | +0.35 | -0.32 | -0.11 | +0.05 | +0.19 | +0.26 | +0.24 | +0.34 | +0.24 | +0.41 | -0.15 | +0.32 |
| **BLM** | +0.27 | — | +0.71 | +0.52 | +0.22 | +0.45 | +0.41 | +0.48 | +0.43 | -0.61 | +0.05 | +0.20 | +0.55 | +0.74 | +0.62 | +0.54 | +0.40 | +0.33 | +0.04 | +0.60 |
| **RPA1** | +0.40 | +0.71 | — | +0.52 | +0.29 | +0.37 | +0.40 | +0.60 | +0.45 | -0.54 | +0.00 | +0.21 | +0.46 | +0.58 | +0.50 | +0.51 | +0.43 | +0.30 | +0.13 | +0.59 |
| **RPA2** | +0.44 | +0.52 | +0.52 | — | +0.53 | +0.35 | +0.50 | +0.78 | +0.22 | -0.57 | +0.32 | -0.01 | +0.09 | +0.53 | +0.27 | +0.21 | +0.40 | +0.02 | -0.40 | +0.33 |
| **RPA3** | +0.32 | +0.22 | +0.29 | +0.53 | — | +0.20 | +0.22 | +0.58 | +0.33 | -0.06 | +0.17 | -0.15 | -0.14 | +0.36 | +0.16 | +0.05 | +0.42 | +0.02 | -0.22 | +0.19 |
| **POLA1** | +0.04 | +0.45 | +0.37 | +0.35 | +0.20 | — | +0.39 | +0.38 | +0.38 | -0.14 | +0.05 | +0.42 | +0.37 | +0.42 | +0.39 | +0.27 | +0.39 | +0.26 | +0.21 | +0.13 |
| **TOP3A** | +0.26 | +0.41 | +0.40 | +0.50 | +0.22 | +0.39 | — | +0.59 | +0.14 | -0.19 | -0.07 | +0.08 | +0.30 | +0.31 | +0.36 | +0.20 | +0.19 | -0.10 | -0.23 | +0.27 |
| **RMI1** | +0.39 | +0.48 | +0.60 | +0.78 | +0.58 | +0.38 | +0.59 | — | +0.19 | -0.30 | +0.13 | -0.08 | +0.19 | +0.47 | +0.17 | +0.26 | +0.41 | -0.04 | -0.17 | +0.40 |
| **RMI2** | +0.35 | +0.43 | +0.45 | +0.22 | +0.33 | +0.38 | +0.14 | +0.19 | — | -0.22 | +0.01 | +0.19 | +0.29 | +0.31 | +0.57 | +0.48 | +0.47 | +0.47 | +0.10 | +0.45 |
| **PML** | -0.32 | -0.61 | -0.54 | -0.57 | -0.06 | -0.14 | -0.19 | -0.30 | -0.22 | — | -0.16 | -0.08 | -0.31 | -0.48 | -0.38 | -0.26 | -0.09 | -0.16 | +0.19 | -0.38 |
| **EXD2** | -0.11 | +0.05 | +0.00 | +0.32 | +0.17 | +0.05 | -0.07 | +0.13 | +0.01 | -0.16 | — | -0.14 | -0.22 | +0.07 | -0.16 | -0.27 | +0.06 | -0.32 | -0.22 | -0.07 |
| **MUS81** | +0.05 | +0.20 | +0.21 | -0.01 | -0.15 | +0.42 | +0.08 | -0.08 | +0.19 | -0.08 | -0.14 | — | +0.39 | +0.48 | +0.23 | +0.20 | +0.30 | +0.24 | +0.35 | +0.04 |
| **SLX4** | +0.19 | +0.55 | +0.46 | +0.09 | -0.14 | +0.37 | +0.30 | +0.19 | +0.29 | -0.31 | -0.22 | +0.39 | — | +0.60 | +0.57 | +0.55 | +0.16 | +0.34 | +0.36 | +0.54 |
| **FANCM** | +0.26 | +0.74 | +0.58 | +0.53 | +0.36 | +0.42 | +0.31 | +0.47 | +0.31 | -0.48 | +0.07 | +0.48 | +0.60 | — | +0.54 | +0.48 | +0.41 | +0.25 | +0.04 | +0.40 |
| **WRN** | +0.24 | +0.62 | +0.50 | +0.27 | +0.16 | +0.39 | +0.36 | +0.17 | +0.57 | -0.38 | -0.16 | +0.23 | +0.57 | +0.54 | — | +0.62 | +0.22 | +0.51 | +0.16 | +0.56 |
| **ATRX** | +0.34 | +0.54 | +0.51 | +0.21 | +0.05 | +0.27 | +0.20 | +0.26 | +0.48 | -0.26 | -0.27 | +0.20 | +0.55 | +0.48 | +0.62 | — | +0.33 | +0.41 | +0.27 | +0.42 |
| **DAXX** | +0.24 | +0.40 | +0.43 | +0.40 | +0.42 | +0.39 | +0.19 | +0.41 | +0.47 | -0.09 | +0.06 | +0.30 | +0.16 | +0.41 | +0.22 | +0.33 | — | +0.23 | +0.09 | +0.13 |
| **TERT** | +0.41 | +0.33 | +0.30 | +0.02 | +0.02 | +0.26 | -0.10 | -0.04 | +0.47 | -0.16 | -0.32 | +0.24 | +0.34 | +0.25 | +0.51 | +0.41 | +0.23 | — | +0.19 | +0.36 |
| **SMARCAL1** | -0.15 | +0.04 | +0.13 | -0.40 | -0.22 | +0.21 | -0.23 | -0.17 | +0.10 | +0.19 | -0.22 | +0.35 | +0.36 | +0.04 | +0.16 | +0.27 | +0.09 | +0.19 | — | -0.01 |
| **ZRANB3** | +0.32 | +0.60 | +0.59 | +0.33 | +0.19 | +0.13 | +0.27 | +0.40 | +0.45 | -0.38 | -0.07 | +0.04 | +0.54 | +0.40 | +0.56 | +0.42 | +0.13 | +0.36 | -0.01 | — |

### Pairs above |rho| = 0.7

| Pair | rho (34 lines) | p | rho (16 ALT+ lines) |
|---|---:|---:|---:|
| `RPA2` – `RMI1` | **+0.7818** | 4.8e-08 | +0.468 |
| `BLM` – `FANCM` | **+0.7372** | 6.6e-07 | +0.871 |
| `BLM` – `RPA1` | **+0.7051** | 3.2e-06 | +0.626 |

**3 of 190 pairs exceed |rho| = 0.7**, and they form two clusters:

- **BLM – RPA1 – FANCM** (BLM–FANCM +0.737, BLM–RPA1 +0.705, RPA1–FANCM +0.578). In the
  16 ALT+ lines BLM–FANCM tightens to **+0.871** and a third pair crosses the line:
  **SLX4 – FANCM +0.750** (it is +0.604 over the 34). Treat `BLM`, `FANCM`, `SLX4` and
  `RPA1` as one axis; picking two of them would pre-register the same variable twice.
- **RPA2 – RMI1** (+0.7818 over the 34, but only +0.468 over the 16 — this pair does
  not hold up in the analysis set).

Checked in the ALT+16 as well, 2 pairs exceed 0.7 there: BLM–FANCM
(+0.871) and SLX4–FANCM (+0.750). The union of flagged pairs across both sets is
therefore: **BLM–FANCM, BLM–RPA1, SLX4–FANCM, RPA2–RMI1**.

### The negative control's collinearity with the primary and its alternate

Reported separately because it is the question that decides whether `ZRANB3` can serve
as a predicted-null control at all:

| Pair | rho (34 lines) | p | rho (16 ALT+ lines) | above \|rho\| = 0.7? |
|---|---:|---:|---:|---|
| `ZRANB3` – `BLM` | +0.5997 | 0.000178 | +0.379 | no |
| `ZRANB3` – `PRIMPOL` | +0.3219 | 0.0634 | +0.185 | no |
| `ZRANB3` – `SMARCAL1` | -0.0063 | 0.972 | -0.094 | no |

**Against `BLM` (the nominated primary): rho = +0.5997 over the 34 lines
(p = 1.8e-04), falling to +0.379 over the 16.** That is below the 0.7 flag, so
the two are not the same variable — but it is also `ZRANB3`'s single strongest partner
in the whole panel, and ~36% shared rank variance in the reference set is not nothing.
A prespec pairing a BLM primary with a ZRANB3 null control should state plainly that the
control is positively correlated with the exposure it is controlling for, so a true
BLM effect would be expected to bleed partially into the control. The bleed is smaller
in the analysis set, where the correlation is only +0.379.

**Against `PRIMPOL` (the first alternate): rho = +0.3219, p = 0.063** — not
significant at 0.05, and only +0.185 over the 16. **`ZRANB3` is a materially cleaner
negative control against a `PRIMPOL` primary than against a `BLM` primary.** If the
prespec wants the control to be independent of the exposure, that pairing is the one the
data supports.

For completeness, `ZRANB3` vs `SMARCAL1` expression is rho = -0.0063 (p = 0.97) —
essentially zero, despite the two being functional paralogues. `ZRANB3`'s other
partners, in descending order, are RPA1 +0.594, WRN +0.561, SLX4 +0.536, RMI2 +0.448;
none reaches 0.7.

### Least collinear candidates

Maximum |rho| against any other candidate, over the 34 lines — low is good, because it
means the gene carries information the rest of the panel does not:

| Gene | max \|rho\| | with |
|---|---:|---|
| **EXD2** | 0.318 | `RPA2` |
| **SMARCAL1** | 0.400 | `RPA2` |
| **PRIMPOL** | 0.437 | `RPA2` |
| **POLA1** | 0.455 | `BLM` |
| **DAXX** | 0.473 | `RMI2` |
| **MUS81** | 0.477 | `FANCM` |
| **TERT** | 0.509 | `WRN` |
| **RMI2** | 0.574 | `WRN` |
| **RPA3** | 0.581 | `RMI1` |
| **TOP3A** | 0.588 | `RMI1` |
| **ZRANB3** | 0.600 | `BLM` |
| **SLX4** | 0.604 | `FANCM` |
| **PML** | 0.606 | `BLM` |
| **WRN** | 0.619 | `ATRX` |
| **ATRX** | 0.619 | `WRN` |
| **RPA1** | 0.705 | `BLM` |
| **BLM** | 0.737 | `FANCM` |
| **FANCM** | 0.737 | `BLM` |
| **RPA2** | 0.782 | `RMI1` |
| **RMI1** | 0.782 | `RPA2` |

`EXD2` (0.318), `SMARCAL1` (0.400) and `PRIMPOL` (0.437) are the three most independent
axes in the panel. PML is the only candidate that correlates *negatively* with the
cluster (PML–BLM −0.61, PML–RPA2 −0.57), so it too carries distinct information.

---

## Usable continuous exposures

Ranked on distributional suitability **alone** — realised spread in the 16, absence of
zero-inflation and floor effects, unimodality or clean bimodality, no single point
carrying the fit, and non-collinearity with another candidate. **No hypothesis is
recommended and no mechanism is weighed here.** Literature grounds and mechanism are the
literature scan's and the prespec's job; this ranking says only what each variable can
statistically carry.

| Rank | Gene | SD / range over the 16 | shape | max \|rho\| | why it ranks here |
|---:|---|---|---|---:|---|
| 1 | **TOP3A** | SD 0.989 / 3.30 | skew +0.27, Shapiro p 0.64 | 0.588 | Widest clean spread in the panel. Unimodal in both sets (Shapiro p = 0.73 / 0.64), no zeros, no dominant point, no collinear partner. It also carries the widest interior ≥4/≥4 gap in the 16 (0.631 units at 2.68× the next eligible, a 4/12 split), so a binary form exists if wanted. Caveat: most dosage-driven of the panel (expr–CN rho +0.722). |
| 2 | **PRIMPOL** | SD 0.761 / 2.17 | skew +0.32, Shapiro p 0.12 | 0.437 | Cleanest distribution in the panel: no zeros anywhere, skew +0.07 over 34, Shapiro p = 0.76, flattest tails, and the third-lowest collinearity (0.437). Spread is good though not the widest. Nothing about it is pathological in either set. |
| 3 | **PML** | SD 0.835 / 3.05 | skew -0.12, Shapiro p 0.91 | 0.606 | Wide realised spread, unimodal, no zeros, no dominant point, and the only candidate anti-correlated with the BLM/FANCM cluster, so it is informationally distinct. Weakest expression–CN coupling in the panel (+0.158), i.e. not a dosage proxy. |
| 4 | **MUS81** | SD 0.837 / 3.05 | skew +0.51, Shapiro p 0.38 | 0.477 | Good spread, no zeros, low collinearity (0.477). Mild right skew over the 34 (+0.91, Shapiro p = 0.034) and one high line in the 16 (gap 0.832 at k = 15) — rank-based methods absorb this, a least-squares slope would feel it. |
| 5 | **RMI2** | SD 0.760 / 2.93 | skew +0.59, Shapiro p 0.58 | 0.574 | Good spread in the 16, no zeros there, low collinearity (0.574, with WRN). One low line at 0.96 over the 34 sits outside the ALT+ set, so the 16-line distribution is clean. |
| 6 | **WRN** | SD 0.662 / 2.57 | skew +0.12, Shapiro p 0.65 | 0.619 | Most symmetric distribution in the panel after PRIMPOL (skew +0.17 / +0.12, Shapiro p = 0.62 / 0.65), no zeros, no outliers. Spread is middling and it correlates 0.619 with ATRX and 0.62 with BLM — below the 0.7 flag but not independent. |
| 7 | **SMARCAL1** | SD 0.607 / 2.13 | skew +0.11, Shapiro p 0.84 | 0.400 | Best-behaved variable *inside the analysis set*: range 2.69–4.83, SD 0.607, skew +0.11, excess kurtosis −0.52, Shapiro p = 0.84, no interior gap above 0.30 with 4 lines either side, and the second-lowest collinearity in the panel (0.400). Two caveats: the single near-zero reference line (0.014) means the 34-line distribution is heavily left-skewed even though the 16-line one is not; and expression tracks its own copy number at rho +0.716. |
| 8 | **ATRX** | SD 1.204 / 5.16 | skew -0.56, Shapiro p 0.68 | 0.619 | Largest raw spread of any candidate (SD 1.204, range 5.16 over the 16) — but the spread is manufactured by one line at each end (1.72 and 6.88, gaps of 1.228 and 1.276 at k = 1 and k = 33 over the 34), it is X-linked, and it is the axis iteration 02 closed. Usable as a continuous variable; not clean. |
| 9 | **ZRANB3** | SD 0.427 / 1.58 | skew +0.17, Shapiro p 0.94 | 0.600 | **Nominated negative control; usable as one, with one caveat each way.** Distributionally it is the most featureless variable in the panel — no zeros, no ties (34 of 34 and 16 of 16 distinct values), skew +0.17 in the 16, Shapiro p = 0.95 over the 34, no interior gap above 0.21 with 4 lines either side, and the weakest expression–CN coupling of any non-degenerate candidate except PML and SLX4 (+0.256, p = 0.15), so it is not a dosage proxy. Two caveats: (i) it has the **smallest SD in the panel** apart from the degenerate TERT (0.427 over the 16, 1.58-unit range), which costs nothing for a rank-based test but makes any least-squares slope poorly determined; (ii) it correlates +0.600 with BLM over the 34 lines, its strongest partner anywhere — fine as a control against a PRIMPOL primary, compromised as a control against a BLM primary (§6). |
| 10 | **FANCM** | SD 0.792 / 2.43 | skew -0.54, Shapiro p 0.18 | 0.737 | Flat, well-spread, no zeros, Shapiro p = 0.27 / 0.18. Demoted only for collinearity: +0.737 with BLM over the 34, +0.871 over the 16, +0.750 with SLX4 over the 16. Pick at most one of BLM / FANCM / SLX4 / RPA1. |
| 11 | **RPA2** | SD 0.652 / 2.59 | skew +0.99, Shapiro p 0.33 | 0.782 | The most Gaussian distribution over the 34 lines in the whole panel (Shapiro p = 1.00, skew +0.11). Demoted because the relative spread is small (CV 0.102), one high line (8.06) dominates the 16-line range, and it is collinear with RMI1 (+0.782 over the 34). |
| 12 | **RPA1** | SD 0.651 / 2.01 | skew +0.23, Shapiro p 0.25 | 0.705 | Clean and unimodal (Shapiro p = 0.50 / 0.25, no zeros) with one of the smallest relative spreads in the panel (CV 0.106 over the 16; only DAXX 0.084 and RPA2 0.102 are lower) and collinear with BLM (+0.705). |
| 13 | **RMI1** | SD 0.572 / 1.99 | skew +0.74, Shapiro p 0.25 | 0.782 | Narrow in the analysis set (SD 0.572, range 1.99, IQR 0.538) and collinear with RPA2 (+0.782). One low line at 1.62 over the 34 is outside the ALT+ set. |
| 14 | **RPA3** | SD 0.722 / 3.30 | skew -1.58, Shapiro p 0.02 | 0.581 | Unimodal over the 34 but left-skewed (−0.73); in the 16 its 3.30-unit range is produced almost entirely by one low line (4.19, next 5.54, gap 1.348 at k = 1) — strip that and the IQR is 0.573. |
| 15 | **POLA1** | SD 0.503 / 1.52 | skew -0.35, Shapiro p 0.17 | 0.455 | Narrowest non-degenerate range in the 16 after DAXX (1.52 units, SD 0.503); left-skewed over the 34 (−1.51, Shapiro p = 1.1e-4) with a 4-line low tail that is **not** in the analysis set; and X-linked, with expression differing by sex at p = 0.017. Usable only with sex adjusted for. |
| 16 | **BLM** | SD 0.895 / 3.83 | skew -1.87, Shapiro p 0.01 | 0.737 | The one genuinely bimodal candidate over the reference set (natural low group of 6, gap 0.815 at 0.64 SD) — but only 2 of those lines are ALT+ screened, so inside the 16 it is 2 low outliers plus 14 lines in a 1.49-unit band (excess kurtosis +4.54, Shapiro p = 0.006). Also the hub of the collinearity cluster. |
| 17 | **SLX4** | SD 0.678 / 2.37 | skew +0.92, Shapiro p 0.16 | 0.604 | Lowest-expressed non-degenerate candidate: median 1.47 in the 16, 4 of 16 below 1.0, 7 of 34 below 1.0. Spread exists but sits where log2(TPM+1) is least reliable, and it is collinear with FANCM inside the 16 (+0.750). |

**Collinearity constraint that cuts across the ranking.** `BLM`, `FANCM`, `SLX4` and
`RPA1` are one axis; `RPA2` and `RMI1` are a second (over the 34 only). At most one gene
from each group can be pre-registered without duplicating a variable. Ranks 1–9 contain
no flagged pair, so any single choice from the top nine is unconstrained.

---

## Not usable, and why

| Gene | reason | the numbers |
|---|---|---|
| **TERT** | **Degenerate in the analysis set.** Zero-inflated at the reference level and effectively constant at the ALT+ level. Cannot serve as a continuous exposure under any transformation. | 15/34 and 8/16 at exactly 0.000000; over the 16: max 0.168, IQR 0.022, SD 0.054, CV 1.79; Shapiro p = 2.9e-09 over the 34 |
| **EXD2** | **A 1-vs-15 indicator wearing a continuous variable's clothes.** One line inside the ALT+ set sits 2.57 units below the next; any slope, correlation or regression would be that line's. Would need pre-specified outlier handling to be usable at all. | 16-line values 0.096, then 2.663 … ; excess kurtosis 7.89, Shapiro p = 4.2e-04; gap 2.567 = 2.27 SD at k = 1 |
| **DAXX** | **Too little dynamic range.** Not degenerate, but the smallest range and the smallest CV in the panel; a continuous slope would be estimated across a 1.3-unit window at n = 16. (ZRANB3 has a marginally smaller SD, 0.427, but over a wider range and on a lower mean, so its CV is more than twice DAXX's.) Its one redeeming feature — the cleanest interior split in the 16 (0.468 at 4.09× the next eligible gap, 11/5) — is a split, not a continuum. | range 1.34, SD 0.468, CV 0.084 over the 16 (TOP3A: 3.30 / 0.989 / 0.225) |
| **Any copy-number exposure** | **No validated event clears the minimum cell size.** The only segment-corroborated focal deletion in the entire 20-gene panel has n = 1 (ATRX in one line); 3 of 4 apparent ATRX losses are refuted, reproducing iteration 02. Every other CN axis that would reach 4 lines is a threshold through a near-diploid continuum, which iteration 03 §4.1 already ruled out. | 19 of 20 genes never fall below 0.3 at gene level in any of 33 lines; only 1 spanning segment below 0.1 exists across 20 loci × 33 lines |
| **Second gene from a flagged pair** | Not a property of the gene but of the pair: pre-registering two collinear exposures registers one variable twice. | BLM–FANCM +0.737 (34) / +0.871 (16); BLM–RPA1 +0.705; SLX4–FANCM +0.750 (16); RPA2–RMI1 +0.782 (34) |

### Stated plainly, for the prespec

- **Exactly one gene in this panel has a natural low group of size ≥ 4**, and it is BLM
  (n = 6 over the 34-line reference set, behind a 0.815-unit gap at 0.64 SD). Only 2 of
  those 6 lines are ALT+ CRISPR-screened, so **that low group cannot be used as a
  stratifier in the analysis set.** No gene has a ≥ 4-line natural low group *within*
  the 16.
- **Nothing in this panel reproduces the WRAP53 pattern** that iteration 01 found: no
  gene has a discrete spike at the quantification floor with a gap at ≥ 1.5 SD and ≥ 3×
  the next, except TERT, whose spike holds half the cohort and is therefore useless.
  SMARCAL1 and EXD2 each have a single near-floor line (n = 1, below MIN_CELL_SIZE).
- **Therefore an exposure in this panel is usable as a continuous variable or not at
  all.** Every candidate binary split would be a percentile cut through a dense body —
  the exact failure mode iteration 01 §6 documented. A rank-based (Spearman) primary on
  a continuous exposure avoids both that arbitrariness and the single-point leverage
  that affects EXD2, RPA3, BLM, ATRX, RPA2 and MUS81 in the 16.
- **The negative control can do its job, but it is narrow.** `ZRANB3`'s 34-line spread
  (SD 0.444, range 1.89, CV 0.185) is 35% of BLM's SD and 63% of PRIMPOL's; over the 16
  it is the smallest SD in the panel after TERT. **It is not, however, too small to show
  a relationship in principle**: all 34 and all 16 values are distinct, so a Spearman
  statistic has full rank resolution and exactly the same null distribution it would
  have for any untied exposure — a rank test cannot be weakened by a narrow scale. What
  the narrow window does weaken is any *variance-scaled* estimate: a Pearson r, an OLS
  slope or a standardised effect size on ZRANB3 would carry a wider confidence interval
  than the same quantity on BLM or PRIMPOL. **Recommendation for the prespec, on
  distributional grounds only: register the control on the same rank-based statistic as
  the primary, and do not interpret a null ZRANB3 slope as evidence of specificity if
  the primary is analysed by OLS.** A second, independent caveat is in §6: ZRANB3 is
  correlated +0.600 with BLM, so against a BLM primary it is a partially dependent
  control.
- **Copy number contributes nothing as an exposure here**, and for the dosage-coupled
  genes (TOP3A +0.722, SMARCAL1 +0.716, RPA3 +0.659, RMI1 +0.652) it is not even an
  independent corroboration channel for expression.

---

## Limitations of this scan

1. **n = 34, and n = 16 where it matters.** Every SD, quartile and rho rests on 34
   observations, and the realised spreads in §4 on 16. No confidence intervals are given
   for any statistic and none should be inferred. A single line added or removed moves
   quartiles and can move a gap across a rank boundary.
2. **No formal modality test.** The unimodal/bimodal verdicts in §2 are descriptive
   readings of sorted values, gaps, skew, kurtosis and Shapiro–Wilk. No Hartigan dip
   test, KDE or mixture fit was run, and at n = 34 any of them would have low power. The
   BLM bimodality claim in particular rests on 6 lines.
3. **One modality, one level, one release.** Gene-level log2(TPM+1) and gene-level WGS CN
   from DepMap 26Q1 only. No isoform, protein, methylation, mutation or fusion data, and
   no second cohort. A line with normal mRNA and a damaging allele is invisible here —
   which matters most for SMARCAL1 itself, where the literature motivation is *loss of
   function*, not low transcript.
4. **The 20 genes are not a complete mechanism panel.** They are the 19 the brief
   specified plus `ZRANB3`, added mid-scan by the literature check. Genes absent from it were not screened, so "least collinear in this panel"
   does not mean least collinear in the transcriptome.
5. **Gene intervals in §5.2 were hand-entered from memory**, not read from an annotation
   file, because the repository contains none (iteration 03 §4.2). They are validated
   only by the internal consistency check in §5.2 (rho ≥ +0.85 against gene-level CN) —
   strong evidence the intervals are right to within a segment, but not a substitute for
   a GTF. A systematically shifted interval that still fell inside the same broad segment
   would pass that check. Before any CN-based claim enters an analysis, the coordinates
   should be re-read from an authoritative annotation and hard-coded in the iteration's
   own notebook, as iteration 03 recommended.
6. **Focality was assessed only at the 19 candidate loci.** No genome-wide segment
   inventory was recomputed here; the cohort-wide focal-zero prior is taken from
   iteration 03 §4.3 rather than re-derived.
7. **No confounder adjustment.** Values are raw log2(TPM+1). Sex was checked only for the
   two X-linked candidates. Nothing was adjusted for growth medium, passage, purity,
   depositor, or the OS-subtype heterogeneity in `Model.csv`, and the 171-pair Spearman
   matrix is unadjusted for multiplicity (the 0.7 flag is a descriptive screen, not a
   test).
8. **Blindness constrains what can be checked.** By construction this scan cannot say
   whether any exposure is associated with SMARCAL1 dependency, whether the 16-line
   spread is adequate *power* for any particular effect size, or whether a line that is
   an exposure outlier is also an outcome outlier. Those are outcome-side questions and
   are out of scope.
9. **Ranking is distributional, not scientific.** The order in "Usable continuous
   exposures" weighs spread, shape and independence. It says nothing about which
   exposure the mechanism literature supports, and it must not be read as a
   recommendation. Two variables can be equally well distributed and differ completely
   in prior plausibility.

---

## Blindness attestation

In producing this file I did **not**:

- import `src.dependency_stats` — **not imported**
- read, stream, parse or grep any gene-effect value from
  `data/raw/depmap/CRISPRGeneEffect.csv`. The file was opened once, by
  `src.cohort.crispr_model_ids()`, which reads the header line and then the first
  comma-delimited field of each row — the ModelID — and discards the rest. That is
  **membership, not outcome**.
- compute, load or display any quantity derived from Chronos gene effect, for SMARCAL1
  or any other gene
- open `autoresearch/iteration_03_bir_module_coherence/data/` or any
  `h4_modules.csv` / `per_gene_results.csv` / `primary_result.json` /
  `stage1_verdict.json` / `composite_scores.csv` / `robustness_table.csv`
- print to a terminal any table pairing ModelIDs with ALT labels. Every ALT+ statistic
  above was computed and written to this file; stdout during the scan carried only
  gene-level aggregates and sorted value vectors with identifiers withheld.

`references/master_manifest.csv` **was** read, via `src.cohort.build_cohort()`, solely
to construct the 16-line ALT+ membership that item 4 of the brief requires. ALT labels
are a frozen input, not an outcome, and no ALT label appears in any table in this file.

No existing repository file was modified. This document is the only file written.
