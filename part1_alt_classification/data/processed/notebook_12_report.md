# Notebook 12 — Aim 2 ALT Signature: Sex-Adjusted Sensitivity Analysis

Run timestamp: `2026-06-05T17:10:00Z`

Sensitivity analysis for sex confounding in the notebook 11 primary top-20 (USP9Y, DDX3Y). Three designs x two definitions = 6 analyses.

## Configuration

- Gene variance filter: top `500` genes (D1, D2); bypassed for D3.
- Model: `LogisticRegression(penalty="elasticnet", solver="saga", l1_ratio=0.5, class_weight="balanced", random_state=42)`.
- C grid: `10**linspace(-2,2,9)`.  Inner CV folds: 5.
- Bootstrap: 1000 resamples.  Permutations: 100 (fixed median C).

## Results

| Definition | Design | AUROC | 95% CI | AUPRC | Sens@90Spec | Perm p |
|---|---|---|---|---|---|---|
| primary | nb11-baseline | 0.4688 | [0.2031, 0.7500] | 0.6422 | 0.0833 | N/A |
| primary | D1-sex-cov | 0.5000 | [0.1904, 0.8001] | 0.5818 | 0.0000 | 0.530 |
| primary | D2-no-Y-no-XCT | 0.3958 | [0.1098, 0.7144] | 0.5386 | 0.0000 | 0.660 |
| primary | D3-curated-ALT | 0.9479 | [0.8125, 1.0000] | 0.9755 | 0.9167 | 0.000 |
| pragmatic | nb11-baseline | 0.6880 | [0.4473, 0.8950] | 0.7024 | 0.1111 | N/A |
| pragmatic | D1-sex-cov (SKIPPED) | SKIPPED |  | SKIPPED | SKIPPED | SKIPPED |
| pragmatic | D2-no-Y-no-XCT | 0.7051 | [0.4772, 0.8955] | 0.7443 | 0.2222 | 0.060 |
| pragmatic | D3-curated-ALT | 0.8291 | [0.6681, 0.9625] | 0.9166 | 0.7778 | 0.000 |

## Sex distribution per arm

### Primary (N=20)

```
alt_call_direction  negative  positive
Sex                                   
Female                     4         7
Male                       4         5
```

### Pragmatic (N=31)

```
alt_call_direction  negative  positive
Sex                                   
Female                     5         8
Male                       5         9
Unknown                    3         1
```

## D1 — Sex covariate notes

**Primary**: Analysis completed. All 20 samples have Female or Male annotation.

**Pragmatic**: SKIPPED — 4 sample(s) with Unknown/missing sex: ['ACH-001813', 'ACH-001817', 'ACH-002779', 'ACH-003573']

## D2 — Excluded genes

Total excluded: 45 columns (41 Y-chr + 4 X CT antigens)

```
AMELY
BPY2
BPY2B
BPY2C
CDY1
CDY1B
CDY2A
CDY2B
DAZ1
DAZ2
DAZ3
DAZ4
DDX3Y
EIF1AY
HSFY1
HSFY2
KDM5D
MAGEA10
MAGEB2
NLGN4Y
PCDH11Y
PRY
PRY2
RBMY1A1
RBMY1B
RBMY1D
RBMY1E
RBMY1F
RBMY1J
RPS4Y1
RPS4Y2
SRY
TBL1Y
TMSB4Y
TSPY1
TSPY2
TSPY3
TSPY4
USP9Y
UTY
VCY
VCY1B
XAGE1A
XAGE1B
ZFY
```

## D3 — Curated ALT gene availability

Found: 19/21 genes

Not found in matrix: ['TERC', 'MRE11A']

## D3 — Curated ALT gene coefficients

### Primary (C=0.3162)

| Rank | Gene | Coefficient |
|---|---|---|
| 1 | TERT | -0.689151 |
| 2 | FANCM | +0.683066 |
| 3 | BRCA2 | +0.196771 |
| 4 | WRN | -0.152990 |
| 5 | NBN | +0.105463 |
| 6 | TP53 | -0.023244 |
| 7 | ATRX | +0.000000 |
| 8 | RAD52 | +0.000000 |
| 9 | GEN1 | +0.000000 |
| 10 | SLX4 | +0.000000 |
| 11 | BRCA1 | +0.000000 |
| 12 | FANCD2 | +0.000000 |
| 13 | WRAP53 | +0.000000 |
| 14 | RAD51AP1 | +0.000000 |
| 15 | DAXX | +0.000000 |
| 16 | RAD50 | +0.000000 |
| 17 | BLM | +0.000000 |
| 18 | SMARCAL1 | +0.000000 |
| 19 | MUS81 | +0.000000 |

### Pragmatic (C=0.3162)

| Rank | Gene | Coefficient |
|---|---|---|
| 1 | FANCM | +0.767177 |
| 2 | TERT | -0.754628 |
| 3 | TP53 | -0.336636 |
| 4 | NBN | +0.237621 |
| 5 | WRN | -0.219992 |
| 6 | ATRX | -0.124507 |
| 7 | BRCA2 | +0.116268 |
| 8 | SMARCAL1 | -0.088260 |
| 9 | RAD52 | +0.000000 |
| 10 | GEN1 | +0.000000 |
| 11 | SLX4 | +0.000000 |
| 12 | BRCA1 | +0.000000 |
| 13 | FANCD2 | +0.000000 |
| 14 | WRAP53 | +0.000000 |
| 15 | RAD51AP1 | +0.000000 |
| 16 | DAXX | +0.000000 |
| 17 | RAD50 | +0.000000 |
| 18 | BLM | +0.000000 |
| 19 | MUS81 | +0.000000 |

