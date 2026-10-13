# Notebook 13 — Aim 2 D3 Curated ALT: Robustness Analysis

Run timestamp: `2026-06-12T18:02:13Z`

Four extensions to the notebook 12 D3 (curated ALT-only) analysis. See notebook 12 for the baseline D3 results.

## Configuration

- Curated gene list: 19/21 found  (not found: ['TERC', 'MRE11A'])
- 1000-perm test: fixed median C, precomputed fold scalers (same results as re-running full pipeline per shuffle, faster).
- Bootstrap stability: 200 resamples, NONZERO_THRESH=1e-10.
- Random seed: 42.  CI bootstrap resamples: 1000.

## Analysis 1: D3 with 1000-permutation test

| Cohort | N | AUROC | 95% CI | AUPRC | Sens@90Spec | Perm p (1000) |
|---|---|---|---|---|---|---|
| primary | 20 | 0.9479 | [0.8125,1.0000] | 0.9755 | 0.9167 | 0.0000 |
| pragmatic | 31 | 0.8291 | [0.6681,0.9625] | 0.9166 | 0.7778 | 0.0100 |
| pragmatic_ks (N=27) | 27 | 0.8765 | [0.7176,0.9835] | 0.9369 | 0.8235 | 0.0030 |

## Analysis 2: Known-sex-only pragmatic sensitivity

Dropped samples: ['ACH-001813', 'ACH-001817', 'ACH-002779', 'ACH-003573']

N=31 AUROC=0.8291  N=27 AUROC=0.8765  delta=+0.0474

### Coefficient comparison (full vs known-sex, sorted by |coef| N=31)

| Rank_N31 | Gene | Coef N=31 | Coef N=27 | Delta |
|---|---|---|---|---|
| 1 | FANCM | +0.767177 | +0.632474 | -0.134702 |
| 2 | TERT | -0.754628 | -0.882270 | -0.127642 |
| 3 | TP53 | -0.336636 | -0.100740 | +0.235896 |
| 4 | NBN | +0.237621 | +0.269259 | +0.031638 |
| 5 | WRN | -0.219992 | -0.300306 | -0.080314 |
| 6 | ATRX | -0.124507 | -0.169728 | -0.045221 |
| 7 | BRCA2 | +0.116268 | +0.000000 | -0.116268 |
| 8 | SMARCAL1 | -0.088260 | +0.000000 | +0.088260 |
| 9 | RAD52 | +0.000000 | +0.000000 | +0.000000 |
| 10 | GEN1 | +0.000000 | +0.000000 | +0.000000 |
| 11 | SLX4 | +0.000000 | +0.000000 | +0.000000 |
| 12 | BRCA1 | +0.000000 | +0.000000 | +0.000000 |
| 13 | FANCD2 | +0.000000 | +0.000000 | +0.000000 |
| 14 | WRAP53 | +0.000000 | +0.000000 | +0.000000 |
| 15 | RAD51AP1 | +0.000000 | +0.000000 | +0.000000 |
| 16 | DAXX | +0.000000 | +0.000000 | +0.000000 |
| 17 | RAD50 | +0.000000 | +0.000000 | +0.000000 |
| 18 | BLM | +0.000000 | +0.000000 | +0.000000 |
| 19 | MUS81 | +0.000000 | +0.000000 | +0.000000 |

## Analysis 3: Bootstrap coefficient stability

| Cohort | N | stable (>80%) | borderline (50-80%) | unstable (<50%) |
|---|---|---|---|---|
| primary | 20 | 3 genes | 2 genes | 14 genes |
| pragmatic | 31 | 3 genes | 5 genes | 11 genes |
| pragmatic_ks | 27 | 3 genes | 3 genes | 13 genes |

Genes stable in both primary and pragmatic (2): ['FANCM', 'TERT']

### Per-gene bootstrap non-zero frequency

| Gene | Primary boot_freq | Primary stability | Pragmatic boot_freq | Pragmatic stability |
|---|---|---|---|---|
| ATRX | 0.235 | unstable | 0.555 | borderline |
| DAXX | 0.165 | unstable | 0.21 | unstable |
| SMARCAL1 | 0.16 | unstable | 0.575 | borderline |
| TP53 | 0.445 | unstable | 0.935 | stable |
| TERT | 1.0 | stable | 1.0 | stable |
| WRN | 0.57 | borderline | 0.705 | borderline |
| BLM | 0.12 | unstable | 0.365 | unstable |
| RAD50 | 0.08 | unstable | 0.155 | unstable |
| NBN | 0.625 | borderline | 0.785 | borderline |
| FANCD2 | 0.03 | unstable | 0.215 | unstable |
| RAD51AP1 | 0.28 | unstable | 0.335 | unstable |
| WRAP53 | 0.17 | unstable | 0.12 | unstable |
| RAD52 | 0.08 | unstable | 0.12 | unstable |
| BRCA1 | 0.21 | unstable | 0.1 | unstable |
| BRCA2 | 0.81 | stable | 0.655 | borderline |
| FANCM | 1.0 | stable | 1.0 | stable |
| SLX4 | 0.11 | unstable | 0.21 | unstable |
| GEN1 | 0.075 | unstable | 0.08 | unstable |
| MUS81 | 0.11 | unstable | 0.16 | unstable |

## Analysis 4: LOOCV fold coefficient frequency

| Gene | Primary fold_freq | Pragmatic fold_freq | Pragmatic_ks fold_freq |
|---|---|---|---|
| ATRX | 0.1000 | 0.6452 | 0.8889 |
| DAXX | 0.0000 | 0.2581 | 0.1481 |
| SMARCAL1 | 0.0000 | 0.6774 | 0.3704 |
| TP53 | 0.5000 | 0.9677 | 0.8889 |
| TERT | 1.0000 | 1.0000 | 1.0000 |
| WRN | 0.9000 | 0.6774 | 0.9259 |
| BLM | 0.0000 | 0.2903 | 0.1852 |
| RAD50 | 0.0000 | 0.1290 | 0.4074 |
| NBN | 0.9000 | 0.7097 | 0.9259 |
| FANCD2 | 0.0000 | 0.2258 | 0.1481 |
| RAD51AP1 | 0.0000 | 0.0323 | 0.0370 |
| WRAP53 | 0.0000 | 0.0968 | 0.1111 |
| RAD52 | 0.0000 | 0.1290 | 0.1852 |
| BRCA1 | 0.0000 | 0.1935 | 0.4074 |
| BRCA2 | 0.9500 | 0.6774 | 0.5185 |
| FANCM | 1.0000 | 1.0000 | 0.9259 |
| SLX4 | 0.0000 | 0.0323 | 0.0741 |
| GEN1 | 0.0000 | 0.1613 | 0.0741 |
| MUS81 | 0.0000 | 0.0968 | 0.0000 |

