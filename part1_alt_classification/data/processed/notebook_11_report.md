# Notebook 11 — Aim 2 ALT Signature Score (DepMap OS RNA-seq)

Regularized ALT signature score on DepMap_26Q1 osteosarcoma RNA-seq using direction-aware labels (`alt_call_direction`), evaluated by leave-one-out cross-validation per `label_spec_v1.0` §2. Signature score, **not** a calibrated classifier and **not** a fixed-threshold predictor.

## Configuration
- Model: `LogisticRegression(penalty="elasticnet", solver="saga", l1_ratio=0.5, max_iter=10000, class_weight="balanced", random_state=42)` in `Pipeline(VarianceGeneFilter -> StandardScaler -> LogisticRegression)`.
- Fold-local gene-variance filter keeps top `500` high-variance genes (training-fold variances only); see notebook §6 / `src.aim2_signature` for rationale.
- Inner 5-fold CV selects `C` from `10**linspace(-2,2,9)` = [0.01, 0.0316, 0.1, 0.3162, 1.0, 3.1623, 10.0, 31.6228, 100.0].
- Bootstrap resamples: 1000; permutation shuffles: 100; seed: 42.

## Filter trace
### Primary (admitted alt_status: ('alt_gold_positive', 'alt_negative_high_conf'))
```
manifest rows (all cohorts)                             N = 115
cohort == 'DepMap_26Q1'                                 N = 35
sample_id present in DepMap RNA ModelID index           N = 34
alt_call_direction in ('positive', 'negative')          N = 31
alt_status in ('alt_gold_positive', 'alt_negative_high_conf') N = 20
direction counts: {'positive': 12, 'negative': 8}
```
### Pragmatic (admitted alt_status: ('alt_gold_positive', 'alt_negative_high_conf', 'alt_silver'))
```
manifest rows (all cohorts)                             N = 115
cohort == 'DepMap_26Q1'                                 N = 35
sample_id present in DepMap RNA ModelID index           N = 34
alt_call_direction in ('positive', 'negative')          N = 31
alt_status in ('alt_gold_positive', 'alt_negative_high_conf', 'alt_silver') N = 31
direction counts: {'positive': 18, 'negative': 13}
```
Count assertions (expected primary 12/8, pragmatic 18/13): **PASS**.
Ground-truth: ACH-000364 (U2OS) = positive, ACH-000748 (SJSA1) = negative: **PASS**. Arms derived from `alt_call_direction`, never `alt_status`.

## Cross-validated performance
| Definition | N (pos/neg) | AUROC [95% CI] | AUPRC [95% CI] | Sens@90%Spec [95% CI] | Permutation p (AUROC) |
|---|---|---|---|---|---|
| primary | 20 (12/8) | 0.4688 [0.2031, 0.7500] | 0.6422 [0.4300, 0.8879] | 0.0833 [0.0000, 0.5000] | 0.6000 |
| pragmatic | 31 (18/13) | 0.6880 [0.4473, 0.8950] | 0.7024 [0.4856, 0.9225] | 0.1111 [0.0000, 0.5789] | 0.0600 |

- primary permutation null: mean AUROC=0.5114, 95th pct=0.7917, fixed C=100 (median observed), n=100.
- pragmatic permutation null: mean AUROC=0.4537, 95th pct=0.7105, fixed C=31.62 (median observed), n=100.

## Top-20 genes by |coefficient| (signed; refit on full N at median C)
### Primary
| rank | gene | coefficient | in ALT gene list |
|---|---|---:|---|
| 1 | GPX3 | +0.2375 |  |
| 2 | RPL39L | -0.2230 |  |
| 3 | CRYAB | +0.2133 |  |
| 4 | XAGE1B | +0.2002 |  |
| 5 | MAGEB2 | +0.1839 |  |
| 6 | AMIGO2 | +0.1828 |  |
| 7 | COLEC12 | -0.1818 |  |
| 8 | LY6E | -0.1772 |  |
| 9 | XAGE1A | +0.1767 |  |
| 10 | IFITM2 | -0.1748 |  |
| 11 | MAGEA10 | +0.1747 |  |
| 12 | KRBOX1 | +0.1731 |  |
| 13 | BMP4 | +0.1709 |  |
| 14 | TDO2 | +0.1701 |  |
| 15 | PLIN2 | -0.1657 |  |
| 16 | GAL | -0.1652 |  |
| 17 | USP9Y | -0.1590 |  |
| 18 | DDX3Y | -0.1586 |  |
| 19 | ACAN | +0.1585 |  |
| 20 | KCNN4 | +0.1534 |  |

ALT-gene-list hits in top-20: none.

### Pragmatic
| rank | gene | coefficient | in ALT gene list |
|---|---|---:|---|
| 1 | GAL | -0.3529 |  |
| 2 | GPX3 | +0.3286 |  |
| 3 | ALDH1A3 | -0.3069 |  |
| 4 | CRYAB | +0.3011 |  |
| 5 | H1-2 | +0.2838 |  |
| 6 | CPM | +0.2760 |  |
| 7 | CDC42EP5 | -0.2687 |  |
| 8 | AMIGO2 | +0.2655 |  |
| 9 | RPL39L | -0.2581 |  |
| 10 | FABP3 | +0.2478 |  |
| 11 | MT1F | +0.2468 |  |
| 12 | LRRN1 | +0.2403 |  |
| 13 | LRRC15 | +0.2401 |  |
| 14 | HLA-DRA | +0.2390 |  |
| 15 | NINJ2 | +0.2380 |  |
| 16 | INA | +0.2363 |  |
| 17 | TSPYL5 | +0.2355 |  |
| 18 | SFRP2 | -0.2238 |  |
| 19 | PHLDA2 | -0.2234 |  |
| 20 | HSPB7 | +0.2175 |  |

ALT-gene-list hits in top-20: none.

## ALT gene-list cross-check (full ranked list)
Curated ALT-related genes: ATRX, DAXX, SMARCAL1, TP53, TERT, TERC, WRN, BLM, MRE11A, RAD50, NBN, FANCD2, RAD51AP1, WRAP53, RAD52, BRCA1, BRCA2, FANCM, SLX4, GEN1, MUS81.

- primary: 1 ALT-list gene(s) with non-zero coefficient: TP53(rank 213, -0.0521).
- pragmatic: 1 ALT-list gene(s) with non-zero coefficient: TP53(rank 59, -0.1417).

## Primary vs pragmatic comparison
- AUROC delta (primary - pragmatic): -0.2193
- Top-20 gene overlap: 5 (Jaccard 0.143); genes: AMIGO2, CRYAB, GAL, GPX3, RPL39L
- Coefficient Pearson r over 416 shared genes: 0.7661

## Methods notes & limitations
1. Signature score by LOOCV (label_spec_v1.0 §2): no calibrated probabilities, no fixed thresholds.
2. All preprocessing (variance filter, scaler) is fit inside the CV training fold only; nested inner CV uses its own per-fold scaler.
3. The gene-variance filter is implemented as a fold-local top-K high-variance selector (K=500) for tractability over ~19k genes at N=20-31.
4. The permutation null fixes C to the median observed value (re-selecting C in every fold of every shuffle is intractable); a documented, slightly conservative concession.
5. DepMap_26Q1 only; not applied to other cohorts (notebook 12) and never pooled with de Nonneville microarray.
6. Very small N (20 strict / 31 pragmatic): treat all estimates as exploratory; the permutation p-value is the primary honesty check.
