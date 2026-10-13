# Notebook 05 — Primary essentiality comparison

First notebook to read the outcome. Pre-specification frozen at `81745ed`
(v1.0) and amended at `a63d56c` (v1.1), both before this notebook ran.

## Endpoint

```
composite = mean(POLD3, RAD52, SMARCAL1) - mean(BRCA1, BRCA2, RAD51)
```

Chronos gene effect is more negative when a gene is more essential, so the
hypothesis predicts a **more negative** composite in TCAB1-low lines, i.e.
Cliff's delta < 0. Tests are two-sided exact Mann-Whitney.

## Primary and parallel sensitivity contrast

### primary_median_8v8 (PRIMARY)

- n: **8** TCAB1-low vs **8** TCAB1-intact
- median composite: low **+0.2528**, intact **-0.0093**
- Mann-Whitney U = 41.0, two-sided **p = 0.3823** (exact)
- Cliff's delta = **+0.2812**, 95% CI [-0.3125, +0.8438] (small)
- Hodges-Lehmann shift = **+0.1865** Chronos units
- direction as predicted (delta < 0): **False**

### sensitivity_tertile_3v13 (pre-specified sensitivity)

- n: **3** TCAB1-low vs **13** TCAB1-intact
- median composite: low **+0.4605**, intact **+0.0384**
- Mann-Whitney U = 30.0, two-sided **p = 0.1893** (exact)
- Cliff's delta = **+0.5385**, 95% CI [-0.2308, +1.0000] (large)
- Hodges-Lehmann shift = **+0.2704** Chronos units
- direction as predicted (delta < 0): **False**

## Minimum detectable effect (prespec §5a)

| Contrast | Observed delta | min \|delta\| for p<0.05 | Significance reachable |
|---|---|---|---|
| primary_median_8v8 | +0.2812 | 0.594 | no |
| sensitivity_tertile_3v13 | +0.5385 | 0.795 | no |

## R1 gold-only sensitivity

| Contrast | n low | n intact | p | delta | Direction as predicted |
|---|---|---|---|---|---|
| primary_median_8v8 | 3 | 7 | 0.1833 | -0.6190 | True |
| sensitivity_tertile_3v13 | — | — | not estimable | — | — |

## Six per-gene comparisons, BH-FDR across the six (§6.1)

### primary_median_8v8

| Gene | Set | median low | median intact | p | q (BH) | delta | 95% CI | HL | Magnitude |
|---|---|---|---|---|---|---|---|---|---|
| POLD3 | BIR | -1.6079 | -1.7463 | 0.9591 | 0.9591 | -0.0312 | [-0.625, +0.562] | -0.0401 | negligible |
| RAD52 | BIR | -0.0133 | +0.0468 | 0.2786 | 0.9591 | -0.3438 | [-0.875, +0.250] | -0.0660 | medium |
| SMARCAL1 | BIR | -0.4020 | -0.4628 | 0.7209 | 0.9591 | +0.1250 | [-0.500, +0.719] | +0.0794 | negligible |
| BRCA1 | HR | -0.5548 | -0.5252 | 0.9591 | 0.9591 | -0.0312 | [-0.594, +0.594] | -0.0185 | negligible |
| BRCA2 | HR | -0.5485 | -0.3582 | 0.4418 | 0.9591 | -0.2500 | [-0.781, +0.344] | -0.1388 | small |
| RAD51 | HR | -1.3193 | -1.0830 | 0.5737 | 0.9591 | -0.1875 | [-0.750, +0.469] | -0.2177 | small |

### sensitivity_tertile_3v13

| Gene | Set | median low | median intact | p | q (BH) | delta | 95% CI | HL | Magnitude |
|---|---|---|---|---|---|---|---|---|---|
| POLD3 | BIR | -1.5598 | -1.7136 | 0.3643 | 1 | +0.3846 | [-0.081, +0.846] | +0.1710 | medium |
| RAD52 | BIR | -0.0142 | -0.0091 | 0.7036 | 1 | -0.1795 | [-0.897, +0.692] | -0.0460 | small |
| SMARCAL1 | BIR | -0.0390 | -0.5547 | 0.2964 | 1 | +0.4359 | [-0.692, +1.000] | +0.3618 | medium |
| BRCA1 | HR | -0.5461 | -0.5635 | 0.7036 | 1 | -0.1795 | [-1.000, +0.590] | -0.1813 | small |
| BRCA2 | HR | -0.3931 | -0.4198 | 1 | 1 | -0.0256 | [-0.641, +0.641] | -0.0178 | negligible |
| RAD51 | HR | -1.1220 | -1.1328 | 0.9 | 1 | +0.0769 | [-0.692, +0.846] | +0.0148 | negligible |

## §6.4 concentration comparison (pre-registered, descriptive, non-gating)

- delta, 8v8 median primary: **+0.2812**
- delta, 3v13 tertile sensitivity: **+0.5385**
- \|delta\| excess of tertile over median: **+0.2572** (threshold 0.3)
- concentrates in extreme-loss lines: **False**

No evidence of concentration in the extreme-loss lines: the tertile contrast exceeds the median primary by +0.257 in |delta|, within the pre-registered 0.3 threshold.

## Figures

- `figures/strip_primary_median_8v8.png`
- `figures/forest_primary_median_8v8.png`
- `figures/strip_sensitivity_tertile_3v13.png`
- `figures/forest_sensitivity_tertile_3v13.png`

## Not yet run

Robustness (R2-R8), the §6.2 ALT- selectivity control, the §6.3 within-line
standardised sensitivity, and the §8/§9/§10 verdict are notebooks 06 and 07.
No verdict is drawn here.
