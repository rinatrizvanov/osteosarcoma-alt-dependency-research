# Notebook 04 — Stratification by TCAB1/WRAP53

**Pre-specification freeze commit:** `81745ed` (v1.0), amended at `a63d56c` (v1.1, before any gene-effect value was read).

Stratifier: direct WRAP53 log2(TPM+1) (notebook 02 branch `direct`, prespec §4c).
Reads no gene-effect values.

## §4b minimum-cell guard: FIRED

- §4a tertile cut (0.3333 quantile over 34 osteosarcoma lines with RNA) = **3.335580**
- ALT+ lines below the cut: **3 of 16**
- Guard threshold: fewer than 4 per arm
- Outcome: tertile primary declared NOT ESTIMABLE; primary falls back to the median split

> TCAB1-low arm n=3 (intact n=13); minimum cell size is 4. Tertile primary declared NOT ESTIMABLE; falling back to the median split within the ALT+ set, per the pre-specified Section 4 guard.

This matches prespec §4d, which recorded the guard outcome before any gene-effect
value was read. The guard was committed in v1.0 before the WRAP53 distribution was
examined, so the fallback is pre-specified, not post-hoc.

## Contrasts carried into notebook 05

| Contrast | Rule | Split | Role |
|---|---|---|---|
| `primary_median_8v8` | median within ALT+ | 8 v 8 | **primary** |
| `sensitivity_tertile_3v13` | §4a OS-reference tertile (superseded) | 3 v 13 | pre-specified sensitivity |
| `alt_negative_median` | median within ALT- | 4 v 5 | §6.2 negative control |

Both ALT+ contrasts are run in parallel in notebook 05 with the same test and the
same output format, and both enter the final report (prespec §4d).

## §4d assertion: WRAP53-null lines

Both lines are required to be in the TCAB1-low arm of both contrasts. **PASSED.**

| ModelID | Cell line | ALT tier | WRAP53 | TP53 | WRAP53 CN | Arm (8v8) | Arm (3v13) |
|---|---|---|---|---|---|---|---|
| ACH-002834 | OS384 | alt_gold_positive | 0.000 | 0.000 | 0.00009 | TCAB1-low | TCAB1-low |
| ACH-003185 | NRHOS4 | alt_silver | 0.000 | 0.106 | 0.00001 | TCAB1-low | TCAB1-low |

Recorded as an exposure-side measurement fact with copy-number corroboration.
**No mechanistic claim is drawn from it here**; interpretation waits for notebooks
05 and 07 (prespec §4d).

## Copy-number concordance

Spearman(WRAP53 expression, WRAP53 CN) within the ALT+ lines: **+0.7638** (n=16). Independent corroboration of the `direct` branch: antisense
cross-assignment would degrade this concordance.

## Stratifier table (ALT+ lines, ascending WRAP53)

| ModelID | Cell line | WRAP53 | TP53 | WRAP53 CN | Arm (8v8) | Arm (3v13) | ALT tier |
|---|---|---|---|---|---|---|---|
| ACH-002834 | OS384 | 0.000 | 0.000 | 0.0001 | TCAB1-low | TCAB1-low | alt_gold_positive |
| ACH-003185 | NRHOS4 | 0.000 | 0.106 | 0.0000 | TCAB1-low | TCAB1-low | alt_silver |
| ACH-003171 | NRHOS1 | 2.748 | 2.531 | 0.7034 | TCAB1-low | TCAB1-low | alt_silver |
| ACH-002433 | CCLFPEDS0019T | 3.336 | 6.541 | 1.3115 | TCAB1-low | TCAB1-intact | alt_silver |
| ACH-003178 | NRHOS2 | 3.427 | 1.057 | 0.6612 | TCAB1-low | TCAB1-intact | alt_silver |
| ACH-001715 | CAL72 | 3.459 | 0.290 | 0.6867 | TCAB1-low | TCAB1-intact | alt_gold_positive |
| ACH-003573 | C337 | 3.484 | 0.458 | 0.7845 | TCAB1-low | TCAB1-intact | alt_silver |
| ACH-002688 | OS457 | 3.489 | 4.089 | 0.7111 | TCAB1-low | TCAB1-intact | alt_gold_positive |
| ACH-000364 | U2OS | 4.100 | 4.994 | 0.7006 | TCAB1-intact | TCAB1-intact | alt_gold_positive |
| ACH-000410 | SAOS2 | 4.286 | 5.008 | 0.8429 | TCAB1-intact | TCAB1-intact | alt_gold_positive |
| ACH-001526 | HUO9 | 4.301 | 4.533 | 1.0901 | TCAB1-intact | TCAB1-intact | alt_gold_positive |
| ACH-003587 | 531MII | 4.415 | 0.537 | 1.3933 | TCAB1-intact | TCAB1-intact | alt_silver |
| ACH-002067 | NOS1 | 4.607 | 5.930 | 1.4940 | TCAB1-intact | TCAB1-intact | alt_gold_positive |
| ACH-002690 | OS526 | 4.619 | 1.228 | 1.0269 | TCAB1-intact | TCAB1-intact | alt_gold_positive |
| ACH-002471 | OS052 | 4.628 | 5.427 | 0.9067 | TCAB1-intact | TCAB1-intact | alt_gold_positive |
| ACH-000082 | G292CLONEA141B1 | 5.888 | 2.475 | 2.3397 | TCAB1-intact | TCAB1-intact | alt_gold_positive |
