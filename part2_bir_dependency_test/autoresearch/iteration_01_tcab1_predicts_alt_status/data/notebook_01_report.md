# Iteration 01 — H1: WRAP53 expression vs ALT status

Pre-specification: `prespec.md` v1.0.1, frozen `8e6ed51`, erratum `2d68239`.
Outcome variable: the frozen ALT label. **No gene-effect value was read** — `load_gene_effects` was replaced with a raiser and never fired (0 attempts).

## Cohort (prespec §4)

| Step | N |
|---|---|
| master_manifest rows | 115 |
| cohort == DepMap_26Q1 | 35 |
| alt_call_direction in {positive, negative} | 31 |
| present in the expression matrix | 31 |
| arms: 18 ALT+ / 13 ALT- | 31 |

Primary is patient-collapsed: **18 ALT+ / 12 ALT−**, N=30. `ACH-000613` dropped as a same-patient derivative of `ACH-001001`, by the rule fixed in §4.

## Primary result (prespec §6)

| Quantity | Value |
|---|---|
| `n_pos` | 18 |
| `n_neg` | 12 |
| `median_pos` | 4.193 |
| `median_neg` | 3.483 |
| `cliffs_delta` | 0.2407 |
| `ci_low` | -0.1852 |
| `ci_high` | 0.6667 |
| `hodges_lehmann` | 0.3878 |
| `auroc` | 0.3796 |
| `u_statistic` | 134 |
| `p_value` | 0.2803 |
| `mwu_method` | asymptotic |
| `min_detectable_delta` | 0.4352 |

## Robustness (prespec §7)

| name | n_pos | n_neg | cliffs_delta | ci_low | ci_high | p_value | auroc | min_detectable_delta | mwu_method |
|---|---|---|---|---|---|---|---|---|---|
| primary_18v12 | 18 | 12 | 0.2407 | -0.1852 | 0.6667 | 0.2803 | 0.3796 | 0.4352 | asymptotic |
| R1_gold_only | 12 | 7 | 0.2857 | -0.3095 | 0.7857 | 0.3402 | 0.3571 | 0.5714 | exact |
| R2_stage1_lines | 16 | 9 | 0.06944 | -0.4167 | 0.5556 | 0.7989 | 0.4653 | 0.4861 | asymptotic |
| R3_drop_D1 | 16 | 12 | 0.3958 | -0.04167 | 0.7917 | 0.08164 | 0.3021 | 0.4479 | exact |
| R3b_drop_all_known | 15 | 11 | 0.4182 | -0.05455 | 0.8182 | 0.0774 | 0.2909 | 0.4667 | exact |
| R4_wrap53_copy_number | 18 | 12 | -0.1667 | -0.5648 | 0.2593 | 0.4648 | 0.5833 | 0.4352 | exact |
| R5_tp53_exposure_descriptive | 18 | 12 | -0.4722 | -0.7963 | -0.08333 | 0.03095 | 0.7361 | 0.4352 | exact |
| R7_uncollapsed_18v13 | 18 | 13 | 0.1709 | -0.2566 | 0.5812 | 0.435 | 0.4145 | 0.4274 | asymptotic |

**R6 leave-one-out:** 30 refits, δ range [+0.192, +0.343], no sign flip = True.

**R8 leave-one-depositor-out:** 12 depositors, no sign flip = True.

| depositor | n dropped | δ | p | note |
|---|---|---|---|---|
| ATCC | 6 | +0.304 | 0.2329 |  |
| CCLF | 1 | +0.255 | 0.2588 |  |
| Clinica Universitaria de Navarra | 1 | +0.216 | 0.341 |  |
| DSMZ | 1 | +0.255 | 0.2588 |  |
| JCRB | 3 | +0.125 | 0.6043 |  |
| MD Anderson Cancer Center | 1 | +0.255 | 0.2588 |  |
| Montefiore Einstein | 3 | +0.086 | 0.7381 |  |
| National Cancer Center Research Institute | 1 | +0.253 | 0.2708 |  |
| Oslo University Hospital, The Norwegian Radium Hospital | 4 | +0.358 | 0.1341 |  |
| RIKEN | 2 | +0.316 | 0.1725 |  |
| The Broad Institute | 1 | +0.263 | 0.2517 |  |
| UCSF Helen Diller Family Comprehensive Cancer Center | 6 | +0.203 | 0.4244 |  |

**R9 bottom-tertile enrichment** (cut 3.335580, strict `<`, descriptive only): low arm [3, 6] (ALT+/ALT−), high arm [15, 6], OR 0.200, Fisher two-sided p = 0.1019.

## Verdict (prespec §8)

**`counter_signal`** — rule `counter-signal`.

Counter-signal: non-significant but substantial effect opposite to the prediction (aligned delta = -0.241, |delta| >= 0.2, p = 0.2803). Not a clean null — the estimate runs the wrong way. This is the cell Stage 1 landed in and that criteria.evaluate left unnamed.

| Condition | Met | Detail |
|---|---|---|
| p<0.05 and delta <= -0.330 | NO | p=0.2803, delta=+0.241 |
| R1 same direction | NO | delta=+0.286 |
| R3 same direction | NO | delta=+0.396 |
| R3b same direction | NO | delta=+0.418 |
| no sign flip on leave-one-out | yes | range [+0.192, +0.343] |
| no sign flip on leave-one-depositor-out | yes | 12 evaluable depositors |
| AUROC >= 0.70 (delta <= -0.400) | NO | AUROC=0.380 |

## Power, as pre-registered (prespec §9)

Minimum detectable |δ| at 18v12 is **0.435**. At Keegan's own 35% lesion prevalence with perfect transcript separation this design has 27% power; under 13% at the realistic attenuated magnitude (|δ| < 0.2). Per §8, a non-significant result in the predicted direction is uninformative and is not evidence against H1.
