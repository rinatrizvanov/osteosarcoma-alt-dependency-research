# Notebook 03 — ALT classification and cohort construction

Labels taken from `references/master_manifest.csv` (frozen). Never re-derived.
Reads no gene-effect values.

## Inclusion filter

| Step | N |
|---|---|
| Model.csv rows matching 'steosarcoma' | 35 |
| present in CRISPRGeneEffect.csv (screened) | 26 |
| present in expression matrix | 26 |
| joined to frozen manifest (cohort == 'DepMap_26Q1') | 26 |
| alt_status in ('alt_gold_positive', 'alt_silver', 'alt_negative_high_conf') | 25 |
| alt_call_direction in ('positive', 'negative') | 25 |
| ALT+ (primary contrast set) | 16 |
| ALT- (negative-control set) | 9 |

## Analysis sets

- ALT+ primary contrast: **N=16** (10 gold + 6 silver)
- ALT- negative control: **N=9**
- Excluded (`alt_unknown`): **N=1**

Total 26 = every CRISPR-screened osteosarcoma line.

## Pre-analysis assertions

All PASSED: ground-truth anchors, set sizes, tier x direction breakdown.

| Anchor | Expected | Observed |
|---|---|---|
| ACH-000364 (U2OS) | positive | positive |
| ACH-000748 (SJSA1) | negative | negative |

| ALT tier | Direction | N |
|---|---|---|
| alt_gold_positive | positive | 10 |
| alt_negative_high_conf | negative | 7 |
| alt_silver | negative | 2 |
| alt_silver | positive | 6 |

## ALT+ primary contrast set

| ModelID | Cell line | ALT tier | ATRX | Telomerase |
|---|---|---|---|---|
| ACH-000082 | G292CLONEA141B1 | alt_gold_positive | not_altered | indeterminate |
| ACH-000364 | U2OS | alt_gold_positive | not_altered | indeterminate |
| ACH-000410 | SAOS2 | alt_gold_positive | not_altered | indeterminate |
| ACH-001526 | HUO9 | alt_gold_positive | not_altered | indeterminate |
| ACH-001715 | CAL72 | alt_gold_positive | not_altered | indeterminate |
| ACH-002067 | NOS1 | alt_gold_positive | not_altered | indeterminate |
| ACH-002433 | CCLFPEDS0019T | alt_silver | not_altered | indeterminate |
| ACH-002471 | OS052 | alt_gold_positive | not_altered | indeterminate |
| ACH-002688 | OS457 | alt_gold_positive | not_altered | indeterminate |
| ACH-002690 | OS526 | alt_gold_positive | not_altered | indeterminate |
| ACH-002834 | OS384 | alt_gold_positive | not_altered | indeterminate |
| ACH-003171 | NRHOS1 | alt_silver | altered | indeterminate |
| ACH-003178 | NRHOS2 | alt_silver | not_altered | indeterminate |
| ACH-003185 | NRHOS4 | alt_silver | not_altered | indeterminate |
| ACH-003573 | C337 | alt_silver | not_altered | indeterminate |
| ACH-003587 | 531MII | alt_silver | not_altered | indeterminate |

## ALT- negative control set

| ModelID | Cell line | ALT tier | ATRX | Telomerase |
|---|---|---|---|---|
| ACH-000359 | MG63 | alt_negative_high_conf | not_altered | positive |
| ACH-000748 | SJSA1 | alt_negative_high_conf | not_altered | positive |
| ACH-001001 | 143B | alt_negative_high_conf | not_altered | positive |
| ACH-001275 | SYO1 | alt_silver | not_altered | positive |
| ACH-001814 | OS252 | alt_negative_high_conf | not_altered | positive |
| ACH-002069 | HSOS1 | alt_negative_high_conf | not_altered | positive |
| ACH-002142 | HUO3N1 | alt_negative_high_conf | not_altered | indeterminate |
| ACH-002470 | OS186 | alt_negative_high_conf | not_altered | indeterminate |
| ACH-003176 | NRHGCT1 | alt_silver | not_altered | positive |

## R5 (ATRX control) feasibility

ATRX-altered ALT+ lines: **1 of 16**. The ATRX-wildtype subset is 15 of 16, so R5 cannot meaningfully control for ATRX. Prespec §7 recorded this expectation in advance; R5 will be reported as uninformative rather than as evidence of no ATRX confounding.

Likely cause: Aim 2 derived `atrx_status` from damaging point mutations only
(`OmicsSomaticMutationsMatrixDamaging`), which misses ATRX deletions and loss of
expression -- the dominant ATRX-inactivation routes in ALT.
