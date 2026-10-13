# Notebook 06 — Robustness and sensitivity

The pre-specified robustness set, run **in full** as frozen, despite the notebook 05
primary being non-significant and opposite in direction
(delta = +0.2812, p = 0.3823).

Sign convention: Cliff's delta **< 0** is the predicted direction (TCAB1-low more
BIR-dependent relative to HR).

## Summary table

| Check | Description | n | delta | 95% CI | p | HL | Predicted direction? |
|---|---|---|---|---|---|---|---|
| **R1** | gold-only ALT+ within primary arms | 3v7 | -0.6190 | [-1.000, +0.143] | 0.1833 | -0.3711 | yes |
| **R2** | median split within ALT+ (DEGENERATE) | 8v8 | +0.2812 | [-0.312, +0.844] | 0.3823 | +0.1865 | no |
| **R3** | tertile within ALT+ only (5v11) | 5v11 | +0.5636 | [+0.055, +0.964] | 0.0897 | +0.3726 | no |
| **R5** | ATRX-wildtype ALT+ only (15 of 16) | 7v8 | +0.1786 | [-0.500, +0.786] | 0.6126 | +0.1317 | no |
| **R6** | ALT- negative control (lineage control per §7) | 4v5 | -0.4000 | [-1.000, +0.500] | 0.4127 | -0.1520 | yes |
| **R8** | BIR set without SMARCAL1 (POLD3, RAD52) | 8v8 | +0.1562 | [-0.438, +0.750] | 0.6454 | +0.1953 | no |
| **6.3** | within-line z-scored composite | 8v8 | +0.2812 | [-0.312, +0.844] | 0.3823 | +0.4484 | no |
| **R4** | Spearman(WRAP53 expression, composite) across ALT+ | 16 | rho -0.4415 | — | 0.0869 | — | no |
| **R7** | leave-one-out range over 16 refits | — | [+0.1786, +0.4643] | — | — | no refit reaches it |

## Notes attached to individual checks

**R1.** n=3 low; at 3v7 the exact null median |delta| is 0.33, so unstable

**R2.** Identical to the primary: the §4b guard made the median split the primary rule. Carries no independent information.

**R5.** Near-vacuous as prespec §7 anticipated: only 1 ALT+ line is ATRX-altered, so this drops 1 line. Reported as UNINFORMATIVE, not as evidence of no ATRX confounding.

**R6.** At 4v5 the exact null median |delta| is 0.30 and P(|delta|>=0.474)=0.286, so this is too small to gate on -- prespec §6.2.

**R8.** Guards against the composite restating the published broad SMARCAL1 ALT dependency.

**6.3.** Stronger control for line-level essentiality scale than the composite difference alone.

## R7 leave-one-line-out detail

Arm assignments held fixed; one line dropped per refit. Observed primary delta +0.2812; refit range **[+0.1786, +0.4643]**.

| Dropped | Cell line | Arm | delta without it |
|---|---|---|---|
| ACH-003171 | NRHOS1 | TCAB1-low | +0.1786 |
| ACH-003185 | NRHOS4 | TCAB1-low | +0.1786 |
| ACH-003573 | C337 | TCAB1-low | +0.1786 |
| ACH-000082 | G292CLONEA141B1 | TCAB1-intact | +0.2143 |
| ACH-002433 | CCLFPEDS0019T | TCAB1-low | +0.2143 |
| ACH-000410 | SAOS2 | TCAB1-intact | +0.2500 |
| ACH-002471 | OS052 | TCAB1-intact | +0.2500 |
| ACH-003178 | NRHOS2 | TCAB1-low | +0.2500 |
| ACH-001526 | HUO9 | TCAB1-intact | +0.2857 |
| ACH-002067 | NOS1 | TCAB1-intact | +0.2857 |
| ACH-003587 | 531MII | TCAB1-intact | +0.2857 |
| ACH-002690 | OS526 | TCAB1-intact | +0.3214 |
| ACH-000364 | U2OS | TCAB1-intact | +0.3571 |
| ACH-002834 | OS384 | TCAB1-low | +0.3571 |
| ACH-001715 | CAL72 | TCAB1-low | +0.4286 |
| ACH-002688 | OS457 | TCAB1-low | +0.4643 |

## §6.2 ALT- selectivity and interaction

- Cliff's delta, ALT+ stratum: **+0.2812**
- Cliff's delta, ALT- stratum: **-0.4000**
- ALT- stratum's own two-sided p: **0.4127** (exact)
- delta difference (ALT+ minus ALT-): **+0.6813**
- permutation p for the interaction (5000 shuffles within stratum): **0.1762**
- §8.5 ALT- contradicts (significant same-direction shift): **False**

At 4v5 the exact null median |delta| is 0.30 with P(|delta| >= 0.474) = 0.286, so the
ALT- stratum is too small to gate on. This is why §8.5 tests for a *significant*
same-direction shift rather than comparing deltas (prespec §6.2, §8.5).

## Files

- `robustness_table.csv`
- `loo_deltas.csv`
- `robustness_extras.json`

The §8/§9/§10 verdict is rendered in notebook 07, not here.
