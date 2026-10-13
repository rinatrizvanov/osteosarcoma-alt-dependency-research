# Notebook 07 — Verdict

Computed by `src.criteria.evaluate` from the frozen §8/§9/§10 rules. No new
statistics are performed here.

## Outcome: `indeterminate`

**Section:** none of Sections 8/9/10

The result matches no pre-specified outcome bucket: p=0.3823, delta=+0.281. Failed support conditions: ['primary', 'bir_genes', 'r7_loo']. Report this as-is; do not round it toward the nearest verdict.

## Plain-language reading

**The hypothesis was not supported.** The primary point estimate ran *opposite* to
the prediction and nothing approached significance:

- primary (8v8): Cliff's delta **+0.2812** (predicted < 0), p = **0.3823**
- pre-specified tertile sensitivity (3v13): delta **+0.5385**, p = **0.1893**
- minimum BH q across the six primary genes: **0.9591**
- R7: no refit of 16 reaches the predicted direction; range **[+0.1786, +0.4643]**
- R4 continuous: Spearman rho **-0.4415** (predicted > 0), p = 0.0869

The two checks that do favour the hypothesis (R1 gold-only, delta −0.619 at 3v7; and
R6, the ALT− control, delta −0.400 at 4v5) are the two **smallest** comparisons, where
the exact null median |delta| is 0.30–0.33 and P(|delta| ≥ 0.474) is 27–29%. R6 is
moreover the negative control, where the hypothesis predicts *no* effect, so a
"predicted-direction" result there counts against a TCAB1-specific account, not for it.

## Verdict conditions

| Condition | Met | Detail |
|---|---|---|
| S8.1 primary p<0.05 and delta <= -0.33 | NO | p=0.3823, delta=+0.281 |
| S8.2 gold-only same direction, |delta| >= 0.30 | yes | gold-only delta=-0.619 |
| S8.3 >=2/3 BIR genes at q<0.10, no HR gene stronger | NO | 0/3 BIR genes at q<0.1 in predicted direction (none); strongest BIR delta -0.344; strongest same-direction HR delta -0.250 |
| S8.4 no sign flip on leave-one-line-out | NO | 16 refits, delta range [+0.179, +0.464] |
| S8.5 ALT- stratum shows no significant same-direction shift | yes | ALT- delta=-0.400, p=0.4127; ALT+ delta=+0.281 (delta difference reported as a diagnostic, not a gate) |

**Section 6.4 concentration comparison** (descriptive, non-gating):

- delta, 8v8 median primary: +0.281
- delta, 3v13 OS-reference tertile: +0.538
- |delta| excess of tertile over median: +0.257 (pre-registered threshold 0.3)

No evidence of concentration in the extreme-loss lines: the tertile contrast exceeds the median primary by +0.257 in |delta|, within the pre-registered 0.3 threshold.

## Why not the other buckets

| Bucket | Rule | Met |
|---|---|---|
| §8  supported | `p<0.05 and delta<=-0.33` | False |
| §9  inconclusive | `p>=0.05 and delta<=-0.4` | False |
| §10 clean null | `p>=0.05 and |delta|<0.2` | False |
| §10 wrong-dir sig | `p<0.05 and delta>0` | False |

Observed: p = 0.3823, delta = +0.2812, |delta| = 0.2812.

`indeterminate` is the honest label: |delta| = 0.281 exceeds the
0.2 clean-null ceiling frozen in §10, yet the sign is opposite to
the prediction and p is far from 0.05. **§11 (negative-result fallback) is therefore
NOT authorized** — §10 did not trigger.

## Methods limitation found by this result

§10's clean-null ceiling (|delta| < 0.20) was calibrated for a null *in the predicted
direction*. It did not anticipate a non-significant estimate of moderate size pointing
the **wrong** way, which is what occurred. The rule is left exactly as frozen and this
is recorded as a lesson for any v2 pre-specification. It was **not** retrofitted after
seeing the result.

## Every pre-specified check

| Check | n | delta | p | Predicted direction? |
|---|---|---|---|---|
| R1 | 3v7 | -0.6190 | 0.1833 | yes |
| R2 | 8v8 | +0.2812 | 0.3823 | no |
| R3 | 5v11 | +0.5636 | 0.0897 | no |
| R5 | 7v8 | +0.1786 | 0.6126 | no |
| R6 | 4v5 | -0.4000 | 0.4127 | yes |
| R8 | 8v8 | +0.1562 | 0.6454 | no |
| 6.3 | 8v8 | +0.2812 | 0.3823 | no |
| R4 | 16 | rho -0.4415 | 0.0869 | no |

In the predicted direction: **2 of 7** checks.

## Next step

See `docs/stage1_outcome_memo.md`.
