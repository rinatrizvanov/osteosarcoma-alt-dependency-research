# Iteration 04 outcome — H6: does PRIMPOL expression predict SMARCAL1 dependency?

**Verdict: `null`** (`rule = clean-null`). Spearman rho between PRIMPOL log2(TPM+1) and
SMARCAL1 Chronos gene effect across the 16 ALT+ lines is **exactly +0.0000**, p = 1.000.

Prespec frozen `5f2f039`; analysis code committed `6eb0908`, before execution.
No post-hoc rule change. Primary exposure had already been switched from BLM to PRIMPOL on
exposure-side grounds before the freeze (prespec §2).

---

## 1. The primary is an exact zero, and that was checked rather than assumed

| | value |
|---|---|
| Spearman rho (PRIMPOL vs SMARCAL1 Chronos) | **+0.0000** |
| two-sided permutation p | **1.0000** |
| n | 16 |
| detection floor at this n | 0.4943 |
| aligned delta (rho x predicted_sign, predicted −1) | −0.0000 |

An exact zero is unusual enough to be a plausible bug, so it was verified against an
independent implementation before being written up. `scipy.stats.spearmanr` returns
statistic 0.0, p = 1.0; Pearson on the ranks returns 6.9e-18; and the rank-difference
identity is exact: **sum of d² = 680**, against n(n² − 1)/6 = **680** for n = 16. There are
no ties in either variable. The zero is real, not a degenerate array or a misalignment.

**Post-hoc descriptive bound, labelled as such.** A 10,000-draw bootstrap over lines — **not
pre-registered**, reported because a null estimate without an interval is uninformative —
gives a 95% interval of **[−0.502, +0.522]** with median −0.007. So the data exclude a very
large effect in the predicted direction (P(rho < −0.5) = 0.025 under the bootstrap) and
exclude nothing smaller: P(rho < −0.2) = 0.231. This is a clean null at the point estimate,
not a demonstration that the mechanism is absent.

## 2. Three of four support conditions failed, and all three failed vacuously

This is the most important thing to state correctly, because the condition table looks much
worse than the evidence is.

| Condition | Result | Is it informative? |
|---|---|---|
| `exposure_variance_gate` | **PASS** — SD 0.761, Shapiro p 0.117, guarded gap 0.454 | yes |
| `negative_control_null` | FAIL — ZRANB3 rho +0.165, p 0.540 | **no — vacuous** |
| `not_screen_confounded` | FAIL — screen-amplitude rho −0.279 vs primary +0.000 | **no — vacuous** |
| `loo_no_sign_flip` | FAIL — range [−0.093, +0.186] | **no — vacuous** |

All three failures are arithmetic consequences of the primary being exactly zero:

- `negative_control_null` required the control's |rho| to be **smaller than the primary's**.
  With the primary at 0.000, *any* non-zero control fails it. ZRANB3's own p is 0.540 —
  i.e. the negative control is behaving exactly as a negative control should.
- `not_screen_confounded` required |rho(exposure, screen amplitude)| < |primary rho|, again
  unachievable against zero.
- `loo_no_sign_flip` required a consistent sign, but the sign of a zero estimate is
  arbitrary, so the leave-one-out range straddles zero by construction.

**These conditions were mis-specified**: each presupposed a non-zero primary and compared
against its magnitude. Written correctly they would have been absolute tests ("the control's
p >= 0.05", "|rho(exposure, screen amplitude)| < 0.30") rather than relative ones. The
failures are therefore **not** evidence against the hypothesis and must not be counted as
such. Recorded as a prespec defect found by the data — the second such defect in two
iterations, after iteration 03's single-control assay gate — and **not** retroactively
rewritten.

## 3. The negative control itself failed the exposure gate

ZRANB3 has SD **0.427** over the 16 lines against the §3 gate's `>= 0.50`, so it is not a
gate-qualified exposure. The prespec **pre-declared this weakness** (§5 S1: "ZRANB3's own
spread is narrow... so it is a **weak** control"), so this is a predicted limitation rather
than a surprise, but it does mean S1 carries reduced weight.

**Methods note for future prespecs, not a waiver now.** `scan_exposure.md` argues the SD
criterion is the wrong instrument for a rank-based control: all 16 ZRANB3 values are
distinct, so Spearman has full rank resolution, and only variance-scaled estimators
(Pearson, OLS slope) lose precision from a narrow spread. That is correct, and a future gate
should apply the SD criterion only to exposures feeding variance-scaled estimates. The gate
is **not** waived here — it was frozen and it applies as written.

## 4. Secondaries

| Analysis | rho | p | reading |
|---|---|---|---|
| S1 ZRANB3 negative control | +0.165 | 0.540 | no relationship, as predicted |
| S2 screen amplitude | −0.279 | 0.297 | PRIMPOL weakly tracks screen amplitude; not significant, but the largest |rho| in the table |
| S2b standardised outcome | +0.012 | 0.970 | the exact zero is not an artefact of line-level scale |
| S3 ALT− specificity (n = 9) | −0.133 | 0.750 | nothing, as the mechanism predicts (floor there is 0.683) |
| S4 BLM — **UNINTERPRETABLE** | +0.188 | 0.481 | fails the §3 gate (Shapiro p 0.0059); 2 outliers + 14 lines in a 1.496-unit band. Not evidence either way |

S2b matters most: standardising the outcome within each line across all ~18.5k genes moves
rho from +0.000 to +0.012. The null is robust to the screen-scale transform that Stage 1
used for the same purpose.

S2 is the one number worth watching. PRIMPOL expression correlates −0.279 with each line's
genome-wide mean gene effect, which is larger in magnitude than the primary itself. It is
not significant and the prespec's confounding rule is keyed to a comparison that a zero
primary makes meaningless, but it is a reminder that expression variables in this cohort
carry some screen-amplitude signal.

## 5. What this iteration establishes

**It does not refute Taglialatela's model, and the prespec pre-committed to not claiming
that** (§4: "A non-significant result in the predicted direction is uninformative about the
mechanism and will not be reported as evidence against it"). Three reasons stand:

1. **Power.** The floor at n = 16 is |rho| >= 0.500. The bootstrap interval spans
   [−0.502, +0.522]. A moderate true effect is entirely compatible with this result.
2. **The exposure is a transcript, not an activity.** PRIMPOL repriming is licensed by fork
   stalling, not by abundance. The prespec declared this gap rather than patching it with
   proteomics that is not in the frozen raw set.
3. **The mechanism's stronger arm was never testable.** Taglialatela's BLM evidence is *full*
   suppression of SMARCAL1-dependent growth defects; PRIMPOL gives only partial rescue and
   needed supraphysiological overexpression to exacerbate. BLM's distribution in these 16
   lines (2 low outliers, 14 lines inside 1.496 units) forbids a rank correlation, and its
   only natural low group has n = 2, below the minimum-cell guard. The prespec said so in
   advance (§2) and called the resulting prior "correspondingly lower".

**What it does establish** is narrower and worth recording: among ALT+ osteosarcoma lines,
PRIMPOL transcript level carries **no monotone information whatsoever** about SMARCAL1
dependency — not a weak signal, an exact zero, robust to within-line standardisation. Given
that SMARCAL1 dependency genuinely varies here (8 of 16 lines below −0.5, 3 below −1.0),
something else explains that variation, and PRIMPOL mRNA is not a trace of it.

Combined with iterations 02 and 03 and with Bennett's published ATRX null, the honest
summary is that **this cohort has now failed to explain within-ALT+ SMARCAL1 dependency from
every axis it can measure**: ATRX status (not estimable, and published null at pan-cancer
scale), BLM (no usable variance), PRIMPOL (exact zero).

## 6. Campaign accounting

**Tests run: 22**, as the prespec enumerated — primary 1, S1 1, S2 1, S2b 1, S3 1, S4 1,
leave-one-out 16. The post-hoc bootstrap interval (§1) is an interval estimate, not a test,
and is not counted. **Campaign total: 101.**

**Campaign BH across all three iterations carrying a primary** (`src.verdict.campaign_fdr`):

| Iteration | raw p | campaign q | survives q < 0.1 |
|---|---|---|---|
| 03 | 0.2285 | 0.4200 | no |
| 01 | 0.2800 | 0.4200 | no |
| 04 | 1.0000 | 1.0000 | no |

**No campaign positive.** Note that adding this iteration's primary *raised* q for
iterations 01 and 03 from 0.280 to 0.420, which is the multiplicity rule working as intended:
each additional pre-registered primary makes the family-wide bar stricter for every member.

## 7. Streak and stop conditions

**Fourth consecutive non-positive iteration** (01 `counter_signal`, 02 `not-estimable`,
03 `not-estimable`/`null`, 04 `null`). Stop condition (b) fires at five consecutive
inconclusive-or-null results **with no new mechanistic angle surfacing**.

**No new mechanistic angle surfaced in this iteration.** That is a deliberate statement
rather than an omission: the result is an exact zero with no secondary pointing anywhere,
and the one quantity that moved (S2's −0.279 screen-amplitude correlation) is a nuisance
signal, not a mechanism. Iteration 03's HR/BRCA2 angle was contaminated by the iteration 03
scan and withdrawn in `DISCLOSURES.md`, so it is not available either.

**Consequence, stated plainly for the next session.** If iteration 05 returns
inconclusive-or-null without surfacing a new mechanistic angle, **stop condition (b) fires
and `FINAL.md` must be written.** The loop's own rule forbids fabricating further hypotheses
to stay alive, so iteration 05 should be the best remaining literature-grounded candidate,
chosen from `Literature/INDEX.md`'s Prespec Conflicts table and the untried rows of
`next_hypotheses.md` — not a variant of anything already in `NULL_SPACE.md`.
