# Iteration 04 pre-specification — FROZEN

**Hypothesis H6.** Among ALT+ osteosarcoma lines, SMARCAL1 dependency is higher in lines
with more of the replication-stress substrate that SMARCAL1 is required to resolve.

**Primary exposure: PRIMPOL mRNA.** Changed from BLM, on exposure-side grounds only, before
any outcome value was read — see §2, which is the most important section of this document.

**Status:** frozen before any exposure–outcome quantity was computed. No Chronos value has
been read in this iteration.
**Supersedes:** `prespec.DRAFT.md` (committed `6d63efe`, erratum `7cec0ea`).
**Prior information declared:** `DISCLOSURES.md`, committed `d184b0e` before any design
decision.
**Scans informing it:** `scan_literature.md`, `scan_exposure.md`.

---

## 1. Hypothesis and directional prediction

> **Prediction.** Spearman rho between PRIMPOL log2(TPM+1) and SMARCAL1 Chronos gene effect
> across the 16 ALT+ CRISPR-screened osteosarcoma lines is **negative** — more PRIMPOL, more
> negative gene effect, more SMARCAL1-dependent. `predicted_sign = -1`.

**Mechanistic basis.** Taglialatela et al. (2026) *Genes Dev*,
https://doi.org/10.1101/gad.353377.125: SMARCAL1 loss causes telomeric single-stranded DNA
to accumulate, and that accumulation is **PRIMPOL-repriming-dependent**. PRIMPOL repriming
is therefore a named source of the lesion SMARCAL1 is needed to resolve, so lines with more
of it should depend on SMARCAL1 more.

**Why the endpoint is unusually well matched.** Chronos gene effect **is** fitness under
CRISPR knockout of SMARCAL1, and Taglialatela's readout is fitness in SMARCAL1-depleted
cells as a function of the substrate-generating activity. The published experiment and this
endpoint measure the same axis; they differ in that the paper varied the substrate source
experimentally while this iteration exploits its natural variation across lines. Stage 1, by
contrast, joined two literatures nobody had joined (`Literature/INDEX.md` Gaps §4), and
iteration 03 showed its endpoint was not internally coherent.

**Not published.** Taglialatela's own DepMap analysis contrasts telomerase-negative against
telomerase-positive lines, not variation within ALT+; their PRIMPOL and BLM work is isogenic
within single lines. Bennett et al. 2026 (ALTitude,
https://doi.org/10.64898/2026.03.04.709438) tested the leading within-ALT+ candidate, ATRX,
and **published a null** (Cohen's d = −0.248, P = 0.5, n = 20 vs 14), while explicitly
flagging within-ALT+ heterogeneity as mechanistically unexplained. **Bennett's ATRX null and
their reciprocal ATRX–SMARCAL1 compensation model must not be re-derived here**; both are
recorded in `NULL_SPACE.md`.

## 2. The exposure change, and why it weakens this iteration's prior

The draft named **BLM** as primary, because Taglialatela's BLM evidence is the stronger of
the two: *"BLM deficiency fully suppressed the SMARCAL1-dependent growth defects"* and
*"overexpression of WT BLM... markedly exacerbated... fitness defects in SMARCAL1-depleted
cells"* — full suppression, against PRIMPOL's merely partial rescue, and PRIMPOL's
exacerbation required supraphysiological overexpression.

**BLM cannot carry a correlation in this cohort.** `scan_exposure.md` ranks it **15th of
16** candidates on distributional grounds, and I verified its numbers independently rather
than taking them on trust. Over the 16 ALT+ lines BLM log2(TPM+1) is:

```
0.96  2.32  |  3.30 3.31 3.35 3.47 3.68 3.77 3.84 3.94 3.95 4.15 4.16 4.19 4.26 4.79
```

Two low lines, then **14 lines inside a 1.496-unit band**. Shapiro p = 0.0059, excess
kurtosis +2.88, skew −1.69, and the largest interior gap respecting
`stratify.MIN_CELL_SIZE = 4` on both sides is only **0.220**. A rank correlation on that
distribution is decided by where those two lines fall, so it is a **2-versus-14 contrast
wearing the clothes of a correlation** — and 2 is below the minimum-cell guard, so the
honest categorical form of the same comparison is barred outright. BLM is also the hub of
the panel's collinearity cluster (|rho| +0.737 with FANCM over the 34 lines, +0.871 over the
16).

**PRIMPOL carries it.** Over the same 16 lines: 2.68 to 4.85, SD 0.761, skew +0.29, excess
kurtosis **−1.38** (flatter than normal, i.e. genuinely spread rather than clumped), Shapiro
p = 0.117, no zeros, largest guarded interior gap 0.454. Over the 34-line reference set it is
the cleanest candidate in the panel (skew +0.07, Shapiro p = 0.761). It is also acceptably
independent of the negative control (§5).

**This swap is legitimate, and it is also bad news.** It is legitimate because it rests
**only** on exposure-side distributions computed on the 34-line reference superset, with no
outcome value read — exactly the decision §6 of the draft reserved for the scan, and exactly
the discipline anti-bias rule (b) requires. It is bad news because the mechanism's
**stronger** arm is the unusable one: this iteration now tests the weaker of two
predictions, on a proxy, at n = 16. The prior that it returns anything is correspondingly
lower, and saying so here is cheaper than discovering it afterwards.

**BLM is not dropped silently.** It is reported in §5 S4 as explicitly descriptive and
**uninterpretable**, with its 2-outlier structure restated at the point of reporting, so a
large value cannot be mistaken for a finding. It is also logged to `NULL_SPACE.md` as a
binding fact about this cohort, alongside the parallel Stage 1 finding that WRAP53 has no
natural low group here.

## 3. Exposure-variance gate — the rule that disqualified BLM

Pre-committed, and applied to **every** exposure this iteration touches, before its outcome
is read. An exposure is usable only if, over the 16 analysed lines, all three hold:

| Criterion | Threshold | PRIMPOL | BLM |
|---|---|---|---|
| SD of log2(TPM+1) | >= 0.50 | **0.761** pass | 0.895 pass |
| lines at exactly 0 | none | **0** pass | 0 pass |
| largest guard-respecting interior gap | < 1.00 unit | **0.454** pass | 0.220 pass |
| not outlier-dominated: Shapiro-Wilk p | >= 0.05 | **0.117** pass | **0.0059 FAIL** |

The Shapiro criterion is what BLM fails, and it is the criterion that encodes the actual
defect: a distribution whose rank order is set by two extreme points cannot support a rank
correlation at n = 16. Thresholds are stated numerically here, before the outcome, and were
chosen against the **reference** distributions rather than against any result.

**Stated plainly to avoid overclaiming:** PRIMPOL was verified to pass this gate *before*
freezing, using exposure data only. So the gate does not function as a live risk to this
iteration — its work was done in §2, in disqualifying BLM. It is written down so that
disqualification is a rule applied, not a judgement made.

## 4. Cohort, outcome, test and power

**Outcome.** SMARCAL1 Chronos gene effect per line, via `src.dependency_stats`. The only
Stage 1 BIR gene with usable dynamic range here: ALT+ median −0.462, SD 0.395, 8 of 16 lines
below −0.5 and 3 below −1.0 (iteration 03 §2.1). POLD3 is essential in all 16 lines and
RAD52 in none, which is why neither can carry this analysis. **SMARCAL1 mRNA appears nowhere
in this design** — not as exposure, covariate or control — because Bennett documents an
mRNA-to-protein dissociation for SMARCAL1 itself (line SKNFI).

**Cohort.** The 16 ALT+ CRISPR-screened osteosarcoma lines from `src.cohort.build_cohort()`,
`assert_all()` first; ALT labels frozen, never re-derived. The 9 ALT− lines are used only for
the §5 S3 specificity check: SMARCAL1 is near-inert there (median +0.056, SD 0.137, 0 of 9
below −0.5), so pooling would add outcome rows with almost no variance while importing the
ALT main effect. **Minimum-cell guard: n >= 12** lines with both exposure and outcome, below
which the primary is not estimable and the iteration is logged as a feasibility failure.

**Test.** Two-sided permutation test on Spearman rho, 10,000 permutations,
`SEED = 20261012`, reusing `iteration_03_bir_module_coherence/src/coherence.py`, whose null
is distribution-free given (n, k = 2) and whose implementation is validated against scipy by
58 tests. Two-sided p with a one-sided success rule, on the Stage 1 §5 precedent: a
two-sided p stays interpretable while requiring the predicted direction stops a significant
*opposite* result being read as confirmation.

**Power — the binding constraint, stated before the result.** From the permutation null at
n = 16, k = 2 (200,000 draws):

| Requirement | smallest \|rho\| attaining it |
|---|---|
| two-sided p < 0.05 | **0.500** |
| **two-sided p < 0.033 — the campaign bar (§6)** | **0.535** |
| one-sided p < 0.05 | 0.426 |

**Pre-committed consequence.** A non-significant result in the predicted direction is
**uninformative about the mechanism** and will not be reported as evidence against
Taglialatela's model. The honest expected outcome of this iteration is inconclusive. It is
run because it is the best-grounded within-ALT+ hypothesis available on held data, and
because iteration 03 established that this cohort *can* carry a real effect when one exists
— BTR reached mean pairwise rho +0.526 at p = 0.0001 on these same 16 lines. What is
unavailable is sensitivity to a moderate effect.

## 5. Controls and secondaries — pre-specified

Reported with their own effects and p; **none can change the primary verdict**; all counted
in the campaign test total.

- **S1 — negative control exposure: ZRANB3 mRNA, predicted null.** A SNF2-family fork
  remodeler functionally paralogous to SMARCAL1 but named nowhere in the motivating
  mechanism. Adequately independent of the primary exposure: rho(ZRANB3, PRIMPOL) = +0.322
  over the 34 lines and +0.185 over the 16, well below the |rho| = 0.7 collinearity bar.
  **Declared weakness:** ZRANB3's own spread is narrow (SD 0.427 over the 16, range 1.58,
  the smallest of any candidate carried forward), so it is a **weak** control — it is better
  at catching a confound that acts on all fork remodelers than at proving specificity. Said
  now rather than after. Note it is *not* an adequate control for BLM (rho +0.600 over the
  34), which is a further reason BLM is not the primary.
- **S2 — screen-amplitude confound.** Spearman rho between the exposure and each line's
  **genome-wide mean gene effect** across all ~18.5k genes. A line with a stronger screen
  looks more dependent on everything, so an exposure tracking screen amplitude would
  manufacture a primary. If |rho| here exceeds the primary's |rho|, the primary is declared
  confounded and reported as such.
- **S2b — standardised repeat.** The primary recomputed on SMARCAL1 effect z-scored within
  each line across all genes (Stage 1 §6.3's transform), removing line-level location and
  scale.
- **S3 — ALT specificity.** The same correlation in the 9 ALT− lines, where the mechanism
  predicts no relationship because SMARCAL1 is near-inert. Descriptive: the floor at n = 9
  is 0.683.
- **S4 — BLM, descriptive and uninterpretable.** Reported for completeness with its
  2-outlier structure restated in the output itself. It **fails the §3 gate** and therefore
  cannot support or undermine anything; a large value is not evidence. Registered this way
  so the failure is visible rather than absent.
- **S5 — leave-one-out.** 16 refits; report the rho range and whether any flips the sign.
- **S6 — rule (c) not triggered, limitation restated.** The outcome is a **single gene**, so
  no composite endpoint is used and rule (c) does not fire. The underlying limitation stands
  as recorded in `NULL_SPACE.md`: no defensible type I/II ALT marker exists in DepMap 26Q1,
  and none is invented here.

**Not used: TOP3A**, despite ranking 1st distributionally. de Nonneville et al.
(https://doi.org/10.15252/emmm.202215859) report TOP3A **amplification** driving ALT in
pediatric osteosarcoma, mutually exclusive with ATRX inactivation. Excluded because the
published claim is a **copy-number** claim, and gene-level copy number is the measurement
class iteration 02 found unreliable here without segment validation — and `scan_exposure.md`
independently flags TOP3A as the most dosage-driven candidate in the panel
(expression–CN rho +0.722). Recorded so the exclusion is a decision, not an oversight.

## 6. Verdict mapping and campaign accounting

`evaluate_iteration(delta=rho, p_value=p_perm, predicted_sign=-1, conditions=...)`. No new
thresholds and no iteration-local evaluator.

Support conditions: `exposure_variance_gate` (§3), `negative_control_null` (S1 shows no
relationship and a smaller |rho| than the primary), `not_screen_confounded` (S2),
`loo_no_sign_flip` (S5).

**Tests this iteration contributes: 22.** Primary 1; S1 1; S2 1; S2b 1; S3 1; S4 1; LOO 16.
Only the primary p enters the campaign BH family.

**The campaign criterion, restated because it is the one that matters.** A result counts as
positive only at **q < 0.1 across every iteration's primary endpoint**, not p < 0.05 within
this iteration. With three iterations carrying a primary that requires **p < 0.033**, hence
|rho| >= 0.535 per §4. Any positive is reported with the number of tests that preceded it —
**79** before this iteration begins — and labelled exploratory.

## 7. What this iteration cannot do

It cannot measure PRIMPOL repriming activity, only PRIMPOL transcript. It cannot establish
causation from 16 observational cell lines. It cannot test the mechanism's stronger arm at
all, because BLM's distribution in this cohort forbids it (§2). It cannot rescue Stage 1 or
iteration 03. And at n = 16 it cannot distinguish a moderate true effect from no effect —
§4 is explicit that the likely outcome is inconclusive, and that outcome will not be dressed
up as either support or refutation.
