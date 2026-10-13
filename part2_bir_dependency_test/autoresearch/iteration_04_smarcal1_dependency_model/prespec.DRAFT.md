# Iteration 04 pre-specification — DRAFT

**Hypothesis H6.** Among ALT+ osteosarcoma lines, SMARCAL1 dependency is higher in lines
with higher BLM expression.

**Status: DRAFT.** The literature scan **has** been read and is cited throughout. The
exposure scan has **not** been read, and everything depending on the realised exposure
distributions is left explicitly open in §6.

> **Erratum, recorded immediately and not silently edited.** The sentence originally here
> claimed this draft was written *before `scan_exposure.md` existed on disk*, and cited the
> `ls` in its own commit as proof. **That claim was false.** The file was present in that
> `ls` output: the scanning agent had finished writing it in the interval between the last
> time its absence was checked — when it genuinely did not exist — and this commit.
>
> The accurate and narrower guarantee is this: **no tool call in this session has read any
> part of `scan_exposure.md` at the time of the draft commit.** That is a claim about what
> was read, not about what existed, and it is the only one the record supports. It is
> weaker than the original, because file absence is checkable by a third party from the
> commit while "I did not read it" rests on the session transcript.
>
> Corrected here rather than rewritten, on the Stage 1 v1.2 precedent: a methods document
> that overstates its own guarantees is worse than one that states a narrower guarantee
> accurately. Nothing analytical changes — the §6 items remain open and are still to be
> fixed from the scan.

**Prior information declared:** `DISCLOSURES.md`, committed `d184b0e` before any design
decision for this iteration.

---

## 1. Hypothesis, prediction, and why the endpoint is unusually well matched

**H6.** Within the 16 ALT+ CRISPR-screened osteosarcoma lines, higher BLM mRNA predicts
greater SMARCAL1 essentiality.

> **Prediction.** Spearman rho between BLM log2(TPM+1) and SMARCAL1 Chronos gene effect
> across the 16 ALT+ lines is **negative** — higher BLM, more negative gene effect, more
> dependent. `predicted_sign = -1`.

**Mechanistic basis.** Taglialatela et al. (2026) *Genes Dev*,
https://doi.org/10.1101/gad.353377.125. Two findings, in the paper's own words as extracted
in `scan_literature.md`:

- *"overexpression of WT BLM, but not the KA mutant, markedly exacerbated DNA damage,
  telomeric RPA accumulation, APB formation and fitness defects in SMARCAL1-depleted
  cells"*
- *"BLM deficiency fully suppressed the SMARCAL1-dependent growth defects"*

SMARCAL1's annealing activity counteracts BLM unwinding, so BLM activity is the source of
the lesion SMARCAL1 is needed to resolve. Raising BLM is close to phenotypically silent on
its own and becomes costly only when SMARCAL1 is absent.

**Why this endpoint is better matched than anything the project has tried.** Chronos gene
effect **is** fitness under CRISPR knockout of SMARCAL1. Taglialatela's readout is fitness
in SMARCAL1-depleted cells as a function of BLM level. The published experiment and this
endpoint are measurements on the same axis, differing in that they varied BLM
experimentally while we exploit its natural variation across lines. Stage 1, by contrast,
joined two literatures that had never been joined (`Literature/INDEX.md` Gaps §4), and
iteration 03 showed its endpoint was not even internally coherent.

Three mechanistic channels — substrate supply, SMARCAL1 recruitment, and ALT activity —
all give the same sign, and `scan_literature.md` found none running the other way.

**Not published.** Taglialatela's own DepMap analysis contrasts telomerase-negative against
telomerase-positive lines, not variation within ALT+; their BLM and PRIMPOL work is
isogenic within single lines. Bennett et al. 2026 (ALTitude) tested the leading
within-ALT+ candidate, ATRX, and **published a null** (d = −0.248, P = 0.5, n = 20 vs 14),
and explicitly flags within-ALT+ heterogeneity as mechanistically unexplained. **Bennett's
ATRX null and their reciprocal ATRX–SMARCAL1 compensation model must not be re-derived
here** — both are now recorded in `NULL_SPACE.md`.

## 2. Why the exposure is mRNA, and why that is not the project's old mistake twice over

This project has been burned by two proxies, so the burden is on this one.

- **WRAP53 mRNA** (Stage 1) stood in for a *structural variant* — a different kind of
  object from the measurement, with no supporting corpus anywhere
  (`Literature/INDEX.md` Gaps §3).
- **Gene-level copy number** (iteration 02) was a measurement error: three of four apparent
  homozygous ATRX deletions were refuted by segment data.

BLM mRNA is a different case. The perturbation in the motivating experiment **is a level
change in the same molecule** — Taglialatela overexpress and deplete BLM and read fitness —
so mRNA is a proxy for the quantity the mechanism actually varies, not for a structural
event. BLM is also a named ALT-up gene in Bennett's DepMap RNA-seq signature.

**It remains a proxy, and the limitation is stated now rather than later.** BLM activity is
licensed by recruitment and fork stalling, not by abundance alone, so mRNA can only track
the mechanism's input imperfectly. The right measurement is BLM protein; CCLE proteomics is
**not** in this project's frozen raw set, and per §1 of `PROGRAM.md` a hypothesis needing
data outside that set is skipped rather than proxied — so this iteration proceeds on mRNA
with the gap declared, and no proteomic substitute is fetched.

**SMARCAL1 mRNA is excluded from this design entirely.** `scan_literature.md` found that
Bennett documents an mRNA-to-protein dissociation for SMARCAL1 itself (line SKNFI), so
SMARCAL1 transcript level cannot appear as exposure, covariate or control.

## 3. Cohort and outcome

**Outcome.** SMARCAL1 Chronos gene effect, per line, from `src.dependency_stats`. It is the
only one of the Stage 1 BIR genes with usable dynamic range in this cohort — ALT+ median
−0.462, SD 0.395, 8 of 16 lines below −0.5 and 3 below −1.0 (iteration 03 §2.1,
pre-registered there). POLD3 is essential in all 16 lines and RAD52 in none, which is why
neither can carry this analysis.

**Cohort.** The 16 ALT+ CRISPR-screened osteosarcoma lines from `src.cohort.build_cohort()`,
`assert_all()` first. ALT labels frozen from `references/master_manifest.csv`, never
re-derived.

**The ALT− arm is not pooled in.** SMARCAL1 is near-inert there (ALT− median +0.056, SD
0.137, 0 of 9 lines below −0.5), so adding those lines would contribute outcome rows with
almost no variance to explain while importing the ALT main effect. The 9 ALT− lines are used
only as a pre-specified specificity check (§7 S3).

**Minimum-cell guard:** n >= 12 lines with both exposure and outcome present. Below that the
primary is not estimable and the iteration is logged to `NULL_SPACE.md` as a feasibility
failure, as iteration 02 was.

## 4. Test and effect size

**Statistic.** Spearman rho between BLM log2(TPM+1) and SMARCAL1 Chronos gene effect across
the analysed lines. Bounded on [−1, +1], the same scale as Cliff's delta, so it feeds
`src.verdict.evaluate_iteration` unmodified with `predicted_sign = -1` and inherits the
campaign's `NULL_CEILING = 0.200`.

**Test.** Two-sided permutation test on rho, 10,000 permutations, `SEED = 20261012`, reusing
the iteration 03 machinery (`iteration_03_bir_module_coherence/src/coherence.py`), whose null
is distribution-free given (n, k=2) and whose implementation is validated against scipy by
58 tests.

**Two-sided p with a one-sided success rule**, following the Stage 1 precedent (prespec §5):
reporting a two-sided p keeps the p-value interpretable, while requiring the predicted
direction for support prevents a significant *opposite* result being read as confirmation.

## 5. Power — the binding constraint, stated before the result

Computed from the permutation null at n = 16, k = 2 (200,000 draws):

| Requirement | smallest \|rho\| that attains it |
|---|---|
| two-sided p < 0.05 | **0.500** |
| **two-sided p < 0.033 — the campaign bar (§8)** | **0.535** |
| two-sided p < 0.025 | 0.562 |
| one-sided p < 0.05 | 0.426 |

**Read this before interpreting any result.** At n = 16 a single correlation can reach
nominal significance only at |rho| >= 0.500, and the campaign-level criterion needs
|rho| >= 0.535. This design is powered **only for a very large effect**, which is the same
structural limitation Stage 1 had and is a property of having 16 screened ALT+ lines, not a
flaw in the hypothesis.

**Consequence, pre-committed.** A non-significant result in the predicted direction is
**uninformative about the mechanism** and must not be reported as evidence against
Taglialatela's model. It will resolve to `inconclusive_underpowered` or `null` by the
partition and be reported as exactly that. The honest expected outcome of this iteration is
inconclusive; it is run because it is the best-grounded within-ALT+ hypothesis available on
held data, and because iteration 03 established that this cohort *can* carry a real effect
when one exists.

Iteration 03's BTR result (mean pairwise rho +0.526, p = 0.0001 at n = 16) is the relevant
calibration: effects of the required magnitude are detectable here. What is not available is
sensitivity to a moderate one.

## 6. Left open for the exposure scan

To be fixed in the frozen prespec from `scan_exposure.md`, and **only** from it:

1. **The exposure-variance gate.** `scan_literature.md` recommends, and this draft adopts,
   a gate read **before the outcome**: if BLM log2(TPM+1) across the analysed lines has
   insufficient spread to support a correlation, the iteration is a feasibility failure and
   no outcome value is read. The numeric threshold is set in the frozen prespec from the
   34-line reference distribution, never from the 16-line analysis subset (anti-bias rule b),
   and is committed before the outcome is touched.
2. **Zero-inflation and shape.** If BLM is zero-inflated or has no usable dynamic range, the
   primary is not estimable and the iteration stops there.
3. **Collinearity of the negative control.** ZRANB3 is registered below as a predicted-null
   control; if the scan finds it collinear with BLM above |rho| = 0.7 it is **not** a
   control and must be replaced or dropped, with the decision recorded.
4. **Whether any candidate has a natural low group of size >= 4**, which would offer a
   categorical sensitivity analysis alongside the continuous primary.

## 7. Controls and secondaries — pre-specified

Each is reported with its own effect and p; **none can change the primary verdict**; all are
counted in the campaign test total.

- **S1 — negative control exposure: ZRANB3 mRNA, predicted null.** A SNF2-family fork
  remodeler functionally paralogous to SMARCAL1 but named nowhere in the motivating
  mechanism. Prediction: no relationship. If ZRANB3 correlates with SMARCAL1 dependency as
  strongly as BLM does, the primary is not specific and that is reported as undermining it.
  Subject to §6 item 3.
- **S2 — screen-amplitude confound.** Spearman rho between the exposure and each line's
  **genome-wide mean gene effect**. A line with a stronger screen looks more dependent on
  everything, so an exposure correlated with screen amplitude would produce a spurious
  primary. Reported as a diagnostic; if |rho| here exceeds the primary's |rho|, the primary
  is declared confounded and reported as such.
- **S2b — standardised repeat.** The primary recomputed on SMARCAL1 effect z-scored within
  each line across all ~18.5k genes (Stage 1 §6.3's transform), which removes line-level
  location and scale.
- **S3 — ALT specificity.** The same correlation in the 9 ALT− lines, where the mechanism
  predicts **no** relationship because SMARCAL1 is near-inert there. Descriptive: at n = 9
  the floor is 0.683.
- **S4 — first alternate exposure: PRIMPOL mRNA**, same predicted sign. The mechanism's only
  other named source gene, and weaker: endogenous PRIMPOL depletion gave only partial
  rescue, and exacerbation required supraphysiological overexpression. Registered as
  secondary so it cannot be promoted to primary after the fact.
- **S5 — leave-one-out.** 16 refits, reporting the rho range and whether any flips the sign.
- **S6 — subtype secondary not triggered, and the standing limitation restated.** Rule (c)
  attaches to composite endpoints; this iteration's outcome is a **single gene**, so no
  composite is used and rule (c) does not fire. The underlying limitation stands and is
  restated as required: no defensible type I/II ALT marker exists in DepMap 26Q1, per
  `NULL_SPACE.md`, and none is invented here.

**Not used: TOP3A.** `scan_literature.md` raises TOP3A as a second alternate on de
Nonneville et al. (https://doi.org/10.15252/emmm.202215859), where TOP3A **amplification**
drives ALT in pediatric osteosarcoma and is mutually exclusive with ATRX inactivation. It is
excluded because the published claim is a **copy-number** claim, and gene-level copy number
is the measurement class iteration 02 found unreliable in this cohort without segment
validation. Recorded so the exclusion is a decision, not an oversight.

## 8. Verdict mapping and campaign accounting

`evaluate_iteration(delta=rho, p_value=p_perm, predicted_sign=-1, conditions=...)`. No new
thresholds, no iteration-local evaluator.

Support conditions: `exposure_variance_gate` (§6 item 1), `negative_control_null` (S1 shows
no relationship), `not_screen_confounded` (S2), `loo_no_sign_flip` (S5).

**The campaign criterion, restated because it is the one that matters.** A result counts as
positive only at **q < 0.1 across every iteration's primary endpoint**, not p < 0.05 within
this iteration. With three iterations carrying a primary, that requires **p < 0.033**, hence
|rho| >= 0.535 per §5. Any positive is reported with the number of tests that preceded it —
79 before this iteration begins — and labelled exploratory.

## 9. What this iteration cannot do

It cannot measure BLM activity, only BLM transcript (§2). It cannot establish causation
from 16 observational cell lines. It cannot rescue Stage 1 or iteration 03. And at n = 16 it
cannot distinguish a moderate true effect from no effect — §5 is explicit that the likely
outcome is inconclusive, and that outcome will not be dressed up as either support or
refutation.
