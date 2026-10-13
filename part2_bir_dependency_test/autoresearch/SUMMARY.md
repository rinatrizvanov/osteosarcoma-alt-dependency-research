# SUMMARY.md — autoresearch campaign standing

**Status: CLOSED 2026-10-12 on stop condition (c) — out of literature-grounded candidates.** 4 of at most 10 iterations ran; there is no iteration 05. **No campaign positive**: final BH q = 0.420, 0.420, 1.000 over the three primaries. The conclusions, the three mandated answers and the honest next step are in `FINAL.md`; the screen that closed the loop is in `candidate_screen_iteration05.md`.

## Formal success criterion

> A result counts as positive only if its pre-specified success criterion
> fires, it is directionally consistent and robust, **and its primary endpoint
> survives Benjamini-Hochberg FDR at q < 0.1 computed across every iteration's
> primary endpoint in the campaign** — not merely p < 0.05 within its own
> iteration.

This is stated here verbatim because it is the criterion that matters. Each
iteration is individually pre-registered, but the hypotheses are chosen in
sequence against one fixed cohort of 16 ALT+ / 9 ALT− osteosarcoma lines, and
each is chosen partly *because* the previous one failed. Per-iteration p-values
therefore understate the campaign error rate. With ten iterations at roughly
seven tests each, three to four nominally significant results are expected
under a complete null. Any positive is reported with the number of tests that
preceded it and labelled exploratory.

## Ranking

| Rank | Iteration | Hypothesis | Verdict | aligned δ | raw p | campaign q | survives q<0.1 |
|---|---|---|---|---|---|---|---|
| 1 | 03 | POLD3/RAD52/SMARCAL1 form a coherent module | `not-estimable` (§6 gate) / mechanical `null` | −0.182 | 0.2285 | 0.420 | **no** |
| 2 | 01 | TCAB1/WRAP53 mRNA predicts ALT status | `counter_signal` | −0.241 | 0.280 | 0.420 | **no** |
| 3 | 04 | PRIMPOL expression predicts SMARCAL1 dependency | `null` | −0.000 | 1.000 | 1.000 | **no** |
| — | 02 | ATRX-gated SMARCAL1 dependency | `not-estimable` | — | — | — | no test run |

Campaign BH is recomputed over all primaries at completion. Across the three iterations
that ran a primary, BH gives **q = 0.420, 0.420 and 1.000**. Adding iteration 04's
primary *raised* q for 01 and 03 from 0.280 to 0.420 — the multiplicity rule working as
intended, since each additional pre-registered primary makes the family-wide bar stricter
for every member. Statistical tests run so far: **101**.

## Best-supported model found

**None.** No iteration has produced a result in its predicted direction at any level of
significance, and nothing approaches the campaign criterion.

## What iteration 03 established about the Stage 1 endpoint

Not a model, but the most substantive result in the campaign so far, and it holds
**deterministically** — independent of any p-value, gate or verdict:

1. **Two of the three Stage 1 BIR genes have no usable dynamic range** in these 16 lines.
   POLD3 is essential in **16/16** (median −1.666), RAD52 in **0/16** (median −0.010), and
   SMARCAL1 is the only member with mid-range spread (median −0.462, 8/16).
2. **The composite measured mostly its own comparator.** The three HR genes supplied
   **68.2%** of the composite's variance against the BIR arm's 31.8%; SMARCAL1 alone was
   23.4%, and **RAD52 contributed −2.7%**. The decomposition is exact.
3. **The comparator arm was the coherent one.** HR (BRCA1/BRCA2/RAD51) is a real module
   here, T = +0.414, p = 0.0039. The BIR triplet is not. Stage 1 subtracted a module from a
   non-module.
4. **The cohort is not too small for this kind of analysis.** BTR (BLM/TOP3A/RMI1/RMI2)
   reaches T = +0.526, p = 0.0001 at n = 16 — while containing both a near-saturated member
   and two members essential in no line. "16 lines is too few" is therefore not available
   as an excuse for any future null here.

Iteration 03's own `not-estimable` verdict of record follows from its frozen §6
assay-sensitivity gate, which fired because the single designated positive control (Pol
delta) failed. The memo records that the gate was mis-specified — a correct gate names
several control complexes and passes if any clears — and deliberately does **not** re-run
the identical test under a corrected label.

## What the campaign has ruled out about SMARCAL1

Within-ALT+ SMARCAL1 dependency genuinely varies in this cohort — 8 of 16 lines below −0.5,
3 below −1.0 — and the campaign has now failed to explain that variation from **every axis
the data can measure**:

| Axis | Outcome | Where |
|---|---|---|
| ATRX status | not estimable here; and **published null** pan-cancer (d = −0.248, P = 0.5, n = 20 v 14) | iteration 02; Bennett 2026 |
| BLM expression | no usable variance — 2 outliers, 14 lines in 1.496 units | iteration 04 §2 |
| PRIMPOL expression | **exactly zero** monotone association, robust to standardisation | iteration 04 |
| Module membership with POLD3/RAD52 | the triplet is not a module; POLD3 and RAD52 have no dynamic range | iteration 03 |

Bennett et al. themselves flag within-ALT+ heterogeneity as mechanistically unexplained, so
this is a real open problem rather than a local failure — but it is not one this cohort's
measurable variables appear able to solve.

## Why the campaign closed

All thirteen remaining candidates — five `next_hypotheses.md` rows and the eight Prespec
Conflicts rows — failed at least one of the four admission filters. The binding constraint
is **not** sample size: iteration 03 proved co-essentiality is estimable at n = 16 here
(BTR T = +0.526, p = 0.0001). It is that every remaining literature-grounded hypothesis
needs either a variable this cohort does not contain — an ALT-activity readout, a WGS SV
call, a protein level, a branch label — or a gene whose dependency does not vary across
these 16 lines. The two strongest survivors died on measurement: Batenburg's cross-arm
antagonism (LIG3 0/16 and PRIMPOL 0/16 below −0.5, so the restart arm is two-thirds screen
noise) and Jiang's CCT/TRiC route (chaperonin essential in 15–16 of 16 lines, and no
proteomics exists). Per the loop's own rule, no fifth hypothesis was fabricated to keep it
alive.

## Candidates carried forward to a different cohort, not a later iteration

**(a) An HR/BRCA2-centred endpoint — WITHDRAWN as contaminated.** Iteration 03 nominated
this: the HR arm is the coherent module in these lines and BRCA2 its largest variance
contributor (25.8%), meeting Lee et al. 2024 (https://doi.org/10.1093/nar/gkae251) where
BRCA2 loss drives ALT/BIR. **It cannot be tested blind.** Iteration 03's essentiality scan
reported ALT+ and ALT− marginals separately for 38 genes, so the ALT+ vs ALT− contrast on
BRCA2 is already known to this session (means −0.537 vs −0.509, i.e. no difference). Recorded
in `iteration_04_smarcal1_dependency_model/DISCLOSURES.md` and withdrawn there rather than
run as though blind.

**(b) TP53 expression, not WRAP53 expression, discriminates ALT status.** Carried from
iteration 01's R5 — δ = −0.472, p = 0.031, AUROC 0.736 at 18v12, ALT+ lines significantly
*lower* in TP53. Pre-registered there as descriptive and non-gating, so **not a result**: 1
of 8 robustness contrasts in a single iteration, and it has not faced campaign FDR.
Mechanistically coherent, because Keegan et al.'s SVs land in TP53 intron 1 and inactivate
TP53 outright, clipping the TCAB1 overlap only incidentally. Note the prior is not high:
TP53 disruption is near-universal in osteosarcoma, so this may be an expression artefact of
the SV architecture rather than a discriminator.

## Baseline for comparison

Stage 1 (`stage1-v1.0-null`), for calibration — not part of this campaign's
FDR family:

- Primary 8v8 median split on WRAP53 log2(TPM+1) within 16 ALT+ lines.
- Endpoint mean(POLD3, RAD52, SMARCAL1) − mean(BRCA1, BRCA2, RAD51), Chronos.
- Cliff's δ = +0.281 (95% CI −0.313 to +0.844), p = 0.382 — opposite to
  prediction. Minimum BH q across the six genes 0.959.
- Under the Stage 1 evaluator: `indeterminate`. Under the exhaustive partition
  in `src/verdict.py`: **`counter_signal`**.
