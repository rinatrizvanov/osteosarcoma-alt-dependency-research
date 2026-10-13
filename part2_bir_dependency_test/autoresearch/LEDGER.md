# LEDGER.md — append-only record of the autoresearch campaign

One row per iteration. Never edit a committed row; corrections go in a new
row with a note. `aligned delta` is `delta * predicted_sign`, so positive is
always "as predicted". Raw p is the per-iteration primary p; the
experiment-level q is computed across all primaries at completion
(`src.verdict.campaign_fdr`) and lives in `SUMMARY.md`.

**Campaign test count: 101** — every statistical test run in any iteration,
cumulative. Anti-bias rule (a) depends on this number being honest.

**CAMPAIGN CLOSED 2026-10-12 on stop condition (c) — out of literature-grounded candidates.** Four iterations ran; there is no iteration 05. See `FINAL.md` and `candidate_screen_iteration05.md`.

**Consecutive non-positive streak: 4** — stop condition (b) fires at 5 and would have fired on the same evidence (no new mechanistic angle surfaced in iteration 04). (c) fired first because the screen of all thirteen remaining candidates found nothing legitimate left to run.

| # | Slug | Hypothesis (one line) | Exposure | Endpoint | n per arm | aligned δ | raw p | Verdict | Tests |
|---|---|---|---|---|---|---|---|---|---|
| 01 | `iteration_01_tcab1_predicts_alt_status` | TCAB1/WRAP53 loss predicts ALT *status*, not dependency within ALT (Keegan 2026) | WRAP53 log2(TPM+1), continuous | frozen ALT label, ALT+ vs ALT− | 18 v 12 | **−0.241** | 0.280 | `counter_signal` | 51 |
| 02 | `iteration_02_atrx_gated_smarcal1` | SMARCAL1 dependency in ALT+ OS is gated by ATRX, not TCAB1 (Brown 2026) | re-derived ATRX status | SMARCAL1 Chronos gene effect | — | — | — | `not-estimable` | 0 |
| 03 | `iteration_03_bir_module_coherence` | POLD3, RAD52 and SMARCAL1 do not behave as a coherent module; a post-mortem on the Stage 1 endpoint (Broderick 2023) | none — no exposure variable; cross-gene correlation structure | mean pairwise Spearman rho over the triplet, Chronos | 16 lines (k=3) | **−0.182** | 0.229 | `not-estimable` (§6 gate) / mechanical `null` | 79 |
| 04 | `iteration_04_smarcal1_dependency_model` | PRIMPOL-generated replication stress predicts SMARCAL1 dependency within ALT+ (Taglialatela 2026) | PRIMPOL log2(TPM+1), continuous (switched from BLM pre-freeze) | SMARCAL1 Chronos gene effect | 16 lines | **−0.000** | 1.000 | `null` | 101 |

## Harness

| Item | Status |
|---|---|
| `src/verdict.py` exhaustive partition | built; 48 tests |
| Exact-null `min_detectable_delta` | validated against all 6 frozen Stage 1 entries |
| `benjamini_hochberg` | validated against `scipy.stats.false_discovery_control` |
| `src/criteria.py` | **unmodified** — Stage 1 stays reproducible at `stage1-v1.0-null` |
| Full suite | **185 passed** at closure (127 + 58 from iteration 03's `src/coherence.py`) |

## Notes per iteration

**01.** Prespec v1.0.1 frozen `8e6ed51` (+ erratum `2d68239`) before any
exposure–outcome join; notebook idempotent. Observed δ = **+0.241** on the raw
scale (ALT+ lines *higher* in WRAP53), so aligned δ = −0.241 against a
predicted-negative direction. Counter-direction held in 30/30 leave-one-out and
12/12 leave-one-depositor-out refits, and **strengthened** under both
known-information guards (R3 +0.396 p=0.082; R3b +0.418 p=0.077). Below the
0.435 detection floor, so suggestive rather than established.

New mechanistic angle surfaced: **R5 TP53 expression** gives δ = −0.472,
p = 0.031, AUROC 0.736 — pre-registered as descriptive and non-gating, not a
result of this iteration, carried to the candidate list. Stop condition (b)
therefore does not advance on the "no new mechanistic angle" clause.

**02.** Killed at the pre-specified feasibility gate in ~15 min; no prespec
frozen, **no Chronos value read**, 0 tests entering campaign BH. Four candidate
ATRX rules give n ∈ {1, 2, 4} altered of 16 ALT+ lines, and the only one
clearing `MIN_CELL_SIZE = 4` (gene-level CN) is refuted for 3 of its 4 calls by
`OmicsCNSegmentsWGS.csv`. The two rules giving n = 1 disagree on which line
(manifest point-mutation → NRHOS1; segment-validated focal deletion → CAL72).
MDE at the unreliable 4v12 would be 0.708 regardless. Counted toward the
stop-condition (b) streak, conservatively, although a feasibility failure is
neither `null` nor `inconclusive_underpowered`.

**03.** Prespec frozen `b24a8ed` before any cross-gene quantity was computed; statistic and
its pre-registered detection floor committed earlier at `cf73508` in a module with no path
to the outcome (asserted against the source text in the test suite). Notebook idempotent —
manifest byte-identical on re-run.

Primary T = **−0.182**, p = 0.229 at n = 16 against a floor of 0.289, so mechanically
`null`. **No pair is positive** (POLD3–RAD52 −0.268, POLD3–SMARCAL1 +0.024,
RAD52–SMARCAL1 −0.303) and the triplet sits at the 12.5th percentile of a mean/SD-matched
random-triplet null — *less* coherent than arbitrary genes with the same marginals.

**The §6 assay-sensitivity gate fired**, so the verdict of record is `not-estimable`
(design-insensitive) and the primary may not be read as evidence against module coherence.
Pol delta reached −0.164 raw / −0.196 standardised against a +0.40 threshold.

**The gate was mis-specified, and the data show it.** BTR (BLM/TOP3A/RMI1/RMI2) gives
T = **+0.526, p = 0.0001**, and the HR module T = **+0.414, p = 0.0039** — both clear their
floors, and BTR does so while containing one near-saturated gene (TOP3A, 16/16 below −0.5)
and two essential in no line (BLM, RMI2, 0/16). So the method is demonstrably sensitive at
n = 16, and "16 lines is too few for co-dependency analysis" is not available as an excuse
for any future null on this cohort. The gate's single-designated-control design conflated
"method insensitive" with "POLD3 carries no signal"; a correct gate names several complexes
and passes if any clears. Recorded as a defect found by the data, not fixed silently, and
**not** re-run under a new number, since T and p would be identical and only the label
would change.

Established regardless of the gate, deterministically: two of the three Stage 1 BIR genes
have no usable dynamic range (POLD3 essential in 16/16, RAD52 in 0/16, SMARCAL1 the only
mid-range member); the Stage 1 composite drew **68.2%** of its variance from the HR
comparator arm, with RAD52 contributing **−2.7%**; and the comparator arm was the coherent
module while the BIR arm was not. Stage 1 subtracted a real module from a non-module.

Secondaries: S1b corroborates Broderick's prediction in sign (rho(RAD52,SMARCAL1) = −0.303
<= 0) but not significantly (k=2 floor 0.494). S4 Batenburg's own MUS81/RAD52/POLD3 module
is also anti-coherent (−0.274, MUS81–RAD52 = −0.512). S2b standardisation changes nothing
(−0.188). S3 the ALT− arm runs the *opposite* sign (+0.378, p = 0.058, n = 9). S5 invokes
rule (c)'s explicit-statement branch: no defensible type I/II marker exists and none was
invented. **S6 was vacuous** — the depositor regex matched 0 of 16 lines so the "refit"
reproduced the primary; counted out of the test total and reported as a failed check.

Test count 28, not the 29 enumerated in the prespec, for that reason. New mechanistic angle
surfaced (HR/BRCA2-centred endpoint), so the stop-condition (b) "no new angle" clause does
not advance.

**04.** Prespec frozen `5f2f039`, analysis code committed `6eb0908` before execution;
notebook idempotent, manifest byte-identical on re-run.

Primary Spearman rho(PRIMPOL expression, SMARCAL1 Chronos) = **exactly +0.0000**, p = 1.000
at n = 16 against a floor of 0.494. Verified against scipy before write-up because an exact
zero is a plausible bug: `spearmanr` returns 0.0, Pearson-on-ranks 6.9e-18, and sum d² = 680
equals n(n²−1)/6 = 680 exactly, with no ties. Robust to within-line standardisation
(+0.012). Post-hoc bootstrap interval (not pre-registered, reported because a null without
an interval is uninformative): [−0.502, +0.522].

**Primary exposure was switched BLM → PRIMPOL before the freeze**, on exposure-side
distributions only. BLM — the mechanism's *stronger* arm, full suppression of
SMARCAL1-dependent growth defects in Taglialatela — is 2 low outliers plus 14 lines inside a
1.496-unit band in these 16 lines (Shapiro p = 0.0059), and its only natural low group has
n = 2, below the minimum-cell guard. The prespec recorded the swap as bad news and lowered
the stated prior accordingly.

**Three of four support conditions failed vacuously and must not be read as evidence.**
Each was written as a comparison against the primary's magnitude, so a primary of exactly
zero fails them arithmetically: any non-zero control beats zero, and the sign of a zero
estimate is arbitrary so leave-one-out straddles it by construction. The negative control's
own p is 0.540 — behaving exactly as a negative control should. Second prespec defect found
by the data in two iterations, after iteration 03's single-control assay gate; recorded, not
rewritten. ZRANB3 also failed the exposure gate itself (SD 0.427 < 0.50), a weakness the
prespec had pre-declared.

Campaign BH now spans three primaries and **raises q for iterations 01 and 03 from 0.280 to
0.420** — the multiplicity rule working as intended. Iteration 04 q = 1.000.

Together with iteration 02 and Bennett's published pan-cancer ATRX null, within-ALT+
SMARCAL1 dependency has now resisted explanation from **every axis this cohort can
measure**: ATRX (not estimable, plus published null), BLM (no usable variance), PRIMPOL
(exact zero).


## Closure — no iteration 05

**Stop condition (c) fired: out of literature-grounded candidates.** Recorded here rather
than as a table row because no hypothesis was run and no test was added; the campaign test
count stays at **101**.

Thirteen candidates were screened against the four admission filters — sourced from the
Prespec Conflicts table or an untried `next_hypotheses.md` row; testable on DepMap 26Q1 as
shipped; backed by a specific indexed paper; not a variant of anything in `NULL_SPACE.md`.
**Zero passed all four.** Full screen in `candidate_screen_iteration05.md`.

One new computation was made for the screen, deliberately marginals-only with no
correlation of any kind: Chronos dependency regimes over the 16 ALT+ lines for the gene
panels the surviving candidates would have needed. It killed the two strongest survivors on
measurement rather than on taste:

- **Batenburg's cross-arm antagonism (conflict row 4)** — the restart arm has no dependency
  variation. LIG3 0/16 and PRIMPOL 0/16 below −0.5; only POLQ (5/16) is usable. Correlating
  one coherent module against a two-thirds-noise set has expectation ≈ 0 regardless of the
  claim's truth.
- **CCT/TRiC control of TCAB1 protein (conflict row 5)** — CCT5 16/16, CCT7 15/16,
  CCT2 16/16, TCP1 16/16 below −0.5. Housekeeping chaperonin dependency, plus the claim is
  about protein and there is no proteomics.

Stop condition (b) was one iteration from firing and rests on the same evidence. (c) is the
more accurate label: the loop did not merely fail five times, it ran out of hypotheses this
cohort can carry.

`FINAL.md` answers the three mandated questions. Final campaign BH over the three primaries
that ran: **q = 0.420, 0.420, 1.000**. No positive.
