# FINAL.md — the autoresearch campaign, closed

**Closed 2026-10-12 on stop condition (c): out of literature-grounded candidates.**
Four iterations ran. No fifth iteration exists, and the screen that reached that
conclusion is `autoresearch/candidate_screen_iteration05.md` — thirteen candidates, the
filter each failed, and the numbers that killed it.

Stop condition (b) was one iteration away and would have fired on the same evidence: the
streak of non-positive iterations stood at 4 of 5 with no new mechanistic angle surfaced in
iteration 04. (c) fired first because the screen found nothing legitimate left to run, and
the loop's rule is explicit — *"Do not fabricate further hypotheses to keep the loop
alive."*

This file answers the three questions the mandate specifies, in order.

---

## (a) Did any iteration survive experiment-wide FDR?

**No. Nothing came close, and nothing was in the predicted direction.**

The formal criterion, stated verbatim as required:

> **The formal success criterion is q < 0.1 at the experiment level, not p < 0.05 at the
> iteration level.**

### Final campaign ledger

Benjamini-Hochberg over every iteration's primary endpoint
(`src.verdict.campaign_fdr`, recomputed at closure):

| # | Hypothesis | Exposure → endpoint | n | aligned δ | raw p | **final q** | survives q<0.1 | Verdict |
|---|---|---|---|---|---|---|---|---|
| 01 | TCAB1/WRAP53 loss predicts ALT *status* (Keegan 2026) | WRAP53 log2(TPM+1) → frozen ALT label | 18 v 12 | **−0.241** | 0.280 | **0.420** | **no** | `counter_signal` |
| 02 | SMARCAL1 dependency is ATRX-gated (Brown 2026) | re-derived ATRX status → SMARCAL1 Chronos | — | — | — | — | no test run | `not-estimable` |
| 03 | POLD3/RAD52/SMARCAL1 are not a module (Broderick 2023) | none → mean pairwise Spearman ρ over the triplet | 16 (k=3) | **−0.182** | 0.229 | **0.420** | **no** | `not-estimable` (§6 gate) / mechanical `null` |
| 04 | PRIMPOL stress predicts SMARCAL1 dependency (Taglialatela 2026) | PRIMPOL log2(TPM+1) → SMARCAL1 Chronos | 16 | **−0.000** | 1.000 | **1.000** | **no** | `null` |

`aligned δ = δ × predicted_sign`, so a positive number is "as predicted". **All three
estimates are negative.** Iteration 01's is a stable counter-signal — it points the
opposite way to the hypothesis and strengthens when the four pre-known label–exposure pairs
are dropped (δ = +0.418 raw direction, p = 0.077), surviving 30/30 leave-one-out and 12/12
leave-one-depositor-out refits.

**Campaign test count: 101.** Test suite at closure: **185 passed**.

Two mechanics worth stating because they are easy to misread:

- Adding iteration 04's primary **raised** q for iterations 01 and 03 from 0.280 to 0.420.
  That is the multiplicity rule working: each additional pre-registered primary makes the
  family-wide bar stricter for every member. It is not a degradation of those results.
- Iteration 02 contributes no q because it ran no test. It was killed at the feasibility
  gate before any outcome was read, which is the correct disposition, not a missing result.

### Baseline, for calibration

Stage 1 (`stage1-v1.0-null`, not part of this FDR family): Cliff's δ = +0.281
(95% CI −0.313 to +0.844), p = 0.382, opposite to prediction; minimum BH q across the six
genes 0.959. Under `src/verdict.py`'s exhaustive partition, Stage 1 is **`counter_signal`**,
not `indeterminate`. The campaign did not rescue Stage 1 and did not overturn it. It
explained a large part of *why* Stage 1 could not have worked — see (b) item 1.

---

## (b) What the pattern of nulls says about the hypothesis space

The four nulls are not four instances of the same thing. They separate cleanly into
**three genuinely ruled out**, **one real finding that is not a hypothesis test**, and
**five untestable-here**. Conflating those categories would be the single most misleading
thing this campaign could do.

### Ruled out on this cohort — do not retry as specified

1. **WRAP53/TCAB1 mRNA as a proxy for 17p13.1 functional status.** Not merely
   unsupported — unusable. WRAP53 log2(TPM+1) over the 34 RNA-backed osteosarcoma lines is
   one continuous mode from 2.041 to 5.888 plus a discrete 2-line spike at exactly 0.0, with
   no low-but-nonzero values. The only naturally separable low group has **n = 2**, below
   `stratify.MIN_CELL_SIZE = 4`, and the `WRAP53 == 0` cut **cannot reach significance at
   any configuration** — the most extreme attainable Fisher two-sided p is 0.497. This
   retrospectively explains a Stage 1 mechanic: the §4b minimum-cell guard fired and forced
   the median fallback because the exposure has no natural low group, not by bad luck.
   Independently, Jiang 2025 and Bergstrand 2020 show TCAB1 protein level is set
   post-translationally by the CCT/TRiC chaperonin, so the transcript was never the right
   measurement — and DepMap 26Q1 ships no proteomics.
2. **The Stage 1 BIR−HR composite endpoint.** Three deterministic facts, independent of any
   p-value, gate or verdict, from iteration 03:
   - Two of the three BIR genes have no usable dynamic range in these 16 lines. POLD3 is
     essential in **16/16** (median −1.666); RAD52 in **0/16** (median −0.010, most negative
     value −0.416). SMARCAL1 is the only member with mid-range spread.
   - The composite measured mostly its own comparator: the HR arm supplied **68.2%** of its
     variance, and **RAD52 contributed −2.7%**. The decomposition is exact.
   - The comparator arm was the coherent module (HR T = +0.414, p = 0.0039) and the BIR
     triplet was not (T = −0.182, no pair positive). **Stage 1 subtracted a module from a
     non-module.**
3. **PRIMPOL transcript as a predictor of SMARCAL1 dependency.** ρ = exactly +0.0000,
   p = 1.000, verified against scipy before write-up because an exact zero is a plausible
   bug (Σd² = 680 = n(n²−1)/6 exactly, no ties). Robust to within-line standardisation
   (+0.012). The association is absent, not weak. This is **not** a refutation of
   Taglialatela's model, and the prespec pre-committed to not claiming one: the detection
   floor is |ρ| ≥ 0.500, the exposure is a transcript rather than repriming activity, and
   the mechanism's stronger arm was never testable here — BLM, whose loss *fully* suppresses
   SMARCAL1-dependent growth defects, is 2 outliers plus 14 lines inside a 1.496-unit band
   in these lines (Shapiro p = 0.0059).

### A real finding that is not a hypothesis test

**Co-essentiality is estimable on these 16 lines.** BTR (BLM/TOP3A/RMI1/RMI2) reaches
T = +0.526, p = 0.0001 and the HR module T = +0.414, p = 0.0039, both clearing their
pre-registered floors — BTR while containing one near-saturated member and two members
essential in no line. **"Sixteen lines is too few" is therefore not available as an
explanation for any null in this campaign.** That is what makes the nulls informative
rather than merely underpowered, and it is the most useful single result the loop produced.

### Untestable here — the hypothesis is open, the cohort is the obstacle

Not evidence against any of these:

| Hypothesis | Missing ingredient |
|---|---|
| ALT branch identity determines which repair genes are essential (Broderick 2023) | Any ALT-activity readout. `alt_assay_evidence` is `unknown` for all 35 manifest rows; no C-circle, APB, TRAP or telomere-length column; `telomerase_status` constant within ALT+; TERT Chronos flat (0/16, SD 0.114). No candidate marker gene is bimodal over the 34 lines (gap/SD 0.22–0.90 for all 16 tested), so any cut is a researcher choice nothing in the cohort could falsify. |
| WRAP53 expression misclassifies 17p13.1 structural disruption | WGS SV calls. DepMap 26Q1 ships none. |
| SMARCAL1 LOF promotes ALT, so the Stage 1 sign was backwards (Liu 2023; Brosnan-Cashman 2021; Sobieszek 2026) | A usable SMARCAL1 exposure. Its mRNA is barred by Bennett's documented transcript↔protein dissociation for SMARCAL1 itself; its CN is confounded with the CRISPR copy-number artefact when the outcome is the same gene; and the direct ALT+/ALT− dependency contrast is **contaminated** — already read in iteration 03's scan. |
| WRAP53β impairs HR and NHEJ non-selectively (Henriksson 2014) | Dependency variation in WRAP53 (1 of 16 lines below −0.5) — and the claim is about perturbation, not co-dependency. |
| BIR and HR arms are coupled with opposite sign (Batenburg 2023) | Dependency variation in the restart arm. LIG3 (0/16 below −0.5) and PRIMPOL (0/16) are two of Batenburg's three named restart genes and both sit in the "not dependency signal" class; only POLQ varies. The statistic would correlate one real module against a two-thirds-noise set, where the expectation is ~0 whether or not the claim is true. |
| ALT is not binary; TMM-DN and TMM-DP lines exist (Wu 2025) | Would require re-deriving the frozen ALT labels, which the loop forbids, or Wu's external 976-line resource. **This is a live limitation of every result above, including Stage 1's.** |

### The standing observation the campaign could not explain

**Within-ALT+ SMARCAL1 dependency genuinely varies in this cohort — 8 of 16 lines below
−0.5, 3 below −1.0, median −0.462, SD 0.395 — and four tested axes do not explain it:**

| Axis | Outcome |
|---|---|
| ATRX status | not estimable here (best rule 4v12 with 3 of 4 calls unsupported by segment data); and a **published pan-cancer null**, Bennett 2026, d = −0.248, P = 0.5, n = 20 v 14 |
| BLM expression | no usable variance — 2 outliers, 14 lines inside 1.496 units |
| PRIMPOL expression | **exactly zero** monotone association, robust to standardisation |
| Module membership with POLD3/RAD52 | the triplet is not a module; POLD3 and RAD52 have no dynamic range |

This is a real open problem, not a local failure: Bennett et al. independently flag
within-ALT+ heterogeneity as mechanistically unexplained and offer only a single-line
account. But it is not one this cohort's measurable variables appear able to solve.

One descriptive observation belongs here, labelled as what it is. The ALT+ versus ALT−
difference in SMARCAL1 dependency is large in these lines — ALT+ mean **−0.4907** (SD 0.395,
8/16 below −0.5) versus ALT− mean **+0.0557** (SD 0.137, 0/9). **This is not a campaign
result.** It was never pre-registered; the values were read during iteration 03's
feasibility scan, which contaminates any primary built on them; and it is Taglialatela's
and Bennett's already-published headline. It is reported as corroboration that this cohort
behaves as the literature says ALT+ cells behave, nothing more.

### What the campaign's own method produced

Three process findings, all found **by the data**, recorded rather than retroactively
rewritten. They are the transferable part:

1. **Iteration 03's assay-sensitivity gate was mis-specified.** It named a single positive
   control (Pol δ, chosen because it contains POLD3); the control failed while BTR and HR
   passed decisively, proving the design *was* sensitive. The frozen rule was applied anyway
   — verdict of record `not-estimable` — and the identical test was deliberately **not**
   re-run under a corrected label, because T and p would be unchanged and only the label
   would differ. A correct gate names several control complexes and passes if any clears.
2. **Iteration 04's support conditions were written relative to the primary's magnitude**,
   so a primary of exactly zero failed three of four arithmetically while carrying no
   information — the negative control's own p was 0.540, i.e. behaving exactly as intended.
   Write support conditions as absolute tests.
3. **Gene-level copy number from `OmicsCNGeneWGS.csv` must be segment-validated.**
   Gene-level ATRX values disagree with `OmicsCNSegmentsWGS.csv` in 3 of 4 apparent
   homozygous deletions. A naive re-derivation would have produced a 4v12 ATRX-gated
   analysis with three of four exposed lines misclassified, undetectably. This is the third
   defect of this kind in the project, after Stage 1's `load_wrap53_cn` layout bug and the
   manifest's artefactual `sv_architecture_auto`.

---

## (c) The honest next step

**Write up the preprint with an expanded discussion, and do not continue searching this
cohort.** Of the three options the mandate names, this is the one the evidence supports.

Why not the alternatives:

- **Continuing the search here is not an option the evidence leaves open.** Thirteen
  literature-grounded candidates were screened; none passes the four filters. The binding
  constraint is not imagination and not sample size — iteration 03 settled the sample-size
  question — it is that every remaining hypothesis needs a variable this cohort does not
  contain: an ALT-activity readout, an SV call, a protein level, or a branch label. Running
  a fifth iteration anyway would mean inventing a proxy, and the clearest single lesson of
  Stage 1 and of this campaign is that invented proxies are how the project went wrong in
  the first place.
- **A different cohort is where the live hypotheses actually point**, but it is a new
  project, not a continuation. The two best-motivated moves are a pan-lineage ALT+ cohort
  with real ALT-activity calls (Bennett's ALTitude set roughly doubles the known ALT+
  models) for the branch-identity and SMARCAL1-polarity questions, and a cohort with paired
  WGS SV calls and RNA at 17p13.1 for the proxy-validity question. Both require relaxing
  this loop's frozen-label and fixed-data rules, which is the user's scope decision and not
  the loop's.
- **Stopping the project is not warranted.** The campaign produced a defensible,
  publishable negative characterisation plus three methodological results. That is a
  smaller claim than the one Stage 1 set out to make, and it is a real one.

### What the preprint can now say that it could not before

1. Stage 1's null is **explained, not merely reported**. The exposure had no natural low
   group, the endpoint subtracted a coherent module from a non-module, and two of its three
   BIR genes carry no dependency signal in this cohort. The design could not have detected
   the effect it was built to detect.
2. The detection claim is **calibrated**. Co-essentiality is demonstrably estimable at
   n = 16 in these very lines (BTR p = 0.0001), so the nulls are informative within their
   stated floors rather than uninformative for want of power.
3. The literature constraints are **documented against the design**: eight prespec
   conflicts, each sourced to a specific paper, with the three that are fatal to the
   original framing identified (SMARCAL1's polarity is inverted; WRAP53β's repair defect is
   non-selective; the BIR and HR arms are coupled rather than orthogonal).
4. The within-ALT+ SMARCAL1 heterogeneity is posed as **an open problem with four axes
   excluded**, which is a contribution in its own right and converges with Bennett's
   independent observation.
5. The loop's outputs belong in **supplementary characterisation**, not in the results as
   positive findings. Every iteration is pre-registered, the campaign-level FDR is reported,
   and no q is below 0.420.

### What must not happen next

- **No reinterpretation of these nulls as support.** All three aligned estimates are
  negative; the campaign has no positive at any level.
- **No re-running of a closed hypothesis under a new label.** `NULL_SPACE.md` names each
  one and why.
- **No new ALT labels derived from this cohort**, and no proxy invented for a variable the
  data does not contain.
- **The preprint remains manual.** This loop does not write it.

---

## Where everything is

| What | Where |
|---|---|
| Standing rules the loop ran under | `autoresearch/PROGRAM.md` |
| Append-only per-iteration record | `autoresearch/LEDGER.md` |
| Campaign standing, ranking, FDR | `autoresearch/SUMMARY.md` |
| Every hypothesis tried and why it is closed | `autoresearch/NULL_SPACE.md` |
| Why there is no iteration 05 | `autoresearch/candidate_screen_iteration05.md` |
| Per-iteration prespec, notebook, data, memo | `autoresearch/iteration_0N_*/` |
| Contamination declarations | `autoresearch/iteration_0{1,4}_*/DISCLOSURES.md` |
| Stage 1, untouched and reproducible | tag `stage1-v1.0-null` |

---

# ERRATUM, 2026-10-12 — the conflict-row-4 kill was overstated

Added the same day, after an exploratory check prompted by a question from the user. Written
as a visible erratum rather than a rewrite, following the Stage 1 v1.2 and iteration-04
draft precedents. **The campaign's four iterations and every q-value above are unaffected.**
What changes is the disposition of one *unrun* candidate.

## What I claimed

In `candidate_screen_iteration05.md` and in (b) above, Batenburg's cross-arm antagonism
(Prespec conflict row 4) was closed as `not-estimable` on the argument that LIG3 (0/16
lines below −0.5) and PRIMPOL (0/16) carry "no dependency signal", so a cross-set statistic
would correlate one real module against "two-thirds screen noise" and have **expectation
≈ 0 whether or not the claim is true.**

## Why that is wrong

The campaign's own positive control refutes it. BTR — the module that reached T = +0.526,
p = 0.0001 — contains two members in exactly that regime, and they are not noise:

```
BLM    median -0.0150  SD 0.1649   0/16 below -0.5
RMI2   median -0.2145  SD 0.1707   0/16 below -0.5
TOP3A  median -1.0780  SD 0.2196  16/16 below -0.5
RMI1   median -0.4045  SD 0.3703   4/16 below -0.5

BLM   x TOP3A  rho +0.6824  p 0.0036   <- non-essential member, strongly coupled
RMI1  x RMI2   rho +0.8000  p 0.0002   <- non-essential member, strongly coupled
BLM   x RMI2   rho +0.4441  p 0.0848   <- BOTH non-essential in all 16 lines
TOP3A x RMI2   rho +0.4559  p 0.0759
BLM   x RMI1   rho +0.3824  p 0.1439
TOP3A x RMI1   rho +0.3882  p 0.1373
```

**Sub-threshold dependency variation is recoverable biology in these 16 lines, not screen
noise.** A gene essential in no line can still have a graded profile that tracks its
obligate complex partners at ρ +0.68. The classification I leaned on — iteration 03 §7.3,
"spread that is not dependency signal" — describes whether a line *depends* on the gene. It
does not license the inference that the gene's across-line ranking is uninformative, and I
used it as though it did.

## What the correct disposition is

**Conflict row 4 is `untested`, not `not-estimable`.** It should move out of the
not-estimable class in `NULL_SPACE.md`. Two of the three reasons I gave do survive, and
they are reasons for care rather than closure:

1. **BRCA2 is still barred** (§6, plus the contamination declared in iteration 04's
   DISCLOSURES), so the HR arm is BRCA1 + RAD51 and the statistic is a 2 × 3 cross-set
   mean, whose null must be derived rather than reused from `src/coherence.py`'s within-set
   floors.
2. **A negative predicted sign is a harder ask than positive coherence**, and the
   construct-validity worry is real even though the estimability worry was not: if no line
   depends on LIG3 or PRIMPOL, then "this line relies on the alternative-restart route" is
   not something the Chronos score directly reports, whatever the correlation comes out at.
   That needs stating in a prespec, not assuming away — in either direction.

What does **not** survive is the claim that the test is uninformative by construction. It
is not, and closing it on that basis was wrong.

## Consequence for section (c)

Section (c)'s recommendation — preprint with an expanded discussion rather than continued
search — rested on *thirteen of thirteen* candidates failing. It is now **twelve of
thirteen**, with one live, literature-grounded, estimable candidate on data already held.
That is a material change and the recommendation should be read accordingly: the preprint
is still the right primary move, but "there is nothing left to test on this cohort" is no
longer an accurate statement and must not be written into the preprint. Separately, the
untouched Pol ζ axis (Wu 2023, INDEX #16) has usable dynamic range here — REV3L median
−0.578, SD 0.367, 9/16 below −0.5 — and failed only the sourcing filter, not a
feasibility one.

Both are new mandates, not resumptions: `PROGRAM.md` is closed and a reopened loop would
need its own prespec, written before any of the values above are used as an outcome.
