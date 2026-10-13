# Draft outline — Stage 1 plus the autoresearch campaign

**Scaffolding only.** Section headers and 2–4 sentence summaries. **All prose is the
author's to write.** This file was produced at the user's instruction after the campaign
closed (`autoresearch-complete-v1.0`); it contains no new analysis and every number is
quoted from a committed artifact.

**Relationship to `docs/preprint_outline.md`.** That file is the Stage 1 short-report
outline written at Stage 1 closeout and it is **not superseded** — §§1–2 and 4 below still
follow it closely. What changes is scope: the campaign supplies a mechanistic explanation
for the Stage 1 null, a calibration that makes the null informative, and an expanded
Discussion. Where the two files disagree, the older one is the record of what was planned
at Stage 1 and this one is the plan for the paper as it now stands.

**Target.** bioRxiv (Cancer Biology / Genomics). Short report, ~3,000–3,500 words — longer
than the Stage 1 plan because §3.3–3.5 are new. 3 main figures, 1 main table, supplementary
pre-specifications.

**Framing, unchanged from Stage 1 and now load-bearing.** The contribution is the
pre-registration and the characterisation, not a discovery. Every section should read as
though the result had come out the other way. The campaign strengthens this rather than
diluting it: five frozen prespecs, 101 tests, campaign-wide FDR, and no positive — reported
as such.

**One framing decision the author must make before drafting.** Either (a) Stage 1 is the
paper and the campaign is supplementary characterisation plus Discussion, or (b) the paper
is about the campaign's method and Stage 1 is its first iteration. This outline assumes
**(a)**, which is the more defensible claim and matches `FINAL.md`'s recommendation. Option
(b) needs a different Introduction and should not be assembled from these sections.

---

## Title and abstract

Working title: *Pre-registered non-confirmation of a TCAB1-gated break-induced replication
dependency in ALT+ osteosarcoma, and why the endpoint could not have detected it.*
Abstract draft in `preprint/abstract_draft.md` (~150 words).

## 1. Introduction

**1.1 ALT and the telomerase-independent escape route.** ALT in osteosarcoma; SMARCAL1 as
an established broad ALT dependency (Taglialatela 2026; Bennett 2026) not previously
stratified by upstream genotype.

**1.2 The 17p13.1 / TCAB1 prediction.** Keegan 2026: SVs in TP53 intron 1 inactivate TP53
and TCAB1 together via the 1,419-nt head-to-head overlap. With telomerase blocked and
WRAP53β lost, the cell should lean harder on BIR. State the prediction explicitly.

**1.3 Why pre-register this one.** Attractive prediction, small N, large analytic degrees
of freedom. Plan, success criteria *and* failure criteria committed to version control
before any outcome value was read.

**1.4 New — what this paper adds beyond the null.** Flag early that the paper also reports
why the endpoint was unusable and a calibration showing the cohort is not too small. A
reader should know by the end of the Introduction that this is not a bare null.

## 2. Results — the pre-registered test

**2.1 Cohort and stratification.** 35 osteosarcoma lines → 26 CRISPR-screened → 16 ALT+,
9 ALT−, 1 excluded. The bottom-tertile rule placed 3 of 16 in the low arm, triggering the
pre-registered minimum-cell guard; the primary became a median split with the tertile kept
as a parallel sensitivity. The fallback was written before the WRAP53 distribution was
examined.

**2.2 The stratifier is a valid measurement.** The strand check: WRAP53 exon 1 overlaps
TP53 antisense by ~1,419 bp, DepMap does not document strandedness, and a pre-committed
rule was resolved empirically — WRAP53 tracks its own copy number (ρ = +0.554) more
strongly than TP53 (ρ = +0.244). Belongs in Results; it is a result about the data.

**2.3 Primary comparison: no support, wrong direction.** **Figure 1.** δ = +0.281,
p = 0.382; tertile sensitivity δ = +0.538, p = 0.189. The point estimate is on the opposite
side of zero from the prediction.

**2.4 No gene carries a signal.** **Figure 2.** All six BH q = 0.959; every 95% interval
spans zero. Note honestly that RAD52 (−0.344) and BRCA2 (−0.250) lean the predicted way and
that nothing is distinguishable from noise.

**2.5 Robustness.** **Table 1.** R1–R8. No leave-one-out refit of 16 reaches the predicted
direction; the two favourable checks are the two smallest comparisons, one of them the
negative control — where a predicted-direction result argues *against* a TCAB1-specific
account.

**2.6 Verdict.** Frozen criteria return `indeterminate`; under the later exhaustive
partition, `counter_signal`. Criteria evaluated mechanically by code so the verdict could
not be shaded.

**2.7 Incidental: OS384 carries the predicted architecture.** Zero expression of both
WRAP53 and TP53 with WRAP53 CN ≈ 0 — the 17p13.1 co-deletion architecture in a public cell
line, recorded in the pre-specification before the outcome was read. Explicitly decoupled
from the dependency result.

## 3. Results — why the endpoint could not have worked (new)

This is the section that distinguishes the paper from a generic null. All three claims are
deterministic: they do not depend on any p-value, gate or verdict.

**3.1 Two of the three BIR genes have no dependency variation.** **Figure 3a.** POLD3
essential in 16/16 ALT+ lines (median −1.666), RAD52 in 0/16 (median −0.010), SMARCAL1 the
only mid-range member (median −0.462, 8/16 below −0.5). Three genes, three regimes.

**3.2 The composite measured mostly its own comparator.** **Figure 3b.** Exact variance
decomposition: HR arm 68.2%, BIR arm 31.8%, SMARCAL1 alone 23.4%, **RAD52 −2.7%**.

**3.3 Stage 1 subtracted a module from a non-module.** **Figure 3c.** Mean pairwise
Spearman over gene-effect profiles: HR (BRCA1/BRCA2/RAD51) T = +0.414, p = 0.0039 — a real
module. BIR triplet T = −0.182, p = 0.229, no pair positive. Batenburg's own MUS81/RAD52/
POLD3 is also anti-coherent (T = −0.274).

**3.4 The cohort is not too small — a calibration.** **Figure 3c**, same panel. BTR
(BLM/TOP3A/RMI1/RMI2) reaches T = +0.526, p = 0.0001 at n = 16, while containing one
near-saturated member and two members essential in no line. This is the paper's licence to
treat its nulls as informative within their stated floors. Emphasise that the statistic's
null depends only on (n, k), because it depends on values only through within-gene ranks —
so the detection floor was pre-registerable and no gene-effect value could move it.

**3.5 The exposure had no natural low group.** WRAP53 log2(TPM+1) over 34 RNA-backed lines:
one continuous mode 2.041–5.888 plus a discrete 2-line spike at exactly 0.0. The only
naturally separable low group has n = 2, and the `WRAP53 == 0` cut cannot reach significance
at any configuration (most extreme attainable Fisher two-sided p = 0.497). This
retrospectively explains why the §4b guard fired in 2.1 — the exposure has no low group,
not bad luck.

## 4. Methods

**4.1 Data and cohort.** DepMap 26Q1 specifically (the ALT labels were derived against it);
five files sha256-verified at load; inclusion filter applied in order with N logged at each
step; ground-truth anchors asserted before analysis.

**4.2 Endpoint and tests.** BIR−HR composite and its sign convention; two-sided exact
Mann–Whitney; Cliff's δ with 2,000-resample bootstrap CIs validated against
δ = 2U/(nm) − 1; Hodges–Lehmann shifts; BH-FDR across exactly the six genes.

**4.3 Module coherence (new).** Mean pairwise Spearman over gene-effect profiles; 10,000-
permutation null; the rank-invariance property and the resulting (n, k)-only null;
mean/variance-matched random-triplet calibration; within-line standardisation as a
screen-amplitude control.

**4.4 Pre-specification and enforcement.** Freeze commits; the amendment made before any
outcome was read; structural enforcement — one module loads gene effects, notebooks assert
it is absent from `sys.modules`, criteria evaluated by code. For the campaign: five frozen
prespecs, exposure-derivation isolation, and campaign-wide BH across all primaries.

**4.5 Power.** Exact minimum detectable |δ| per split, computed exhaustively over rank
configurations and committed in advance.

**4.6 Code and data availability.** Releases `stage1-v1.0-null` and
`autoresearch-complete-v1.0`; DepMap not redistributed, fetch policy given.

## 5. Discussion

**5.1 What the result excludes, and what it does not.** CI (−0.313 to +0.844) does not
exclude a moderate predicted-direction effect; §5a committed in advance that only
|δ| ≥ 0.594 was detectable at this N. A non-confirmation, not a refutation.

**5.2 Low power does not explain a wrong-direction estimate.** The primary, both tertile
variants, the continuous correlation and all 16 leave-one-out refits agree on the opposite
sign. Keep this short and hard.

**5.3 The construct, not the power, is the binding limitation (new).** Carry §3 into the
Discussion: the endpoint averaged a pan-essential gene, a never-essential gene and a
mid-range gene, then subtracted a coherent module. No pre-registration downstream of that
choice could have detected it. This is the paper's main methodological claim and it
generalises beyond ALT.

**5.4 The proxy is the second binding limitation.** DepMap 26Q1 ships no SV calls, so
WRAP53 expression stood in for 17p13.1 structural disruption throughout. Jiang 2025 and
Bergstrand 2020 make this worse than a power problem: TCAB1 protein level is set
post-translationally by CCT/TRiC, so the transcript was never the right measurement. A
cohort with real SV calls, or with protein, tests the prediction as stated.

**5.5 What the campaign ruled out, and what it could not test (new).** Four further
pre-registered iterations; campaign-wide q = 0.420, 0.420, 1.000; no positive. Keep the
three categories separate and say so explicitly: *ruled out on this cohort* (WRAP53 mRNA as
a proxy; the BIR−HR composite; PRIMPOL transcript as a predictor of SMARCAL1 dependency,
ρ = exactly 0.000) versus *untestable here for want of a variable the data does not
contain* (ALT branch identity — no ALT-activity readout of any kind; 17p13.1 architecture —
no SV calls; SMARCAL1's inverted polarity — no usable exposure). Conflating these would be
the most misleading thing the paper could do.

**5.6 An open problem: within-ALT+ SMARCAL1 heterogeneity (new).** SMARCAL1 dependency
genuinely varies among ALT+ lines — 8 of 16 below −0.5, 3 below −1.0 — and four tested axes
do not explain it: ATRX (not estimable here, and a published pan-cancer null, Bennett 2026,
d = −0.248, P = 0.5), BLM expression (no usable variance), PRIMPOL expression (exactly
zero), and module membership with POLD3/RAD52. Bennett et al. independently flag this
heterogeneity as unexplained, so this is a real open problem, posed with four axes excluded.

**5.7 ALT is probably not binary — a limitation on everything above (new).** Wu 2025
reports TMM-double-negative and TMM-double-positive lines, including mixtures of ALT and
telomerase subclones. Every result here, Stage 1's included, is conditioned on a binary ALT
label that the current literature does not fully support. State plainly that this was not
resolvable with the frozen labels and is the first thing a replication should change.

**5.8 What pre-registration bought, including where it bit.** Two criteria repaired before
the freeze on exact-null grounds (the selectivity rule would have blocked a true positive
20–30% of the time). Two mis-calibrated and **not** repaired afterwards: §10's clean-null
ceiling assumed a null in the predicted direction, so a wrong-direction estimate fell into
`indeterminate`; and the campaign's assay-sensitivity gate named a single positive control
that failed while two others passed decisively. Both left frozen deliberately. Design
lessons, reported as such.

**5.9 Reporting nulls in dependency screens.** DepMap makes hypothesis-driven dependency
tests cheap, which makes selective reporting cheap too. A pre-registered null with published
criteria, code and campaign-wide FDR is the low-cost correction. Note that the campaign's
own closing judgement required a published erratum
(`autoresearch/FINAL.md`, 2026-10-12) — the loop wrongly declared one hypothesis
not-estimable, and the correction is in the record. That belongs in the paper: it is the
clearest available evidence that the discipline was real rather than decorative.

---

## Drafting notes

- Lead the abstract with the prediction and the pre-registration, not with the null.
- Give the exact pre-registered threshold wherever a result is quoted.
- Do not write "trending against the hypothesis": p = 0.382 does not support a trend claim.
  "Opposite in sign, not distinguishable from noise" is accurate.
- Name the two favourable robustness checks explicitly rather than leaving them in the
  table unmentioned.
- §3 is the new centre of gravity. If length pressure forces cuts, cut from §5, not §3.
- **Do not write that nothing further is testable on this cohort.** The erratum in
  `FINAL.md` establishes that Batenburg's cross-arm antagonism (Int J Mol Sci 2023) stands
  `untested` on data already held, and the Pol ζ axis (Wu 2023) is untouched and estimable.
  Both are legitimate future work and the paper should name them as such.
- Every number in this outline traces to a committed artifact; `key_figures.md` gives the
  path for each. Re-check each against its source before it goes into prose.
