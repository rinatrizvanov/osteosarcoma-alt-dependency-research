# NULL_SPACE.md — hypotheses tried and found unsupportable

The purpose of this file is that **no future session repeats work recorded
here.** Every entry carries the reason and the numbers, so a reader can tell a
genuine null from a feasibility failure from a missing-data skip.

Three kinds of entry:

- **tested-and-null** — pre-registered, run, did not support. The numbers are
  the finding.
- **not-estimable** — the cohort cannot support the test (group sizes, no
  usable exposure variable, contradictory measurements). Killed at the
  feasibility gate, before any outcome was read.
- **skipped-missing-data** — requires data DepMap 26Q1 does not ship. Not an
  argument against the hypothesis, only against testing it here.

---

## skipped-missing-data

### H5 — WRAP53 expression as a proxy for 17p13.1 structural disruption

Requires a cohort with both WGS SV calls and RNA at 17p13.1. DepMap 26Q1 ships
no SV calls (`docs/analytical_prespec.md` §12; `references/label_spec_v1.0.md`
§6.1). External data (TARGET-OS, PCAWG) or a collaborator-held PDX cohort would be
needed; the latter is stop condition (d).

Recorded before the loop began, from `Literature/openalex_scan/next_hypotheses.md`
H5 and `Literature/INDEX.md` Gaps §3 — axis G of the literature scan found **no
corpus at all** on when expression proxies fail for SV-driven phenotypes. The
assumption is the weakest-evidenced part of the Stage 1 design and it remains
unresolved rather than refuted.

---

## not-estimable

### H2 — ATRX-gated SMARCAL1 dependency in ALT+ osteosarcoma

Iteration 02, killed at the feasibility gate. **Do not retry on DepMap 26Q1
osteosarcoma.** Full numbers in `iteration_02_atrx_gated_smarcal1/outcome_memo.md`.

Four candidate ATRX rules over the 16 ALT+ CRISPR-screened lines: manifest
point-mutation **1** (NRHOS1), expression `TPM < 3.0` **2**, gene-level
`CN < 0.01` **4**, segment-validated focal homozygous loss **1** (CAL72). Only
the gene-level CN rule clears `MIN_CELL_SIZE = 4`, and it is wrong: for NOS1,
OS457 and U2OS the ATRX-spanning segment reads 0.41–0.67 with no focal event,
so 3 of its 4 calls are unsupported. The two n = 1 rules disagree on the
identity of the line. MDE at the unreliable 4v12 is 0.708.

**Still open, but needs a user decision.** A pan-cancer version using
ALTitude's expanded ALT+ set across lineages would have the group sizes, but
requires ALT labels for non-osteosarcoma lines. The frozen manifest covers
osteosarcoma only, and the loop's hard rule forbids deriving new ALT labels.
Relaxing that is a scope decision for the user, not the loop.

**Independently corroborated, and now closed on published evidence too** (found by
iteration 04's literature scan, `iteration_04_smarcal1_dependency_model/scan_literature.md`).
Bennett et al. 2026 (ALTitude, https://doi.org/10.64898/2026.03.04.709438) ran exactly this
contrast at pan-cancer scale and **published a null**: SMARCAL1 dependency in ATRX-loss
ALT+ lines (n = 20) versus ATRX-wildtype ALT+ lines (n = 14), Cohen's d = **−0.248**,
**P = 0.5**.

This matters twice over. It corroborates iteration 02's feasibility kill on independent
data rather than on our group sizes alone, and it removes the "pan-cancer version with
proper group sizes" escape hatch described in the paragraph above — that study has been
done, at n = 34, and found nothing. **Do not spend a user decision on relaxing the ALT-label
rule for this purpose.** Bennett also flags within-ALT+ heterogeneity as mechanistically
unexplained and offers only a single-line account, so the question of what *does* predict
SMARCAL1 dependency inside ALT+ remains genuinely open — which is what iteration 04 takes
up. Bennett's reciprocal ATRX-SMARCAL1 compensation model must not be re-derived here.

### Gene-level copy number from OmicsCNGeneWGS.csv must be segment-validated

Binding on every future iteration. Gene-level ATRX values disagree with
`OmicsCNSegmentsWGS.csv` in 3 of 4 apparent homozygous deletions. A naive
"CN loss OR expression silencing" re-derivation would have produced a 4v12
ATRX-gated analysis with three of four exposed lines misclassified, undetectably.

Third defect of this kind in the project: Stage 1's `load_wrap53_cn` layout bug;
the manifest's artefactual `sv_architecture_auto` (33/35 rows from a single CN
segment spanning both WRAP53 and TP53); and now gene-level ATRX.

---

### H4 — POLD3, RAD52 and SMARCAL1 as a coherent module in ALT+ osteosarcoma

Iteration 03. **Verdict of record `not-estimable`** by the frozen §6 assay-sensitivity
gate; mechanical verdict `null` (T = −0.182, p = 0.229, n = 16, floor 0.289). Full numbers
in `iteration_03_bir_module_coherence/outcome_memo.md`.

**Do not retest the triplet's coherence on this cohort under a relabelled gate.** T and p
would be identical; only the interpretive label would change, which is laundering rather
than a new test. The corrected gate design (several control complexes, pass if any clears)
belongs in any *future* prespec that uses this statistic, applied to a different question.

What was observed: no pair positive (POLD3–RAD52 −0.268, POLD3–SMARCAL1 +0.024,
RAD52–SMARCAL1 −0.303); the triplet at the 12.5th percentile of a mean/SD-matched
random-triplet null; Batenburg's own MUS81/RAD52/POLD3 module also anti-coherent (−0.274,
MUS81–RAD52 −0.512); within-line standardisation changes nothing (−0.188); the ALT− arm
runs the opposite sign (+0.378, p = 0.058, n = 9).

**Binding on every future iteration — three facts, deterministic, not verdict-dependent:**

1. **Chronos co-essentiality is estimable on these 16 lines.** BTR reaches T = +0.526,
   p = 0.0001 and the HR module T = +0.414, p = 0.0039, both clearing their floors, with
   BTR doing so while containing a near-saturated gene (TOP3A) and two essential in no line
   (BLM, RMI2). Low n is **not** an available explanation for a null on this cohort.
2. **POLD3 and RAD52 are unusable as module members here.** POLD3 is essential in 16/16
   lines, RAD52 in 0/16. Any future endpoint averaging either with a mid-range gene inherits
   the same construct defect. Check marginal regimes before composing any multi-gene score.
3. **The Stage 1 BIR−HR composite drew 68.2% of its variance from the HR arm**, with RAD52
   at −2.7%. Do not reuse that endpoint.

### Type I vs type II ALT subtyping on DepMap 26Q1 — no defensible marker

Established in iteration 03 (prespec §7 S5) and binding on rule (c) for every later
iteration: the explicit-statement branch is the correct one for this cohort, and no proxy
should be invented.

Type I/II ALT is defined in **mouse** ESCs via the ~7 kb mTALT template, and Lee et al.
2026's own supplementary material (https://doi.org/10.1016/j.mocell.2026.100405) reports
that **no compact human genomic segment resembling it was identified in GRCh38**. The three
human two-mode models (Broderick 2023 `10.1038/s41467-023-38029-z`; Zhang 2026
`10.1002/bies.70187`; Folini 2026 `10.21203/rs.3.rs-11116403/v1`) all require activity,
imaging or C-circle readouts.

The cohort carries **no ALT-activity measure of any kind**: `alt_assay_evidence` is
`unknown` for all 35 rows, there is no C-circle, APB, TRAP or telomere-length column, and
`telomerase_status` is constant (`indeterminate`) within the ALT+ set. Over the 34
RNA-backed lines **no candidate marker gene is bimodal** — the largest guard-respecting
interior gap is below 1 SD for all 16 tested (gap/SD 0.22–0.90) — so any cut is a researcher
choice producing a label that nothing in the cohort could falsify. The only
expression/CN-obtainable partition is ATRX/DAXX-driven vs SMARCAL1-driven, which is not
type I/II and which iteration 02 already found non-estimable.

---

## tested-and-null

### H1 — WRAP53 **mRNA** does not discriminate ALT status in DepMap 26Q1 osteosarcoma lines

Iteration 01, verdict `counter_signal`. 18 ALT+ v 12 ALT− (patient-collapsed),
δ = +0.241 against a predicted-negative direction, p = 0.280, AUROC 0.380.
ALT+ lines are **higher** in WRAP53 (median 4.193 v 3.483).

**Do not retest this as specified.** The counter-direction is stable: 30/30
leave-one-out refits, 12/12 leave-one-depositor-out refits, and it *strengthens*
when the four pre-known label–exposure pairs are dropped (R3b δ = +0.418,
p = 0.077). Gold-only labels agree (δ = +0.286). It is below the 0.435 detection
floor, so it is suggestive rather than established — but no variant that keeps
WRAP53 mRNA as the exposure and binary ALT status as the outcome is worth
another iteration on this cohort.

**What remains open.** Keegan's tumour-level SV claim is untouched — it concerns
structural variants in tumour WGS assayed by protein IHC. What is closed is the
**mRNA proxy in cell lines**. Two pointers: WRAP53 copy number runs the
*opposite* way to its own mRNA (R4, δ = −0.167), and TP53 expression
discriminates where WRAP53 does not (R5, δ = −0.472, p = 0.031, descriptive).

### WRAP53 expression has no natural low group in this cohort — binding on any future WRAP53 stratum

From `iteration_01/scan_exposure.md`, established blind to the ALT label. Over
the 34 RNA-backed osteosarcoma lines WRAP53 log2(TPM+1) is **not bimodal**: one
continuous mode spanning 2.041–5.888 plus a discrete **2-line spike at exactly
0.0**, with no low-but-nonzero values. Largest interior gap 2.041 units (3.3× the
next largest, 1.58 SD).

Consequences for any future iteration: the only naturally separable low group
has **n = 2**, below `stratify.MIN_CELL_SIZE = 4`; every percentile cut
(0.25/0.33/0.50 → 3.078/3.336/3.603) sits inside the dense body and is
arbitrary; and the natural `WRAP53 == 0` cut **cannot reach significance at any
configuration** (most extreme attainable Fisher two-sided p = 0.497, one-sided
0.329 = C(18,2)/C(31,2)).

This also retrospectively explains a Stage 1 mechanic: the §4b minimum-cell
guard fired and forced the median fallback because the exposure has no natural
low group, not by bad luck.

### H6 — PRIMPOL expression as a predictor of SMARCAL1 dependency within ALT+ osteosarcoma

Iteration 04, verdict `null`. Spearman rho = **exactly +0.0000**, p = 1.000, n = 16, floor
0.494. Verified against scipy (sum d² = 680 = n(n²−1)/6 exactly, no ties) because an exact
zero is a plausible bug. Robust to within-line standardisation (+0.012). Post-hoc bootstrap
interval [−0.502, +0.522]. Full numbers in
`iteration_04_smarcal1_dependency_model/outcome_memo.md`.

**Do not retest PRIMPOL mRNA against SMARCAL1 dependency on this cohort.** The association
is not weak, it is absent, and standardisation does not recover it.

**Not a refutation of Taglialatela's model**, and the prespec pre-committed to not claiming
one: the floor is |rho| >= 0.500, the exposure is a transcript rather than repriming
activity, and the mechanism's stronger arm was never testable here (below).

### BLM expression cannot serve as a continuous exposure in this cohort — binding

Established in iteration 04 §2 on exposure data alone, before any outcome was read, and
verified independently of the scan that first reported it. Over the 16 ALT+ CRISPR-screened
lines BLM log2(TPM+1) is:

```
0.96  2.32  |  3.30 3.31 3.35 3.47 3.68 3.77 3.84 3.94 3.95 4.15 4.16 4.19 4.26 4.79
```

Two low lines, then **14 inside a 1.496-unit band**. Shapiro p = 0.0059, excess kurtosis
+2.88, skew −1.69, largest guard-respecting interior gap **0.220**. A rank correlation is
decided by where those two lines fall — a 2-versus-14 contrast wearing the clothes of a
correlation — and n = 2 is below `stratify.MIN_CELL_SIZE = 4`, so the honest categorical
form is barred too. BLM is also the hub of the panel's collinearity cluster (+0.737 with
FANCM over the 34 lines, +0.871 over the 16).

**This matters more than it looks**, because BLM is the *stronger* arm of Taglialatela's
mechanism: *"BLM deficiency fully suppressed the SMARCAL1-dependent growth defects"*, against
PRIMPOL's merely partial rescue. So the one prediction this cohort could test was the weak
one. Any future attempt on BLM needs a cohort with real BLM variance, not this one. The
parallel Stage 1 fact — WRAP53 has no natural low group here — is recorded above.

### A general lesson, binding on every later prespec: write support conditions as absolute tests

Iteration 04's conditions compared the control's and the confounder's |rho| **against the
primary's magnitude**. The primary came back at exactly zero, so three of four failed
arithmetically while carrying no information — the negative control's own p was 0.540, i.e.
it behaved exactly as intended. A condition phrased relative to the primary is undefined
when the primary is null, which is the most likely outcome in an underpowered campaign.

Write them absolutely instead: "control p >= 0.05", "|rho(exposure, screen amplitude)| <
0.30", "leave-one-out range excludes zero **when the point estimate does**". This is the
second prespec defect the data have exposed in two iterations, after iteration 03's
single-designated-control assay gate, and both were recorded rather than retroactively
rewritten.

---

## Closing additions — the iteration-05 candidate screen (2026-10-12)

Added when the campaign closed on stop condition (c). Full screen with all thirteen
candidates in `autoresearch/candidate_screen_iteration05.md`; the four entries below are
the ones newly killed, each with the numbers that killed it. All marginals are Chronos gene
effect over the 16 ALT+ CRISPR-screened lines, computed **marginals-only, no correlation of
any kind**, per binding fact #2 above.

### not-estimable — BIR and HR arms coupled with opposite sign (Batenburg 2023)

Prespec-conflict row 4, the strongest survivor of the screen. **Not** a variant of
iteration 03: cross-set rather than within-set statistic, negative rather than positive
predicted sign, and no cross-arm pairwise value exists anywhere in the campaign's outputs.
It dies on the restart arm's dynamic range.

Of Batenburg's three named restart genes, two are in iteration 03 §7.3, the
"non-essential throughout — not dependency signal" class:

```
LIG3     median -0.1432  SD 0.1358  min -0.3357   0/16 below -0.5
PRIMPOL  median +0.1323  SD 0.1488  min -0.2362   0/16 below -0.5
POLQ     median -0.4089  SD 0.1780  min -0.6501   5/16 below -0.5   <- the only usable one
```

The HR arm is usable and is a real module here (iteration 03: internal T = +0.414,
p = 0.0039), but BRCA2 is barred by §6 and by iteration 04's DISCLOSURES, leaving
BRCA1 (SD 0.440, 9/16) and RAD51 (near-saturated, 16/16). A cross-set statistic between one
coherent module and a two-thirds-noise set **has expectation ≈ 0 whether or not the claim
is true**, so a null would be uninformative by construction. The one clean cross-pair,
BRCA1 × POLQ, has floor |rho| = 0.494 at n = 16 and would be one cell chosen from a 3 × 3
grid after seeing which cells were usable.

**Needs a cohort where the alternative-restart routes are actually depended upon.** The
claim itself is untouched.

### not-estimable — CCT/TRiC sets TCAB1 protein level (Jiang 2025; Bergstrand 2020)

Prespec-conflict row 5. The claim is about protein and DepMap 26Q1 ships **no proteomics**,
so the direct form is `skipped-missing-data`. Recast as a dependency model it dies on
dynamic range — the chaperonin is essential in every line, as a chaperonin should be:

```
CCT5  median -1.2029  SD 0.2533  16/16 below -0.5  12/16 below -1.0
CCT7  median -0.8339  SD 0.2003  15/16 below -0.5   4/16 below -1.0
CCT2  median -1.4525  SD 0.3000  16/16 below -0.5  14/16 below -1.0
TCP1  median -1.9882  SD 0.4411  16/16 below -0.5  16/16 below -1.0
```

The row's real content — WRAP53 mRNA is the wrong measurement for TCAB1 function — is a
measurement-validity verdict on the Stage 1 exposure, already recorded above.

### not-estimable — WRAP53β impairs HR and NHEJ non-selectively (Henriksson 2014)

Prespec-conflict row 2. The *selective* version is removed outright by the loop's §6 hard
constraint. The non-selective version would have to be recast as a co-dependency
prediction, translating a perturbation-biology claim into a correlational one — proxy
invention. Independently the channel is absent: `WRAP53` median −0.1726, SD 0.2473,
min −0.5695, **1 of 16** lines below −0.5. Its dominant role in these cells is
Cajal-body/scaRNA trafficking (Mahmoudi 2010), not DSB repair, so even a positive
correlation would not be attributable to the repair claim.

### not-estimable — SMARCAL1 LOF promotes ALT, so the Stage 1 sign was inverted

Prespec-conflict row 1 (Liu 2023; Brosnan-Cashman 2021; Sobieszek 2026), the most
attractive row in the table and osteosarcoma-specific. All three operationalisations are
closed:

- **SMARCAL1 mRNA as exposure** — barred by iteration 04's frozen prespec on Bennett's
  documented transcript↔protein dissociation for SMARCAL1 itself (line SKNFI):
  "SMARCAL1 transcript level cannot appear as exposure, covariate or control." Written
  pre-freeze and against interest; reviving it to keep the loop alive would be the post-hoc
  rule change §3 forbids.
- **SMARCAL1 copy number as exposure** — a copy-number claim class the campaign already
  excluded once on stated grounds (TOP3A), confounded with the known CRISPR copy-number
  artefact when the outcome is the same gene's Chronos score, and rho +0.716 with the barred
  expression variable.
- **The direct ALT+ vs ALT− dependency contrast** — **contaminated.** The values were read
  during iteration 03's feasibility scan (`scan_essentiality.md` §3/§4: ALT+ mean −0.4907,
  SD 0.395, 8/16 below −0.5; ALT− mean +0.0557, SD 0.137, 0/9) and the contamination is
  declared in `iteration_04/DISCLOSURES.md`. It is also Taglialatela's and Bennett's
  published headline. Reported in `FINAL.md` as a descriptive corroboration, explicitly not
  as a campaign result.

### The binding generalisation — why the campaign closed

Every remaining literature-grounded hypothesis needs **either** a variable this cohort does
not contain (an ALT-activity readout, a WGS SV call, a protein level, a branch label)
**or** a gene whose dependency does not vary across these 16 lines. Sample size is not the
constraint: iteration 03 proved co-essentiality is estimable here (BTR T = +0.526,
p = 0.0001). Any future work on these hypotheses needs a different cohort, not a cleverer
statistic on this one.

---

## ERRATUM 2026-10-12 — row 4 is `untested`, not `not-estimable`

**The "BIR and HR arms coupled with opposite sign (Batenburg 2023)" entry above is
withdrawn from the not-estimable class.** Reason, numbers and the surviving cautions are in
`FINAL.md`'s erratum section. In short: I argued that LIG3 (0/16 below −0.5) and PRIMPOL
(0/16) carry no dependency signal, so the cross-set statistic would have expectation ≈ 0
regardless of the claim's truth. The campaign's own positive control refutes that — BTR
contains BLM (0/16) and RMI2 (0/16), and BLM × TOP3A is ρ +0.6824 (p 0.0036), RMI1 × RMI2
ρ +0.8000 (p 0.0002), BLM × RMI2 ρ +0.4441 with both members essential in no line.
Sub-threshold variation is recoverable biology here.

A future prespec on row 4 must still handle: BRCA2 barred (§6 + iteration 04 DISCLOSURES),
so a 2 × 3 cross-set null has to be derived rather than reused; and the construct-validity
question of whether a Chronos score reports *graded reliance* on a route no line depends on
— which is a question to pre-specify, not to assume either way.

**The two entries added above that remain correct:** CCT/TRiC (no proteomics; chaperonin
essential in 15–16 of 16 lines) and the SMARCAL1-polarity operationalisations (mRNA barred
on Bennett's transcript↔protein dissociation; CN confounded; the direct contrast
contaminated). The WRAP53/Henriksson entry rests partly on the same overstated
noise argument and should be re-read as "requires a co-dependency recast that is proxy
invention", which is the part of it that stands.

**Also noted, not a candidate the loop screened:** the Pol ζ axis (Wu 2023, INDEX #16 —
MiDAS requires Pol ζ then Pol δ) has usable dynamic range in these 16 lines — REV3L median
−0.5777, SD 0.3674, 9/16 below −0.5; REV1 median −0.0341, SD 0.2861, 1/16; MAD2L2
pan-essential 16/16. It failed only filter 1 (not in the Prespec Conflicts table), never a
feasibility test. The RAD51 paralog complex is likewise untouched and estimable: RAD51C
12/16 (SD 0.263), RAD51D 8/16 (SD 0.294), XRCC2 4/16, RAD51B 2/16.
