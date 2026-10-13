# Preprint outline — short report

**Working title.** Pre-registered test of the Recombination Paradox prediction in DepMap
osteosarcoma cell lines: non-confirmation

**Target.** bioRxiv (Cancer Biology / Genomics). Short report format: ~2,000–2,500 words,
2 main figures, 1 main table, supplementary pre-specification.

**Status.** Outline only — section headers with 2–3 sentence summaries. Prose to be
written by the author.

**Framing note.** The paper's contribution is the pre-registration and the clean
non-confirmation, not a discovery. Every section should read as though the result had
come out the other way: the method is the claim. Avoid both over-claiming the
counter-trend and apologising for a null.

---

## Title and abstract

**Abstract (~200 words, structured).** State the mechanistic prediction, that it was
pre-registered before any outcome value was read, the cohort (16 ALT+ DepMap
osteosarcoma lines), and the headline: Cliff's δ = +0.281 (95% CI −0.313 to +0.844),
p = 0.382, opposite to the predicted direction, minimum BH q = 0.959. Close on what the
result does and does not exclude, and on the CI containing moderate predicted-direction
effects.

## 1. Introduction

**1.1 ALT and the telomerase-independent escape route.** Introduce ALT as the mechanism
used by a large fraction of osteosarcomas, and SMARCAL1 as an established broad ALT
dependency (Bennett et al. 2026) that has not been stratified by upstream genotype.

**1.2 The 17p13.1 / TCAB1 prediction.** Keegan et al. 2026 describe structural variants
disrupting TCAB1 and WRAP53β together; TCAB1 loss blocks telomerase and WRAP53β loss
cripples classical HR (Henriksson & Farnebo 2015). State the resulting prediction
explicitly: a cell denied both routes should lean harder on non-classical BIR.

**1.3 Why pre-register this one.** The prediction is attractive, the available N is
small, and the analytic degrees of freedom are large — exactly the combination that
produces irreproducible positives. State that the full analysis plan, success criteria
*and* failure criteria were committed to version control before the outcome was read.

## 2. Results

**2.1 Cohort and stratification.** 35 osteosarcoma lines → 26 CRISPR-screened → 16 ALT+,
9 ALT− control, 1 excluded. Report that the pre-specified bottom-tertile rule placed only
3 of 16 ALT+ lines in the low arm, triggering the pre-registered minimum-cell guard, so
the primary became a median split with the tertile retained as a parallel sensitivity.
Emphasise that the fallback was written before the WRAP53 distribution was examined.

**2.2 The stratifier is a valid measurement.** Summarise the strand check: WRAP53 exon 1
overlaps TP53 antisense by ~1,419 bp, DepMap does not document strandedness, and a
pre-committed rule was resolved empirically — WRAP53 tracks its own copy number
(ρ = +0.554) more strongly than TP53 (ρ = +0.244), the opposite of cross-assignment.
This belongs in Results, not Methods, because it is a result about the data.

**2.3 Primary comparison: no support, wrong direction.** **Figure 1.** Report δ = +0.281,
p = 0.382 and the parallel tertile sensitivity δ = +0.538, p = 0.189. State plainly that
the point estimate is on the opposite side of zero from the prediction.

**2.4 No gene carries a signal.** **Figure 2.** All six BH q = 0.959; every 95% interval
spans zero. Note honestly that RAD52 (−0.344) and BRCA2 (−0.250) have negative point
estimates — some genes do lean the predicted way — but nothing is distinguishable from
noise.

**2.5 Robustness: eight of ten checks agree.** **Table 1.** Walk the pre-specified set
(R1–R8, §6.2, §6.3). Highlight that no leave-one-out refit of 16 reaches the predicted
direction, and that the two favourable checks are the two smallest comparisons — one of
them the negative control, where a predicted-direction result argues *against* a
TCAB1-specific account.

**2.6 Verdict.** The frozen criteria return `indeterminate`: not `falsified` only because
|δ| exceeds the 0.20 clean-null ceiling while the sign is wrong. Explain that the criteria
were evaluated mechanically by code so the verdict could not be shaded by reading the
numbers.

**2.7 Incidental: OS384 carries the predicted architecture.** OS384 shows zero expression
of both WRAP53 and TP53 with WRAP53 copy number ≈ 0 — the 17p13.1 co-deletion
architecture, visible in a public cell line, recorded in the pre-specification before the
outcome was read. Present as an observation of interest to others working on the locus,
explicitly decoupled from the dependency result.

## 3. Discussion

**3.1 What the result excludes, and what it does not.** The CI (−0.313 to +0.844) does
not exclude a moderate effect in the predicted direction, and §5a committed in advance
that only |δ| ≥ 0.594 was detectable at this N. Be precise: this is a non-confirmation,
not a refutation.

**3.2 Low power does not explain a wrong-direction estimate.** Underpowering explains
failure to detect. It does not explain why the primary, both tertile variants, the
continuous correlation and all 16 leave-one-out refits agree on the opposite sign. This
is the paragraph that distinguishes the paper from a generic null.

**3.3 The proxy is the binding limitation.** DepMap 26Q1 ships no SV calls, so WRAP53
expression stood in for 17p13.1 structural disruption throughout. A cohort with real SV
calls would test the prediction as stated; this tests a proxy for it, and that gap is
the most likely place the prediction survives.

**3.4 The SMARCAL1 counter-observation.** SMARCAL1 ran backwards (δ +0.436 at 3v13) and
dropping it moved the composite from +0.281 to +0.156. Flag as post-hoc, underpowered and
non-significant, and resist building an account on it — note only that it is the
opposite of what a TCAB1-gated BIR dependency predicts.

**3.5 What pre-registration bought, including where it bit.** Two criteria were repaired
before the freeze on exact-null grounds (the selectivity rule would have blocked a true
positive 20–30% of the time). One was mis-calibrated and was *not* repaired afterwards:
§10's clean-null ceiling assumed a null in the predicted direction, so a wrong-direction
estimate fell into `indeterminate`. Report this as a design lesson, and say that the rule
was left frozen deliberately.

**3.6 Reporting nulls in dependency screens.** Short closing argument: DepMap makes
hypothesis-driven dependency tests cheap, which makes selective reporting cheap too. A
pre-registered null with published criteria and code is the low-cost correction.

## 4. Methods

**4.1 Data and cohort.** DepMap 26Q1 specifically (the ALT labels were derived against
it), five files sha256-verified at load. Inclusion filter applied explicitly in order with
N logged at each step; ground-truth anchors asserted before analysis.

**4.2 Endpoint and tests.** The BIR−HR composite with its sign convention; two-sided exact
Mann-Whitney; Cliff's δ with 2,000-resample bootstrap CIs validated against
δ = 2U/(nm) − 1; Hodges-Lehmann shifts; BH-FDR across exactly the six genes.

**4.3 Pre-specification and enforcement.** The freeze commit, the amendment made before
any outcome was read, and the structural enforcement: one module loads gene effects,
notebooks 02–04 assert it is absent from `sys.modules`, criteria evaluated by code.

**4.4 Power.** The exact minimum detectable |δ| per split, computed exhaustively over
rank configurations and committed in advance.

**4.5 Code and data availability.** Repository at release `stage1-v1.0-null`; DepMap not
redistributed, fetch policy given.

---

## Figures and tables

| Item | Content | Source |
|---|---|---|
| **Figure 1** | BIR−HR composite by TCAB1 stratum; panel a primary 8v8, panel b tertile 3v13; OS384 and NRHOS4 labelled | `figures/figure1_composite_by_stratum.pdf` |
| **Figure 2** | Per-gene Cliff's δ with 95% bootstrap CIs, both contrasts, predicted-direction region shaded | `figures/figure2_per_gene_forest.pdf` |
| **Table 1** | All pre-specified checks: n, δ, CI, p, direction | `data/processed/robustness_table.csv` |
| **Supplementary** | The frozen pre-specification, verbatim | `docs/analytical_prespec.md` |
| **Supplementary** | Strand-check diagnostics and null distribution | `data/processed/notebook_02_report.md`, `figures/strand_check_null.png` |

## Drafting notes

- Lead the abstract with the prediction and the pre-registration, not with the null.
  A reader deciding whether to cite needs to know what was tested and that the plan was
  fixed in advance.
- Give the exact pre-registered thresholds wherever a result is quoted, so no number
  looks chosen after the fact.
- Keep §3.2 short and hard. It is the one genuinely informative asymmetry in the data.
- Do not describe the result as "trending against the hypothesis" — p = 0.382 does not
  support a trend claim. "Opposite in sign, not distinguishable from noise" is accurate.
- Name the two favourable checks explicitly rather than burying them; a reader who finds
  them in the table unmentioned will trust the rest less.
