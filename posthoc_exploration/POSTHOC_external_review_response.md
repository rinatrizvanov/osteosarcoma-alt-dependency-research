# POST-HOC EXTERNAL-REVIEW RESPONSE — NOT A CAMPAIGN ITERATION

> **NO ANALYSIS RESULTS ARE REPORTED IN THIS DOCUMENT.** Nothing here was run, tested or
> estimated. One proposed test was declined on substance, one was scoped but not executed,
> and one observation requires re-verification before it can be cited.
>
> **The campaign is closed and unchanged.** `part2_bir_dependency_test/autoresearch/FINAL.md`
> remains the conclusion of record; the campaign's Benjamini-Hochberg family still contains
> exactly three primary endpoints at q = 0.420, 0.420 and 1.000; the cumulative test count
> is still 101. No campaign file — `FINAL.md`, `LEDGER.md`, `NULL_SPACE.md`, `SUMMARY.md`,
> `PROGRAM.md` — was modified in producing this memo, and nothing here enters the campaign
> FDR family.

**Date:** 2026-10-14. **Trigger:** external review of the campaign results after closure
(tag `autoresearch-complete-v1.0`). Two additional tests were proposed by the reviewer and
assessed here. Any future execution of either would require its own pre-registration and
would carry weaker scientific weight than the pre-registered iterations, being
hypothesis-generating only.

---

## 1. DECLINED ON SUBSTANCE — Pol ζ / REV3L in ALT+ versus ALT− lines

**Proposed test.** That REV3L essentiality is more negative in the 16 ALT+ lines than in the
ALT− comparison lines, on the basis of Wu et al. 2023.

**Source.** Wu W, Barwacz SA, Bhowmick R, Lundgaard K, Gonçalves Dinis MM, Clausen M,
Kanemaki MT, Liu Y. *Mitotic DNA synthesis in response to replication stress requires the
sequential action of DNA polymerases zeta and delta in human cells.* Nat Commun 2023.
doi:10.1038/s41467-023-35992-5

**Declined.** Four independent reasons, in order of weight.

**(a) The paper makes no claim about ALT or telomeres.** A full-text search of the
Supplementary Information (`41467_2023_35992_MOESM1_ESM.pdf`, 13.0 MB) and the complete peer
review file (`41467_2023_35992_MOESM2_ESM.pdf`, 5.4 MB) returns **zero occurrences** of
`alternative lengthening`, `ALT-positive` or `telomere`. The source-data workbook is 38
sheets of mitotic-DNA-synthesis and common-fragile-site quantification. The paper's subject
is MiDAS at under-replicated regions and CFSs; its models are U2OS, HCT116 and HeLa. The
proposed test extrapolates from MiDAS at fragile sites to ALT telomere synthesis — two
different processes — and the paper does not address the second.

**(b) The paper's own cited prior work predicts the null.** Wu 2023 reference 77 is
Ozer O, Bhowmick R, Liu Y, Hickson ID. *Human cancer cells utilize mitotic DNA synthesis to
resist replication stress at telomeres **regardless of their telomere maintenance
mechanism**.* Oncotarget 9:15836–15846 (2018) — sharing a senior author with Wu 2023. If
telomeric MiDAS is TMM-independent, Pol ζ dependency should not separate ALT+ from ALT−.
The prediction's **direction is not supported**; the nearest relevant published claim
anticipates no difference.

**(c) The model systems argue for TMM-independence.** Wu's principal model, U2OS, is an
ALT-positive osteosarcoma line and is one of this project's 16 ALT+ lines. But HCT116 and
HeLa are both telomerase-positive and display the same sequential Pol ζ → Pol δ requirement.
A mechanism operating in one ALT+ and two TEL+ lines is evidence of generality, not of
ALT-specificity.

**(d) The literature chain implicates POLD3, which this cohort cannot test.** The
ALT-relevant polymerase citation in Wu 2023 is reference 42: Lydeard JR, Jain S, Yamaguchi M,
Haber JE. *Break-induced replication and telomerase-independent telomere maintenance require
Pol32.* Nature 448:820–823 (2007) — Pol32 being the yeast POLD3 ortholog. POLD3 is essential
in **16/16** ALT+ lines here (median −1.666) and 9/9 ALT− lines (median −1.797), with no
usable dynamic range. This is the construct defect already recorded as binding fact #2 in
`NULL_SPACE.md`.

**Additionally, the contrast was not blind.** REV3L's marginals in both arms were reported
during iteration 03's feasibility scan
(`part2_bir_dependency_test/autoresearch/iteration_03_bir_module_coherence/scan_essentiality.md`,
§3 and §4) and have been read in this project:

| Arm | n | mean | median | SD | n < −0.5 | n < −1.0 |
|---|---|---|---|---|---|---|
| ALT+ | 16 | −0.678 | −0.578 | 0.367 | 9/16 | 4/16 |
| ALT− | 9 | −0.602 | −0.632 | 0.127 | 7/9 | 0/9 |

The means differ by 0.076 in the predicted direction while the **medians run the opposite
way**, and the arms differ far more in spread than in location. Running a formal test on
values already read would have produced the appearance of a test without information gain —
the same contamination that caused iteration 04 to withdraw its own BRCA2 candidate
(`iteration_04_smarcal1_dependency_model/DISCLOSURES.md`). These numbers are reproduced here
as the reason for declining, **not as a result**.

---

## 2. SCOPED BUT NOT EXECUTED — cross-arm coupling between HR and BIR

**Status: designed, never run. No prespec was frozen, no code written, no estimate
produced.** The design is recorded so that it is preserved as future work.

**Hypothesis.** That dependency on the HR arm and on the BIR arm are negatively coupled
across ALT+ lines — the arms being antagonistic rather than orthogonal, so that the Stage 1
difference endpoint *added* two anticorrelated arms rather than contrasting two independent
ones.

**Source.** Batenburg NL, Walker JR, Zhu X-D. *CSB regulates pathway choice in response to
DNA replication stress induced by camptothecin.* Int J Mol Sci 24(15):12419 (2023).
doi:10.3390/ijms241512419 — BIR proceeds via MUS81–RAD52–POLD3, while RAD51, BRCA1 and
BRCA2 suppress the POLQ/LIG3/PRIMPOL restart routes. Supported by Richards et al. 2023
(doi:10.1038/s41467-023-40071-w) and Shah et al. 2024 (doi:10.1038/s41467-024-52916-z).

This hypothesis **stands untested rather than closed**, per the published erratum in
`FINAL.md`. An earlier draft of the campaign's closing screen declared it not-estimable on
the argument that genes essential in no line carry no signal; the campaign's own positive
control refutes that (BLM × TOP3A ρ = +0.682, p = 0.0036 and RMI1 × RMI2 ρ = +0.800,
p = 0.0002, with BLM and RMI2 each essential in 0/16 lines). Sub-threshold dependency
variation is recoverable biology in this cohort.

### Design as scoped

| Element | Specification |
|---|---|
| Cohort | the 16 ALT+ CRISPR-screened osteosarcoma lines, unchanged |
| HR arm | BRCA1, BRCA2, RAD51 |
| BIR arm | POLD3, RAD52, SMARCAL1 |
| **Primary endpoint** | mean of the **nine** cross-arm pairwise Spearman ρ, predicted sign **negative** |
| Secondary 1 | the same statistic on within-line-standardised gene effects |
| Secondary 2 | the same statistic with BRCA2 dropped from the HR arm |
| Multiplicity | Bonferroni α = 0.025 as originally proposed for two tests; with only one test in scope this requires restating |
| Detection floor | **must be derived by permutation before the data are read.** `src/coherence.py` implements a within-set null only; a cross-set statistic needs its own. For reference, a single Spearman at n = 16 requires \|ρ\| ≥ 0.494 at p < 0.05 and roughly 0.55 at α = 0.025 |

### Design notes that must survive into any future prespec

**Operationalisation is not arbitrary.** "Correlation between HR and BIR essentiality" can
mean the mean of the nine cross-arm pairwise correlations, or a single correlation between
the two arm means. These are different statistics and can disagree. The mean-of-nine was
chosen as primary for two reasons: it matches Batenburg's claim about pathway arms, and the
arm-means version is **partially contaminated** — iteration 03 published the exact variance
decomposition of `mean(BIR) − mean(HR)`, which together with the per-gene SDs constrains the
covariance between the arms. No cross-arm *pairwise* value exists anywhere in the project's
outputs (`iteration_03/data/h4_pairwise.csv` holds only the three within-triplet pairs), so
the mean-of-nine is clean.

**The screen-amplitude confound has a known direction, and it works against the
hypothesis.** All six genes share a per-line screen-quality component, which induces
*positive* correlation between any two gene-effect profiles. A negative finding would
therefore be credible despite the confound, while a null would be ambiguous. This matters
more than usual because POLD3 is essential in 16/16 lines, so its variation is variation
among already-essential lines — the most confound-prone kind. Reporting both raw and
within-line-standardised is not optional.

**BRCA2 may be on the wrong side of the contrast.** Lee et al. 2024
(doi:10.1093/nar/gkae251) show BRCA2 loss *driving* ALT/BIR via G4 stabilisation and TERRA
R-loops, so placing BRCA2 in the HR arm whose anticorrelation with BIR is the hypothesis
muddles the predicted sign. The campaign's own hard constraints also bar BRCA2 from any
HR-neutral control set, and iteration 04's disclosures withdrew it as a candidate on
contamination grounds. Hence secondary 2.

### Known information, disclosed

Already read and therefore not blind: all six genes' marginals in both arms; the within-BIR
pairwise values (POLD3–RAD52 −0.268, POLD3–SMARCAL1 +0.024, RAD52–SMARCAL1 −0.303); HR
module internal coherence (T = +0.414, p = 0.0039); BIR triplet coherence (T = −0.182,
p = 0.229); and the full variance decomposition of the Stage 1 composite. Not read: any
cross-arm pairwise correlation, which is the primary endpoint.

---

## 3. OPEN CROSS-CHECK — NOT PERFORMED, NO FINDINGS RECORDED

**Nothing is reported in this section, by design.** During the same session an external
telomere-maintenance resource was briefly opened and raised one question worth answering.
The source file was then replaced on disk before any figure could be re-checked, so **no
numbers from it are recorded here.** A value read once and not re-verified does not belong in
this repository, flagged or otherwise — the reader cannot check it, and a caveat is not a
substitute for verification.

**The question, which stands open.** Whether the ALT status of one or more lines in the
**ALT− comparison arm** is contradicted by direct assay. The frozen labels in
`references/master_manifest.csv` are computational calls sourced to `ALTitude_2024`; a
measured C-circle or telomerase readout would be stronger evidence, and the two could
disagree. If any ALT− line is in fact ALT-positive, it would bear on iteration 01's ALT+
versus ALT− contrast and on the n = 9 ALT− comparison arm used in iteration 03 — and it
would turn the standing "ALT status is not binary" limitation (prespec conflict 7 in
`Literature/INDEX.md`) from a general caveat into a specific, named instance.

**A second open question.** Whether any external resource supplies enough assay-confirmed
ALT+ lines that are *also* CRISPR-screened in DepMap 26Q1 to improve on the 16 ALT+ lines
already in this cohort. This matters because several hypotheses in `NULL_SPACE.md` are
recorded as untestable for want of an ALT-activity measure, and the limiting factor may be
sample size rather than measurement.

**How to answer both.** Obtain the supplementary tables of Wu et al. 2025,
doi:10.1038/s41467-025-67190-w, which reports telomere-maintenance metrics across cancer
cell lines. Join its cell-line table to this cohort on `DepMap_ID`, compare its telomere
maintenance mechanism calls against `master_manifest.csv` line by line, and count how many
of its ALT-classified lines appear in `CRISPRGeneEffect.csv`. Record the result only after
the join has been run and checked against the file.

---

## What this memo does not do

- It does not report a test, an effect size, a p-value or a q-value.
- It does not alter the campaign record, its FDR family or its test count.
- It does not reopen the campaign. `PROGRAM.md` remains closed; either test above would need
  a fresh mandate and its own pre-registration.
- It does not claim the biological hypotheses are false. Item 1 is declined as *not supported
  by its cited source*; item 2 is **untested**; item 3 is unverified.
