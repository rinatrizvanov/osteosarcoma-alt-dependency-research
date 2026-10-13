# scan_literature.md — Iteration 03, H4 (BIR module coherence)

**Hypothesis under evaluation.** H4: POLD3, RAD52 and SMARCAL1 do not behave as a
coherent module in ALT+ osteosarcoma cells — their CRISPR gene-effect profiles across
ALT+ lines are mutually uncorrelated, so no composite over them is a legitimate
endpoint. Pre-specified in the **module-premise** direction (premise: positive
co-dependency; predicted sign **+1**), so that premise failure reads as a verdict
against the Stage 1 endpoint rather than as confirmation of a null.

**Sources.** `Literature/INDEX.md` (87 readable PDFs, Priority Reading, Prespec
Conflicts table) used as the primary interface. Two PDFs pulled via `pdftotext`:
Broderick 2023 (`10.1038/s41467-023-38029-z`) and Taglialatela 2026
(`10.1101/gad.353377.125`). Supplementary text of the otherwise-unretrievable
Lee 2026 (`10.1016/j.mocell.2026.100405`) read from its SI `.docx`. Two WebSearch
queries for Q3. Scan date 2026-10-12.

---

## 1. Is the module premise supported?

**Mechanistic coupling: partially, for POLD3+RAD52 only. Correlated gene-effect across
cell lines: no paper claims it, for any pair in the triplet.**

### 1a. What the papers actually assert

- **Batenburg et al. 2023** (`10.3390/ijms241512419`) is the strongest module claim in
  the corpus. At high camptothecin, CSB promotes replication restart through
  **MUS81–RAD52–POLD3-dependent BIR**. This is a genuine one-pathway statement for
  RAD52 and POLD3. Note the third member is **MUS81, not SMARCAL1**.
- **Richards et al. 2023** (`10.1038/s41467-023-40071-w`) couples the same two genes
  one layer deeper: Mre11/ATM-induced PARP1/PARP2 activity under replication stress
  promotes BIR via RAD52 assembly at stalled forks, and **site-specific
  ADP-ribosylation of PolD3 is required for BIR activity and fork recovery** —
  explicitly in the regime "when HR is dysfunctional". RAD52 and POLD3 are therefore
  mechanistically coupled, and the coupling is the closest mechanistic match in the
  literature to the Stage 1 intuition.
- **Broderick et al. 2023** (`10.1038/s41467-023-38029-z`) states that *both* ALT BIR
  branches "rely on BLM and **PCNA-RFC-Polδ**" — i.e. POLD3 sits downstream of the
  branch point, shared by both routes.
- **No paper in the 87-PDF corpus names POLD3, RAD52 and SMARCAL1 together as one
  module.** SMARCAL1 appears in neither Batenburg's nor Richards's BIR module. The
  Stage 1 triplet substituted SMARCAL1 for MUS81.
- **Bennett et al. 2026** (`10.64898/2026.03.04.709438`), the source of the frozen
  `alt_silver` labels and of the SMARCAL1-dependency claim Stage 1 was built on,
  mentions POLD3, RAD52 and "BIR" **zero times** (INDEX.md entry 5). The composite was
  testing a link that paper never asserts.

### 1b. The distinction that decides this iteration

Mechanistic coupling within a pathway and correlated gene-effect across cell lines are
different claims, and the first does not imply the second.

- Coupling in a **linear** pathway predicts **epistasis**: knocking out one member
  masks the phenotype of knocking out another. In a CRISPR gene-effect matrix that
  produces *redundancy and non-additivity*, not positive correlation.
- Correlated gene-effect across lines requires three further conditions: (i) the
  *degree of reliance* on the pathway varies across the lines in the panel; (ii) each
  gene is rate-limiting in proportion to that reliance; (iii) no member carries
  substantial pathway-external function or redundancy.
- The co-essentiality literature is explicit that correlation recovers **obligate,
  stoichiometric complexes** far more reliably than it recovers *pathways*:
  Pan et al. 2018 (`10.1016/j.cels.2018.04.011`) derives complex structure, function
  and membership from fitness-screen correlation; Wainberg et al. 2021
  (`10.1038/s41588-021-00840-z`) notes that ubiquitous false positives are the central
  obstacle and that co-essentiality yields **pan-cell-type**, not context-specific,
  interactions. POLD3/RAD52/SMARCAL1 is not a complex.
- The triplet violates all three conditions: a documented **branch point** (§2,
  Broderick 2023), **POLD3 pathway-external function** (§2, Wu 2023), and
  **post-translational gating** of the switch (Richards 2023,
  `10.1038/s41467-023-40071-w`: the ADP-ribosylation switch is invisible to a CRISPR
  knockout of POLD3).

**Verdict on Q1.** The module premise as a *mechanistic* claim is supported for the
POLD3–RAD52 pair and unsupported for SMARCAL1. The module premise as a *co-dependency*
claim — the one the Stage 1 composite actually required — is asserted by nobody.

---

## 2. Is it contradicted?

**Yes, and specifically enough to predict the sign of each pairwise correlation.**

### 2a. Broderick 2023 — RAD52 and SMARCAL1 are on opposite sides of a switch

Read from the PDF, not the index. Broderick et al. 2023
(`10.1038/s41467-023-38029-z`) establish that ALT "relies upon **two distinct** ...
break-induced-replication (BIR) processes", that EXD2 fork-processing activity
"triggers a switch between **RAD52-dependent and -independent** ALT-associated BIR",
and that the RAD52-independent branch "is suppressed by EXD2 but **depends
specifically on the fork remodeler SMARCAL1 and the MUS81 nuclease**".

The architecture this implies for the Stage 1 triplet is:

| Gene | Position in Broderick's model | Predicted gene-effect behaviour across ALT+ lines |
|---|---|---|
| RAD52 | marks branch A (dominant telomere-lengthening route) | varies with branch usage |
| SMARCAL1 | required specifically for branch B (RAD52-**independent**) | varies **inversely** with branch usage |
| POLD3 | shared downstream polymerase (PCNA-RFC-Polδ, both branches) | near-constant; low variance |

A pathway **switch** predicts ρ(RAD52, SMARCAL1) ≤ 0, not > 0. A shared downstream node
predicts ρ(POLD3, ·) ≈ 0 because low variance leaves nothing to correlate. The paper the
loop cited for "BIR" therefore predicts the *opposite* of the module premise's +1 sign
for the one pair that carries branch information. Broderick also reports that the
RAD52-independent route "is purported to maintain telomeres independently of RAD51 and
MRE11, while the ablation of these HR factors is associated with C-circle
overproduction" — relevant to Q4.

### 2b. Carson & Flynn 2023 — RAD51 and RAD52 run distinct ALT branches

Carson & Flynn 2023 (`10.1016/j.coph.2023.102380`) describe RAD51 and RAD52 as
regulating **distinct branches** of ALT. Per INDEX.md, "BRCA" appears **0 times** in
that review. The BIR-versus-HR contrast Stage 1 built is not how the ALT-vulnerability
literature carves the problem.

### 2c. Wu 2023 — POLD3 is a general replication-stress reporter

Wu et al. 2023 (`10.1038/s41467-023-35992-5`) show MiDAS requires both gap-filling
(REV1, Pol ζ) and replicative (Pol δ) polymerases with a polymerase switch, and that
oncogene-activated cancer-cell survival depends on it. POLD3 essentiality therefore
reports a **general oncogenic replication-stress requirement**, not an ALT/BIR-specific
one — so within a uniformly replication-stressed ALT+ panel POLD3 gene effect should be
flat, which is the low-variance condition that mechanically suppresses correlation.

### 2d. The SMARCAL1-polarity problem

- Liu et al. 2023 (`10.1093/neuonc/noad022`): inducible **rescue** of SMARCAL1 in a
  natively SMARCAL1-deficient ALT+ line suppressed ALT indicators and blocked de novo
  telomere synthesis in GBM **and osteosarcoma** cells. SMARCAL1 presence *suppresses*
  ALT — the opposite polarity to the BIR module's assumption.
- Brosnan-Cashman et al. 2021 (`10.1038/s41379-021-00841-7`): of 19 ALT+ giant-cell GBM
  cases, **6/19 (31.5%)** showed loss or mislocalisation of SMARCAL1. ALT+ cells can
  lack SMARCAL1 entirely.
- Sobieszek & Idzik 2026 (`10.36922/GPD025480078`): **germline biallelic** pathogenic
  SMARCAL1 variants are viable in humans (Schimke immuno-osseous dysplasia). SMARCAL1
  is not constitutively essential, so a near-zero SMARCAL1 gene effect in some lines is
  expected rather than informative.
- Cox, Maréchal & Flynn 2016 (`10.1016/j.celrep.2016.01.011`) supports the *other*
  sign: SMARCAL1 resolves replication stress at ALT telomeres and protects them. Both
  signs are supported; the polarity is context-dependent, which alone makes a
  single-sign prediction fragile.
- Taglialatela et al. 2026 (`10.1101/gad.353377.125`) adds a further decoupling: the
  SMARCAL1–FANCM synthetic-lethal interaction is also present in **non-transformed,
  ALT-negative** epithelial cells, pointing to "a redundant essential function of these
  enzymes independent of ALT". SMARCAL1 gene effect is not a clean BIR-reliance readout.

### 2e. Collateral contradictions already logged

Prespec Conflicts #1, #3 and #4 in `Literature/INDEX.md` all bear on this iteration:
BRCA2 loss itself drives ALT/BIR (Lee 2024, `10.1093/nar/gkae251`); and RAD51, BRCA1
and BRCA2 actively *suppress* the POLQ/LIG3/PRIMPOL restart routes while BIR itself is
MUS81–RAD52–POLD3 (Batenburg 2023, `10.3390/ijms241512419`) — so the two arms of the
Stage 1 composite are coupled with **opposite sign**, not orthogonal.

**Verdict on Q2.** The premise is contradicted, and the contradiction is directional:
the literature predicts ρ(RAD52, SMARCAL1) ≤ 0 and ρ(POLD3, ·) ≈ 0.

---

## 3. Has it already been done?

**The method is published and standard. This specific application is not. The iteration
falls on the permitted side of the no-duplication rule.**

### 3a. The method is established

- Pan et al. 2018, Cell Systems (`10.1016/j.cels.2018.04.011`) — "Interrogation of
  Mammalian Protein Complex Structure, Function, and Membership Using Genome-Scale
  Fitness Screens". The founding demonstration that correlated gene-effect profiles
  recover complex membership.
- Wainberg et al. 2021, Nat Genet (`10.1038/s41588-021-00840-z`) — genome-wide atlas of
  co-essential modules; recapitulates diverse pathways and protein complexes and
  assigns function to >100 uncharacterised genes. Also the canonical source for the
  false-positive problem in naive Pearson co-essentiality.
- Shimada et al. 2021, eLife (shinyDepMap) (`10.7554/eLife.57116`) — shinyDepMap;
  clusters genes with strong efficacy scores into functional units, many of which are
  complexes or pathways.
- Chen et al. 2023, Nat Commun (`10.1038/s41467-023-39990-5`) — "Context-defined cancer
  co-dependency mapping". This is the closest published *methodological* precedent to
  what iteration 03 proposes, i.e. co-dependency restricted to a biological context.
  It is also the paper Taglialatela 2026 cites for its DepMap dependency-score
  convention.

Implication: the loop cannot claim methodological novelty, and should not try to. It
should cite Pan 2018 / Wainberg 2021 for the estimator and Chen 2023 for the
context-restricted variant.

### 3b. The specific result is not published

Two WebSearch queries returned no source reporting the co-dependency structure of
POLD3, RAD52 and SMARCAL1 — pan-cancer or in ALT+ lines. Nothing in the 87-PDF corpus
does either. Negative evidence from the corpus is strong: Bennett 2026
(`10.64898/2026.03.04.709438`) mentions POLD3/RAD52/BIR 0 times, and INDEX.md Gap #4
records that **across all 87 readable PDFs no paper substantively mentions both
TCAB1/WRAP53 and break-induced replication**.

### 3c. The nearest published work, and why it is not the same analysis

Taglialatela et al. 2026 (`10.1101/gad.353377.125`) is the closest. Read from the PDF:
they used DepMap 23Q4 gene-effect scores for 1009 pan-cancer lines, ran a differential
dependency analysis stratified by hTERT expression, repeated it **within osteosarcoma
(13 lines)**, and then "performed hierarchical clustering on gene dependencies of cell
lines" (R `pheatmap`, `ward.D2`, euclidean, `cutree_cols = 2`) over the **top 20
differential dependency genes**, which "successfully stratified osteosarcoma cell lines
into two groups with largely distinct ALT/telomerase status".

Three reasons this is a different analysis:

1. **Wrong axis.** They cluster *cell lines* by their dependency profile over a gene
   set. Iteration 03 proposes to correlate *genes* by their gene-effect profile over a
   line set. Co-dependency is the gene-by-gene axis.
2. **Wrong gene set.** Their top-20 set is driven by SMARCAL1 and FANCM, with HIRA,
   RMI1 and TRIM28 named as corroborating ALT dependencies. POLD3 and RAD52 do not
   appear among the reported differential dependencies in osteosarcoma.
3. **Wrong stratification.** Theirs is ALT+ vs ALT−/telomerase+. Iteration 03 is
   *within* ALT+.

Prior ALT-dependency analyses that Taglialatela 2026 cites — Lu & Pickett 2022,
Salgado et al. 2024, Dubois et al. 2025, Wu et al. 2025
(`10.1038/s41467-025-67190-w`) — are all ALT+ vs ALT− differential-dependency analyses,
not within-ALT+ co-dependency analyses.

### 3d. Which side of the line, and the osteosarcoma-specific justification

**Permitted side.** The iteration is (i) a methodological post-mortem on *the loop's own
Stage 1 endpoint*, which no external paper has any reason to have tested, and (ii) an
osteosarcoma-restricted application. The osteosarcoma-specific justification is
explicit and threefold:

- Osteosarcoma is the ALT-enriched tumour type at issue: >50% of osteosarcomas are ALT+
  (Taglialatela 2026, `10.1101/gad.353377.125`); 30–50% per Broderick 2023
  (`10.1038/s41467-023-38029-z`).
- The endpoint being post-mortemed was defined on 16 ALT+ DepMap **osteosarcoma** lines.
  The question "was this composite legitimate in this panel" is only answerable in this
  panel.
- INDEX.md Gap #2 records that the corpus contains **no osteosarcoma-specific ALT
  heterogeneity source at all** — every readable ALT-heterogeneity paper is glioma
  (Diplas 2018 `10.1038/s41467-018-04448-6`; Liu 2023 `10.1093/neuonc/noad022`;
  Brown 2026 `10.1093/neuonc/noaf300`; Brosnan-Cashman 2021
  `10.1038/s41379-021-00841-7`) or pan-cancer, and
  `10.1016/j.cancergen.2026.06.001` ("Beyond ATRX: the heterogeneous basis of ALT in
  osteosarcoma") could not be retrieved.

---

## 4. What is the right comparator?

The coherence statistic must be judged against **two** references, and the second is
not optional: with n ≈ 13–16 ALT+ osteosarcoma lines, a null result for the triplet is
uninterpretable unless a known-coherent set is shown to be *detectable in the same n*.
Otherwise "the premise fails" and "the design is underpowered" are indistinguishable.

### 4a. Negative reference — the genome-wide null

**Random gene triplets drawn from the same gene-effect matrix, matched on per-gene
gene-effect variance and mean within the same ALT+ OS lines.** Variance matching is
required, not cosmetic: correlation magnitude scales with dynamic range, and Wainberg
et al. 2021 (`10.1038/s41588-021-00840-z`) identifies exactly this class of artefact —
ubiquitous false positives driven by profile shape rather than function — as the
central obstacle in naive co-essentiality. An unmatched random null would make a
low-variance triplet look non-coherent for purely mechanical reasons, which would
contaminate the post-mortem's conclusion.

### 4b. Positive reference — recommended primary: the BTR dissolvasome

**BLM, TOP3A, RMI1, RMI2.** This is the best-justified positive control available:

- It is an obligate stoichiometric complex, which is the regime where correlation-based
  co-essentiality actually works (Pan et al. 2018, `10.1016/j.cels.2018.04.011`).
- It is ALT-relevant in the *same* biology, so it is not an off-axis control:
  Broderick 2023 (`10.1038/s41467-023-38029-z`) reports that **both** ALT BIR branches
  rely on BLM, and that EXD2 loss is synthetic sick with BLM.
- Its essentiality is *selective* in ALT, so it should have usable variance in an ALT+
  panel: Taglialatela 2026 (`10.1101/gad.353377.125`) names **RMI1** among the top
  significant differential dependencies in telomerase-negative cells.

### 4c. Positive reference — recommended secondary: the Pol δ holoenzyme

**POLD1, POLD2, POLD3 (± POLD4).** This is the sharpest possible test of "if anything
co-depends, these do", and it has the specific virtue of containing POLD3 itself: if
POLD3 does not correlate with POLD1/POLD2, then POLD3's gene-effect profile carries no
recoverable pathway signal in this panel and the Stage 1 composite was unsalvageable
for a reason that has nothing to do with ALT. Justified by Wu et al. 2023
(`10.1038/s41467-023-35992-5`), which treats Pol δ as acting as a unit in MiDAS, and by
Broderick 2023's "PCNA-RFC-Polδ" framing (`10.1038/s41467-023-38029-z`).

**Caveat to pre-specify:** POLD1 and POLD2 are near pan-essential, so their within-panel
variance may be too small for correlation to be estimable. If both the Pol δ and BTR
controls fail, the correct reading is *underpowered*, not *premise refuted*. State this
as an explicit interpretive rule in the prespec.

### 4d. The HR set — re-purpose it, do not use it as a neutral comparator

The hard constraint is wider than BRCA2. **No member of {BRCA1, BRCA2, RAD51} is
ALT-neutral:**

- **BRCA2 must be excluded** from any HR-neutral control set: BRCA2 depletion
  stabilises telomeric G4s, accumulates TERRA R-loops, drives LLPS and APB assembly,
  and a subset of BRCA2-null cells bypasses senescence into BIR-mediated telomere
  synthesis (Lee et al. 2024, `10.1093/nar/gkae251`). BRCA2 loss can phenocopy ALT
  induction. It was also one of only two genes with a negative δ in the Stage 1 result.
- **RAD51 is not neutral either**: Broderick 2023 (`10.1038/s41467-023-38029-z`) reports
  the RAD52-independent ALT route maintains telomeres independently of RAD51 and MRE11
  and that ablating these HR factors is associated with **C-circle overproduction**;
  Carson & Flynn 2023 (`10.1016/j.coph.2023.102380`) has RAD51 and RAD52 regulating
  distinct ALT branches.
- **BRCA1 is not neutral either**: Batenburg 2023 (`10.3390/ijms241512419`) has BRCA1
  actively suppressing an alternative restart route.

**Recommendation.** Do not use BRCA1/BRCA2/RAD51 as an HR-neutral comparator — that was
Stage 1's error (Prespec Conflict #3) and repeating it would reproduce the defect under
investigation. Instead use the HR triplet as a **second module whose internal coherence
is tested in parallel**: BRCA1/BRCA2/RAD51 are a genuine functional HR unit, so if the
estimator has power in this panel they should co-depend. This converts the HR arm from
an invalid neutral baseline into a third informative positive control, at no extra
multiplicity cost beyond one test. If the HR triplet is coherent and the BIR triplet is
not, the post-mortem conclusion is unusually clean.

### 4e. One composition note

If a BIR module is ever re-specified rather than post-mortemed, the literature's module
is **MUS81–RAD52–POLD3** (Batenburg 2023, `10.3390/ijms241512419`), with SMARCAL1
belonging to the RAD52-*independent* branch alongside MUS81 (Broderick 2023,
`10.1038/s41467-023-38029-z`). MUS81 should be reported descriptively in this iteration
so the record shows what the correctly-composed module would have looked like. It must
not enter the primary endpoint, which is fixed on the three Stage 1 genes.

---

## 5. Is there a defensible type I vs type II ALT subtype marker from DepMap 26Q1 expression or copy number alone?

**No. There is no defensible type I vs type II marker, and the nomenclature does not
transfer to human cells.** This is a straight answer, and the prespec should take the
explicit-statement option under PROGRAM.md rule (c).

### 5a. The type I/II distinction is not a human distinction

- The only paper in scope that directly profiles the subtypes is Lee et al. 2026,
  Mol Cells (`10.1016/j.mocell.2026.100405`), flagged in INDEX.md as the corpus'
  highest-value gap (main text unretrievable; five SI files present). Its SI, read this
  scan, settles the question: the type I/type II comparison is performed in **mouse
  embryonic stem cells** ("type I and type II ALT mESCs"), type I ALT being defined by
  an ~7 kb **mTALT** mouse telomere-associated template. Supplementary Table 5 reports
  a survey of human GRCh38 subtelomeres for mTALT-like regions and concludes: "**No
  compact (~7 kb) human genomic segment with an organization closely resembling the
  mTALT template was identified.**" The paper's own data say type I ALT as defined has
  no established human counterpart.
- The type I/II nomenclature originates in *S. cerevisiae* (type I = Y′ subtelomeric
  amplification; type II = long heterogeneous terminal tracts). Broderick 2023
  (`10.1038/s41467-023-38029-z`) warns directly against importing it: yeast
  "mechanisms may not be exactly analogous to RAD52-dependent and -independent ALT in
  humans".

### 5b. The three human "two-mode" models carve ALT on different axes, none measurable from expression or CN

- Broderick 2023 (`10.1038/s41467-023-38029-z`): RAD52-dependent vs RAD52-independent
  BIR, switched by **EXD2 fork-processing activity** and read out by RAD52 *telomere
  association*, C-circles, T-SCE and segregated CO-FISH. None of these is an mRNA or
  copy-number quantity; EXD2 transcript level is not the switch.
- Zhang 2026, BioEssays (`10.1002/bies.70187`): two kinetically distinct APB formation
  pathways modulating ALT strength. Requires condensate imaging — excluded by the data
  constraint.
- Folini 2026 (`10.21203/rs.3.rs-11116403/v1`): two complementary replication-coupled
  C-circle biogenesis pathways. Requires the C-circle assay — excluded by the data
  constraint.

### 5c. What expression/CN can and cannot give

- Published expression-based ALT classifiers exist but classify **ALT status, not
  subtype**: the 28-gene signature of Armendáriz-Castillo et al. 2026
  (`10.3390/biology15171531`) and the ALTitude tabular model of Bennett et al. 2026
  (`10.64898/2026.03.04.709438`, WGS-based). Wu et al. 2025
  (`10.1038/s41467-025-67190-w`) partitions telomere-maintenance classes
  (ALT / TEL / TMM-double-negative / TMM-double-positive) — a different and coarser
  axis than type I/II, and also not derivable from expression or CN alone.
- The only literature-grounded partition that *is* obtainable from expression or CN is
  **ATRX/DAXX-driven vs SMARCAL1-driven ALT** (Diplas et al. 2018,
  `10.1038/s41467-018-04448-6`, establishing SMARCAL1 inactivation as an
  ATRX-independent route to ALT; Brown et al. 2026, `10.1093/neuonc/noaf300`,
  SMARCAL1 synthetic lethality specifically in **ATRX-deficient** ALT+ glioma;
  Brosnan-Cashman 2021, `10.1038/s41379-021-00841-7`, only 9/19 ALT+ cases ATRX-lost
  and 6/19 SMARCAL1-lost). This is **not** type I/II, it is glioma-derived, and the
  loop has already established it is not estimable here: Prespec Conflict #6 records
  that only **1 of 16** ALT+ OS lines is ATRX-altered, and iteration 02 closed with
  "H2 not estimable — gene-level ATRX copy number is unreliable".
- Richards 2023 (`10.1038/s41467-023-40071-w`) supplies the general reason this is not
  a retrievable omission: the BIR switch it identifies is **post-translational**
  (site-specific ADP-ribosylation of PolD3). A transcript or copy-number marker cannot
  see it in principle.

### 5d. Recommended prespec language

Declare **no defensible type I vs type II ALT subtype marker exists** in DepMap 26Q1
expression or copy number, citing Lee 2026 (`10.1016/j.mocell.2026.100405`, mouse
system, no human mTALT analogue), Broderick 2023 (`10.1038/s41467-023-38029-z`, yeast
analogy warning; switch is an activity not a transcript) and Richards 2023
(`10.1038/s41467-023-40071-w`, post-translational switch). If a secondary split is
wanted for completeness, report ATRX/DAXX-low vs intact (Diplas 2018,
`10.1038/s41467-018-04448-6`; Brown 2026, `10.1093/neuonc/noaf300`) as a **descriptive,
non-inferential** secondary with an explicit non-estimability note carried forward from
iteration 02. Do not label it a subtype marker.

---

## Bottom line for the prespec

**(i) Worth an iteration: yes.** The module premise has no co-dependency support
anywhere — Bennett 2026 (`10.64898/2026.03.04.709438`) never mentions POLD3/RAD52/BIR,
and no paper names the triplet as one module. It is actively contradicted with a sign:
Broderick 2023 (`10.1038/s41467-023-38029-z`) puts SMARCAL1 in the RAD52-*independent*
BIR branch, so RAD52 and SMARCAL1 mark a switch (ρ ≤ 0), while POLD3 is the shared
downstream polymerase and a general replication-stress reporter (Wu 2023,
`10.1038/s41467-023-35992-5`), so ρ(POLD3,·) ≈ 0. Method is established (Pan 2018
`10.1016/j.cels.2018.04.011`; Wainberg 2021 `10.1038/s41588-021-00840-z`; Chen 2023
`10.1038/s41467-023-39990-5`); this application is not. Osteosarcoma-specific and
self-referential, so not duplication.

**(ii) Controls.** Negative: variance- and mean-matched random triplets. Positive:
BTR (BLM/TOP3A/RMI1/RMI2) primary; Pol δ (POLD1/POLD2/POLD3) secondary; HR
(BRCA1/BRCA2/RAD51) tested as a *module*, never as a neutral comparator — BRCA2 drives
ALT (Lee 2024, `10.1093/nar/gkae251`) and RAD51 is branch-specific.

**(iii) Subtype marker: none defensible** (Lee 2026, `10.1016/j.mocell.2026.100405` is
mouse; no human mTALT analogue). Invoke rule (c)'s explicit-statement option.

**(iv) Constraints.** Keep predicted sign **+1** on positive co-dependency. Pre-commit
the rule that if all positive controls also fail, the verdict is underpowered, not
refuted.
