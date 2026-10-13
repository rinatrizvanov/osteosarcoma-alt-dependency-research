# Iteration 04 — literature scan: what should predict SMARCAL1 dependency *within* ALT+ lines?

**Scope discipline.** Started from `Literature/INDEX.md`. Pulled **two** PDFs: Taglialatela
et al. 2026 Genes Dev (`10.1101/gad.353377.125`) and Bennett et al. 2026 ALTitude
(`10.64898/2026.03.04.709438`), both via `pdftotext`. Two WebSearches, used only for
question 4. Brown 2026 (`10.1093/neuonc/noaf300`) was assessed from `INDEX.md` entry 7,
its abstract, and search results, not from its PDF.

Nothing proposed here appears in `autoresearch/NULL_SPACE.md`. The ATRX stratifier,
gene-level copy number without segment validation, type I/II ALT subtyping, WRAP53 mRNA,
the BIR triplet's coherence, and any ALT+ vs ALT− group contrast on the 38 disclosed genes
(`DISCLOSURES.md` §1) are all treated as closed.

---

## 1. Taglialatela et al. 2026 — the mechanism, stated literally

**Taglialatela A, Lee J, Azeroglu B, et al. "SMARCAL1 is a candidate therapeutic target for
ALT-positive tumors." Genes Dev 2026;40(15–16):1283–1302. `10.1101/gad.353377.125`**
(preprint `10.64898/2026.02.15.704070`).

### 1.1 What the paper claims generates the dependency

The abstract states the mechanism in two clauses:

> "Mechanistically, we showed that SMARCAL1 loss leads to telomeric ssDNA accumulation in
> ALT-positive cells, dependent in part on DNA repriming mediated by the DNA
> primase/polymerase PRIMPOL. Moreover, SMARCAL1's ssDNA annealing activity counteracts
> DNA unwinding by the BLM helicase, limiting telomeric ssDNA accumulation and DNA damage
> in ALT-positive cells."

The summary sentence of the results section is the one the prespec's sign should rest on:

> "Taken together, these findings suggest that BLM facilitates the formation of ssDNA
> intermediates at telomeres of ALT-positive cells, which in turn promote SMARCAL1
> recruitment. By inhibiting PRIMPOL-mediated repriming and counteracting BLM-dependent DNA
> unwinding, SMARCAL1 could prevent excessive formation of ssDNA intermediates that would
> otherwise aggravate telomere instability and activate aberrant DDR signaling, resulting in
> an impairment of cancer cell fitness."

And the discussion's causal reading:

> "The above observations indicate that unscheduled BLM unwinding activity, licensed at
> telomeres of ALT-positive cells by frequent DNA polymerase stalling, may generate the DNA
> intermediates that trigger ALT. These same DNA structures could also serve as substrates
> for SMARCAL1, which in turn opposes BLM to restrain ALT."
>
> "In this context, SMARCAL1 may act as a 'brake' to restrain ALT activity and maintain
> tolerable levels of telomeric stress."

**So the proximate cause of the dependency is the rate at which RPA-coated telomeric ssDNA
is generated.** SMARCAL1 is the annealing activity that removes it. A cell dies when
SMARCAL1 is lost because the ssDNA supply exceeds what the cell can tolerate without
SMARCAL1.

### 1.2 Therefore: what must a line have **more of** to be **more** SMARCAL1-dependent?

**More BLM-driven telomeric unwinding** (primary), and secondarily **more PRIMPOL repriming**.
The direction is not inferred — the paper manipulates it in both directions:

> "loss of BLM, but not PRIMPOL, significantly impaired SMARCAL1 enrichment at telomeres in
> ALT-positive cells"
>
> "PRIMPOL loss partially mitigated the cell fitness defects induced by SMARCAL1 loss in
> SAOS2 and G292 cells, while **BLM deficiency fully suppressed the SMARCAL1-dependent
> growth defects** in those cells"
>
> "overexpression of WT BLM, but not the KA mutant, **markedly exacerbated** DNA damage,
> telomeric RPA accumulation, APB formation and fitness defects in SMARCAL1-depleted cells"

And, critically for the sign and for the interaction structure:

> "consistent with a protective role for SMARCAL1 in counteracting BLM-helicase-driven
> telomeric stress, **overexpression of WT BLM in control G292 cells did not induce
> appreciable increases** in telomeric RPA localization, APB formation, or cell fitness
> impairment"

**Direction, stated so the prespec cannot get the sign wrong.** Raising BLM is harmless on
its own but worsens fitness *specifically when SMARCAL1 is absent*. Removing BLM abolishes
the fitness cost of losing SMARCAL1 entirely. A Chronos gene effect for SMARCAL1 **is**
"fitness under SMARCAL1 loss", so Taglialatela's BLM experiment and the Chronos endpoint
are measuring the same quantity on the same axis. More BLM activity → larger fitness cost
of SMARCAL1 loss → **more negative SMARCAL1 Chronos**.

Three independent channels in the paper all push the same way, which is the main reason to
trust this sign:

1. **Substrate supply.** BLM unwinding generates the RPA-coated ssDNA that SMARCAL1 anneals
   (in vitro, topoisomerase-coupled plasmid assay: "addition of SMARCAL1 to the same
   reaction suppressed the formation of negatively supercoiled DNA by BLM, revealing
   SMARCAL1's ability to counteract BLM's DNA unwinding activity … through the annealing of
   RPA-coated ssDNA bubbles"; mechanism established in Halder et al. 2022,
   `10.1093/nar/gkac583`).
2. **Recruitment.** BLM loss abolishes SMARCAL1's telomeric enrichment, and BLM
   overexpression "significantly enhanced SMARCAL1 localization at ALT telomeres". More BLM
   → more SMARCAL1 at telomeres → SMARCAL1 is more load-bearing.
3. **ALT activity.** BLM is required for ALT itself, and "suppressing ALT abrogates their
   dependency on SMARCAL1". More BLM → more ALT → more dependency.

There is no channel in the paper that runs the other way. (The only conceivable
counter-channel — a line with so little BLM that ALT is weak — also predicts *less*
dependency, i.e. the same sign.)

### 1.3 What the paper does **not** claim

- It does **not** claim SMARCAL1 acts in break-induced replication. It argues the opposite:
  "SMARCAL1 is not enriched at telomeres in either ALT-positive or ALT-negative cells
  following induction of site-directed DNA double-strand breaks (DSBs), suggesting that its
  primary role is preventing telomeric replication stress rather than facilitating the
  repair of DSBs through BIR or other DSB repair pathways" (citing Zhang et al. 2023). This
  is a further reason the Stage 1 BIR framing was mis-aimed, independent of iteration 03.
- It does **not** claim PRIMPOL is sufficient: "depletion of endogenous PRIMPOL only
  partially suppressed RPA2 phosphorylation and telomeric ssDNA accumulation … suggesting
  that additional PRIMPOL-independent mechanisms contribute".
- It does **not** report any cross-line analysis of what makes one ALT+ line more
  SMARCAL1-dependent than another (see §4).

---

## 2. The genes the mechanism names, and the predicted direction for each

Direction convention: "**more**" = higher activity of that gene makes the cell **more**
SMARCAL1-dependent = **more negative** SMARCAL1 Chronos.

| Gene | Named by Taglialatela `10.1101/gad.353377.125`? | Higher activity → | Confidence | Basis |
|---|---|---|---|---|
| **BLM** | Yes, centrally; in abstract, results heading, Fig. 6, discussion | **MORE dependent** | **High** | Bidirectional perturbation: Cas12a/siRNA KO *fully* suppresses the SMARCAL1-loss growth defect; WT-cDNA overexpression (but not the K695A helicase-dead mutant) markedly exacerbates it; plus direct biochemistry and a recruitment relationship. BLM OE alone is phenotypically silent, so the effect is specifically an interaction with SMARCAL1 loss. |
| **PRIMPOL** | Yes; abstract, results heading, Fig. 5 | **MORE dependent** | **Medium** | WT-PRIMPOL overexpression exacerbates fitness defects and RPA2 phosphorylation; repriming-dead (CH), catalytically dead (AxA) and RPA-binding (RBD) mutants do not. But endogenous depletion gives only partial rescue, the paper invokes "additional PRIMPOL-independent mechanisms", and explicitly cannot exclude "nontelomeric functions of PRIMPOL". Mechanism from Tirman et al. 2021 `10.1080/10409238.2020.1841089`; González-Acosta et al. 2021 `10.15252/embj.2020106355`. |
| **RPA1 / RPA2 / RPA3** | Yes, but as SMARCAL1's partner and as the *readout* (phospho-RPA2 S4/S8, telomeric RPA32 foci), never as an exposure | **Direction undetermined** | **Low — do not use** | Two opposing channels. RPA is required for SMARCAL1 function (ΔN1–115 and ΔN1–30 RPA-binding mutants are impaired), so more RPA → more effective SMARCAL1 → *more* dependency. But RPA also coats and stabilises the ssDNA bubbles, i.e. it buffers the lesion → *less* dependency. The paper gives no direction. Reviewed in Dueva & Iliakis 2020 `10.1093/narcan/zcaa022`. Also pan-essential, so expression variance across 16 lines will be minimal. |
| **POLA1 / Polα-primase** | **No — never named** | **No direction available** | **n/a — do not use** | Taglialatela's only primase is PRIMPOL. The Polα-primase fill-in mechanism belongs to the BIR-at-DSBs literature (Shah et al. 2024 `10.1038/s41467-024-52916-z`; Wu et al. 2025 `10.1016/j.celrep.2025.116386`), which Taglialatela explicitly argues is *not* SMARCAL1's role (§1.3). Using POLA1 would import a mechanism this paper rejects. |
| **TOP3A** | Yes, once — only as an "ALT inducer" whose hyperactivation can establish ALT | **MORE dependent** (by BTR analogy) | **Medium-low** | "Less commonly, ALT-positive tumors arise from hyperactivation of ALT inducers, such as TOP3A", citing de Nonneville et al. 2022 `10.15252/emmm.202215859` (TOP3A amplification in pediatric osteosarcoma, mutually exclusive with ATRX inactivation). As a BTR subunit it should act with BLM, but Taglialatela never tests it, and de Nonneville's claim is at the **copy-number** level — the exact measurement class `NULL_SPACE.md` requires segment validation for. |
| **RMI1** | Yes, once — as a corroborating *top differential dependency gene* in their own screen (Loe et al. 2020 `10.1101/gad.333963.119`) | MORE dependent (by BTR analogy only) | **Low** | Named as a validated dependency, not as a source of the SMARCAL1 dependency. No direction tested. |
| **RMI2** | No | — | — | Not named. |
| **PML** | Yes, but only as the APB marker in foci quantification | MORE dependent (indirect) | **Low** | PML is not a mechanism source here. Loe et al. 2020 `10.1101/gad.333963.119` makes PML the determinant of BTR *localisation* to telomeres, which would make PML an activity-licensing proxy for BLM — mechanistically attractive but untested by Taglialatela and one inferential step further out. |
| **FANCM** | Yes, repeatedly | **Sign ambiguous** | **Medium on mechanism, low on sign** | FANCM is an ALT suppressor whose depletion "further increases telomere damage upon SMARCAL1 loss", so higher FANCM → less ALT → *less* dependency. But the paper also records "a strong synthetic lethal interaction between SMARCAL1 and FANCM in nontransformed, ALT-negative epithelial cells … point[ing] to a redundant essential function of these enzymes independent of ALT that remains to be fully defined" — which predicts the opposite sign on a fitness readout. Two channels, opposite signs, unresolved. Unsuitable as a signed primary. |
| **ATRX / DAXX** | Yes, extensively | (claimed synthetic lethal) | **Closed** | Non-estimable on this cohort (iteration 02), and published as a within-ALT+ null by Bennett (§4). Do not use. |
| **ZRANB3** | Yes — as the **negative** control | **No effect** | **High (for the null)** | ZRANB3 loss causes no telomeric damage or fitness defect in ALT+ cells, and "ZRANB3 was unable to counteract BLM's activity … reinforcing the unique antagonistic functional relationship between BLM and SMARCAL1". This makes ZRANB3 a pre-specifiable **negative-control exposure** (see §3.4). |
| **HLTF** | Yes, briefly | lower HLTF → more dependent (speculative) | **Low** | HLTF loss increases telomeric damage and ALT in U2OS, but Taglialatela states "it remains to be determined whether their loss can exacerbate the dependency of ALT-positive cancer cells on SMARCAL1". Explicitly open. |
| **TERT (hTERT) expression** | Yes — "a robust inverse correlation between hTERT expression and SMARCAL1 dependency" | higher hTERT → **LESS** dependent | **High but useless here** | This *is* the ALT-status axis, already published as the headline finding, and `NULL_SPACE.md` records that `telomerase_status` is constant (`indeterminate`) within the project's ALT+ set. |

**The mechanism names exactly two genes as the source of the dependency: BLM and PRIMPOL.**
Everything else in the candidate list is either a readout (RPA, PML), a dependency rather
than an exposure (RMI1), imported from a different mechanism (POLA1), sign-ambiguous
(FANCM), or closed (ATRX/DAXX). Of the two, **BLM is the stronger by the paper's own
arithmetic**: BLM loss *fully* suppresses the phenotype where PRIMPOL loss only *partially*
does, and only BLM loss affects SMARCAL1's telomeric recruitment.

---

## 3. Is expression the right measurement? — the honest answer

Verdict: **for BLM, mRNA is a defensible proxy — materially better grounded than the
project's two previous proxy failures — but it is a proxy, it should be labelled one, and
its adequacy should be gated before the outcome is touched.** It is not the same mistake,
but it is adjacent to it, and the prespec should carry a falsifiable guard rather than an
assumption.

### 3.1 Why it is *not* the WRAP53 mistake

The two previous failures were of a different kind.

- **WRAP53 mRNA → 17p13.1 structural disruption** (`NULL_SPACE.md` H1, H5) asked a transcript
  level to report a **genomic rearrangement**, a quantity with no monotone relation to
  transcript abundance. `NULL_SPACE.md` records that axis G of the OpenAlex scan found **no
  corpus at all** supporting that substitution. That was a cross-modality leap.
- **Gene-level copy number** (`NULL_SPACE.md`, ATRX) was not a proxy failure at all; it was a
  **measurement error** — a derived gene-level value contradicted by the underlying segments
  in 3 of 4 calls.

Here the mechanism's own causal perturbation is a **level change in the same molecule**.
Taglialatela raised BLM by WT-cDNA transduction and the dependency got worse; lowered it by
Cas12a and siRNA and the dependency disappeared. Abundance → phenotype is the paper's
experimental axis, not an inference the iteration adds on top. That is a categorically
different epistemic position from asking mRNA to stand in for a structural variant.

### 3.2 Empirical support that BLM mRNA carries signal *in this exact data modality*

Bennett et al. 2026 (`10.64898/2026.03.04.709438`) derived an ALT transcriptional signature
from DepMap RNA-seq and **BLM is one of the named ALT_up genes** — upregulated in ALT+
lines:

> "These 2 resulting gene sets included well-established ALT-associated factors, known to
> also play a role in replicative stress and DNA damage response, including **BLM**, BRCA2,
> RAD51, XRCC3, and FEN1 in ALT_up, as well as genes expected to be depleted in ALT contexts
> such as ATRX and TERT, as well as TP53 and RB1 in ALT_down."

And BLM transcript moves ALT-specifically after SMARCAL1 knockout: "Among the genes showing
the most ALT-specific changes were **BLM** (G292: +0.33; Huh7: −0.52) … (all P < 0.05)". So
BLM mRNA in DepMap RNA-seq already tracks ALT biology at group level and responds to
SMARCAL1 loss in the ALT-specific direction. That is exactly the kind of modality-level
validation axis G could never find for WRAP53.

### 3.3 Why it is still a proxy — stated plainly, because this is the main threat

1. **The mechanism requires *unscheduled helicase activity at telomeres*, not abundance.**
   Taglialatela: "unscheduled BLM unwinding activity, **licensed** at telomeres of
   ALT-positive cells by frequent DNA polymerase stalling". Activity is a product of at
   least three terms — abundance, telomeric recruitment (BTR assembly; PML-dependent
   localisation, Loe et al. 2020 `10.1101/gad.333963.119`), and the frequency of polymerase
   stalling that licenses it. mRNA measures one of the three.
2. **The K695A result shows the relevant quantity is specifically catalytic.** No transcript
   measurement reports catalytic competence.
3. **The dose–response was demonstrated by overexpression, far outside the physiological
   range.** Whether the narrower endogenous range that distinguishes 16 cell lines produces
   any measurable effect is untested, by anyone.
4. **A published mRNA↔protein dissociation exists for the central gene in this question.**
   Bennett's one within-ALT+ explanation is exactly this failure mode, for SMARCAL1 itself:
   > "we observed that while the **transcript level of SMARCAL1 is comparable by RNA-seq to
   > other DepMap cell lines, it expresses very low SMARCAL1 protein** … Thus, we hypothesize
   > that this low protein level causes ALT in SKNFI, explaining why these cell lines are
   > resistant to SMARCAL1 KO, despite being ALT+."

   This is decisive for one design choice: **SMARCAL1 mRNA must not appear anywhere in this
   iteration's exposure**, as primary, alternate or covariate. It is also a genuine (n = 1)
   warning about transcript proxies for SNF2-family translocases and helicases generally.

### 3.4 What the prespec should do about it instead of assuming

Three cheap, pre-registerable guards turn §3.3 from an unexamined assumption into something
this iteration partly answers:

- **An exposure feasibility gate, evaluated before the outcome is read.** Report BLM
  log2(TPM+1) range, SD, and bimodality over the 16 ALT+ lines *first*. If the spread is not
  materially above measurement noise, declare non-estimable at the gate rather than running
  a correlation on a near-constant variable. This is the WRAP53 lesson in its correct form:
  `NULL_SPACE.md` records that iteration 01's minimum-cell guard fired because the exposure
  had no natural low group, and that was discoverable blind.
- **A pre-registered negative-control exposure: ZRANB3 mRNA.** Taglialatela establishes that
  ZRANB3 cannot substitute for SMARCAL1 and cannot counteract BLM, so ZRANB3 expression is
  predicted to carry **no** association with SMARCAL1 Chronos. If ZRANB3 mRNA predicts
  SMARCAL1 dependency as well as BLM mRNA does, the result is an artefact of expression-wide
  structure, not mechanism.
- **A proliferation/screen-amplitude confounder control.** BLM is a replication gene;
  fast-cycling, replication-stressed lines plausibly have both higher BLM mRNA and larger
  absolute Chronos effects across many genes. Pre-register (a) the same test on within-line
  standardised Chronos — iteration 03 already has this machinery — and (b) a small set of
  negative-control *outcomes* (mid-range genes unrelated to telomere metabolism), to show
  BLM mRNA does not predict everything.

And one prohibition carried forward from `NULL_SPACE.md` lesson 2: **do not compose a
multi-gene "ssDNA-generation score"** from BLM + PRIMPOL + TOP3A. Marginal regimes must be
checked before anything is composed. One gene, one exposure.

### 3.5 The measurement that would actually be right, and why it is not proposed

The defensible measurement is BLM **protein** (or, better, telomeric BLM foci). DepMap ships
CCLE proteomics, but it is not in the frozen raw set
(`data/raw/depmap/` holds `OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv`,
`CRISPRGeneEffect.csv`, `OmicsCNGeneWGS.csv`, `OmicsCNSegmentsWGS.csv`, `Model.csv` only),
osteosarcoma coverage in CCLE proteomics is partial, and adding a data source is a scope
decision for the user rather than for the loop. **Flagging it, not silently substituting
mRNA for it, is the point of this section.** If the user wants the stronger measurement,
that is the redirection this scan supports.

---

## 4. Has this specific question been published?

**No. The cross-line, within-ALT+ form of the question is unpublished, and one paper
explicitly flags it as mechanistically unexplained.**

### 4.1 Taglialatela 2026 (`10.1101/gad.353377.125`) — no within-ALT+ predictor

Their DepMap work is entirely a **between-group** contrast: "we initially classified >1300
cell lines … as telomerase-positive or -negative, based on hTERT expression", then "a
differential analysis of gene dependencies within a single tumor type (osteosarcoma)
stratified by hTERT expression" — Fig. 1A pan-cancer (1009 lines), Fig. 1B osteosarcoma (13
lines). The only continuous relationship they report is "a robust inverse correlation
between hTERT expression and SMARCAL1 dependency", i.e. the ALT-status axis itself.

Their BLM and PRIMPOL results are **isogenic perturbations within single lines** (G292,
SAOS2, G292iDAXX, OS384) — never a correlation of endogenous BLM or PRIMPOL level with
dependency magnitude across lines. The cross-line test is a natural corollary of their
mechanism that they did not run.

### 4.2 Bennett 2026 / ALTitude (`10.64898/2026.03.04.709438`) — publishes a within-ALT+ *null*, and flags the gap

They tested the field's leading within-ALT+ candidate and got nothing:

> "We next revisited the relationship between ATRX status, ALT, and dependence on SMARCAL1,
> hypothesizing that ATRX loss may potentiate SMARCAL1 dependency in ALT. Surprisingly
> however, **SMARCAL1 dependency did not differ detectably between ATRX loss (n = 20) and
> ATRX WT (n = 14) ALT+ groups** (Fig. 3j; Cohen's d = −0.248, P = 0.5). Thus, SMARCAL1
> dependency appears to be a general feature of ALT, regardless of whether ALT is driven by
> ATRX loss or some other event(s)."

This independently corroborates iteration 02's `not-estimable` closure of the ATRX
stratifier, and it means the axis is doubly closed: non-estimable here *and* published null
at n = 34.

They also establish that within-ALT+ heterogeneity exists and is unexplained — "there are a
small number of ALT+ cell lines that are not SMARCAL1-dependent … though the mechanisms are
not known (see Discussion)" (48 DepMap lines at Chronos < −0.35; 28 of them ALT+). Their
single offered explanation is **n = 1** (SKNFI's low SMARCAL1 protein, §3.3) and is framed
as a hypothesis.

Two Bennett results that **must not be re-derived** by this iteration: the ATRX within-ALT+
null above, and "ATRX itself was a strong and selective dependency in ALT+ cell lines, when
ALT was not driven by ATRX loss" (the reciprocal ATRX↔SMARCAL1 compensation model, validated
in NY and OS384 osteosarcoma lines).

### 4.3 Brown 2026 (`10.1093/neuonc/noaf300`) — wrong stratifier, wrong lineage

Brown et al., Neuro-Oncology 2026;28(4):895–910 identify SMARCAL1 as synthetic lethal
specifically in **ATRX-deficient, ALT+ glioma**, using DepMap plus patient-derived lines and
xenografts; the effect holds in both IDH-mutant astrocytoma and IDH-WT GBM. The stratifier
is ATRX — the closed axis — in a different lineage, and the comparison is ALT+ vs non-ALT
rather than a within-ALT+ gradient. Bennett's n = 34 null (§4.2) is the better-powered test
of the same idea.

### 4.4 Web search (two queries, question 4 only)

Neither query returned any study relating a continuous exposure — BLM expression or
otherwise — to SMARCAL1 dependency *magnitude*. The returned corpus (Taglialatela Genes Dev
and preprint; Bennett; Brown; Liu 2023 `10.1093/neuonc/noad022`; Diplas-lineage SMARCAL1-LOF
papers) establishes the BLM–SMARCAL1 antagonism qualitatively and the ALT+/ALT− dependency
split quantitatively, and stops there. **BLM expression as a biomarker of SMARCAL1
dependency strength appears nowhere.**

### 4.5 Conclusion for question 4

The iteration does **not** duplicate published work. It must avoid two things that would:
(i) any ATRX-based within-ALT+ split (Bennett, published null; iteration 02, non-estimable);
(ii) any ALT+ vs ALT− group contrast on the 38 genes, which is both published and
contaminated per `DISCLOSURES.md` §1.

Note also that this is the first iteration of the campaign whose primary is a
**correlation between an exposure variable and a line-level gene effect** — which
`DISCLOSURES.md` §2 states is uncontaminated, because no line-level gene-effect value and no
exposure–outcome relationship has been seen. The SMARCAL1 Chronos *marginals* are known
(median −0.462, SD 0.395, 8/16 below −0.5) and are already declared; they bear on the MDE
calculation, not on the predicted direction, which is fixed by §1.2 independently of any
data in this cohort.

---

## 5. Recommended exposure, alternates, and how each fails

### Primary — **BLM mRNA**, predicted **negative**

- **Exposure:** BLM log2(TPM+1) from `OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv`,
  over the 16 ALT+ osteosarcoma lines.
- **Outcome:** SMARCAL1 Chronos gene effect per line (`CRISPRGeneEffect.csv`).
- **Signed prediction:** **Spearman rho(BLM expression, SMARCAL1 Chronos) < 0**, one-sided.
  Higher BLM → more negative Chronos → more SMARCAL1-dependent.
- **Strongest single reason it fails:** BLM transcript level may not track *telomeric
  helicase activity* at all, because that activity is licensed by recruitment and polymerase
  stalling (BTR assembly, PML-dependent localisation) rather than set by abundance — in which
  case a null is uninformative about the mechanism and is the WRAP53 proxy failure in a new
  costume. §3.4's exposure gate and ZRANB3 negative control are what make that outcome
  distinguishable from a real null instead of indistinguishable from one.

### Alternate 1 — **PRIMPOL mRNA**, predicted **negative**

Same sign, same construction. Mechanistically named in the abstract, with a clean
mutant series (CH/AxA/RBD) establishing that it is specifically repriming activity.

- **Strongest single reason it fails:** the paper's own data say PRIMPOL is a *partial*
  contributor — endogenous depletion "only partially suppressed" the phenotype, and
  exacerbation required *over*expression above physiological levels, so the endogenous
  cross-line range may produce no effect at all. Compounded by the authors' own caveat about
  nontelomeric PRIMPOL functions.

### Alternate 2 — **TOP3A mRNA**, predicted **negative**

Grounded in de Nonneville et al. 2022 `10.15252/emmm.202215859`: TOP3A amplification is an
ALT *inducer* in pediatric osteosarcoma and is **mutually exclusive with ATRX inactivation**
— and this cohort is almost entirely ATRX-intact (1/16 altered by point mutation), which is
precisely the TOP3A-driven class. Taglialatela cites it as the alternative route to ALT.
Same sign as BLM, via the BTR complex.

- **Strongest single reason it fails:** the published claim is about **copy-number
  amplification**, not expression. Converting it to an mRNA exposure is an extra unvalidated
  step, while running it as copy number lands in the measurement class `NULL_SPACE.md`
  requires segment validation for — and a segment-validated focal amplification is unlikely
  to reach `MIN_CELL_SIZE = 4`, i.e. probably non-estimable.

### Explicitly rejected as the primary

RPA1/2/3 (direction undetermined by the mechanism; pan-essential, minimal variance);
POLA1 (mechanism not named by this paper and explicitly contradicted by it);
FANCM (two mechanistic channels with opposite signs); ATRX/DAXX (closed twice over);
SMARCAL1's own mRNA (Bennett's published transcript↔protein dissociation);
any multi-gene composite (`NULL_SPACE.md` lesson 2).

---

## Recommendation for the prespec

**Primary exposure: BLM mRNA** (log2(TPM+1), DepMap 26Q1) across the 16 ALT+ osteosarcoma
lines. **Outcome: SMARCAL1 Chronos.** **Predicted sign: negative** — higher BLM, more
negative Chronos, more SMARCAL1-dependent. One-sided.

It rests on this sentence from Taglialatela et al. 2026 (`10.1101/gad.353377.125`):
*"overexpression of WT BLM, but not the KA mutant, markedly exacerbated DNA damage,
telomeric RPA accumulation, APB formation and fitness defects in SMARCAL1-depleted cells"* —
paired with *"BLM deficiency fully suppressed the SMARCAL1-dependent growth defects"*.
Raising BLM is harmless alone and lethal only when SMARCAL1 is gone; Chronos **is** fitness
under SMARCAL1 loss, so the experiment and the endpoint share an axis.

**Main threat to validity:** mRNA measures BLM abundance, while the mechanism requires
unscheduled BLM helicase *activity at telomeres*, which is licensed by recruitment and
polymerase stalling. Bennett (`10.64898/2026.03.04.709438`) documents an mRNA↔protein
dissociation for SMARCAL1 itself. Mitigate by pre-registering, in order: an exposure
variance/bimodality gate read before the outcome; **ZRANB3 mRNA** as a predicted-null
negative-control exposure (Taglialatela: ZRANB3 cannot counteract BLM); and a
proliferation/screen-amplitude control via within-line standardised Chronos plus
negative-control outcomes. No composite endpoints. Nobody has published this test.
