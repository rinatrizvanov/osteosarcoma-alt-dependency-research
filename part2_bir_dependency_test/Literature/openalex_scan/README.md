# OpenAlex literature scan — scoping after the Stage 1 null

Run 2026-10-11 against the OpenAlex API (polite pool, authenticated key). Seven
pre-specified axes plus one extra term query. Each axis file records its exact query
URLs with the key redacted.

**Scope.** Reading recommendations and hypothesis candidates only. No new analysis was
run and the Stage 1 result (`stage1-v1.0-null`) is untouched.

| Axis | Topic | Matching works | File |
|---|---|---|---|
| A | WRAP53β DDR post-2015 | 35 | `axis_a_wrap53b_ddr.md` |
| B | TCAB1/WRAP53 cancer dependency | 41 | `axis_b_tcab1_dependency.md` |
| C | 17p13.1 SV architecture | 115 | `axis_c_17p13_architecture.md` |
| D | BIR module composition | 38 | `axis_d_bir_module.md` |
| E | SMARCAL1 context-dependency in ALT | 50 | `axis_e_smarcal1_context.md` |
| F | ALT subtype heterogeneity | 42 | `axis_f_alt_subtypes.md` |
| G | Expression vs SV proxies | 6 → reframed, 548 | `axis_g_expression_vs_sv_proxies.md` |
| Z | "Recombination Paradox" as a term | 3 | `axis_z_recombination_paradox_term.md` |

---

## Read this first: three reasons for the null are in papers you already own

Found by reading `Literature/*.pdf` before searching, not by the scan.

**1. Keegan frames TCAB1 loss as ALT-*initiating*, not as a dependency modifier within
ALT+.** The abstract's claim is that TCAB1 SVs are "an early and previously unrecognized
event contributing to the **activation** of the ALT pathway," in ~40% of ALT+
osteosarcoma. The preprint contains **zero** occurrences of BIR, break-induced, POLD3,
RAD52, WRAP53β or synthetic lethality. Stage 1 stratified *within* ALT+ — conditioning on
the very state TCAB1 loss is proposed to cause. An initiation effect is largely absorbed
by that conditioning. This is the single most likely reason for the null.

**2. The "BIR module" is not a construct either source paper supports.** Bennett/ALTitude
mentions SMARCAL1 264 times and POLD3, RAD52, BIR, WRAP53 and TCAB1 **zero** times — it is
a standalone SMARCAL1 dependency paper. Carson & Flynn 2023 never mentions POLD3 at all.
The three-gene composite was assembled from general BIR biology, not from the literature
being tested.

**3. Flynn's own review says RAD51 and RAD52 are "not functionally redundant and in fact
regulate distinct branches of the ALT pathway."** Stage 1 placed RAD52 in the BIR arm and
RAD51 in the HR comparator arm — opposite sides of one difference. If both are pro-ALT
branches, the endpoint differenced two arms of the same pathway. "BRCA" appears **zero**
times in that review, so BRCA1/2 as ALT comparators also lack support. This matches the
observed data exactly: POLD3 behaved as a common essential (≈−1.6 in both arms), RAD52 as
non-essential (≈0), SMARCAL1 backwards. Three genes behaving three ways is what "not a
module" looks like.

---

## Axis A — WRAP53β DNA damage response

1. **https://doi.org/10.1101/gad.246546.114** — Henriksson et al. 2014, *Genes Dev*. The
   primary experimental paper behind your entire HR-scaffold premise. Read it for what it
   actually shows: WRAP53β scaffolds the **RNF8/RNF168 ubiquitin response** at DSBs,
   governing recruitment of 53BP1 and BRCA1. That is DSB *signalling*, upstream of and
   distinct from "the classical HR scaffold." The Stage 1 mechanism chained this to a BIR
   dependency shift, a link this paper does not make.
2. **https://doi.org/10.1016/j.molcel.2016.11.027** — Mahmoudi et al. 2016, *Mol Cell*.
   The WRAP53 **antisense transcript** is required for p53 induction after DNA damage.
   WRAP53 and TP53 are coupled at the RNA level, not merely adjacent on 17p13.1. This
   complicates interpreting WRAP53 expression independently of TP53 and is a deeper
   version of the concern your strand check addressed.
3. **https://doi.org/10.3390/cimb44110372** — Gadelha et al. 2022. Current review of
   WRAP53 in homeostasis and carcinogenesis; the fastest way to see what has accumulated
   since Henriksson & Farnebo 2015.
4. **https://doi.org/10.1038/s41419-020-2421-4** — Biallelic WRAP53 mutations cause
   dysfunctional telomeres, Cajal bodies and DNA repair. The human loss-of-function
   phenotype, useful for calibrating what complete TCAB1 loss does.

## Axis B — TCAB1/WRAP53 cancer dependency

**The headline is a gap, not a paper: there is no published CRISPR dependency or
synthetic-lethality screen for TCAB1/WRAP53.** The query surfaced CRISPR-lncRNA
methodology and telomerase-assembly biology instead. Nobody has done the screen.

1. **https://doi.org/10.1038/cddis.2010.90** — Mahmoudi et al. 2011. WRAP53 promotes
   cancer cell survival and is proposed as a therapeutic target — i.e. WRAP53 loss is
   argued to be *bad* for cancer cells, the opposite polarity to a dependency-conferring
   lesion.
2. **https://doi.org/10.1128/mcb.00379-12** — Stern et al. 2012. Telomerase recruitment
   needs TCAB1 and Cajal bodies independently; relevant to what "TCAB1 loss" means
   mechanistically.
3. **https://doi.org/10.1002/1878-0261.13324** — pan-cancer multi-omic signatures of
   telomerase complex components. Closest existing thing to a systematic TCAB1 analysis.

## Axis C — 17p13.1 structural variant architecture

Mostly pre-genomic cytogenetics. Three worth the time:

1. **https://doi.org/10.3324/haematol.2026.301572** — Mueller et al. 2026, *Haematologica*.
   "Beyond 17p13.1: recurrent 17p11.2 breakpoints in TP53-altered AML." The nearest
   analogue to Keegan outside osteosarcoma, and evidence the locus architecture recurs
   across cancers.
2. **https://doi.org/10.1002/gcc.20579** — transcriptional characterization of 17p loss in
   B-CLL. One of the few papers linking 17p structural loss to downstream expression,
   which is the Stage 1 proxy question in another disease.
3. **https://doi.org/10.1016/j.ajhg.2014.10.006** — dosage changes at 17p13.1 in
   developmental disease; clean germline mapping of the segment.

## Axis D — BIR module composition and context-specificity

1. **https://doi.org/10.1038/s41467-023-38029-z** — "Pathway choice in the alternative
   telomere lengthening in neoplasia is dictated by replication fork processing." The most
   important paper on this axis. ALT is not one pathway; which branch a tumour uses is
   determined upstream. If branch identity sets which genes are essential, then an
   unstratified BIR composite averages across branches — exactly the dilution your R8 and
   per-gene results suggest.
2. **https://doi.org/10.1038/s41467-023-40071-w** — RAD52-dependent fork recovery is
   regulated through ADP-ribosylation of **PolD3**. Direct evidence RAD52 and POLD3 are
   mechanistically coupled, which is the strongest case *for* a module — and worth reading
   against your finding that they behaved nothing alike.
3. **https://doi.org/10.1016/j.celrep.2024.114116** — RAD52-dependent mitotic DNA synthesis
   is required in Cyclin E1-overexpressing cells. A worked example of RAD52 essentiality
   being conditional on an oncogenic background rather than constitutive.
4. **https://doi.org/10.1038/s41467-023-35992-5** — mitotic DNA synthesis requires
   sequential polymerase ζ then δ action. Bears on whether POLD3 is a specific BIR
   component or general replication machinery; your data had it as a common essential.

## Axis E — SMARCAL1 context-dependency in ALT

**This axis answers your question directly, and the answer inverts an assumption.**

1. **https://doi.org/10.1093/neuonc/noad022** — "Cancer-associated SMARCAL1
   **loss-of-function** mutations **promote** alternative lengthening of telomeres and
   tumorigenesis." SMARCAL1 loss *drives* ALT here. That is the opposite polarity to
   "SMARCAL1 is a selective dependency of ALT+ cells."
2. **https://doi.org/10.1038/s41379-021-00841-7** — SMARCAL1 loss and ALT are **enriched**
   in giant cell glioblastoma. Independent support for the same inverted polarity.
3. **https://doi.org/10.1093/neuonc/noaf300** — "SMARCAL1 is a targetable synthetic lethal
   vulnerability in **ATRX-deficient** gliomas that use ALT." **The most actionable paper
   in this scan.** Someone has stratified SMARCAL1 essentiality within ALT+ — and the
   stratifier that works is **ATRX**, not TCAB1. Note the bind this puts you in: only
   1 of your 16 ALT+ lines was ATRX-altered, which is why R5 was vacuous. The field's
   working stratifier is the one your cohort could not test.
4. **https://doi.org/10.1016/j.celrep.2016.01.011** — Cox et al. 2016. The foundational
   "SMARCAL1 resolves replication stress at ALT telomeres" paper.

Taken together, 1–3 give a coherent account of why SMARCAL1 ran *backwards* in your
WRAP53-null lines: if SMARCAL1 loss promotes ALT, then ALT+ lines that arrived via a
SMARCAL1-low route need not depend on it.

## Axis F — ALT subtype heterogeneity

1. **https://doi.org/10.1016/j.cancergen.2026.06.001** — "Beyond ATRX: the heterogeneous
   basis of Alternative Lengthening of Telomeres in **osteosarcoma**." Directly your
   disease, directly the subtype question, 2026. Read first on this axis.
2. **https://doi.org/10.1038/s41467-025-67190-w** — large-scale drug sensitivity, gene
   dependency and proteogenomic analysis of telomere maintenance mechanisms. A ready-made
   dependency resource stratified by TMM; the natural external check on anything you find
   in DepMap.
3. **https://doi.org/10.1093/neuonc/noae016** — telomere maintenance mechanisms and
   associated therapeutic vulnerabilities in gliomas. The best-developed example of
   mapping vulnerabilities onto TMM subtypes.
4. **https://doi.org/10.3390/cells14201613** — "ALT: a prognostic paradox in cancer"
   (2025). Note for your preprint title: this is the closest published use of "paradox"
   in ALT, and it means something different from your usage.

## Axis G — Expression vs structural proxies

**The axis as posed has no corpus.** The pre-registered query returned 6 works, none on
topic; four broadenings drifted into maize, rice, yeast and Parkinson's. "When does an
expression proxy fail to capture an SV-driven phenotype" is not a named research area.
That is a finding, not a failed search: your Stage 1 limitation is under-studied in
general, which is mildly in your favour when writing it up.

The nearest real literature inverts the question — how much expression variance SVs
explain, which bounds how well expression can proxy SV status:

1. **https://doi.org/10.1038/s41576-018-0007-0** — Spielmann et al. 2018, *Nat Rev Genet*,
   structural variation in the 3D genome. Why an SV's expression consequence depends on
   chromatin context and is frequently *not* a simple dosage effect on the disrupted gene.
2. **no DOI available — OpenAlex W2135185276** — Hyman et al. 2002, impact of DNA
   amplification on expression in breast cancer. The classic quantification of how loosely
   copy number predicts expression.
3. **https://doi.org/10.1101/gr.145144.112** — transcriptional landscape and mutational
   profile of lung adenocarcinoma; paired SV/expression at scale.

## Axis Z — "Recombination Paradox"

**Not a published term.** Three matches across all of OpenAlex: a 2005 Ku/telomere review
using the phrase incidentally, a 1972 prophage chapter, and a 2026 perovskite solar-cell
paper. It appears in none of your seven PDFs, including Flynn's own 2023 review.

Treat it as an internal coinage: define it on first use in the preprint, do not cite it as
established, and expect a reviewer to ask where it comes from. If it originated as an advisor's
verbal framing, that is worth saying explicitly in acknowledgements rather than leaving
the term unsourced.
