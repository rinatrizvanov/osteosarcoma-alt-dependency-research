# Candidate follow-on hypotheses

Five candidates emerging from the scan. **Scoping only** — nothing here is designed,
pre-specified or costed as a protocol, and none should be run without the same
pre-registration discipline Stage 1 used.

Ordered by what the literature actually supports, not by appeal.

Effort scale: **S** ≈ days on data already held · **M** ≈ weeks, new public data ·
**L** ≈ months, external/collaborator data.

---

## H1 — TCAB1 loss predicts ALT *status*, not dependency within ALT

**Hypothesis.** 17p13.1/TCAB1 disruption is an ALT-**initiating** lesion: it predicts
whether an osteosarcoma line is ALT+ at all, and carries little or no information about
which recombination machinery an already-ALT+ line depends on.

**Motivation.** Keegan et al. 2026 (`Literature/4 TCAB.pdf`) — TCAB1 SVs are "an early
event contributing to the **activation** of the ALT pathway," present in ~40% of ALT+
tumours, with no claim about differential dependency. The preprint never mentions BIR,
POLD3 or RAD52.

**Data needed.** Already held. DepMap 26Q1 expression + the frozen ALT labels. The
contrast moves from *within* ALT+ to **ALT+ vs ALT−**, using the 16/9 split Stage 1
already built and set aside.

**Effort.** **S.** The cohort, loaders, stratifier and test machinery all exist.

**Relation to Recombination Paradox.** Not a variant — it is the **diagnosis**. It says
Stage 1 conditioned on the mediator. Testing it is how you find out whether the
underlying biology is real and Stage 1 simply looked in the wrong stratum.

**Caveat worth stating up front.** With 16 ALT+ vs 9 ALT− and WRAP53 near-zero in only
2 lines, this is underpowered too, and at N=25 a null would again be uninformative.
Worth doing mainly because it is nearly free and directly interprets the Stage 1 result.

---

## H2 — SMARCAL1 dependency in ALT+ is gated by ATRX, not TCAB1

**Hypothesis.** Within ALT+ lines, SMARCAL1 essentiality is selectively elevated in
**ATRX-deficient** lines, and TCAB1/WRAP53 status adds nothing once ATRX is accounted for.

**Motivation.** https://doi.org/10.1093/neuonc/noaf300 — SMARCAL1 as a synthetic lethal
vulnerability specifically in **ATRX-deficient** ALT gliomas. Someone has already done
the stratification Stage 1 attempted; the stratifier that works is ATRX. Supported by
https://doi.org/10.1093/neuonc/noad022 and
https://doi.org/10.1038/s41379-021-00841-7, in which SMARCAL1 **loss** promotes ALT and is
enriched in ALT+ tumours — the inverted polarity that plausibly explains why SMARCAL1 ran
backwards in your WRAP53-null lines.

**Data needed.** DepMap pan-cancer, not osteosarcoma-only. Your osteosarcoma cohort cannot
test this: **only 1 of 16 ALT+ lines is ATRX-altered**, which is exactly why R5 was
vacuous. Needs ALT+ lines across lineages (Bennett's ALTitude calls roughly double the
known ALT+ models) plus a properly derived ATRX status — deletion and expression loss, not
just damaging point mutations, which is the flaw in the inherited `atrx_status`.

**Effort.** **M.** New cohort construction, new ATRX labelling, cross-lineage confounding.

**Relation.** Adjacent. Same locus-to-dependency logic, different gate, and an
already-published effect to replicate rather than discover — which makes a null
interpretable for once.

---

## H3 — ALT branch identity, not TCAB1, determines which repair genes are essential

**Hypothesis.** ALT+ lines partition into replication-fork-processing branches, and BIR
gene essentiality tracks branch identity; a composite averaged over branches has no
detectable effect, which is what Stage 1 observed.

**Motivation.** https://doi.org/10.1038/s41467-023-38029-z — ALT pathway choice is
dictated by replication fork processing. Reinforced by Carson & Flynn 2023's statement
that RAD51 and RAD52 "regulate distinct branches of the ALT pathway," and by
https://doi.org/10.1016/j.cancergen.2026.06.001 ("Beyond ATRX: the heterogeneous basis of
ALT in osteosarcoma").

**Data needed.** A branch assignment for each ALT+ line, which does not currently exist as
a label — it would have to be derived from fork-processing or APB-related markers, and
that derivation is itself the hard part and the main risk. Then DepMap dependencies.

**Effort.** **M–L**, and the honest reading is that it is closer to L: the exposure
variable has to be invented before anything can be tested, which is how Stage 1's proxy
problem arose in the first place.

**Relation.** Adjacent, and the most mechanistically interesting. Also the easiest to
fool yourself with, since a derived subtype label fitted on the same data that supplies
the outcome is circular. Would need the label frozen before any dependency is read.

---

## H4 — POLD3, RAD52 and SMARCAL1 do not behave as a module in ALT+ cells

**Hypothesis.** The three genes' essentiality profiles across ALT+ lines are mutually
uncorrelated, so no composite over them is a coherent endpoint.

**Motivation.** Bennett/ALTitude treats SMARCAL1 alone and never mentions POLD3 or RAD52;
Carson & Flynn never mention POLD3. Against that,
https://doi.org/10.1038/s41467-023-40071-w shows RAD52 and PolD3 *are* mechanistically
coupled in fork recovery, so the question is genuinely open rather than rhetorical. Your
own data is the third strand: POLD3 common-essential, RAD52 non-essential, SMARCAL1
reversed.

**Data needed.** Already held — DepMap CRISPR across ALT+ lines, pairwise correlation of
gene effects. Extends naturally pan-cancer.

**Effort.** **S** on held data, **M** pan-cancer.

**Relation.** A direct methodological post-mortem on the Stage 1 endpoint rather than a
biological hypothesis. Its value is defensive: if the module construct fails here, it
should not be reused, and that finding belongs in the preprint's discussion whatever else
you do.

---

## H5 — WRAP53 expression is a poor proxy for 17p13.1 structural disruption

**Hypothesis.** Among tumours with known 17p13.1 SV status, WRAP53 expression
misclassifies disruption at a rate high enough to have made Stage 1 untestable regardless
of the underlying biology.

**Motivation.** Axis G found **no literature** on when expression proxies fail for
SV-driven phenotypes, and the nearest work
(https://doi.org/10.1038/s41576-018-0007-0, https://doi.org/10.1101/gr.145144.112) shows
SV-to-expression coupling is weak and context-dependent. Sharpened by
https://doi.org/10.1016/j.molcel.2016.11.027: the WRAP53 antisense transcript is required
for p53 induction, so WRAP53 and TP53 are functionally coupled at the RNA level, not just
neighbours.

**Data needed.** A cohort with **both** WGS SV calls and RNA at 17p13.1. DepMap 26Q1 has
no SV calls, so this requires external data — TARGET-OS, PCAWG, or a collaborator-held PDX cohort.

**Effort.** **M** on public data (TARGET-OS/PCAWG), **L** if it needs the PDX cohort.

**Relation.** A measurement-validity study, not a biological one. It is the honest
prerequisite for *any* expression-proxied version of the Recombination Paradox, including
H1. If the proxy fails here, Stage 1's result is "untested," not "not supported" — the
distinction prespec §11.3 was written for.

---

## Reading order if you only have an afternoon

1. https://doi.org/10.1093/neuonc/noaf300 — the stratification you attempted, done with ATRX
2. https://doi.org/10.1093/neuonc/noad022 — SMARCAL1 loss *promotes* ALT; the polarity problem
3. https://doi.org/10.1038/s41467-023-38029-z — ALT branch choice
4. https://doi.org/10.1016/j.cancergen.2026.06.001 — ALT heterogeneity in osteosarcoma specifically
5. https://doi.org/10.1101/gad.246546.114 — what Henriksson 2014 actually demonstrated

## One scoping judgement

H1 and H4 are nearly free and interpret the result you already have; H2 replicates a
published effect, which is the only candidate here where a null would be informative.
H3 and H5 are more interesting and considerably more expensive, and both require
inventing or validating an exposure variable before any outcome is touched — which is
precisely where Stage 1 went wrong. If only one goes forward, H2 is the one with an
external anchor.
