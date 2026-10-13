# Literature scan — H1: does WRAP53/TCAB1 loss predict ALT *status*?

Iteration 01. Literature checker. Built 2026-10-12 from `Literature/INDEX.md` (read in
full) plus targeted `pdftotext` + grep over all 156 PDFs in `Literature/` (no PDF opened
with the Read tool; no data file loaded).

---

## VERDICT

**Partially supported, and the supported half is already published** — Keegan 2026 has
already established the TCAB1–ALT-status link in osteosarcoma *tumours*, at the
**structural-variant and protein (IHC) level**; the specific prediction under test here
(ALT+ **cell lines** have **lower WRAP53 mRNA**) is **not addressed** by any indexed
paper and is **actively undercut** by four of them, including Keegan's own RNA-seq
control.

---

## 1. Does Keegan et al. claim TCAB1 inactivation is ALT-initiating? Tumours or cell lines?

**Yes, it claims initiation — and the claim is established in TUMOURS and PDX tumour
tissue, not in cell lines.**

`10.1186/s12885-026-16250-5` (BMC Cancer 2026;26:948). Verbatim, from the abstract:

> "Here, using whole-genome sequencing of osteosarcoma samples we identified SVs within
> the TP53 gene that not only compromise TP53, but also inactivate TCAB1. These TCAB1 SVs
> were prevalent in **approximately 40% of ALT positive osteosarcoma tumors** suggesting
> that functional inactivation of the telomerase holoenzyme is an **early** and previously
> unrecognized event contributing to the activation of the ALT pathway."

The in-text figure is slightly lower:

> "Together, our data demonstrate that **approximately 35% of the ALT positive samples and
> none of the ALT negative samples, harbor SVs in TCAB1** (Fig. 4B)"

**Evidence types, cohort by cohort:**

| Cohort | n | Material | ALT assay | TCAB1 assay | ALT− comparator? |
|---|---|---|---|---|---|
| GCCRI | 16 | Early-passage, viably frozen **PDX tumour** samples | TRAP (telomerase-negative) + TEL-FISH (ultrabright foci) on FFPE sections | **IHC** protein; WGS on the TCAB1-deficient subset (5 ALT+ models by WGS) | **No** — "all sixteen samples displayed these specific molecular characteristics indicative of ALT" |
| TARGET-OS | 21 of 23 cases with quality WGS | **Patient tumours** (public WGS, phs000218) | **TelFusDetector computational prediction** (not a tissue assay) — 10 ALT+ / 11 ALT− | WGS SV calling (DELLY + MANTA + SvABA, ≥2 callers) | **Yes** — the only ALT+ vs ALT− contrast in the paper |
| TCGA SARC vs BRCA | — | Patient tumours | — | focal CN at the TP53/TCAB1 locus | **No** — "Although we cannot stratify these data based on ALT status" |

The GCCRI IHC result: "we found that 75% of the tumors in our cohort lacked TP53
expression. Likewise, **37.5% (6/16)** of the cohort demonstrated, marked reduction
(0560-000), or loss (OS33, OS37, OS39, OS52, 1981) of TCAB1 expression."

**Cell lines in Keegan 2026: two, both as controls only.** SJSA1 (telomerase-positive /
TP53-mutant / TCAB1-positive control for TRAP, TEL-FISH and both IHC stains) and HEK293
(PCR control). There is **no cell-line panel, no DepMap analysis, and no cell-line ALT+ vs
ALT− comparison anywhere in the paper.** The preprint (`10.64898/2026.01.02.691864`,
`Literature/4 TCAB.pdf`) is the same study; cite the published version.

**Three further caveats from Keegan's own text that bear directly on H1:**

- TCAB1 loss is framed as **collateral damage from TP53 intron-1 SVs**, not an independent
  event: "all the tumors demonstrating defects in TCAB1 expression also showed loss of TP53
  expression… suggest that **loss of TCAB1 may be a consequence of TP53 SVs**", and "These
  SVs in TCAB1 were only identified in tumors that also retain SVs in TP53."
- TCAB1 loss is **not sufficient** for ALT in Keegan's own model: "TCAB1 deficient tumors
  that acquire additional pathogenic mutations (i.e. ATRX, DAXX, SMARCAL1, SLX4IP) may
  ultimately activate ALT whereas, TCAB1 deficient tumors that do not gain these mutations
  could either progress to an EST phenotype or simply undergo cell death." And: "**In fact,
  there were two ALT negative tumors with SVs in TCAB1.**" (Reconciled with "none of the
  ALT negative samples" because those SVs were judged not functionally inactivating — see §6.)
- The TP53 SV contrast itself does **not** discriminate ALT: "approximately 50% of the ALT
  positive samples and **36% of the ALT negative samples** showed SVs in TP53."

---

## 2. Has anyone tested whether WRAP53/TCAB1 *expression* discriminates ALT+ from ALT−?

**Not at the mRNA level, in any lineage, in this corpus.** The only ALT+ vs ALT− TCAB1
contrast ever published in the indexed literature is **structural-variant based**, in
tumours:

- `10.1186/s12885-026-16250-5` — TARGET-OS WGS, **10 ALT+ vs 11 ALT−**, functionally
  inactivating TCAB1 SVs in ~35% of ALT+ vs **0%** of ALT−. ~35% of n=10 implies 3–4
  samples versus 0 of 11. Effect metric: a prevalence difference, not an expression
  difference; n is in single digits per arm and ALT labels are computationally predicted.
- The same paper's IHC arm (6/16 = 37.5% TCAB1 protein loss/reduction) is **protein-level
  and entirely within ALT+** — all 16 PDX models were ALT+, so it yields no contrast.

A corpus-wide scan for papers containing both a WRAP53/TCAB1 alias and an ALT term returns
only 10 files; of these, only Keegan (both versions) uses both substantively. The near
misses:

- **Wang et al. 2022 Mol Oncol, `10.1002/1878-0261.13324`** is the one paper that contains
  *both* a TCAB1 analysis and an ALT classification — TCGA/CCLE pan-cancer, 11,352 tumours
  partitioned into TERT+ (74.2%), ALT+ (8.4%, defined by ATRX **or** DAXX mutation) and
  TERT−/ALT− (17.4%) — but it **never crosses the two**. Its §3.11 ALT analysis reports only
  telomerase score; TCAB1 appears only in the CNA/expression sections. So even here the
  decisive 2×2 was not run.
- **Gadelha et al. 2026 Cancer Rep, `10.1002/cnr2.70464`** measures low WRAP53 expression in
  a real malignancy — qPCR in **110 AML patients** vs 10 healthy donors, mean WRAP53
  −2.490 vs −0.1111, p < 0.0001, with 79/110 patients at FC ≤ −1 — but AML is not an ALT
  lineage and no ALT status was measured. If anything it points the wrong way: a ~2.4 log2
  unit WRAP53 deficit in a malignancy that is overwhelmingly telomerase-driven.

**Conclusion: H1's expression-level formulation has not been tested and is therefore novel
— but it is novel in the same way the Stage 1 BIR hypothesis was novel (Gap 4 in the
index: no paper joins TCAB1 and BIR). Novelty here coexists with specific published
reasons to expect it to fail.**

---

## 3. Do any indexed papers report WRAP53/TCAB1 levels in ALT+ versus ALT− cells or tumours, even incidentally?

Four incidental data points, three of which point **against** the predicted direction.

1. **Keegan 2026, Supplemental Fig. 3 — the most damaging single observation.** For the
   three TARGET samples with a single TCAB1 SV disrupting the α
   (ENST00000431639.6) and γ (ENST00000457584.6) promoters but not β
   (ENST00000396463.7): "analysis of RNA sequencing data that corresponded to 2 of the
   three samples … **confirmed that TCAB1 RNA expression levels were similar to normal
   controls**." In the only place where Keegan measured TCAB1 mRNA in SV-bearing ALT
   tumours, mRNA was **normal**.
2. **Mahmoudi et al. 2011, `10.1038/cddis.2010.90`** — the Western panel of
   WRAP53-overexpressing cancer lines (Fig. 1a) explicitly includes **U2OS and Saos-2**,
   i.e. two osteosarcoma lines the corpus elsewhere identifies as ALT models (Bennett 2026:
   "the U2OS line, an ALT+ osteosarcoma cell line"; Liu 2023: "ATRX-deficient ALT models
   (e.g. U2OS, SaOS2)"). Verbatim: "WRAP53 protein expression is increased in immortalized
   cells and **up to 20 times higher in cancer cells compared with primary cells**." ALT+ OS
   lines are therefore WRAP53-**high** relative to primary cells. (Caveat: the comparator is
   primary/immortalized non-cancer cells, not ALT− OS lines — this is not a test of H1, but
   it is the wrong sign for the "ALT+ lines are WRAP53-low" intuition.)
3. **Wang et al. 2022, `10.1002/1878-0261.13324`** — TCAB1 homozygous deletion in 6.6% of
   TCGA SARC (also 11.7% PRAD, 6.9% LIHC), and "**40% of blood cancer-derived cell lines
   carry homozygous TCAB1 deletion, and there was a highly significant reduction of TCAB1
   expression in these cell lines**." So deep, expression-visible TCAB1 loss *does* occur in
   CCLE lines — but in telomerase-competent ones (see §6).
4. **Gadelha 2026** (above) — low WRAP53 mRNA is measurable in tumours, decoupled from hTERT
   transcript level.

---

## 4. Is there a published ALT-status classifier that includes WRAP53 as a feature?

**No. Every published ALT classifier or ALT resource in the corpus omits WRAP53 entirely.**

| Source | DOI | WRAP53/TCAB1 mentions in full text | WRAP53 a feature? |
|---|---|---|---|
| Armendáriz-Castillo 2026 Biology, 28-gene ALT signature | `10.3390/biology15171531` | **0** | **No** |
| Bennett 2026 bioRxiv, ALTitude (WGS→ALT tabular foundation model, doubled known ALT+ DepMap lines) | `10.64898/2026.03.04.709438` | **0** | **No** |
| Taglialatela 2026 Genes Dev | `10.1101/gad.353377.125` | **0** | n/a |
| Wu 2025 Nat Commun, 976-line telomere-maintenance resource | `10.1038/s41467-025-67190-w` | **1**, in the reference list only (Freund 2014, TRiC folding of TCAB1) | **No** |

The Armendáriz-Castillo 28 genes verbatim: ARID1A, SMG7, ATF6, SENP5, ATR, RFC4, RAD1,
**CCT5**, CDKN1A, HUS1, MCPH1, PRKDC, TERF1, LRATD2, NBN, RAD54B, UBR5, RAD21, NSMCE2,
RECQL4, FBH1, MRE11, ATM, UBE2I, TERF2, TP53, DNMT1, PCNA. WRAP53 is absent. Note the
irony that **CCT5 is in** — the chaperonin subunit that Jiang 2025
(`10.1016/j.stemcr.2025.102750`) shows sets TCAB1 **protein** abundance. The field's own
ALT signature captures the TCAB1 axis at the chaperonin, not at the transcript.

**Does the external literature agree with `references/notebook_13_report.md`?** Yes,
emphatically. In that prior project's 19-gene signature WRAP53's L1 coefficient was exactly
**+0.000000 in all three cohorts** (primary N=20, pragmatic N=31, pragmatic_ks N=27),
bootstrap non-zero frequency 0.17 / 0.12 (both "unstable"), and LOOCV fold frequency
**0.0000** (primary) / 0.0968 / 0.1111 — the second-weakest gene in the set. The external
literature independently reaches the same place by never selecting WRAP53 at all. Two
independent feature-selection exercises (one ours, one Armendáriz-Castillo's) and one
state-of-the-art WGS-based ALT classifier (Bennett) all exclude it.

---

## 5. What does Mahmoudi 2011 imply about the achievable low-WRAP53 range in a viable screened panel?

`10.1038/cddis.2010.90`: WRAP53 is up to **20× higher** in cancer lines than in primary
cells; overexpression transforms NIH 3T3 (4× more, larger anchorage-independent colonies);
knockdown kills HeLa and **U2OS** within 2–3 days via Bax/Bak mitochondrial apoptosis,
rescuable by Bcl-2; "cancer cells, but not normal cells, depend on WRAP53." Corroborated
by **Niu 2021** (`10.1186/s12935-020-01745-3`: TCAB1 knockdown → rapid senescence in A549,
H1299, HeLa via p21 protein stabilisation) and **Yuan 2014** (`10.3892/ol.2014.2306`:
TCAB1 knockdown antiproliferative, G1 arrest in A549).

**Implications for H1, in order of severity:**

1. **Range compression / survivorship.** Any line in a CRISPR-screened, continuously
   passaged panel has already survived selection *for* WRAP53 sufficiency. Genuinely
   WRAP53-low lines are a selected minority that necessarily carry compensation (Bcl-2-like
   apoptotic buffering, or p21-pathway loss — and OS lines are near-uniformly TP53-null,
   which supplies exactly that). The observable low tail is truncated and non-random.
2. **Direction.** Mahmoudi's baseline is *overexpression*, including in two ALT+ OS lines.
   A panel-wide expectation of "ALT+ ⇒ WRAP53-low" has to overcome a lineage-wide upward
   baseline. The plausible signature of TCAB1 inactivation in a viable panel is therefore a
   **small low-expression subpopulation inside an otherwise high distribution** — a mixture,
   not a location shift.
3. **Analysis-design corollary.** A median split is close to the worst available test for a
   ~35%-prevalence minority subgroup: it puts roughly half the "high" subgroup into the low
   arm. The signal, if any, lives in the lower ~tertile/tail.

---

## 6. Counter-evidence: TCAB1 loss as a consequence of ALT, or independent of it?

Nothing in the corpus argues TCAB1 loss is **downstream of ALT** — that specific claim is
**not established from the index** and I found no paper making it. But there are three
indexed lines of counter-evidence that are stronger than a reverse-causation argument,
because they attack the premise directly.

**(a) Low/absent TCAB1 does not abolish telomerase, so it does not force ALT.** Wang et al.
2022 (`10.1002/1878-0261.13324`), abstract: "Surprisingly, primary tumors and cancer cell
lines frequently displayed a homozygous deletion of the TCAB1 gene… However, **tumors or
cells carrying a TCAB1 deletion still exhibited telomerase activity comparable to or even
higher than their wildtype counterparts.**" In CCLE: 40% of blood-cancer lines are TCAB1
homozygous-deleted with "highly significant reduction of TCAB1 expression", yet EXTEND
telomerase activity shows "no differences between cell lines with and without TCAB1
deletion"; across 33 TCGA cancer types, 21/33 show comparable activity across deletion
classes; and in six hematological lines "TCAB1-deleted cells expressed significantly lower
levels of TCAB1 protein, but their telomerase activity was not decreased or even higher…
There was **no correlation between TCAB1 protein levels and telomerase activity**", with
HL60 expressing "negligible TCAB1" but the **highest** telomerase activity of the six.
This is cell-line-level evidence, in CCLE, that WRAP53-low lines are telomerase-positive —
i.e. the exact population H1 predicts should be ALT+.

**(b) TCAB1 loss is a passenger of TP53 inactivation, and TP53 SVs do not discriminate
ALT.** Keegan 2026: "loss of TCAB1 may be a consequence of TP53 SVs"; TCAB1 SVs occurred
only in TP53-SV tumours; and TP53 SVs were present in 50% of ALT+ **and 36% of ALT−**
samples. In a DepMap OS panel, where TP53 loss is near-universal, the upstream exposure is
non-discriminating by construction.

**(c) TCAB1 deficiency is explicitly insufficient for ALT, and its most likely outcome is a
no-TMM state that a binary classifier calls ALT−.** Keegan's own model: TCAB1-deficient
tumours "immediately lack any telomere length maintenance mechanism" and need ATRX/DAXX/
SMARCAL1/SLX4IP hits to reach ALT; otherwise they become EST (ever-shorter telomeres) or
die. "In fact, there were two ALT negative tumors with SVs in TCAB1." Wu et al. 2025
(`10.1038/s41467-025-67190-w`) independently shows this intermediate class is real and
populated in cell lines: TMM-double-negative, TMM-double-positive, and
EverShorter-Telomere lines. Under a binary silver label, TMM-DN TCAB1-low lines land in the
**ALT− arm** — which would push the measured effect toward zero or invert it.

Also relevant to the premise, though not causal: only 1 of the 16 ALT+ OS lines in the
Stage 1 panel was ATRX-altered (index, Prespec conflict #6), so the cofactor Keegan's model
requires for the TCAB1→ALT transition is essentially absent from this panel.

---

## EXPECTED EFFECT

**Direction.** If H1 is true at the mRNA level, ALT+ lines should be **lower** in WRAP53 —
sign as predicted. But note that two indexed observations point the other way in OS
specifically (Mahmoudi 2011's ALT+ U2OS/Saos-2 in the overexpressing panel; Wang 2022's
WRAP53-low-but-telomerase-high lines), so a null or sign-inverted result is fully
compatible with the corpus.

**Magnitude — the literature supports no published effect size for an expression contrast.**
No paper in the corpus reports a WRAP53 expression difference between ALT+ and ALT− groups.
What is defensible is an **arithmetic ceiling derived from Keegan's prevalence figure**:

- Keegan: functionally inactivating TCAB1 lesions in ~35–40% of ALT+ and ~0% of ALT−.
- Best case for H1 — assume every lesioned ALT+ line is transcript-low and ranks below
  *every* ALT− line, and the remaining 60–65% of ALT+ lines are exchangeable with ALT−.
  Then Cliff's δ = 0.35 × (−1) + 0.65 × 0 ≈ **−0.35** (ceiling, |δ| ≤ 0.35–0.40). This is
  the *maximum* a prevalence of 35% can generate in a rank statistic; it is not a point
  estimate.
- Realistic attenuation, all from indexed sources: Keegan's own RNA-seq found **normal**
  TCAB1 mRNA in the SV-bearing samples it checked (2/2); Jiang 2025 and Bergstrand 2020
  show TCAB1 abundance is set post-translationally by CCT/TRiC with intact transcript;
  Sánchez-Morán 2020 shows the DDR function is set by localisation; the lesions disrupt the
  α/γ promoters while sparing β, and total WRAP53 TPM sums all three isoforms (Gadelha 2022,
  `10.3390/cimb44110372`). Each of these moves the transcript-level δ toward zero
  independently. A realistic expectation is **|δ| < 0.2**, plausibly indistinguishable from
  zero in a panel of this size.
- **Shape, not location.** The literature's prediction is a **mixture**: a high-WRAP53 bulk
  with a minority low tail of ~1/3 of ALT+ lines. A median split on a 16-line arm has
  essentially no power against that alternative. If this iteration proceeds, the
  better-specified endpoint is the **lower tail / tertile** or an explicit mixture
  comparison, not a median contrast.
- **For reference, the one real magnitude in the corpus** is Gadelha 2026's AML result
  (mean WRAP53 −2.490 vs −0.111, p < 0.0001, n = 110 vs 10) — ~2.4 log2 units. It shows
  large WRAP53 deficits *are* measurable by expression when a biological mechanism drives
  them; it says nothing about ALT and is in a telomerase-positive lineage.

---

## RISKS TO THIS HYPOTHESIS

Ordered by how likely each is to produce a null or inverted result, each grounded in an
indexed source.

1. **Wrong measurement layer — the single largest risk.** Keegan established TCAB1
   inactivation by **SV calling and IHC**, and where he measured mRNA he found it
   **normal** (`10.1186/s12885-026-16250-5`, Supplemental Fig. 3). Three further indexed
   papers show TCAB1 function is decoupled from transcript: Jiang 2025 CCT5/CCT7 → TCAB1
   protein loss (`10.1016/j.stemcr.2025.102750`); Bergstrand 2020 HHS point mutants
   destabilised by TRiC with intact transcript (`10.1038/s41419-020-2421-4`);
   Sánchez-Morán 2020 function set by nuclear translocation (`10.1126/sciadv.abc5702`).
   This is index Prespec conflict #5, carried over unchanged into H1.
2. **Isoform conflation.** The lesions Keegan describes disrupt the α and γ transcription
   start sites while sparing β, which is the protein product that matters. Total WRAP53
   TPM cannot resolve this (Gadelha 2022, `10.3390/cimb44110372`). The index's Gap 5 also
   flags that nothing in the corpus addresses whether an unstranded protocol can resolve
   the 1,419-nt head-to-head TP53/WRAP53 overlap at all.
3. **Low WRAP53 does not imply loss of telomerase, and therefore does not imply ALT.**
   Wang 2022 (`10.1002/1878-0261.13324`): TCAB1-homozygous-deleted CCLE lines with
   significantly reduced expression retain equal or higher telomerase activity, with **no
   correlation** between TCAB1 protein and telomerase activity. If that holds in OS, the
   WRAP53-low arm will be a mix of ALT+ and telomerase+ lines.
4. **Insufficiency and the missing cofactor.** In Keegan's own model TCAB1 loss produces a
   *no-TMM* state; ALT requires additional ATRX/DAXX/SMARCAL1/SLX4IP hits. Only 1/16 ALT+
   lines in this panel is ATRX-altered. Two ALT− TARGET tumours carried TCAB1 SVs.
5. **Binary labels mis-sort exactly the relevant lines.** Wu 2025
   (`10.1038/s41467-025-67190-w`) documents TMM-DN, TMM-DP and EverShorter-Telomere lines,
   plus ALT/TEL subclonal mixtures. TCAB1-low/no-TMM lines are predicted to be labelled
   **ALT−**, which would flatten or invert the contrast. The silver labels come from
   Bennett 2026's WGS inference, which never uses WRAP53 — good for independence, but it
   also means the label carries its own ALT-definition error.
6. **Survivorship and range truncation.** Mahmoudi 2011 (`10.1038/cddis.2010.90`): WRAP53
   up to 20× higher in cancer than primary cells, knockdown lethal in U2OS; Niu 2021
   (`10.1186/s12935-020-01745-3`): TCAB1 knockdown → senescence; Yuan 2014
   (`10.3892/ol.2014.2306`): antiproliferative. A screened panel's WRAP53-low arm is a
   compensated, non-random minority.
7. **The real exposure may be TP53, which cannot discriminate.** TCAB1 SVs occur only in
   TP53-SV tumours and TP53 SVs appear in 36% of ALT− samples too. OS cell lines are
   near-uniformly TP53-defective.
8. **Every independent feature-selection exercise has already discarded WRAP53.** The
   Armendáriz-Castillo 28-gene ALT signature omits it (0 mentions,
   `10.3390/biology15171531`), Bennett's ALTitude omits it (0 mentions), and this repo's own
   19-gene model gave it coefficient 0.000 in all three cohorts with bootstrap frequency
   0.17/0.12 (`references/notebook_13_report.md`). Three independent nulls is weak evidence
   of a true null, but it is the only convergent evidence available, and it is all negative.
9. **Power against a mixture alternative.** Even at the arithmetic ceiling (|δ| ≈ 0.35) a
   median split on a panel of this size will not reach significance; the Stage 1 95% CI on
   Cliff's δ spanned [−0.313, +0.844].
10. **No corpus support for expression-as-SV-proxy, at all.** Index Gap 3: axis G of the
    OpenAlex scan returned no on-topic corpus; the pre-registered query drifted into maize,
    rice, yeast and Parkinson's. This remains the weakest-evidenced assumption inherited
    from Stage 1 and H1 does not fix it — it re-uses it.
11. **No OS-specific ALT-heterogeneity source exists in the corpus.** Index Gap 2: "Beyond
    ATRX: the heterogeneous basis of ALT in osteosarcoma" (`10.1016/j.cancergen.2026.06.001`)
    was never retrieved; every readable ALT-heterogeneity paper is glioma or pan-cancer. If
    OS-specific ALT classes exist, this design cannot see them. Also unretrieved and directly
    on point: `10.1080/15384101.2025.2609646` ("TCAB1: a key promoter of tumour growth
    through telomere maintenance and senescence evasion", 2026). Both are worth chasing by
    hand before freezing the prespec.

---

## Provenance

- `Literature/INDEX.md` read in full (846 lines, 85 indexed works + 42 unretrieved DOIs).
- All 156 PDFs under `Literature/` converted with `pdftotext` to a scratchpad cache and
  grepped for WRAP53/TCAB1/WDR79 × ALT co-occurrence; 10 files hit, all triaged.
- Sections read via `sed` on the extracted text, not via the Read tool, for: Keegan 2026
  BMC Cancer (abstract, IHC/ALT-status results, TARGET prevalence, discussion), Wang 2022
  Mol Oncol (abstract, §3.7–3.8, §3.11), Mahmoudi 2011 (results §1, Fig. 1 legend),
  Armendáriz-Castillo 2026 (methods, 28-gene list), Gadelha 2026 (abstract, results).
  Zero PDFs opened with the Read tool.
- `references/notebook_13_report.md` read in full (prior-project 19-gene signature).
- No data file loaded. No existing file modified.
