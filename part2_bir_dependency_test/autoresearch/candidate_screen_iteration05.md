# Candidate screen for iteration 05 — the result is that no iteration 05 exists

**Date:** 2026-10-12 · **Outcome: stop condition (c), hypothesis-space exhaustion.**

This file exists because a negative screening result needs the same auditability as a
negative analysis. It records every candidate considered for iteration 05, the filter it
failed, and the number that killed it. No iteration directory was created, because no
hypothesis was run.

## The four filters

A candidate had to pass all four:

1. Sourced from `Literature/INDEX.md`'s **Prespec conflicts** table or an **untried** row
   of `Literature/openalex_scan/next_hypotheses.md`.
2. Testable on DepMap 26Q1 as shipped — no missing-data workaround, no invented proxy.
3. At least one specific paper indexed in `Literature/INDEX.md` as mechanistic backing.
4. Not a variant of anything already in `NULL_SPACE.md`.

Thirteen candidates: five `next_hypotheses.md` rows and eight Prespec-conflict rows.
**Zero passed all four.**

## Feasibility evidence generated for this screen

One new computation, deliberately **marginals only — no correlation of any kind** —
following `NULL_SPACE.md`'s binding fact #2 ("check marginal regimes before composing any
multi-gene score") and iteration 03's precedent of recording dependency regimes pre-freeze.
Chronos gene effect over the 16 ALT+ CRISPR-screened lines:

| gene | role in a candidate | median | SD | min | n < −0.5 | iteration-03 class |
|---|---|---|---|---|---|---|
| POLQ | Batenburg restart arm | −0.409 | 0.178 | −0.650 | 5/16 | §7.1 usable |
| LIG3 | Batenburg restart arm | −0.143 | 0.136 | −0.336 | **0/16** | §7.3 non-essential throughout |
| PRIMPOL | Batenburg restart arm | +0.132 | 0.149 | −0.236 | **0/16** | §7.3 non-essential throughout |
| BRCA1 | HR arm | −0.555 | 0.440 | −1.841 | 9/16 | §7.1 usable |
| BRCA2 | HR arm | −0.406 | 0.375 | −1.221 | 7/16 | §7.1 usable, but barred (below) |
| RAD51 | HR arm | −1.127 | 0.317 | −1.821 | 16/16 | §7.4 near-saturated |
| WRAP53 | Henriksson candidate | −0.173 | 0.247 | −0.570 | **1/16** | thin |
| CCT5 | Jiang chaperonin candidate | −1.203 | 0.253 | −1.580 | 16/16 | §7.4 near-saturated |
| CCT7 | Jiang chaperonin candidate | −0.834 | 0.200 | −1.272 | 15/16 | §7.4 near-saturated |
| CCT2 | Jiang chaperonin candidate | −1.453 | 0.300 | −2.118 | 16/16 | §7.4 near-saturated |
| TCP1 | Jiang chaperonin candidate | −1.988 | 0.441 | −2.542 | 16/16 | §7.4 near-saturated |
| TERT | telomerase axis | +0.002 | 0.114 | −0.179 | **0/16** | flat |

Also computed and consistent with iteration 03: MUS81 −0.378/0.153/3-of-16,
LIG1 −0.295/0.164/3, PALB2 −0.365/0.234/3, RAD54L −0.354/0.192/3, XRCC3 −0.533/0.357/8,
RNF8 −0.532/0.312/10, RNF168 −0.556/0.261/10, TP53BP1 −0.106/0.164/1, PARP1 −0.296/0.182/2,
PARP2 −0.030/0.118/0, XRCC1 −0.360/0.113/1, POLB +0.115/0.130/0, EXD2 +0.143/0.147/0,
NBN −0.066/0.163/0, RAD51AP1 +0.045/0.128/0, FANCM −1.020/0.388/16, POLD3 −1.666/0.319/16,
RAD52 −0.010/0.173/0, SMARCAL1 −0.462/0.395/8, ATRX −0.000/0.369/2, DAXX −0.388/0.297/4.

Script: `scratchpad/feas05.py`, `scratchpad/feas05b.py` (session scratch, not committed —
the numbers above are the record, and every one reproduces from
`src.cohort.build_cohort()` plus `src.dependency_stats.load_gene_effects`).

## `next_hypotheses.md` — five rows

| H | Filter failed | Reason |
|---|---|---|
| H1 TCAB1 loss predicts ALT status | **4** | Run as iteration 01. `NULL_SPACE` tested-and-null (`counter_signal`, δ = +0.241 against prediction, p = 0.280). |
| H2 SMARCAL1 dependency gated by ATRX | **4** | Run as iteration 02, `not-estimable` (best rule 4v12 with 3 of 4 calls unsupported by segments). Independently closed by Bennett 2026's published pan-cancer null (d = −0.248, P = 0.5, n = 20 v 14). |
| H3 ALT branch identity determines essentiality | **2** | Requires a branch label that does not exist and would have to be invented. `NULL_SPACE` establishes that this cohort carries **no ALT-activity measure of any kind** (`alt_assay_evidence` unknown for all 35 rows; no C-circle, APB, TRAP or telomere-length column; `telomerase_status` constant within ALT+) and that **no candidate marker gene is bimodal** over the 34 RNA-backed lines (largest guard-respecting gap/SD 0.22–0.90, all 16 tested). TERT gene effect is flat here too (0/16, SD 0.114), so there is no dependency-side telomerase axis either. |
| H4 POLD3/RAD52/SMARCAL1 are not a module | **4** | Run as iteration 03. |
| H5 WRAP53 expression is a poor SV proxy | **2** | Needs WGS SV calls. DepMap 26Q1 ships none. Already `NULL_SPACE` skipped-missing-data. |

## Prespec conflicts table — eight rows

### Row 1 — SMARCAL1 LOF *promotes* ALT (Liu 2023; Brosnan-Cashman 2021; Sobieszek 2026)

Passes 1 and 3. The osteosarcoma-specific re-expression result makes this the most
attractive row in the table. **Fails 2 and 4** on all three available operationalisations:

- **SMARCAL1 mRNA as exposure** is barred on literature grounds already established in
  this campaign. Iteration 04's frozen prespec records Bennett's documented
  transcript-to-protein dissociation for SMARCAL1 itself (line SKNFI) and states that
  "SMARCAL1 transcript level cannot appear as exposure, covariate or control." That
  constraint was written pre-freeze, against interest, and reviving the variable now to
  keep the loop alive would be exactly the post-hoc rule change §3 forbids.
- **SMARCAL1 copy number as exposure** is a copy-number claim class the campaign has
  already excluded once on stated grounds (TOP3A, iteration 04 prespec), is confounded with
  the known CRISPR copy-number artefact when the outcome is the same gene's Chronos score,
  and is 0.716 correlated with the barred expression variable anyway
  (`iteration_04/scan_exposure.md`).
- **The direct ALT+ vs ALT− contrast on SMARCAL1 dependency is contaminated.** The values
  are already on disk and have been read this campaign: `iteration_03/scan_essentiality.md`
  §3 and §4 give ALT+ mean −0.4907 (SD 0.395, 8/16 below −0.5) versus ALT− mean +0.0557
  (SD 0.137, 0/9). Declared in `iteration_04/DISCLOSURES.md`. It is also Taglialatela's and
  Bennett's published headline, so it is corroboration, not a campaign test. **Recorded as
  a descriptive observation in `FINAL.md`, not as a result.**

### Row 2 — WRAP53β loss impairs HR *and* NHEJ, non-selectively (Henriksson 2014)

Passes 1 and 3. **Fails 2.** The selective version is removed outright by the loop's §6
hard constraint. The non-selective version would have to be recast as a co-dependency
prediction (WRAP53 co-essential with the RNF8/ubiquitin arm), which translates a
perturbation-biology claim into a correlational one — proxy invention of the kind §6
forbids. Independently, the channel is not there: WRAP53 Chronos is −0.173 median,
SD 0.247, with **1 of 16** lines below −0.5, so its line-to-line variation sits inside the
non-essential band. WRAP53's dominant cellular role here is Cajal-body/scaRNA trafficking
(Mahmoudi 2010), not DSB repair, so even a positive correlation would not be attributable.

### Row 3 — BRCA2 loss *drives* ALT (Lee 2024)

Passes 1 and 3. **Fails 4.** Withdrawn in `iteration_04/DISCLOSURES.md` before iteration 04
was designed, because iteration 03's scan had already reported BRCA2's ALT+/ALT− marginals
(−0.537 v −0.509, no difference). The loop's §6 also bars BRCA2 from any HR-neutral set.

### Row 4 — the BIR and HR arms are coupled with *opposite sign* (Batenburg 2023)

**This was the strongest survivor on paper, and it is the one this screen had to test
rather than assume.** It passes 1 and 3 (Batenburg 2023, INDEX #11; supported by Richards
2023 #17 and Shah 2024 #15), and it is genuinely *not* a variant of iteration 03: the
statistic is cross-set rather than within-set, the predicted sign is negative rather than
positive, and no cross-arm pairwise value exists anywhere in the campaign's outputs
(`iteration_03/data/h4_pairwise.csv` holds only the three within-triplet pairs). The
hypothesis would have been: dependency on the HR arm and on the alternative-restart arm are
negatively coupled across the 16 ALT+ lines, so the Stage 1 difference endpoint *added* two
anticorrelated arms rather than contrasting two orthogonal ones.

**It fails filter 2 on measurement, established above.** Of Batenburg's three named restart
genes, two — LIG3 (median −0.143, SD 0.136, **0/16** below −0.5) and PRIMPOL (+0.132,
0.149, **0/16**) — fall in iteration 03 §7.3, the class that scan defined as spread that
"is not dependency signal." Only POLQ carries dependency variation. The HR arm is usable and
is a real module here (iteration 03: internal T = +0.414, p = 0.0039), but BRCA2 is barred
by row 3, leaving BRCA1 + RAD51.

So the test would correlate one coherent module against a set that is two-thirds screen
noise. **The expectation of that statistic is approximately zero whether or not Batenburg's
claim is true**, which makes a null uninformative by construction — the definition of
not-estimable, and precisely the construct defect `NULL_SPACE` binding fact #2 was written
to prevent repeating. The one clean cross-pair, BRCA1 × POLQ, has a detection floor of
|ρ| = 0.494 at n = 16 and would be a single cell selected from a 3 × 3 grid after seeing
which cells were usable.

### Row 5 — TCAB1 protein level is set post-translationally by CCT/TRiC (Jiang 2025; Bergstrand 2020)

Passes 1 and 3. **Fails 2.** The claim is about protein, and DepMap 26Q1 ships no
proteomics (loop §1). Recast as a dependency model it dies on dynamic range: CCT5 16/16,
CCT7 15/16, CCT2 16/16 and TCP1 16/16 lines below −0.5, all in the near-saturated
housekeeping class — the chaperonin is essential in every line, as expected for a
chaperonin. The row's real content is a measurement-validity verdict on the Stage 1
exposure, which `NULL_SPACE` already records.

### Row 6 — SMARCAL1 dependency stratifies on ATRX (Taglialatela 2026; Brown 2026)

**Fails 4.** Iteration 02 plus Bennett's published pan-cancer null. ATRX Chronos here is
flat anyway (median −0.000, 2/16 below −0.5).

### Row 7 — ALT status is not binary; TMM-DN and TMM-DP lines exist (Wu 2025)

Passes 1 and 3, and it is a real criticism of the design. **Fails 2.** Acting on it means
either re-deriving ALT labels — forbidden: the manifest is frozen and "ALT labels live here
and are never re-derived" — or importing Wu's external 976-line resource, which is outside
the permitted data. The honest treatment is to state the limitation, which `FINAL.md` does.

### Row 8 — TCAB1 SV inactivation is an early, ALT-*initiating* event (Keegan 2026)

**Fails 4.** This is H1, run as iteration 01.

## Conclusion

No candidate passes all four filters. The binding constraint is not imagination and not
sample size — iteration 03 proved co-essentiality is estimable at n = 16 (BTR T = +0.526,
p = 0.0001). It is that **every remaining literature-grounded hypothesis needs either a
variable this cohort does not contain (an ALT-activity readout, an SV call, a protein
level, a branch label) or a gene whose dependency does not vary in these 16 lines.**

Per the loop's own rule — *"Do not fabricate further hypotheses to keep the loop alive"* —
iteration 05 is not run and **`FINAL.md` is written on stop condition (c)**.
