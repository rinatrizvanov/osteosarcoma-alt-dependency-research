# Analytical Pre-Specification v1.0 — Stage 1

**Project:** TCAB1-stratified break-induced replication (BIR) dependency in ALT+ osteosarcoma
**Author:** Rinat Rizvanov
**Status:** FROZEN at commit `81745ed` (v1.0). Amended to **v1.1** before any gene-effect
value was read — see §4d, §6.4 and the changelog. **v1.2** is a post-result factual
correction to §0 with no analytical content; see the changelog.
**Freeze commit:** _recorded in `data/processed/notebook_04_report.md` at run time_

Changes after the freeze require a new version number and a changelog entry. Section 11
bounds what may be run if the hypothesis is falsified; anything beyond that is a new
pre-specification, not an amendment to this one.

---

## 0. Hypothesis and what would change my mind

**Hypothesis.** Within ALT+ osteosarcoma cell lines, those with low TCAB1/WRAP53
expression show selectively elevated CRISPR essentiality for BIR-module genes (POLD3,
RAD52, SMARCAL1) relative to canonical HR genes (BRCA1, BRCA2, RAD51), compared with
ALT+ lines with intact TCAB1.

**Mechanism.** 17p13.1 structural variants (Keegan et al. 2026) disrupt TCAB1 and
WRAP53β together. TCAB1 loss blocks telomerase; WRAP53β loss cripples classical HR
(Henriksson & Farnebo 2015, *Front Genet* 6:91). An ALT+ cell that loses both the
telomerase escape route and the classical HR scaffold is forced onto non-classical BIR
machinery. Bennett et al. 2026 (ALTitude) established SMARCAL1 as a broad ALT dependency
**without** stratifying by TCAB1 or 17p13.1 architecture; this project adds that
stratification.

**Prior work this builds on.** Aim 2 of the adjacent project produced a 19-gene curated
ALT signature on DepMap osteosarcoma lines (AUROC 0.948 primary, N=20, permutation
p<0.001 at 1000 shuffles; TERT and FANCM the only genes stable in bootstrap resampling
under both definitions). Frozen copies of that work are in `references/`.

**What this pre-specification is for.** Every analytical choice below is committed before
any gene-effect *value* is read. `src/dependency_stats.py` is the only module that loads
gene-effect values, and `load_gene_effects` is first called in notebook 05.

Notebooks 02, 03 and 04 carry an explicit gate check asserting that
`src.dependency_stats` is absent from `sys.modules`, so they structurally cannot reach the
outcome; all three passed.

**Correction (v1.2).** An earlier wording of this paragraph claimed that *no* module used
by notebooks 00–04 imports `src.dependency_stats`. That was inaccurate: **notebook 01
imports it**, for `SIX_GENES` and `crispr_gene_columns`, in order to resolve the exact
`SYMBOL (entrez)` column names and count missing values. It calls no loader and retains no
gene-effect value — missingness is counted by scanning for empty fields with the raw `csv`
module. No outcome information therefore entered this pre-specification, and nothing
analytical changes. The claim is corrected rather than quietly dropped because a methods
document that overstates its own guarantees is worse than one that states a narrower
guarantee accurately.

---

## 1. Data

DepMap **26Q1**, local copies under `data/raw/depmap/`, all read-only (`chmod -w`) and
sha256-verified at load against `src.utils.EXPECTED_SHA256`:

| File | Role | sha256 (first 16) |
|---|---|---|
| `CRISPRGeneEffect.csv` | Chronos gene effect — **the outcome** | `e610a4cefb13a82b` |
| `OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv` | WRAP53/TP53 expression — the exposure | `0377be80c525fde9` |
| `Model.csv` | cell line metadata, osteosarcoma inclusion | `ea4e0b2a3bc806f8` |
| `OmicsCNGeneWGS.csv` | WRAP53 copy number — strand-check corroboration | `d2fe65b516802fb2` |
| `OmicsCNSegmentsWGS.csv` | focal 17p13.1 segments | `690e861b3a35db89` |
| `references/master_manifest.csv` | **frozen ALT labels** | `d51c5f180391ffbc` |

**Why 26Q1 and not a newer release.** The ALT labels in the frozen manifest were derived
against 26Q1. A newer release changes which lines are screened and which models exist,
which would silently break the label join. No 26Q2 public release was confirmed at
planning time (a forum thread title surfaced; no release notes or downloads page). If a
file must ever be re-fetched it is fetched from the **26Q1 archive specifically**, never
from "latest".

**Chronos scale.** Gene effect is already Chronos-scaled: 0 ≈ non-essential, −1 ≈ the
median common essential. No further normalisation is applied to the primary endpoint.

---

## 2. Cell line inclusion

Filters applied explicitly and in order, N logged after each step
(`src.cohort.FilterTrace`). The manifest's `include_*` flags are deliberately **not**
trusted — the same convention Aim 2 used.

| Step | N |
|---|---|
| Model.csv rows matching 'steosarcoma' | 35 |
| present in CRISPRGeneEffect.csv (screened) | 26 |
| present in expression matrix | 26 |
| joined to frozen manifest (cohort == 'DepMap_26Q1') | 26 |
| alt_status in ('alt_gold_positive', 'alt_silver', 'alt_negative_high_conf') | 25 |
| alt_call_direction in ('positive', 'negative') | 25 |
| ALT+ (primary contrast set) | 16 |
| ALT- (negative-control set) | 9 |

Yielding three sets, fixed before any analysis:

- **ALT+ primary contrast set, N=16** (10 gold + 6 silver)
- **ALT− negative-control set, N=9** (7 high-confidence negative + 2 silver)
- **Excluded, N=1** (`alt_unknown`: ACH-002778 / NCCOS1X2C1)

16 + 9 + 1 = 26 = every CRISPR-screened osteosarcoma line. Nothing is unaccounted for.

`src.cohort.assert_all` aborts the analysis if any of these counts, the tier breakdown,
or the ground-truth anchors differ at run time. A mismatch means the release or the
manifest changed; the correct response is to investigate, never to adapt the expectation
to the data.

### ALT+ set (N=16)

| ModelID | Cell line | ALT tier | ALT source | ATRX |
|---|---|---|---|---|
| ACH-000082 | G292CLONEA141B1 | alt_gold_positive | ALTitude_2024 | not_altered |
| ACH-000364 | U2OS | alt_gold_positive | ALTitude_2024 | not_altered |
| ACH-000410 | SAOS2 | alt_gold_positive | ALTitude_2024 | not_altered |
| ACH-001526 | HUO9 | alt_gold_positive | ALTitude_2024 | not_altered |
| ACH-001715 | CAL72 | alt_gold_positive | ALTitude_2024 | not_altered |
| ACH-002067 | NOS1 | alt_gold_positive | ALTitude_2024 | not_altered |
| ACH-002471 | OS052 | alt_gold_positive | ALTitude_2024 | not_altered |
| ACH-002688 | OS457 | alt_gold_positive | ALTitude_2024 | not_altered |
| ACH-002690 | OS526 | alt_gold_positive | ALTitude_2024 | not_altered |
| ACH-002834 | OS384 | alt_gold_positive | ALTitude_2024 | not_altered |
| ACH-002433 | CCLFPEDS0019T | alt_silver | ALTitude_2024_predicted | not_altered |
| ACH-003171 | NRHOS1 | alt_silver | ALTitude_2024_predicted | altered |
| ACH-003178 | NRHOS2 | alt_silver | ALTitude_2024_predicted | not_altered |
| ACH-003185 | NRHOS4 | alt_silver | ALTitude_2024_predicted | not_altered |
| ACH-003573 | C337 | alt_silver | ALTitude_2024_predicted | not_altered |
| ACH-003587 | 531MII | alt_silver | ALTitude_2024_predicted | not_altered |
### ALT− negative-control set (N=9)

| ModelID | Cell line | ALT tier | ALT source | ATRX |
|---|---|---|---|---|
| ACH-000359 | MG63 | alt_negative_high_conf | ALTitude_2024 | not_altered |
| ACH-000748 | SJSA1 | alt_negative_high_conf | ALTitude_2024 | not_altered |
| ACH-001001 | 143B | alt_negative_high_conf | ALTitude_2024 | not_altered |
| ACH-001814 | OS252 | alt_negative_high_conf | ALTitude_2024 | not_altered |
| ACH-002069 | HSOS1 | alt_negative_high_conf | ALTitude_2024 | not_altered |
| ACH-002142 | HUO3N1 | alt_negative_high_conf | ALTitude_2024 | not_altered |
| ACH-002470 | OS186 | alt_negative_high_conf | ALTitude_2024 | not_altered |
| ACH-001275 | SYO1 | alt_silver | ALTitude_2024_predicted | not_altered |
| ACH-003176 | NRHGCT1 | alt_silver | ALTitude_2024_predicted | not_altered |
---

## 3. ALT classification

Labels are **taken from the frozen manifest, not re-derived**: `alt_call_direction` gives
the arm, `alt_status` gives the tier. The 19-gene signature was the *product* of Aim 2,
not a classifier for unlabelled lines, and every line here already carries a curated call
from ALTitude plus literature.

Ground-truth anchors, asserted before any analysis, abort on failure:

- ACH-000364 (U2OS) → `positive` — canonical ALT-positive
- ACH-000748 (SJSA1) → `negative` — canonical ALT-negative / telomerase-positive

The 19-gene signature is scored on these lines as an **orthogonal QC overlay only** (does
the signature score separate the label arms?). It never assigns or revises a label.

**Primary set is gold+silver (N=16); gold-only (N=10) is a reported sensitivity analysis
(R1), not the primary.** Rationale: at gold-only N=10 the split is 3v7 or 5v5, where the
smallest attainable two-sided exact Mann-Whitney p is 0.0167 or 0.0079 respectively — so
after BH correction across six genes at most one gene could ever clear FDR 0.05, and only
on perfect separation. N=16 makes the primary estimable. The cost is that 6 of 16 ALT+
lines carry ALTitude prediction-only evidence; R1 is the guard against a silver-driven
result, and this limitation belongs in any write-up.

---

## 4. TCAB1-low vs TCAB1-intact split

### 4a. Primary rule

1. WRAP53 log2(TPM+1) for **all osteosarcoma lines with RNA** (N=34) — not only the 16
   ALT+ lines entering the contrast.
2. `cut` = 33.3rd percentile (bottom tertile) of that distribution.
3. Within the ALT+ set: `TCAB1-low` if WRAP53 < `cut`, else `TCAB1-intact`.
4. Resulting group sizes are an empirical result and are reported, not fixed in advance.

Estimating the cutpoint on 34 lines rather than 16 makes it both better determined and
independent of which lines happen to be ALT+, so it cannot be tuned toward the outcome.
`tests/test_stratify.py::TestCutpointIsOutcomeIndependent` pins that property.

### 4b. Minimum-cell guard

If either arm has fewer than **4** lines, the tertile primary is declared NOT ESTIMABLE
and the analysis falls back to a **median split within the ALT+ set**. The switch and its
reason are recorded in the notebook 04 report.

This rule is written here before the WRAP53 distribution has been examined, which is what
makes the fallback pre-specified rather than post-hoc.

### 4d. Guard outcome (recorded 2026-10-11, v1.1 — before any outcome was read)

**The §4b minimum-cell guard fired.** Exactly **3 of 16** ALT+ lines fall below the
bottom-tertile cut of 3.3356 computed over the 34 osteosarcoma lines with RNA, and the
guard threshold is <4 per arm. Per §4b the tertile primary is therefore declared NOT
ESTIMABLE and the primary becomes the **median split within the ALT+ set (8v8)**.

This is the guard behaving exactly as written. It was committed in v1.0 before the WRAP53
distribution was examined, so the fallback is pre-specified, not post-hoc.

| | Rule | Split | min \|delta\| for p<0.05 |
|---|---|---|---|
| **Primary** | median within ALT+ | **8 v 8** | 0.594 |
| **Pre-specified sensitivity** | §4a OS-reference tertile (superseded) | **3 v 13** | ≈0.90 |

The superseded §4a tertile contrast is **retained as a pre-specified sensitivity
analysis**, not demoted to a secondary: it was the primary rule in v1.0 and it is the
biologically cleanest contrast available. Notebook 05 runs both in parallel, with the same
test and the same output format, and both enter the final report.

**Exposure-side observation, interpretation withheld.** The bottom of the WRAP53
distribution is occupied by two ALT+ lines with homozygous WRAP53 loss:

| Line | ALT tier | WRAP53 log2(TPM+1) | TP53 log2(TPM+1) | WRAP53 CN |
|---|---|---|---|---|
| ACH-002834 / OS384 | gold+ | 0.000 | 0.000 | 0.00009 |
| ACH-003185 / NRHOS4 | silver+ | 0.000 | 0.106 | 0.00001 |

Recorded here as a measurement fact with its copy-number corroboration. **No mechanistic
claim is drawn from it before the outcome is read.** Both lines must appear in the
TCAB1-low arm of both the 8v8 and 3v13 contrasts; notebook 04 asserts this and aborts if
either drops out.

The CN/expression concordance at both extremes (CN≈0 ↔ expression 0; CN≈2.35 ↔ expression
≈5.3–5.9) is independent corroboration of the §4c `direct` branch: antisense
cross-assignment would degrade that concordance.

### 4c. Strandedness branch — decision rule committed in advance

WRAP53 exon 1 carries a ~1,419 bp head-to-head antisense overlap with TP53. On an
unstranded library, reads in that interval cannot be assigned to the correct gene, so
WRAP53 TPM would be contaminated by TP53 — which would invalidate the stratifier.

Notebook 02 runs three diagnostics and the branch is selected by the rule below.
**The rule is committed here; notebook 02 only evaluates it.** Notebook 02 reads no
gene-effect data, so this runs before the freeze without compromising it.

Diagnostics: (i) Spearman ρ(WRAP53, TP53) across all lines and within osteosarcoma lines;
(ii) the same ρ against a null of 1,000 **expression-matched** genes — matching matters
because ρ is attenuated for low-expressed genes and an unmatched null would be biased
downward; (iii) whether WRAP53 expression tracks WRAP53 **copy number** better than it
tracks TP53 expression, copy number being strand-agnostic.

| Condition | Branch | Stratifier |
|---|---|---|
| ρ(WRAP53,TP53) below the null's 99th percentile | `direct` | WRAP53 log2(TPM+1) as in 4a |
| ρ exceeds the null **but** WRAP53 tracks its own CN more strongly than it tracks TP53 | `direct_cn_corroborated` | WRAP53 log2(TPM+1), CN concordance reported alongside |
| ρ exceeds the null **and** CN does not track better | `fallback_cn_conjunction` | conjunction of low WRAP53 expression **and** WRAP53 CN; TP53-residualised split as sensitivity |

**Hard constraint, flagged.** Exon-aware requantification masking WRAP53 exon 1 is
**not possible on public DepMap data** — DepMap distributes no FASTQ/BAM, and CCLE
alignments are controlled-access. The CN-conjunction fallback replaces it. If the
CN-conjunction and TP53-residualised stratifiers disagree on the arm assignment of more
than 3 of 16 lines, Stage 1 cannot test the hypothesis on public data; that is reported
as **"proxy invalid, hypothesis untested"**, never as a falsification (see 11.3).

---

## 5. Primary analysis

**Endpoint.** One composite per cell line:

```
composite_i = mean(POLD3, RAD52, SMARCAL1)_i − mean(BRCA1, BRCA2, RAD51)_i
```

**Directional prediction.** Chronos gene effect is more negative when a gene is more
essential. TCAB1-low lines are predicted to be *more* dependent on BIR relative to HR, so
their composite is **more negative**:

```
median(composite | TCAB1-low)  <  median(composite | TCAB1-intact)
Cliff's delta (low vs intact)  <  0
```

**Test.** Mann-Whitney U, two-sided, exact null
(`scipy.stats.mannwhitneyu(..., method='exact')`). Exact p-values are invalid under ties,
so if ties are present the test falls back to the asymptotic form and records which was
used — the choice is never silent.

**Why two-sided with a one-sided success rule.** Reporting a two-sided p keeps the
p-value interpretable, while requiring the predicted direction for support prevents a
significant *opposite* result from being read as confirmation.

**Effect size.** Cliff's delta with a 2000-resample percentile bootstrap 95% CI, **and**
the Hodges-Lehmann median shift in Chronos units. Cliff's delta is validated in the test
suite against the identity `delta = 2U/(nm) − 1`.

**Multiplicity.** None — this is one test.

**Why a composite rather than six tests.** A line with globally more negative gene
effects shifts BIR and HR together, so the difference is internally controlled for screen
quality and line-level essentiality scale. The test suite pins this: adding an arbitrary
per-line offset leaves the composite unchanged.

### 5a. Minimum detectable effect — read this before interpreting any result

Computed exhaustively over rank configurations, the smallest |Cliff's delta| that can
reach two-sided exact p<0.05:

| Split | min \|delta\| for p<0.05 | min attainable p |
|---|---|---|
| 5v11 (expected tertile) | **0.673** | 0.00046 |
| 8v8 (median fallback) | **0.594** | 0.00016 |
| 3v7 (gold-only tertile) | 0.905 | 0.0167 |
| 5v5 (gold-only median) | 0.840 | 0.0079 |

**Consequence, stated before the analysis runs.** At the expected 5v11 split the p<0.05
requirement alone forces |delta| ≥ 0.673 — roughly 84% of low-versus-intact pairs ordered
as predicted. The Section 8.1 magnitude threshold of 0.33 is therefore **not the binding
constraint** and is retained only for completeness. Stage 1 can return `supported` only
for a very large effect; a genuine moderate effect (|delta| 0.40–0.67) lands in the
Section 9 inconclusive band by construction. This is a property of N=16, not a flaw in
the hypothesis, and it is the honest reason Stage 2 exists.

---

## 6. Secondary analyses (pre-specified)

**6.1 Six per-gene comparisons.** Mann-Whitney per gene with Cliff's delta and
Hodges-Lehmann, BH-FDR across exactly those six. Visualised as per-gene strip plots and a
forest plot of effect sizes with CIs.

**6.2 ALT × TCAB1 selectivity.** The same composite and the same split rule applied to
the 9 ALT− lines, where the hypothesis predicts **no** shift. Reported as: each stratum's
own Cliff's delta and Mann-Whitney p, the delta difference, and a 5000-shuffle
permutation test that shuffles the TCAB1 label **within** each ALT stratum (preserving
stratum sizes and the ALT main effect, so the p-value speaks only to the interaction).

> **Power caveat, measured not assumed.** At the ALT− group sizes (3v6 or 4v5) the exact
> null distribution of |Cliff's delta| places 26–29% of its mass at or above 0.474, and
> its median is ≈0.30. A noise-only ALT− stratum therefore routinely produces a large
> delta. This is why Section 8.5 gates on *significance* in the ALT− stratum rather than
> on a delta point comparison — see 8.5. The delta difference and interaction p-value are
> reported as diagnostics, never as gates.

**6.4 Concentration comparison (v1.1, pre-registered before any outcome was read).**
Compare |Cliff's delta| for the 3v13 tertile contrast against the 8v8 median primary. If
the tertile exceeds the median primary by more than **0.300** in |delta|, in the predicted
direction, the signal concentrates in the extreme-WRAP53-loss lines rather than being
spread across the diluted median split.

**Descriptive and explicitly non-gating.** It does not enter the §8/§9/§10 verdict and
cannot change the outcome bucket. It is reported in notebook 07 alongside the verdict
because it bears on whether Stage 2 is justified: a concentrated signal argues that real
17p13.1 SV calls would sharpen the contrast, which is precisely what the PDX cohort
provides. Implemented as `src.criteria.concentration_comparison`; threshold frozen as
`CONCENTRATION_DELTA_THRESHOLD`.

**6.3 Within-line standardised sensitivity.** Recompute the composite from gene effects
z-scored across all ~18.5k genes within each line, then re-run the primary. A stronger
control for line-level essentiality scale than the composite difference alone.

---

## 7. Robustness and sensitivity (pre-specified)

| # | Check | Pre-specified expectation |
|---|---|---|
| R1 | **Gold-only ALT+ (N=10) — reported in main results** | same direction, \|delta\| ≥ 0.30 |
| R2 | Median split within ALT+ (8v8) | direction agrees |
| R3 | Tertile computed within ALT+ only (5v11) | direction agrees |
| R4 | Continuous: Spearman ρ(WRAP53 expression, composite) across the 16 | ρ > 0 |
| R5 | ATRX control: re-run within ATRX-wildtype ALT+ lines | direction agrees |
| R6 | Lineage: constant (all osteosarcoma), so the lineage control **is** the ALT− contrast in 6.2 | no shift in ALT− |
| R7 | Leave-one-line-out: 16 refits, report the delta range | no single line flips the sign |
| R8 | Drop SMARCAL1 from the BIR set and re-run | direction agrees |

**R5 is expected to be near-vacuous and that is recorded now, not discovered later.**
Only **1 of 16** ALT+ lines carries `atrx_status == altered` in the frozen manifest
(ACH-003171 / NRHOS1), so the ATRX-wildtype subset is 15 of 16 and R5 cannot
meaningfully control for ATRX. The likely cause is that Aim 2 derived `atrx_status` from
damaging point mutations only (`OmicsSomaticMutationsMatrixDamaging`), which misses ATRX
deletions and loss of expression — the dominant ATRX-inactivation routes in ALT. R5 is
retained because it was pre-specified, and its result will be reported as uninformative
rather than as evidence of no confounding. **Flagged for independent verification.**

**R8 matters for novelty.** If the composite is driven by SMARCAL1 alone, the result is a
TCAB1-stratified refinement of an already-published finding (Bennett et al. 2026) rather
than a new BIR-module claim. Both are reportable; they are not the same claim.

---

## 8. Success criterion — hypothesis supported

Evaluated mechanically by `src.criteria.evaluate`, which notebook 07 calls and reports
verbatim. **All** of:

1. **8.1** Primary p < 0.05 **and** Cliff's delta ≤ −0.33, in the predicted direction.
   (Per 5a, p<0.05 is the binding constraint; the magnitude clause is not.)
2. **8.2** R1 gold-only: same direction, |delta| ≥ 0.30.
3. **8.3** ≥2 of 3 BIR genes at BH q < 0.10 in the predicted direction, **and** no HR gene
   showing a stronger same-direction effect than the strongest BIR gene.
4. **8.4** R7: no sign flip on any leave-one-line-out refit.
5. **8.5** The ALT− stratum shows **no significant same-direction shift** (its own
   two-sided Mann-Whitney p ≥ 0.05).

> **8.5 was revised before freezing.** It originally read "ALT− delta closer to zero than
> ALT+ delta". Given the exact null in 6.2, that formulation would have failed roughly
> 20–30% of the time *even when the hypothesis is true*, because a noise-only ALT−
> stratum frequently out-deltas a genuine ALT+ effect at n=4v5. The non-contradiction
> form has a 5% false-failure rate by construction and still catches what the criterion
> exists to catch: a TCAB1 effect that is general rather than ALT-specific.

---

## 9. Inconclusive band — explicitly a third outcome

Primary p ≥ 0.05 **but** Cliff's delta ≤ −0.40 in the predicted direction, with R1, R2 and
R7 directionally consistent → **"directionally consistent, N-limited"**.

Not a positive result. Not falsified. Per 5a this is the *likely* outcome for any genuine
moderate effect, and it is the finding that justifies Stage 2, because the PDX cohort both
adds N and replaces the expression proxy with real 17p13.1 SV calls. Pre-specifying this
band is what stops it being written up afterwards as either a win or a loss.

---

## 10. Failure criterion — hypothesis falsified, stop

Primary p ≥ 0.05 **and** |delta| < 0.20 (a null, not merely an underpowered result),
**or** a significant effect in the **opposite** direction; **and** no BIR gene at q < 0.10.

On falsification: run only the Section 11 set, write `docs/negative_result_memo.md`, stop.
Specifically — **do not** approach the collaborating lab for Stage 2 PDX access on a falsified Stage 1,
and do not expand the gene set, swap the composite, or add cohorts hunting for signal.

A clean negative here is informative: given that SMARCAL1 is an established broad ALT
dependency, it would say TCAB1 status does **not** modulate BIR dependency in ALT+
osteosarcoma lines, which constrains the Keegan model.

**If the result matches none of Sections 8, 9 or 10**, `src.criteria.evaluate` returns
`indeterminate` and names the failed conditions. That is reported as-is and not rounded
toward the nearest verdict.

---

## 11. Negative-result fallback — bounded

Permitted only after Section 10 triggers, and only these three:

1. **Broader BIR/ALT gene set**, declared here so it cannot be assembled post-hoc, BH-FDR
   across the whole set, labelled exploratory: POLD3, POLD1, POLD2, RAD52, SMARCAL1,
   SLX4, BLM, WRN, FANCM, BRCA1, BRCA2, RAD51, RAD51AP1, MUS81, GEN1, NBN, MRE11, RAD50,
   ATRX, DAXX. (Frozen as `src.dependency_stats.BROADER_SET`.)
2. **Telomerase-status re-slice.** Re-run the primary split on `telomerase_status` from
   the frozen manifest instead of WRAP53 expression. Tests whether the proxy, rather than
   the hypothesis, failed.
3. **Expression-proxy validity audit.** If notebook 02 selected
   `fallback_cn_conjunction`, or the 4c disagreement threshold was exceeded, the negative
   is uninterpretable and the memo records **"proxy invalid, hypothesis untested"** rather
   than "hypothesis falsified". This distinction is the single most important line in that
   memo.

No fourth analysis. Anything further is a new pre-specification with a new version number.

---

## 12. Limitations binding on any deliverable

1. **WRAP53 expression is a proxy for 17p13.1 structural disruption, not a measurement of
   it.** Keegan et al. 2026 established the SVs; DepMap 26Q1 ships no SV calls
   (`references/label_spec_v1.0.md` §6.1). This is the central limitation of Stage 1 and
   the reason Stage 2 exists.
2. **N=16.** See 5a for the minimum detectable effect. The study is powered only for a
   large effect.
3. **6 of 16 ALT+ lines carry ALTitude prediction-only (silver) evidence.** R1 is the
   guard.
4. **The ALT− negative control (N=9) is too small to gate on.** See 6.2.
5. **ATRX confounding is not controlled** in practice — see R5.
6. **Cell lines are not tumours.** Any claim is about dependency in cultured osteosarcoma
   lines under a genome-wide CRISPR screen.
7. **DepMap strandedness is not documented.** The authoritative Broad pipeline
   documentation (`depmap_omics/documentation/DepMap_processing_pipeline.md`) and the 26Q1
   release notes do not state RNA-seq library strandedness; the per-quarter pipeline
   configs are not publicly reachable; only non-authoritative third-party data
   dictionaries describe RSEM "unstranded mode", with the batch-corrected variant
   described as stranded-or-unstranded per sample. The `direct` branch therefore rests on
   the §4c **empirical** diagnostics, not on documentation. This must be stated as a
   limitation in any write-up, with the diagnostics given as the evidence base.
8. **The primary split is the 8v8 median fallback, not the tertile.** See §4d. The low arm
   contains lines with mid-range WRAP53 that may not be TCAB1-disrupted, which biases the
   primary toward the null. §6.4 is the pre-registered check on that dilution.

---

## 13. Changelog

| Version | Date | Notes |
|---|---|---|
| v1.2 | 2026-10-11 | **Factual correction only, no analytical change.** §0 previously claimed no module used by notebooks 00–04 imports `src.dependency_stats`; notebook 01 does import it for column-name resolution and missing-value counting, while loading no gene-effect value. The accurate guarantee — notebooks 02/03/04 carry a passing gate check, and `load_gene_effects` is first called in notebook 05 — is now stated. Entered after the result was known and recorded as such; it corrects a description of the code, not a decision rule. |
| v1.1 | 2026-10-11 | Amended **before any gene-effect value was read**. §4d records that the §4b minimum-cell guard fired (3 of 16 ALT+ lines below the tertile cut), making the 8v8 median split the primary and retaining the 3v13 OS-reference tertile as a pre-specified sensitivity run in parallel; records the two WRAP53-null ALT+ lines as an exposure-side observation with interpretation withheld, and requires both to appear in the TCAB1-low arm of both contrasts. §6.4 pre-registers the concentration comparison (tertile vs median \|delta\| excess > 0.300), explicitly descriptive and non-gating. Limitations 7 and 8 added. §8/§9/§10 success, inconclusive and failure criteria are **unchanged**. |
| v1.0 | 2026-10-11 | Initial pre-specification. Primary ALT+ = gold+silver (N=16) with gold-only (N=10) as reported sensitivity; primary endpoint = BIR−HR composite, one two-sided exact Mann-Whitney; TCAB1 split = bottom tertile of WRAP53 over 34 osteosarcoma lines with a <4-per-arm median-split guard; DepMap 26Q1. Section 8.5 revised from a delta point comparison to a non-contradiction test before freezing, on the basis of the exact null at the ALT− group sizes. Section 5a added to state the minimum detectable effect before any result is seen. |
