---
marp: true
theme: default
paginate: true
style: |
  section {
    font-family: 'Helvetica Neue', Arial, sans-serif;
    font-size: 22px;
    padding: 40px 50px;
    color: #1a1a2e;
  }
  h1 {
    font-size: 38px;
    color: #0f3460;
    border-bottom: 3px solid #e94560;
    padding-bottom: 10px;
  }
  h2 {
    font-size: 30px;
    color: #0f3460;
    border-left: 5px solid #e94560;
    padding-left: 12px;
    margin-bottom: 20px;
  }
  h3 {
    font-size: 22px;
    color: #16213e;
  }
  section.title {
    background: linear-gradient(135deg, #0f3460 0%, #16213e 100%);
    color: white;
    display: flex;
    flex-direction: column;
    justify-content: center;
  }
  section.title h1 {
    color: white;
    border-bottom-color: #e94560;
    font-size: 42px;
  }
  section.title h2 {
    color: #a8d8ea;
    border: none;
    padding: 0;
    font-size: 24px;
    font-weight: normal;
    margin-top: 10px;
  }
  section.title p {
    color: #ccc;
    margin-top: 30px;
    font-size: 18px;
  }
  section.section-header {
    background: #0f3460;
    color: white;
    display: flex;
    flex-direction: column;
    justify-content: center;
    align-items: center;
    text-align: center;
  }
  section.section-header h1 {
    color: white;
    border-bottom-color: #e94560;
    font-size: 46px;
  }
  section.section-header p {
    color: #a8d8ea;
    font-size: 22px;
  }
  table {
    width: 100%;
    border-collapse: collapse;
    font-size: 18px;
  }
  th {
    background: #0f3460;
    color: white;
    padding: 8px 12px;
    text-align: left;
  }
  td {
    padding: 7px 12px;
    border-bottom: 1px solid #ddd;
  }
  tr:nth-child(even) td {
    background: #f0f4f8;
  }
  .highlight {
    background: #fff3cd;
    border-left: 4px solid #e94560;
    padding: 10px 16px;
    margin: 10px 0;
    border-radius: 4px;
  }
  code {
    background: #f0f4f8;
    padding: 2px 6px;
    border-radius: 3px;
    font-size: 18px;
    color: #c0392b;
  }
  pre {
    background: #1a1a2e;
    color: #e0e0e0;
    padding: 16px;
    border-radius: 6px;
    font-size: 16px;
  }
  ul li, ol li {
    margin-bottom: 6px;
  }
  .verdict-go {
    color: #27ae60;
    font-weight: bold;
  }
  .verdict-no {
    color: #e94560;
    font-weight: bold;
  }
  footer {
    font-size: 14px;
    color: #888;
  }
---

<!-- _class: title -->

# Multi-Source ALT Classification in Osteosarcoma

## A Computational Feasibility Study

**Rinat Rizvanov**
June 2026 · Independent computational analysis

---

## Agenda

1. Scientific motivation — OS, ALT, TCAB1
2. Research aims and go/no-go framework
3. Pipeline architecture and data sources
4. Labeling framework — four-tier hierarchy
5. Feasibility audit — Pass A and Pass B
6. RNA availability crisis and three recovery routes
7. Aim 2 discovery attempt — a null result
8. Sex confounding diagnosis (Notebook 12)
9. Hypothesis-driven curated gene score — D3
10. Robustness analysis (Notebook 13)
11. Limitations and future directions

---

<!-- _class: section-header -->

# Scientific Background

*Why osteosarcoma, ALT, and TCAB1?*

---

## Osteosarcoma: The Problem

- Most common **primary bone malignancy** in children and adolescents
- 5-year survival: ~65–70% localized / **< 20% metastatic or relapsed**
- Standard of care = MAP chemotherapy (methotrexate, doxorubicin, cisplatin) — **unchanged since 1980s**
- Defined by extreme genomic instability: near-universal TP53 disruption, complex SVs

**Defining molecular feature:**

> ~50% of osteosarcomas use **Alternative Lengthening of Telomeres (ALT)** instead of telomerase

This is the highest ALT prevalence of any tumor type studied systematically — making OS an ideal system to study ALT as a therapeutic vulnerability.

---

## What Is ALT?

ALT is a homologous-recombination-based mechanism for telomere maintenance that replaces telomerase in ~50% of osteosarcomas.

**Molecular hallmarks:**
- ALT-associated PML bodies (APBs)
- Extrachromosomal telomeric repeat (ECTR) DNA — detectable by **C-circle assay**
- Loss of ATRX/DAXX → aberrant histone H3.3 loading at telomeres
- Heterogeneous, ultra-long telomeres

**Why it matters therapeutically:**
ALT cells are dependent on HR machinery and FANCM-mediated break-induced replication at telomeres — potential synthetic lethal targets that telomerase-positive cells do not share.

---

## TCAB1 / WRAP53

- **TCAB1** (Telomerase Cajal Body protein 1), encoded by **WRAP53** at 17p13.1
- Localizes telomerase to Cajal bodies — required for efficient telomerase-mediated telomere extension
- In an **ALT-positive tumor**: TCAB1 function may be irrelevant or actively harmful
- Some evidence: TCAB1-disrupted tumors may be more dependent on ALT-mediated repair

**The two research questions:**

> 1. Is TCAB1 genomic disruption (17p13.1 SV breakpoints) associated with a distinct transcriptional program?
> 2. Can a transcriptional signature distinguish ALT+ from ALT− osteosarcoma cell lines?

---

<!-- _class: section-header -->

# Pipeline Architecture & Data

*What I built and where the data came from*

---

## Repository Architecture

```
tcab1-alt-osteo pipeline/
├── src/                      # All reusable logic
│   ├── manifest_builder.py   # cohort loaders
│   ├── label_rules.py        # deterministic label rules
│   ├── sv_annotator.py       # 17p13.1 locus SV calls
│   ├── audit.py              # Pass A / Pass B logic
│   └── aim2_signature.py     # LOOCV / signature training
├── notebooks/  00–13         # thin orchestration layers
├── tests/
│   └── test_sv_annotator.py  # 26 unit tests, all passing
└── data/
    ├── manifests/            # master_manifest.csv (115 × 46)
    └── processed/            # reports, predictions, sha256 hashes
```

**Design principles:**
- Notebooks import from `src/` — logic lives in one place
- Every output is **sha256-verified** for idempotent re-runs
- Analysis spec frozen in `label_spec_v1.0.md` before any modeling

---

## Five Cohorts, 115 Samples

| Cohort | N | Labels | RNA-seq |
|---|---|---|---|
| **DepMap 26Q1** | 35 cell lines | ALTitude predictions | ✅ 37 with RNA |
| PCGP 2014 | 34 pediatric tumors | SV + TP53 (peer-reviewed) | ✗ public release |
| de Nonneville 2022 | 22 bone tumors | C-circle / APB (peer-reviewed) | ✗ microarray only |
| Lawlor 2026 GCCRI | 17 PDX | SV / IHC (preprint) | ✗ internal only |
| Taglialatela 2026 | 7 cell lines | ALT list (preprint) | ✗ |

**Key tension:** Labels and expression are in **completely disjoint datasets**

The DepMap RNA matrix (1,719 cell lines × 19,177 genes) had no ALT labels. The labeled cohorts had no RNA-seq on disk.

---

<!-- _class: section-header -->

# Labeling Framework

*Four-tier hierarchy, frozen before any modeling*

---

## ALT Status Hierarchy (`label_spec_v1.0.md`)

| Tier | Label | Source | Used in |
|---|---|---|---|
| 1 | `alt_gold_positive` | C-circle / APB / TRAP, peer-reviewed | Primary + Pragmatic |
| 1 | `alt_negative_high_conf` | Same experimental evidence | Primary + Pragmatic |
| 2 | `alt_silver` | ALTitude prediction only | Pragmatic only |
| 3 | `alt_unknown` | No usable evidence | **Excluded** |

**Three parallel definitions** evaluated in audit:
- **Strict** (gold only) → primary analysis
- **Pragmatic** (gold + silver) → sensitivity arm
- Broad (gold + silver + proxy) → no additional samples, not used

**Critical column separation:**
`alt_call_direction` (positive/negative/unknown) is the arm-assignment column.
`alt_status` is provenance metadata only. All model code uses `alt_call_direction`.

---

## Why `alt_call_direction` Matters

The silver labels include both ALT+ and ALT− predictions:
- 6 silver → **positive** direction
- 5 silver → **negative** direction

If arm assignment used `alt_status` directly, all 11 silver samples would go to the positive arm regardless of their actual predicted direction — inflating the positive class and corrupting the analysis.

`alt_call_direction` was added in Notebook 10 after recognizing this during the manifest merge.

---

<!-- _class: section-header -->

# Feasibility Audit

*Two checks before any modeling code was written*

---

## Pass B: Confounding Census

Checks whether potential confounders are imbalanced across TCAB1 and ALT arms.

**Critical finding (Aim 1):**
> PCGP 2014 dominates the TCAB1-intact arm — 100% of intact samples in PCGP are from that cohort alone.

Any Aim 1 DE result would be **completely confounded** with cohort of origin.
Documented with `Halt override: True` and a written explanation — audit runs to completion for documentation.

**Pass A: RNA Availability Census**

Checks how many labeled samples have file-backed RNA-seq on disk.

| Aim | RNA-backed samples (initial) |
|---|---|
| Aim 1 (disrupted arm) | **0** |
| Aim 2 (ALT-labeled) | **0** |

This was not a pipeline failure — Stage 0 diagnostics were designed to catch this *before* any analysis code was written.

---

<!-- _class: section-header -->

# RNA Crisis and Recovery

*Three routes; one succeeded*

---

## The Gap: Labels and RNA Are Disjoint

```
                    [PCGP 2014]   [de Nonneville]   [Lawlor 2026]
                         ↓               ↓                ↓
                    SV labels      C-circle labels    SV labels
                    (no RNA)        (no RNA)           (no RNA)

                    [DepMap 26Q1]       [TARGET-OS]
                         ↓                  ↓
                    RNA-seq            RNA-seq
                    (no ALT labels)    (no labels)
```

**Three recovery routes:**

| Route | Description | Outcome |
|---|---|---|
| **A** | de Nonneville → ArrayExpress (E-MTAB-11743) | Microarray, NOT RNA-seq → external validation only |
| **B** | ALTitude predictions → DepMap cell-line-name join | **✅ Succeeded — rescued Aim 2** |
| **C** | Lawlor PDX internal RNA-seq | Pending PI collaboration — critical for Aim 1 |

---

## Route B: ALTitude → DepMap

Notebook 07 had tried exact sample-ID matching and failed (0 matches).
Notebook 08 tried **cell-line-name fuzzy join**:

- 38 DepMap OS cell lines matched to ALTitude table by name
- 37 had RNA-seq in the DepMap matrix
- Of these: **12 ALT+ / 7 ALT− gold labels** + 6+/5− silver

**The column-binding bug (Notebook 09):**
Initial run returned NOT ESTIMABLE — suspicious given 37 RNA-backed lines.
Investigation: code was reading `Curated.Marker` (ATRX mutation feature) instead of `ALT_pred_class`.
Caught by ground-truth sanity check: U2OS came back as `ambiguous_or_missing`. Fixed immediately.

**Final cohort after Notebook 10 direction-aware merge:**
- Primary (gold): **N = 20, 12 ALT+ / 8 ALT−**
- Pragmatic (gold + silver): **N = 31, 18 ALT+ / 13 ALT−**

---

<!-- _class: section-header -->

# Aim 2 Modeling

*Notebook 11 → 12 → 13*

---

## Notebook 11: Whole-Transcriptome Discovery — Null

**Method:** Nested LOOCV, elastic-net logistic regression, 19,177 genes

Every preprocessing step is **fold-local only** (no leakage):
VarianceGeneFilter (top K=500) → StandardScaler → ElasticNet (inner 5-fold C selection)

**Result:**

| Definition | N | AUROC [95% CI] | Perm p |
|---|---|---|---|
| Primary (gold) | 20 | 0.4688 [0.20, 0.75] | 0.60 |
| Pragmatic (gold+silver) | 31 | 0.6880 [0.45, 0.90] | 0.06 |

**Primary = worse than random.** Pragmatic = borderline, not significant.

**Top features exposed the problem:**

| Rank | Gene | Note |
|---|---|---|
| 1–4 | XAGE1B, MAGEB2, XAGE1A, MAGEA10 | X-linked CT antigens (sex-linked) |
| 17–18 | USP9Y, DDX3Y | Y-chromosome genes |

The model was learning **sex, not ALT biology.** Zero ALT-gene-list genes in the top 20.

---

## Notebook 12: Sex Confounding Diagnosis

Three designs to characterize and rule out the sex confound:

| Design | Method | Primary AUROC | Perm p |
|---|---|---|---|
| **D1** — sex covariate | Append Female/Male after gene filter | 0.5000 | 0.53 |
| **D2** — no Y/XCT genes | Remove 41 Y-chr + 4 X-CT antigen columns | 0.3958 | 0.66 |
| **D3** — curated ALT genes | 19-gene prespecified list, no variance filter | **0.9479** | **< 0.001** |

D1 and D2 confirm: the borderline pragmatic signal in nb11 **depended entirely on sex-linked features**.
After sex correction: no ALT signal survives the unguided approach.

D3 is a different kind of test — covered next.

---

## Notebook 12 D3: Hypothesis-Driven Curated Gene Score

**Not data-driven discovery — a prespecified biological hypothesis test**

21-gene ALT pathway list: TERT, TERC, ATRX, DAXX, SMARCAL1, TP53, WRN, BLM, MRE11A, RAD50, NBN, FANCD2, RAD51AP1, WRAP53, RAD52, BRCA1, BRCA2, FANCM, SLX4, GEN1, MUS81

19 of 21 found in DepMap matrix (TERC = ncRNA; MRE11A = symbol mismatch)

**Results:**

| Definition | N | AUROC [95% CI] | AUPRC | Sens@90%Spec | Perm p |
|---|---|---|---|---|---|
| **Primary** | 20 | **0.9479 [0.81, 1.00]** | 0.9755 | 91.7% | **< 0.001** |
| Pragmatic | 31 | 0.8291 [0.67, 0.96] | 0.9166 | 77.8% | 0.010 |

**Primary nonzero coefficients (full-data refit, C = 0.3162):**

| Gene | Coefficient | Interpretation |
|---|---|---|
| TERT | −0.689 | Telomerase ↓ in ALT cells (mutual exclusivity) |
| FANCM | +0.683 | BIR at ALT telomeres — upregulated in ALT |
| BRCA2 | +0.197 | HR factor — ALT-dependent |
| WRN | −0.153 | RecQ helicase — complex, borderline-stable |
| NBN | +0.105 | HR factor — ALT-dependent |

---

## Gene Interpretation: TERT and FANCM

**TERT (negative coefficient)**
- Encodes the catalytic subunit of telomerase
- Telomerase and ALT are **mutually exclusive** telomere maintenance mechanisms
- Lower TERT expression → ALT-positive
- Most stable predictor across every analysis (bootstrap freq = 1.00, LOOCV fold freq = 1.00)

**FANCM (positive coefficient)**
- Fanconi anemia helicase required for break-induced replication at ALT telomeres
- Bhowmick et al. (2016): FANCM depletion selectively inhibits ALT
- Proposed **synthetic lethal target** in ALT-positive tumors
- Higher FANCM expression → ALT-positive
- Second most stable predictor (bootstrap freq = 1.00)

**WRAP53 (zero coefficient)**
- Despite being central to the original question, WRAP53 mRNA expression adds no information beyond TERT + FANCM in this dataset at N = 20
- Does not mean it's irrelevant — just that expression-level signal is not recoverable here

---

## Notebook 13: Robustness Analysis

**Four stress tests on D3:**

**1. Extended permutation test (1,000 shuffles):**

| Cohort | N | AUROC | Perm p (1,000 shuffles) |
|---|---|---|---|
| Primary | 20 | 0.9479 | **< 0.001** (0/1,000) |
| Pragmatic | 31 | 0.8291 | 0.010 (10/1,000) |
| Pragmatic, known-sex only | 27 | 0.8765 | 0.003 (3/1,000) |

**2. Known-sex-only pragmatic (drop 4 unknown-sex samples):**
AUROC improves from 0.8291 → 0.8765; CI tightens; TP53 coefficient shrinks substantially (+0.101 vs −0.337), suggesting its contribution in the full pragmatic cohort was entangled with label uncertainty from unknown-sex samples.

---

## Bootstrap Coefficient Stability (200 Resamples)

Fraction of resamples with |β| > 10⁻¹⁰ per gene:

| Gene | Primary freq | Primary class | Pragmatic freq | Pragmatic class |
|---|---|---|---|---|
| **TERT** | **1.00** | stable | **1.00** | stable |
| **FANCM** | **1.00** | stable | **1.00** | stable |
| BRCA2 | 0.81 | stable | 0.655 | borderline |
| TP53 | 0.445 | unstable | **0.935** | stable |
| WRN | 0.57 | borderline | 0.705 | borderline |
| NBN | 0.625 | borderline | 0.785 | borderline |
| ATRX | 0.235 | unstable | 0.555 | borderline |
| 12 others | < 0.30 | unstable | — | — |

**TERT and FANCM: the only genes stable in both definitions.** This is the signal worth taking forward.

---

## LOOCV Fold Coefficient Frequency

Fraction of N fold-refits where each gene has nonzero coefficient:
*(bootstrap tests resample variation; this tests sensitivity to individual sample removal)*

| Gene | Primary (N=20) | Pragmatic (N=31) | Pragmatic KS (N=27) |
|---|---|---|---|
| **TERT** | **1.000** | **1.000** | **1.000** |
| **FANCM** | **1.000** | **1.000** | 0.926 |
| BRCA2 | 0.950 | 0.677 | 0.519 |
| NBN | 0.900 | 0.710 | 0.926 |
| WRN | 0.900 | 0.677 | 0.926 |
| TP53 | 0.500 | 0.968 | 0.889 |

TERT and FANCM appear in every held-out fold in the primary definition. The result is not driven by any single sample.

---

<!-- _class: section-header -->

# Limitations & Next Steps

*Honest assessment*

---

## Limitations

**Sample size**
N = 20 primary, N = 31 pragmatic. AUROC CI of [0.81, 1.00] is wide. The expected AUROC in new samples could be substantially lower.

**Label quality**
Gold labels from different groups using different protocols — no guarantee of perfect internal consistency. Silver labels (ALTitude predictions) are unvalidated computationally-derived.

**Single platform**
DepMap 26Q1 log2(TPM+1) only. Cell lines in long-term culture may not represent primary tumor expression profiles.

**This is confirmatory, not discovery**
D3 tests prior biological knowledge about TERT, FANCM, and HR genes. The high AUROC confirms the prior is recoverable in this dataset — it does not discover new biology. Framing it as a novel signature would overstate the finding.

**No probability calibration**
LOOCV scores are not calibrated probabilities. This is a signature score, not a diagnostic classifier.

---

## Future Directions

**1. External validation on de Nonneville microarray (independent, near-term)**
- 22 osteosarcoma tumors with peer-reviewed C-circle ALT calls
- Apply 19-gene curated score after within-platform normalization
- Cross-platform validation substantially strengthens the finding
- Can be pursued independently, no new data needed

**2. Aim 1 via Route C — Lawlor PDX RNA (PI collaboration)**
- 17 PDX models, 8 TCAB1-disrupted vs 32 intact
- This is the critical missing piece for Aim 1
- The manifest infrastructure is ready to integrate the data immediately

**3. CRISPR dependency integration (DepMap Chronos)**
- FANCM's consistent positive coefficient in the ALT score connects to FANCM synthetic lethality literature
- DepMap Chronos CRISPR data: does FANCM dependency score correlate with ALT status?
- Directly testable with the existing pipeline and label set

**4. TCAB1-expression-quantile contrast within DepMap**
- Within 35 DepMap OS lines, stratify by WRAP53 expression quantile
- Alternative Aim 1 reframe that requires no new data (but answers a different question)

---

## Summary

| What was done | Result |
|---|---|
| 5-cohort master manifest (115 samples, 46 columns) | Built and audited |
| Four-tier ALT label hierarchy, frozen specification | Documented in label_spec_v1.0.md |
| Stage 0 RNA availability audit | Both aims initially inestimable |
| Three recovery routes evaluated | Route B rescued Aim 2 |
| Column-binding bug caught and fixed | U2OS sanity check |
| Aim 2 whole-transcriptome LOOCV (nb11) | **Null** — sex confound |
| Sex confounding diagnosis (nb12 D1, D2) | Confirmed: no ALT signal after correction |
| Curated gene score, primary cohort N=20 (nb12 D3) | **AUROC = 0.948, p < 0.001** |
| 1,000-shuffle permutation test (nb13) | Primary 0/1,000; pragmatic 10/1,000 |
| Bootstrap + LOOCV fold stability (nb13) | TERT and FANCM stable in all analyses |
| Aim 1 (TCAB1 DE) | **Not estimable** — pending Route C |

---

<!-- _class: title -->

# Thank You

**Key takeaway:** The transcriptional signal of ALT status in osteosarcoma is real but concentrated in a small number of functionally motivated genes. It cannot be found by unguided gene selection at N = 20, but it is recoverable when you know where to look.

TERT and FANCM are the most robust predictors across all analyses.
WRAP53 mRNA expression does not contribute to the classification signal.

**Code:** `tcab1-alt-osteo pipeline/` · 14 notebooks · 26 unit tests · sha256-verified outputs

Rinat Rizvanov · github.com/rinatrizvanov
