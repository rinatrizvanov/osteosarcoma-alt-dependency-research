# Metadata Scan — Iteration 03 (BIR module coherence)

**Scope.** Which stratification axes can this cohort actually support, with honest
group sizes, over (a) the **16 ALT+ CRISPR-screened lines** and (b) the **34
RNA-backed osteosarcoma lines**. `src.stratify.MIN_CELL_SIZE = 4` is the binding
guard: an axis is usable only if **both** arms have n ≥ 4.

**Outcome untouched.** `src.dependency_stats` was never imported. No Chronos gene
effect was read. `CRISPRGeneEffect.csv` was touched only through
`src.cohort.crispr_model_ids()`, which streams the first field of each line — row
membership, not values. Nothing here is derived from gene-effect data.

**Anti-bias rule b0 observed.** No ModelID↔ALT-label pairing was printed to
stdout at any point. Where a per-line vector is informative it appears below as a
**sorted, de-identified value list**, which carries no ID and therefore no join
key. The per-line CN-fragmentation table (§4.3) has its ModelID column dropped
deliberately.

## What I ran

| Step | Script (scratchpad) | Inputs |
|---|---|---|
| Cohort sets | `s1_sets.py` | `src.cohort.build_cohort()` + `assert_all()`, `osteosarcoma_rna_models()` |
| Manifest column census over 35 `DepMap_26Q1` rows | `s2_manifest.py` | `references/master_manifest.csv` (sha256-verified) |
| ALT+16 / rna34 level counts | `s3_pos16.py` | manifest |
| `Model.csv` census (49 columns) | `s4_model.py` | `Model.csv` |
| Named binary splits, both cohorts | `s5_splits.py` | `Model.csv` + manifest |
| Marker expression distributions | `s6_markers.py` | `OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv` via `src.stratify.load_expression` |
| Marker gene-level CN | `s7_cn.py` | `OmicsCNGeneWGS.csv` |
| Ploidy / segment validation | `s8_segval.py`, `s9_focal.py` | `OmicsCNSegmentsWGS.csv` |

`uv run python`, Python 3.12, pandas 2.x. All five raw inputs sha256-verified
against `src.utils.EXPECTED_SHA256`; the manifest matched.

Cohort reproduced exactly, all `src.cohort` assertions passing:

| Step | N |
|---|---|
| `Model.csv` rows matching `steosarcoma` | 35 |
| present in `CRISPRGeneEffect.csv` (screened) | 26 |
| present in expression matrix | 26 |
| joined to frozen manifest (`cohort == DepMap_26Q1`) | 26 |
| `alt_status` in admitted tiers | 25 |
| **ALT+ (primary set)** | **16** |
| ALT− (negative control) | 9 |
| **RNA-backed OS lines (cut-point reference set)** | **34** |

All 16 ALT+ and all 9 ALT− lines are inside the 34 RNA-backed set.

---

## 1. `references/master_manifest.csv` — every column as a candidate stratifier

Frozen, read-only, authoritative. ALT labels were **not** re-derived. Value
counts are over the 35 `cohort == DepMap_26Q1` rows; the "min arm, ALT+16" column
is the largest *minimum* arm reachable by any grouping of that column's
informative levels inside the 16 ALT+ lines (sentinels `unknown` /
`not_applicable` excluded from both arms).

### 1.1 Columns with no stratifying content at all (constant over 35)

| Column | Single value |
|---|---|
| `cohort` | `DepMap_26Q1` |
| `sample_type` | `cell_line` |
| `has_sc` | `False` |
| `has_clinical` | `True` |
| `sv_architecture_17p13_1` | `unknown` |
| `sv_architecture_evidence` | `pending_sv_review` |
| **`alt_assay_evidence`** | **`unknown` (all 35)** |
| `purity_estimate` | `not_applicable` |
| `purity_quartile` | `not_applicable` |
| `include_telomerase_negative_secondary` | `False` |
| `notes` | empty (all 35) |

### 1.2 Non-constant columns, with ALT+16 partitionability

| Column | distinct / 35 | `unknown`-type / 35 | Value counts over 35 | min arm, ALT+16 | ≥4 both arms? |
|---|---|---|---|---|---|
| `alt_status` | 4 | 0 | gold+ 12, silver 11, neg_high_conf 8, unknown 4 | **6** (gold 10 / silver 6) | **YES** |
| `alt_call_direction` | 3 | 4 | positive 18, negative 13, unknown 4 | 0 (constant in ALT+16) | no |
| `alt_source` | 4 | 0 | ALTitude_2024 19, ALTitude_2024_predicted 11, pending_literature_merge 4, ALTitude+C-circle 1 | **6** | **YES** — identical partition to `alt_status` |
| `source_confidence` | 3 | 4 | peer_reviewed 20, computational 11, unknown 4 | **6** | **YES** — identical partition to `alt_status` |
| **`telomerase_status`** | 3 | 27 | indeterminate 26, positive 8, unknown 1 | **0 — constant (`indeterminate` ×16) inside ALT+** | **no** |
| `tert_expression` | 21 | 1 | continuous log2(TPM+1), 0.0 ×15, max 3.397 | 8 (zero 8 / nonzero 8) | YES, but see §1.4 |
| `telomerase_evidence` | 20 | 0 | free text restating `tert_expression` | 8 | same channel as above |
| **`atrx_status`** | 3 | 1 | not_altered 33, altered 1, unknown 1 | **1** | **no** |
| `daxx_status` | 2 | 1 | not_altered 34, unknown 1 | 0 (constant in ALT+16) | no |
| `smarcal1_status` | 2 | 1 | not_altered 34, unknown 1 | 0 (constant in ALT+16) | no |
| `top3a_status` | 2 | 1 | not_altered 34, unknown 1 | 0 (constant in ALT+16) | no |
| `tp53_status` | 3 | 1 | not_altered 26, altered 8, unknown 1 | 3 (altered 3 / not 13) | no |
| `tp53_evidence` | 4 | 0 | mutation-matrix text; no CN check performed | 3 | no |
| `has_rna` | 2 | 0 | True 34, False 1 | 0 (constant in ALT+16) | no |
| `has_wes` | 2 | 0 | True 34, False 1 | 0 (constant in ALT+16) | no |
| `has_wgs` | 2 | 0 | True 33, False 2 | 0 (constant in ALT+16) | no |
| `tcab1_rna_value` / `_quantile` / `_quantile_within_cohort` | 34 | 1 | continuous WRAP53 log2(TPM+1) 0.000–5.888 | — | iteration-01 exposure, see §1.5 |
| `tcab1_functional_status` | 3 | 1 | intact_supported 25, loss_possible 9, unknown 1 | 3 | no |
| `tcab1_genomic_status` | 2 | 32 | unknown 32, disrupted_high_conf 3 | 3 | no |
| `tcab1_overall_status` | 2 | 26 | unknown 26, candidate_disrupted 9 | 3 | no |
| `tcab1_confidence` | 2 | 0 | 0.0 ×26, 0.5 ×9 | 3 | no |
| `sv_architecture_auto` | 2 | 0 | Both-single-SV 33, Neither 2 | 0 (constant in ALT+16) | no — and artefactual (it-01 §3.3) |
| `include_aim1_primary_strict` / `_pragmatic` | 2 | 0 | False 33, True 2 | 2 | no |
| `include_aim1_primary_broad` | 2 | 0 | True 19, False 15 | 8 (8/8) | YES mechanically — stale Stage-0 flag, see §1.5 |
| `include_aim2_strict` | 2 | 0 | True 20, False 14 | 6 | YES — but = gold+high-conf tier, i.e. §1.2 row 1 again |
| `include_aim2_pragmatic` / `_broad` | 2 | 0 | True 31, False 4 | 0 (constant in ALT+16) | no |
| `include_tp53_sensitivity` | 2 | 0 | False 34, True 1 | 1 | no |

### 1.3 Is there any continuous ALT-activity measure? **No.**

This is the single most important answer in this scan, so it is stated flatly.

**There is no quantitative ALT-activity variable anywhere in this cohort.**

Specifically:

- `alt_assay_evidence` is the literal string `unknown` for **all 35** rows. The
  column that was designed to carry assay evidence is empty.
- **No C-circle column exists.** The only C-circle datum in the entire frozen
  manifest is a free-text fragment inside `alt_source` for **one** line:
  `ALTitude_2024; PMID:33344901; ALT_assay=C-circle; ALT_assay_strength=Negative
  (low C-circle signal)`. It is **qualitative** ("Negative", "low signal" — no
  number), and that line is **in neither the 16 ALT+ nor the 9 ALT−** sets (it is
  inside the 34 RNA-backed set but was not CRISPR-screened). It is therefore
  unusable for anything in iteration 03.
- **No APB column exists.** No PML-body count, no immunofluorescence field, no
  colocalisation score. Nothing.
- **No telomere-length column exists.** No TRF, no qPCR T/S ratio, no Flow-FISH,
  no telomere-FISH intensity, no computational estimate (`TelSeq`/`Telomerecat`
  style) — and no raw reads in this repo from which one could be computed.
- `purity_estimate` and `purity_quartile` are `not_applicable` ×35, so even the
  tumour-purity continuum that might have supported a weighting is absent.
- `references/manifest_builder.py` **does** contain code paths for TRAP assays,
  telomere FISH and ATRX IHC (lines ~180, ~227, ~294–376), but those branches
  serve the *other* upstream cohorts (Chen 2014 PCGP Table S6 etc.). None of them
  populated a single `DepMap_26Q1` row: every one of the 35 rows has
  `alt_source` beginning `ALTitude_2024`.

The only continuous, ALT-adjacent number in the manifest is `tert_expression`,
and it measures the **opposite** mechanism (telomerase), is a thresholded RNA
read-out rather than an assay, and is near-degenerate inside the ALT+ set (§1.4).

**Consequence for rule (c).** A per-subtype secondary cannot be operationalised
on an ALT-activity quantification, because no such quantification exists. If a
subtype secondary is reported at all it must rest on a surrogate marker gene, and
§3 shows what that costs.

### 1.4 `telomerase_status` and `tert_expression` inside the ALT+ set

`telomerase_status` is **constant** across the 16 ALT+ lines: `indeterminate`
×16. Zero ALT+ lines are `positive`. The axis does not exist inside the primary
set. (Over the 34 RNA-backed lines it is 8 positive / 26 indeterminate, which
does clear the floor — but 8 of the 9 telomerase-positive calls sit in the ALT−
arm, so over the 34 it is a near-proxy for the ALT label itself, exactly as
iteration 01 §5.5 found.)

`tert_expression` as a continuous variable:

| Set | n | exact zeros | nonzero | max |
|---|---|---|---|---|
| ALT+16 | 16 | 8 | 8 | **0.168** |
| ALT−9 | 9 | 1 | 8 | 3.397 |
| rna34 | 34 | 15 | 19 | 3.397 |

A zero-vs-nonzero split gives **8/8** inside ALT+16 and formally clears the
guard. It should not be used. Every nonzero ALT+ value is ≤ 0.168 log2(TPM+1),
i.e. TPM ≤ 0.12 — at or below the detection floor of an RSEM TPM matrix. The
8/8 split is a split on **detection noise**, not on telomerase activity: TERT is
uniformly off in all 16 ALT+ lines, which is the biologically expected result and
carries no within-ALT information.

### 1.5 Columns that are not exogenous to anything iteration 03 would test

- `tcab1_rna_value`, `tcab1_rna_quantile`, `tcab1_rna_quantile_within_cohort`,
  `tcab1_functional_status`, `tcab1_overall_status`, `tcab1_confidence`,
  `tcab1_functional_evidence` — all functions of WRAP53 expression, which was
  iteration 01's exposure. Iteration 01 returned `counter_signal`; reusing them
  re-tests a closed question.
- `source_confidence` and `alt_source` are **deterministic in `alt_status`**
  inside ALT+16: gold ⇔ `ALTitude_2024` ⇔ `peer_reviewed` (10 lines), silver ⇔
  `ALTitude_2024_predicted` ⇔ `computational` (6 lines). They are one axis with
  three names, not three axes.
- `include_aim1_primary_broad` gives a tidy 8/8 but is a **stale Stage-0
  inclusion flag** whose construction rule lives in `references/label_rules.py`
  and is not a biological property. `src.cohort` deliberately does not trust the
  `include_*` flags. Do not use it.
- `sv_architecture_auto` is constant inside ALT+16 and is artefactual anyway
  (33/35 calls trace to one CN segment spanning both WRAP53 and TP53).

---

## 2. `data/raw/depmap/Model.csv` — demographic, clinical and provenance axes

49 columns over the 35 osteosarcoma rows. Identifier columns (`ModelID`,
`PatientID`, `CellLineName`, `StrippedCellLineName`, `CCLEName`, `RRID`,
`CatalogNumber`, `SangerModelID`, `COSMICID`, `WTSIMasterCellID`) are **excluded
from the usable list**: a combinatorial search over near-singleton levels will
always report a balanced "split", and that split is meaningless. They are
reported below only in their present/absent form, which is a genuine (if weak)
provenance axis.

### 2.1 Columns unusable before any group-size question

Constant over the 35 OS rows: `OncotreeLineage` (all `Bone`),
**`OncotreePrimaryDisease` (all `Osteosarcoma`)**, `ModelType`,
`PediatricModelType`, `SerumFreeMedia`.

100% missing: `TissueOrigin`, `ModelTreatment`, `PatientTreatmentType`,
`PatientTreatmentDetails`, `Stage`, `StagingSystem`, `PatientTumorGrade`,
`PlateCoating`, `EngineeredModel`, `EngineeredModelDetails`,
`CulturedResistanceDrug`, `PublicComments`, `ModelSubtypeFeatures`,
`ModelIDAlias`.

≥ 74% missing: `PatientSubtypeFeatures` (34/35), `PatientTreatmentResponse`
(34/35), `HCMIID` (34/35), `PatientTreatmentStatus` (30/35),
`ModelDerivationMaterial` (28/35), `PatientRace` (26/35).

**Stage, grade and prior treatment are unavailable, not merely sparse.** No
clinical-severity or treatment-history axis can be built, described or adjusted
for in this cohort.

**`OncotreeSubtype` cannot stratify.** 34 `Osteosarcoma` vs 1 `Osteoblastic
Osteosarcoma`; inside ALT+16 it is 15/1. Same for `OncotreeCode` and
`DepmapModelType` (identical partition).

### 2.2 Named binary splits, both cohorts

| Candidate axis | ALT+16 (excl.) | ≥4 both? | rna34 (excl.) | ≥4 both? |
|---|---|---|---|---|
| `Sex`: Male vs Female | **8 / 7** (1) | **YES** | **17 / 13** (4) | **YES** |
| `AgeCategory`: Pediatric vs Adult | 11 / 3 (2) | no | **21 / 8** (5) | **YES** |
| `Age` ≥ 13 | **9 / 4** (3) | **YES** (at floor) | **20 / 6** (8) | **YES** |
| `Age` ≥ 14 | **6 / 7** (3) | **YES** | **14 / 12** (8) | **YES** |
| `Age` ≥ 15 | **6 / 7** (3) | **YES** | **13 / 13** (8) | **YES** |
| `Age` ≥ 18 | 3 / 10 (3) | no | **8 / 18** (8) | **YES** |
| `PrimaryOrMetastasis`: Primary vs Metastatic | 14 / 1 (1) | no | **25 / 5** (4) | **YES** |
| `SampleCollectionSite`: bone+rib vs soft_tissue | 12 / 3 (1) | no | **26 / 4** (4) | **YES** (at floor) |
| `GrowthPattern`: Adherent vs Mixed | 13 / 3 (0) | no | **30 / 4** (0) | **YES** (at floor) |
| `SourceType`: commercial repository vs academic/hospital | **6 / 10** (0) | **YES** | **13 / 21** (0) | **YES** |
| `SourceType`: UCSF Sweet-Cordero vs all other | **4 / 12** (0) | **YES** (at floor) | **6 / 28** (0) | **YES** |
| `FormulationID`: DMEM-family vs other | **6 / 10** (0) | **YES** | **13 / 21** (0) | **YES** |
| `OnboardedMedia`: `MF-002-001` vs other | **5 / 11** (0) | **YES** | **9 / 25** (0) | **YES** |
| ModelID block `ACH-000/001` vs `ACH-002/003` (onboarding-era proxy) | **5 / 11** (0) | **YES** | **15 / 19** (0) | **YES** |
| `RRID` present vs absent | **8 / 8** (0) | **YES** | **20 / 14** (0) | **YES** |
| `SangerModelID` present vs absent | **6 / 10** (0) | **YES** | **12 / 22** (0) | **YES** |
| `SourceDetail` present vs absent | **9 / 7** (0) | **YES** | **20 / 14** (0) | **YES** |
| `CCLEName` present vs absent | 14 / 2 (0) | no | 31 / 3 (0) | no |
| `ModelAvailableInDbgap`: CCLE-approved vs HMB-MDS | 3 / 3 (10) | no | **7 / 8** (19) | **YES**, 56% unusable |
| `OncotreeSubtype`: OS vs OSOS | 15 / 1 (0) | no | 33 / 1 (0) | no |

Caveat on `SampleCollectionSite`: the only way to reach a 4-line arm inside
ALT+16 is to group the single `rib` line with `soft_tissue`, which is
anatomically wrong (a rib is bone). On the correct grouping it is 12/3 and
fails. It is listed as failing.

### 2.3 Depositor detail and medium

`SourceDetail` (the depositing investigator, the sharpest batch handle) inside
ALT+16: Alejandro Sweet-Cordero 4, Ola Myklebost 3, Richard Gorlick 1, Ana
Patiño-García Lab 1, **7 missing**. "Sweet-Cordero vs other named depositor"
is 4/5 — clears the floor at the floor, with 44% of the set unlabelled.

`OnboardedMedia` / `FormulationID` inside ALT+16: `MF-002-001` (DMEM+10%FBS) 5,
`MF-001-001` (RPMI+10%FBS) 3, `MF-023-001` (McCoy's 5A+10%FBS) 2, then 6
singletons. The two columns are 1:1 aliases.

**Medium and depositor are the same axis.** Iteration 01 §5.2 established that
`SampleCollectionSite`, `GrowthPattern`, `OnboardedMedia` and `SourceType` are
not independent in this cohort — the 3 ALT+ `Mixed`/`soft_tissue` lines are the
same Oslo NRH block. Adjusting for one partly adjusts for all of them, and at
n=16 they cannot be separated.

### 2.4 `PatientID` — depositor/patient collapsing

| Set | patients / lines | multi-line patients |
|---|---|---|
| ALT+16 | 16 / 16 | none |
| rna34 | 33 / 34 | `PT-epgMPc` ×2 |
| OS35 | 33 / 35 | `PT-epgMPc` ×3 |

**No patient-level pseudo-replication exists inside the 16 ALT+ lines.** All 16
are distinct patients, so no collapsing or leave-pair-out sensitivity run is
needed for an ALT+-only analysis. The HOS/143B/HTK triple (`PT-epgMPc`) is
entirely outside the ALT+ set; if a secondary analysis widens to the 34
RNA-backed lines, two of those three enter and must be collapsed or
leave-pair-out tested (iteration 01 §5.4).

---

## 3. Type I vs type II ALT — what, if anything, could stand in

The literature (Broderick 2023 EXD2 `10.1038/s41467-023-38029-z`; Zhang 2026 APB
kinetics `10.1002/bies.70187`; Folini 2026 C-circles) argues ALT runs in two
modes. §1.3 establishes that **no assay-based subtype variable exists**, so the
only candidates are expression or copy-number surrogates.

### 3.1 Availability

`TERC` is **absent** from
`OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv` — it is a non-coding RNA
and the matrix is protein-coding only. There is no telomerase-RNA read-out in
this release. All other requested genes are present in both the expression and
the gene-level CN matrices: EXD2, ATRX, DAXX, TERT, NBN, BLM, FANCM, SMARCAL1,
PML. I added NR2F2 and TERF2 (APB / shelterin), RAD51AP1 and SLX4 (ALT
recombination machinery), and POLD3, RAD52, WRAP53 for reference.

### 3.2 Distribution shape over the 34 RNA-backed OS lines

`max_gap` is the largest interior gap in the sorted 34 values **subject to ≥ 4
values on each side** (so the gap, if taken as a cut, would itself clear
`MIN_CELL_SIZE`). `gap/SD` is the bimodality diagnostic requested: a genuinely
bimodal distribution has its largest interior gap **comparable to or larger
than** the SD. `bimod_coef` is the Pearson bimodality coefficient
`(skew²+1)/kurtosis`; > 0.555 is the conventional flag.

| Gene | min | median | max | SD | max_gap | **gap/SD** | gap cut | gap splits 34 as | bimod_coef |
|---|---|---|---|---|---|---|---|---|---|
| EXD2 | 0.096 | 3.855 | 5.164 | 0.936 | 0.359 | **0.384** | 2.843 | 4 / 30 | 0.508 |
| ATRX | 1.719 | 4.669 | 6.882 | 0.947 | 0.429 | **0.453** | 3.711 | 5 / 29 | 0.304 |
| DAXX | 4.604 | 5.396 | 6.496 | 0.508 | 0.282 | **0.556** | 5.687 | 21 / 13 | 0.551 |
| TERT | 0.000 | 0.013 | 3.397 | 0.846 | 0.757 | **0.895** | 1.302 | 30 / 4 | 0.833 |
| NBN | 3.106 | 4.930 | 6.478 | 0.626 | 0.197 | **0.314** | 5.525 | 29 / 5 | 0.229 |
| BLM | 0.462 | 3.575 | 4.793 | 1.268 | 0.815 | **0.642** | 1.686 | 6 / 28 | 0.731 |
| FANCM | 0.989 | 2.523 | 3.848 | 0.847 | 0.233 | **0.275** | 3.394 | 27 / 7 | 0.506 |
| SMARCAL1 | 0.014 | 3.809 | 4.825 | 0.872 | 0.248 | **0.284** | 3.182 | 7 / 27 | 0.576 |
| PML | 2.618 | 4.477 | 5.762 | 0.748 | 0.220 | **0.294** | 4.045 | 11 / 23 | 0.446 |
| NR2F2 | 2.745 | 4.570 | 8.012 | 1.244 | 0.467 | **0.376** | 6.254 | 29 / 5 | 0.424 |
| TERF2 | 3.327 | 4.450 | 5.407 | 0.490 | 0.143 | **0.292** | 4.086 | 11 / 23 | 0.421 |
| RAD51AP1 | 1.321 | 5.441 | 7.243 | 1.493 | 0.621 | **0.416** | 3.493 | 6 / 28 | 0.591 |
| SLX4 | 0.771 | 1.418 | 3.139 | 0.627 | 0.238 | **0.379** | 2.241 | 28 / 6 | 0.612 |
| POLD3 | 1.521 | 4.070 | 6.452 | 1.023 | 0.229 | **0.224** | 3.142 | 4 / 30 | 0.300 |
| RAD52 | 1.084 | 2.679 | 4.446 | 0.846 | 0.228 | **0.270** | 1.661 | 4 / 30 | 0.373 |
| WRAP53 | 0.000 | 3.603 | 5.889 | 1.295 | 0.360 | **0.278** | 2.297 | 4 / 30 | 0.455 |

**Not one gene is bimodal.** Every `gap/SD` is **below 1.0**, and 13 of 16 are
below 0.5. The distributions are unimodal continua. The sorted value lists make
this visible — e.g. DAXX over the 34 spans only 4.60–6.50 log2(TPM+1), a total
range of 1.9, with the largest step anywhere in it being 0.28:

```
EXD2      0.10 1.94 2.18 2.66 3.02 3.11 3.11 3.12 3.34 3.53 3.56 3.61 3.63 3.69 3.72 3.75 3.85 3.86 3.89 3.90 4.00 4.01 4.01 4.02 4.13 4.19 4.20 4.21 4.34 4.45 4.62 4.66 4.92 5.16
ATRX      1.72 2.95 3.22 3.34 3.50 3.93 4.05 4.09 4.21 4.23 4.25 4.50 4.53 4.56 4.62 4.63 4.67 4.67 4.77 4.79 4.93 4.94 5.18 5.18 5.21 5.28 5.32 5.39 5.42 5.46 5.55 5.59 5.61 6.88
DAXX      4.60 4.85 4.89 4.98 5.02 5.15 5.16 5.17 5.17 5.21 5.26 5.27 5.30 5.32 5.33 5.33 5.36 5.44 5.45 5.46 5.55 5.83 5.85 5.97 6.01 6.01 6.06 6.13 6.17 6.27 6.30 6.32 6.33 6.50
TERT      0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.01 0.01 0.02 0.02 0.02 0.03 0.05 0.09 0.11 0.14 0.17 0.50 0.53 0.77 0.92 1.68 1.70 3.08 3.40
NBN       3.11 4.13 4.14 4.19 4.27 4.42 4.42 4.47 4.53 4.55 4.56 4.67 4.69 4.74 4.78 4.80 4.93 4.93 4.94 4.98 4.99 5.02 5.13 5.16 5.20 5.38 5.39 5.41 5.43 5.62 5.72 5.77 5.80 6.48
BLM       0.46 0.57 0.79 0.81 0.96 1.28 2.09 2.32 3.09 3.11 3.21 3.28 3.30 3.31 3.35 3.39 3.47 3.68 3.77 3.84 3.94 3.95 3.99 4.00 4.04 4.10 4.15 4.16 4.16 4.19 4.26 4.50 4.73 4.79
FANCM     0.99 1.07 1.18 1.38 1.42 1.48 1.64 1.81 1.97 2.00 2.02 2.03 2.08 2.16 2.25 2.32 2.49 2.55 2.59 2.68 2.71 2.86 2.88 3.00 3.01 3.18 3.28 3.51 3.52 3.55 3.66 3.74 3.84 3.85
SMARCAL1  0.01 1.85 2.69 2.74 2.89 3.01 3.06 3.31 3.31 3.41 3.44 3.45 3.53 3.54 3.65 3.73 3.80 3.82 3.83 3.86 3.88 3.91 3.98 3.98 3.98 4.09 4.14 4.17 4.21 4.28 4.35 4.48 4.51 4.83
PML       2.62 2.71 3.18 3.36 3.39 3.49 3.59 3.64 3.81 3.93 3.94 4.15 4.25 4.30 4.32 4.43 4.47 4.48 4.61 4.63 4.64 4.68 4.69 4.73 4.75 4.84 4.87 4.90 4.92 5.12 5.17 5.24 5.37 5.76
NR2F2     2.75 2.87 3.14 3.33 3.34 3.56 3.74 4.04 4.17 4.21 4.26 4.29 4.37 4.41 4.43 4.54 4.57 4.57 4.64 4.71 4.79 5.06 5.13 5.41 5.53 5.79 5.90 5.92 6.02 6.49 6.61 6.63 7.04 8.01
TERF2     3.33 3.50 3.65 3.74 3.77 3.80 3.84 3.91 3.98 3.99 4.01 4.16 4.21 4.27 4.34 4.36 4.45 4.45 4.49 4.53 4.53 4.54 4.60 4.67 4.67 4.73 4.79 4.85 4.86 4.87 4.87 4.93 5.05 5.41
RAD51AP1  1.32 1.89 2.28 2.71 3.09 3.18 3.80 4.03 4.10 4.27 4.36 4.64 4.85 4.89 4.96 5.17 5.19 5.70 5.74 5.83 5.87 5.96 6.06 6.07 6.07 6.13 6.14 6.25 6.25 6.26 6.35 6.50 6.81 7.24
SLX4      0.77 0.83 0.95 0.97 0.97 0.98 1.00 1.12 1.15 1.17 1.18 1.22 1.26 1.27 1.38 1.39 1.41 1.42 1.47 1.52 1.63 1.65 1.73 1.78 1.91 1.95 1.97 2.12 2.36 2.54 2.70 2.72 2.79 3.14
POLD3     1.52 2.29 2.88 3.03 3.26 3.33 3.34 3.40 3.48 3.49 3.49 3.54 3.56 3.67 3.78 3.96 3.99 4.15 4.17 4.18 4.29 4.31 4.48 4.70 4.84 4.85 4.91 5.08 5.12 5.26 5.41 5.62 5.62 6.45
RAD52     1.08 1.24 1.34 1.55 1.77 1.80 1.95 1.99 2.04 2.09 2.28 2.30 2.45 2.57 2.61 2.68 2.68 2.68 2.70 2.75 2.76 2.77 2.85 2.91 3.06 3.10 3.22 3.41 3.55 3.63 3.68 4.17 4.40 4.45
WRAP53    0.00 0.00 2.04 2.12 2.48 2.50 2.75 2.94 3.07 3.11 3.25 3.34 3.43 3.46 3.46 3.48 3.49 3.72 3.85 4.10 4.15 4.28 4.29 4.30 4.42 4.61 4.62 4.63 4.79 4.80 4.88 5.02 5.27 5.89
```

The six genes whose `bimod_coef` exceeds 0.555 (TERT, BLM, SLX4, RAD51AP1,
SMARCAL1, DAXX) are **not** bimodal in any useful sense; the coefficient is being
driven by one or two low outliers in a long left tail (SMARCAL1's 0.014 and
EXD2's 0.096 are single points). Separating one outlier is not finding a subtype.

The one genuinely large gap in the table — **TERT**, gap/SD 0.895 — separates the
30 telomerase-null lines from the 4 telomerase-high lines. That is the **TMM
axis** (ALT vs telomerase), not a subtype axis within ALT, and inside ALT+16 it
puts **all 16 lines on one side (16/0)**, as §1.4 predicted.

### 3.3 Cut values and the ALT+16 arm sizes each cut would produce

All cuts estimated on the **34 RNA-backed lines** (rule b: widest non-outcome
superset), then applied to the 16 ALT+ lines. Arms shown as **low / high**.

| Gene | t33 cut | median cut | t67 cut | gap cut | ALT+16 @ t33 | ALT+16 @ median | ALT+16 @ t67 | ALT+16 @ gap |
|---|---|---|---|---|---|---|---|---|
| **EXD2** | 3.606 | 3.855 | 4.013 | 2.843 | **4 / 12 YES** | **6 / 10 YES** | **8 / 8 YES** | 2 / 14 no |
| ATRX | 4.497 | 4.669 | 5.176 | 3.711 | **7 / 9 YES** | **8 / 8 YES** | **10 / 6 YES** | 3 / 13 no |
| DAXX | 5.270 | 5.396 | 5.855 | 5.687 | **5 / 11 YES** | **7 / 9 YES** | **11 / 5 YES** | **11 / 5 YES** |
| TERT | 0.000 | 0.013 | 0.093 | 1.302 | 0 / 16 no | **10 / 6 YES** | 13 / 3 no | 16 / 0 no |
| NBN | 4.675 | 4.930 | 5.128 | 5.525 | **4 / 12 YES** | **6 / 10 YES** | **9 / 7 YES** | 14 / 2 no |
| BLM | 3.277 | 3.575 | 3.990 | 1.686 | 2 / 14 no | **6 / 10 YES** | **11 / 5 YES** | 1 / 15 no |
| FANCM | 2.027 | 2.523 | 2.881 | 3.394 | 3 / 13 no | **4 / 12 YES** | **7 / 9 YES** | **10 / 6 YES** |
| SMARCAL1 | 3.450 | 3.809 | 3.978 | 3.182 | **6 / 10 YES** | **10 / 6 YES** | **11 / 5 YES** | **4 / 12 YES** |
| **PML** | 4.155 | 4.477 | 4.690 | 4.045 | **5 / 11 YES** | **10 / 6 YES** | **12 / 4 YES** | **5 / 11 YES** |
| NR2F2 | 4.290 | 4.570 | 5.126 | 6.254 | **5 / 11 YES** | **9 / 7 YES** | **11 / 5 YES** | 13 / 3 no |
| TERF2 | 4.158 | 4.450 | 4.602 | 4.086 | 3 / 13 no | **7 / 9 YES** | **9 / 7 YES** | 3 / 13 no |
| RAD51AP1 | 4.637 | 5.441 | 6.056 | 3.493 | **4 / 12 YES** | **6 / 10 YES** | **7 / 9 YES** | 1 / 15 no |
| SLX4 | 1.221 | 1.418 | 1.729 | 2.241 | **5 / 11 YES** | **7 / 9 YES** | **10 / 6 YES** | 13 / 3 no |
| POLD3 | 3.541 | 4.070 | 4.475 | 3.142 | **4 / 12 YES** | **8 / 8 YES** | **9 / 7 YES** | 0 / 16 no |
| RAD52 | 2.298 | 2.679 | 2.850 | 1.661 | **6 / 10 YES** | **9 / 7 YES** | **11 / 5 YES** | 3 / 13 no |
| WRAP53 | 3.336 | 3.603 | 4.286 | 2.297 | 3 / 13 no | **8 / 8 YES** | **9 / 7 YES** | 2 / 14 no |

### 3.4 What this means — read before picking anything

**I am not picking a marker.** Per the brief, here is what the data can and
cannot support.

**Mechanically available.** A **median split** on any of the 16 genes clears
`MIN_CELL_SIZE` inside ALT+16. So does a t33 split on 11 of them. There is no
shortage of splits that pass the guard.

**But the guard is the only thing they pass.** Three problems, in order of
severity:

1. **No threshold is data-determined.** Because nothing is bimodal (§3.2, all
   gap/SD < 1), the cut point is a researcher choice, not a feature of the
   distribution. t33, median and t67 on the same gene give materially different
   ALT+ partitions — EXD2 goes 4/12 → 6/10 → 8/8, PML goes 5/11 → 10/6 → 12/4.
   Choosing among them after seeing this table is exactly the degree of freedom
   rule (b) exists to close.
2. **No way to validate the assignment.** With no C-circle, APB or
   telomere-length measurement in the cohort (§1.3), a marker-based subtype call
   **cannot be checked against the phenotype it claims to proxy**. The label
   would be unfalsifiable within this dataset.
3. **The specific marker rationales do not survive contact with this cohort.**
   - **EXD2** (Broderick 2023) is the brief's strongest literature candidate and
     the one I would expect to be proposed. Its distribution over the 34 is a
     tight unimodal cluster (3.0–5.2 for 32 of 34 lines) plus two low outliers;
     gap/SD = 0.384. A median split at 3.855 would divide a cluster whose
     interquartile span is a few tenths of a log2 unit. Broderick's type I/II
     distinction was drawn from functional ALT assays, not from an expression
     tertile, and nothing here reconstructs it.
   - **ATRX / DAXX** are the field's canonical type-I axis, and iteration 02
     already closed them: `atrx_status` is 1/16, the segment-validated focal ATRX
     deletion count is 1/16, and ATRX **expression** is actively misleading — the
     four lines with gene-level CN ≈ 0 are the four *highest* ATRX expressers.
     An ATRX-expression median split (8/8) would therefore be a split that
     contradicts the genomic data for the lines where genomic data exists.
     `daxx_status` is constant.
   - **PML** (Zhang 2026, APB kinetics) is the second-strongest candidate.
     Problem: APB formation is a **kinetic, imaging** phenotype. PML *mRNA* over
     these 34 lines spans 2.62–5.76 with gap/SD = 0.294 and no structure.
     Steady-state PML transcript is not an APB count, and there is no basis for
     claiming a threshold in it marks a subtype.
   - **TERT** is the TMM axis, degenerate inside ALT+ (§1.4, §3.2).
   - **SMARCAL1, POLD3, RAD52** are iteration 03's *endpoint* genes. Stratifying
     a POLD3/RAD52/SMARCAL1 dependency-coherence endpoint on POLD3/RAD52/SMARCAL1
     *expression* is a semi-circular design even though expression and gene
     effect are different measurements. Flagged so it is not done by accident.
   - **WRAP53** was iteration 01's exposure and that question returned
     `counter_signal`. Not a subtype marker.

**Honest conclusion for rule (c).** There is **no defensible type I / type II
ALT marker in DepMap 26Q1 for this cohort.** Rule (c)'s own escape clause — "if
no literature-grounded subtype marker exists for that hypothesis, say so
explicitly as a limitation in the memo" — is the applicable branch, and §1.3 is
the evidence for invoking it.

If the main agent nevertheless wants a reported secondary, the least indefensible
construction, in decreasing order of defensibility, is:

1. **Report no subtype split, and discharge rule (c) by stating §1.3 as the
   limitation.** The cohort has no ALT-activity measurement of any kind.
2. If a split must be shown: **a single, pre-specified median split on one
   literature-named gene, chosen and committed before this table is consulted**,
   threshold fixed on the 34 RNA-backed lines, declared **exploratory**, and
   described in the memo as "an arbitrary cut on a unimodal expression
   distribution", never as a subtype assignment. EXD2 (median 3.855 → 6/10) and
   PML (median 4.477 → 10/6) are the two with literature names attached.
3. **Not** a gap-based or tertile-based cut selected for its arm sizes. **Not**
   ATRX or DAXX (iteration 02). **Not** SMARCAL1/POLD3/RAD52 (endpoint genes).

---

## 4. Copy-number reliability

Binding rule from iteration 02: gene-level values in `OmicsCNGeneWGS.csv` must be
segment-validated against `OmicsCNSegmentsWGS.csv` before use — three of four
apparent homozygous ATRX deletions were refuted that way.

### 4.1 Which CN-based candidates even exist

Gene-level CN over the 16 ALT+ lines, all 16 genes covered with no missing values:

| Gene | min | median | max | n < 0.1 | n < 0.3 | n < 0.6 | n > 1.5 |
|---|---|---|---|---|---|---|---|
| **ATRX** | **0.0002** | 0.607 | 1.073 | **4** | **4** | 8 | 0 |
| **WRAP53** | **0.0000** | 0.814 | 2.340 | **2** | **2** | 2 | 1 |
| PML | 0.463 | 1.013 | 1.519 | 0 | 0 | 1 | 1 |
| RAD52 | 0.399 | 1.216 | 1.897 | 0 | 0 | 2 | 3 |
| SLX4 | 0.485 | 0.847 | 1.282 | 0 | 0 | 1 | 0 |
| SMARCAL1 | 0.535 | 0.936 | 1.296 | 0 | 0 | 1 | 0 |
| TERF2 | 0.625 | 0.952 | 1.495 | 0 | 0 | 0 | 0 |
| NR2F2 | 0.611 | 1.285 | 3.447 | 0 | 0 | 0 | 5 |
| RAD51AP1 | 0.661 | 1.174 | 2.038 | 0 | 0 | 0 | 3 |
| FANCM | 0.748 | 1.334 | 1.780 | 0 | 0 | 0 | 4 |
| POLD3 | 0.814 | 1.142 | 2.006 | 0 | 0 | 0 | 3 |
| TERT | 0.819 | 1.033 | 1.442 | 0 | 0 | 0 | 0 |
| NBN | 0.824 | 1.173 | 2.670 | 0 | 0 | 0 | 4 |
| EXD2 | 0.839 | 1.145 | 1.825 | 0 | 0 | 0 | 3 |
| BLM | 0.853 | 1.355 | 2.312 | 0 | 0 | 0 | 7 |
| DAXX | 0.907 | 1.151 | 1.761 | 0 | 0 | 0 | 2 |

**Only two genes in the entire candidate set show any apparent deletion**: ATRX
(4 lines, already adjudicated) and WRAP53 (2 lines, Stage 1's exposure). The
ATRX count reproduces iteration 02 exactly — 4 lines below 0.1, 8 below 0.6 — so
that finding replicates and its verdict carries: **1 of 4 segment-validated
(CAL72), 3 of 4 refuted.**

Formally clearing the ≥4 floor inside ALT+16 are these CN splits:

| Split | ALT+16 arms | Segment check corroborates? |
|---|---|---|
| ATRX CN < 0.1 (or < 0.3) | 4 / 12 | **NO — 3 of 4 refuted (iteration 02).** Do not use. |
| ATRX CN < 0.6 | 8 / 8 | **Not applicable** — mid-distribution cut, no focal event to validate |
| WRAP53 CN < 0.9 | 9 / 7 | Not applicable — mid-distribution cut; and WRAP53 is iteration 01's closed exposure |
| SMARCAL1 / SLX4 CN < 0.9 | 7/9, 9/7 | Not applicable — mid-distribution cut |
| PML / RAD51AP1 / NR2F2 / TERT CN < 0.9 | 5/11, 5/11, 4/12, 4/12 | Not applicable — mid-distribution cut |
| TERF2 CN < 0.9 | 6 / 10 | Not applicable — mid-distribution cut |
| BLM / NBN / NR2F2 / FANCM CN > 1.5 | 7/9, 4/12, 5/11, 4/12 | Not validated — see §4.2 |

**This is the key point.** Every CN split that clears the floor other than the
refuted ATRX one is **a cut through the middle of a near-diploid continuum**, not
a call of a discrete event. The ATRX-style segment check has nothing to
corroborate: there is no focal deletion at the locus to confirm or refute,
because the gene-level value is ~0.85, not ~0. Such a split would be an arbitrary
threshold on noisy relative copy ratio, with no biological event behind it.

### 4.2 The gains are not ploidy artefacts, but they are not validated either

Per-line autosomal length-weighted median `RelCopyRatio` over the 16 ALT+ lines
(ploidy proxy), sorted, IDs withheld:

```
0.870 0.877 0.934 0.959 0.971 0.980 0.988 1.028 1.037 1.041 1.045 1.049 1.056 1.063 1.092 1.114
```

The proxy spans only 0.87–1.11, so `OmicsCNGeneWGS.csv` is already
ploidy-normalised and dividing gene CN by this proxy barely moves anything (BLM
> 1.5 goes 7 → 5, NBN 4 → 3, everything else within ±1 line). The apparent gains
are therefore **relative gains, not whole-genome ploidy artefacts**.

That is as far as I can take them. Confirming whether a given `CN > 1.5` call is
a *focal* amplification or a passenger inside a broad arm-level gain requires
locating the gene in `OmicsCNSegmentsWGS.csv`, which requires base-pair
coordinates. **This repository contains no gene annotation** — no GTF, no RefSeq
table, no coordinate column; the only coordinates anywhere in `src/` or
`references/` are `LOCUS_WINDOW = ("chr17", 7_560_000, 7_710_000)` for the
17p13.1 window, and iteration 02's hand-entered ATRX interval, which lives in a
memo rather than in code. So each ATRX-style check is a **manual, hard-coded,
one-gene-at-a-time exercise**, and I have not performed it for genes whose
coordinates I would have to supply from memory. The main agent should budget for
that if it wants any CN-based candidate, and should hard-code the interval in the
iteration's own notebook so it is reviewable.

### 4.3 Why gene-level CN is unreliable cohort-wide, not just at ATRX

Autosomal segment structure of the 16 ALT+ lines (one row per line, ModelID
column dropped per rule b0, sorted by `n_zero`):

| autosomal segments | `RelCopyRatio` < 0.10 segments | of which focal (< 2 Mb) | median focal length (kb) | fraction of autosomal genome at < 0.10 |
|---|---|---|---|---|
| 851 | 3 | 2 | 42.5 | 0.001 |
| 982 | 5 | 4 | 72.0 | 0.001 |
| 881 | 5 | 4 | 780.0 | 0.007 |
| 882 | 8 | 6 | 150.5 | 0.007 |
| 591 | 8 | 6 | 178.5 | 0.007 |
| 1471 | 10 | 8 | 114.0 | 0.007 |
| 981 | 11 | 9 | 100.0 | 0.002 |
| 287 | 13 | 11 | 129.0 | 0.007 |
| 527 | 14 | 12 | 84.8 | 0.008 |
| 546 | 14 | 12 | 98.5 | 0.007 |
| 274 | 15 | 14 | 51.5 | 0.006 |
| 862 | 16 | 15 | 136.0 | 0.002 |
| 515 | 85 | 84 | 7.0 | 0.005 |
| 816 | 89 | 88 | 8.0 | 0.004 |
| 646 | 89 | 89 | 12.0 | 0.002 |
| 586 | 107 | 105 | 8.0 | 0.010 |

**All 16 ALT+ lines carry focal near-zero segments** — 2 to 105 of them, median
11 per line, median length 7–780 kb, together covering ~0.6% of the autosomal
genome. Four lines have 84–105 such segments at a median length of 7–12 kb,
which is the signature of segmentation noise rather than 100 independent
homozygous deletions.

This generalises iteration 02's lesson: **a gene-level CN ≈ 0 in an osteosarcoma
line is not a rare, self-evidently credible event.** These are chromothriptic
genomes with hundreds of segments and dozens of sub-megabase near-zero calls
each. Any single-gene "homozygous deletion" has a material prior probability of
being one of those artefacts, which is exactly what happened to 3 of the 4 ATRX
calls. Gene-level CN from this release should be treated as **unusable for
single-gene event calls without a per-line, per-locus segment check**, and the
repository does not currently make that check cheap (§4.2).

---

## Usable axes

Both arms ≥ 4. Identifier columns excluded. Axes that are redundant restatements
of another listed axis are marked.

### Over the 16 ALT+ CRISPR-screened lines

| Axis | Arms | Note |
|---|---|---|
| `alt_status` tier: gold vs silver | **10 / 6** | Label-quality axis. `alt_source` and `source_confidence` and `include_aim2_strict` give the **identical** partition — one axis, four names. |
| `Sex`: Male vs Female | **8 / 7** | 1 `Unknown` excluded. Cleanest demographic axis. |
| `Age` ≥ 13 | **9 / 4** | 3 missing. At the floor. |
| `Age` ≥ 14 | **6 / 7** | 3 missing. Most balanced age cut. |
| `Age` ≥ 15 | **6 / 7** | 3 missing. Same partition as ≥ 14 in ALT+16. |
| `SourceType`: commercial repository vs academic/hospital | **6 / 10** | Batch axis. |
| `SourceType`: UCSF Sweet-Cordero vs other | **4 / 12** | At the floor. |
| `FormulationID`: DMEM-family vs other | **6 / 10** | Alias of depositor (it-01 §5.2). |
| `OnboardedMedia` = `MF-002-001` vs other | **5 / 11** | Alias of depositor. |
| `SourceDetail` present vs absent | **9 / 7** | Provenance-completeness, not biology. |
| ModelID block `ACH-000/001` vs `ACH-002/003` | **5 / 11** | Unvalidated onboarding-era proxy. |
| `RRID` present vs absent | **8 / 8** | Provenance-completeness. |
| `SangerModelID` present vs absent | **6 / 10** | Provenance-completeness. |
| Expression median split on any of the 16 marker genes | see §3.3 | Clears the guard; **threshold not data-determined** — read §3.4 first. |
| Expression t33 split on EXD2, ATRX, DAXX, NBN, PML, NR2F2, RAD51AP1, SLX4, POLD3, RAD52, SMARCAL1 | see §3.3 | As above. |
| Gene-level CN mid-distribution cuts (< 0.9, > 1.5) | see §4.1 | Clear the guard but are **not event calls**; §4.1–4.3. |

**Of these, exactly one is both exogenous and biologically meaningful at n=16:
`Sex` (8/7).** The `alt_status` gold/silver tier axis (10/6) is the natural
robustness axis — it is already the prespec's R1 sensitivity cohort — but it is a
label-confidence axis, not a biological one. Age (6/7 at ≥14) is the only other
exogenous biological axis, at the cost of 3 missing lines.

### Over the 34 RNA-backed osteosarcoma lines

Everything usable at n=16 remains usable, plus:

| Axis | Arms | Note |
|---|---|---|
| `AgeCategory`: Pediatric vs Adult | **21 / 8** | 5 `Unknown`. Fails at n=16, passes here. |
| `Age` ≥ 18 | **8 / 18** | 8 missing. |
| `PrimaryOrMetastasis`: Primary vs Metastatic | **25 / 5** | 4 excluded. Fails at n=16. |
| `SampleCollectionSite`: bone+rib vs soft_tissue | **26 / 4** | At the floor; the 4 are one depositor block. |
| `GrowthPattern`: Adherent vs Mixed | **30 / 4** | At the floor; same 4 lines as above. |
| `tp53_status`: altered vs not | **8 / 26** | Fails at n=16 (3/13). |
| `telomerase_status`: positive vs indeterminate | **8 / 26** | **Constant inside ALT+.** Near-proxy for the ALT label over the 34 — do not treat as exogenous. |
| `tcab1_functional_status`: loss vs intact | **9 / 25** | Function of iteration 01's exposure. |
| `ModelAvailableInDbgap`: CCLE-approved vs HMB-MDS | **7 / 8** | 19 of 34 excluded. |
| `Sex`, `alt_status` tier, `SourceType`, media, era-proxy, ID-presence axes | 17/13, 12/11, 13/21, 13/21, 15/19, 20/14 | As at n=16. |

Caution for any analysis widening to the 34: the `PT-epgMPc` patient contributes
**2** of the 34 lines, so a collapse or leave-pair-out sensitivity run is needed
(§2.4). Inside the 16 ALT+ lines this problem does not exist.

---

## Not usable, and why

### Nothing here supports it — stated plainly

| Wanted | Verdict |
|---|---|
| **Any continuous ALT-activity measure** | **Does not exist.** `alt_assay_evidence` = `unknown` ×35. No C-circle, APB, telomere-length, TRAP or ALT-score column. `purity_estimate`/`purity_quartile` = `not_applicable` ×35. One qualitative C-circle string, for one line outside both analysis sets. §1.3. |
| **A defensible type I / type II ALT marker** | **Not available from this data.** No candidate gene is bimodal (all gap/SD < 1.0, §3.2); no threshold is data-determined; and with no ALT-activity assay in the cohort, any marker-based subtype call is unfalsifiable here. Rule (c)'s "no such marker exists" branch applies. §3.4. |
| **An ATRX-based axis** | **Closed by iteration 02 and reconfirmed here.** `atrx_status` 1/16; gene-level CN < 0.1 gives 4/16 but **3 of the 4 are refuted by the segment file**; segment-validated focal loss = 1/16 (CAL72). ATRX expression contradicts the CN calls (CN-zero lines are the highest expressers). |
| **A `DAXX`, `SMARCAL1` or `TOP3A` genomic-status axis** | **Constant.** `daxx_status`, `smarcal1_status`, `top3a_status` are `not_altered` for all 16 ALT+ lines (and 34 of 35 overall, the 35th being `unknown`). |
| **A `telomerase_status` axis inside ALT+** | **Constant** (`indeterminate` ×16). §1.4. |
| **A TERT-expression axis inside ALT+** | Formally 8/8 (zero vs nonzero), but every nonzero value is ≤ 0.168 log2(TPM+1), below the TPM detection floor. It is a split on detection noise. §1.4. |
| **`TERC` expression** | **Absent from the matrix** — non-coding RNA, and the file is protein-coding only. §3.1. |
| **Stage, grade, prior treatment, chemo response** | **100% or ≥86% missing.** Cannot be used, adjusted for, or even described. §2.1. |
| **Sequencing / profiling batch** | **Does not exist in this repository.** No `OmicsProfiles.csv` or sequencing-metadata table; platform, library prep, lane and profiling date are unknowable (iteration 01 §5.6, unchanged). |
| **An era or onboarding-date axis** | No date column anywhere in `Model.csv`. The ModelID-block proxy (5/11 in ALT+16) is unvalidated. |
| **Independent ALT-label replication** | All 35 labels trace to a single upstream source, ALTitude 2024. No second assay exists inside the cohort. |

### Clears the floor but should not be used

| Axis | Arms (ALT+16) | Why not |
|---|---|---|
| `include_aim1_primary_broad` | 8 / 8 | Stale Stage-0 inclusion flag; not a biological property; `src.cohort` deliberately distrusts `include_*`. §1.5. |
| `include_aim2_strict` | 10 / 6 | Identical to the `alt_status` gold/silver partition. Redundant. |
| `alt_source`, `source_confidence` | 10 / 6 | Identical to the `alt_status` gold/silver partition. Redundant. |
| `tert_expression` > 0 | 8 / 8 | Detection noise. §1.4. |
| `tcab1_*` columns | 3 / 13 (fail anyway) | Functions of iteration 01's exposure; that question is closed (`counter_signal`). |
| ATRX CN < 0.1 / < 0.3 | 4 / 12 | 3 of 4 calls refuted by the segment file. §4.1. |
| ATRX / WRAP53 / SMARCAL1 / SLX4 / PML / TERF2 / NR2F2 / TERT / RAD51AP1 CN < 0.6–0.9 | 4/12 … 9/7 | Mid-distribution cuts on a near-diploid continuum. No focal event exists to segment-validate; the threshold is arbitrary. §4.1. |
| BLM / NBN / NR2F2 / FANCM CN > 1.5 | 4/12 … 7/9 | Not ploidy artefacts (§4.2), but focality is **unverified** — the repo has no gene annotation, so the ATRX-style check was not performed. Do not use without hard-coding the interval and running it. §4.2. |
| SMARCAL1 / POLD3 / RAD52 expression splits | 4/12 … 10/6 | Semi-circular: these are iteration 03's own endpoint genes. §3.4. |
| `SampleCollectionSite` 4-line arm | 12/3 correctly, 11/4 only by grouping `rib` with `soft_tissue` | The rib line is bone. The passing version is an anatomical error. §2.2. |
| Any identifier column treated as categorical | spuriously "8/8" | A combinatorial split over 16 singleton levels always balances; it means nothing. §2. |

---

## Limitations of this scan

1. **Metadata and exposure-side data only.** No Chronos value was read, and
   `src.dependency_stats` was never imported. `CRISPRGeneEffect.csv` was used
   only as a row-membership set. Nothing here speaks to effect size or power for
   iteration 03's endpoint.
2. **No statistical testing was performed**, by instruction. Every number above
   is a count, a quantile or a distributional summary. Nothing here enters
   campaign BH accounting.
3. **Group sizes are necessary, not sufficient.** An axis passing
   `MIN_CELL_SIZE = 4` says only that a contrast is *estimable*, not that it is
   *powered*. At 4 v 12 the minimum detectable |Cliff's δ| is 0.708 (iteration
   02); at 8 v 8 it is smaller but still large. The main agent should pull
   `src.verdict.power_note(n, m)` for whichever split it chooses.
4. **The segment-level focality check was not run per gene.** §4.2 explains why:
   no gene annotation exists in the repository and each check needs a hand-entered
   interval. I ran the coordinate-free ploidy and fragmentation checks instead,
   which are sufficient to show that no *new* CN candidate deserves one, but a
   main agent that insists on a CN axis must do the per-locus check itself.
5. **`gap/SD` and the bimodality coefficient are screens, not tests.** I did not
   run a Hartigan dip test or fit a mixture model. The conclusion "nothing is
   bimodal" rests on every gene's largest guard-respecting interior gap being
   smaller than its SD, plus visual inspection of the sorted value lists, which
   is adequate for a feasibility call at n=34 but is not a formal test of
   unimodality.
6. **Marker selection is the brief's, plus six additions.** I added NR2F2,
   TERF2, RAD51AP1, SLX4, POLD3 and RAD52 to the required ten. A gene not in that
   list of 16 was not examined, and absence from this report is not evidence
   against it.
7. **Iteration 01's §3 and §5 were not re-derived.** Where this scan's §1–2
   overlaps iteration 01's covariate inventory, the numbers agree; the new
   content is the ALT+16 and rna34 group sizes, which iteration 01 did not
   compute (it worked on a 31-line RNA-only labelled cohort).
