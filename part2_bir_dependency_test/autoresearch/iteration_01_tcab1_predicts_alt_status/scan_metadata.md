# Metadata Scan — Iteration 01

**Role:** metadata scanner. **Date:** 2026-10-12.
**Scope:** stratification / covariate axes available for DepMap 26Q1 osteosarcoma
cell lines, and exact cohort sizes for an ALT+ vs ALT− contrast that requires RNA
but **not** CRISPR screening.

**Sources read (read-only):**

| Path | Use |
|---|---|
| `references/master_manifest.csv` | frozen labels; 115 rows total, 35 with `cohort == "DepMap_26Q1"` |
| `data/raw/depmap/Model.csv` | DepMap model metadata; 2,154 models, 49 columns |
| `references/label_spec_v1.0.md` | §2 and §6 (see §2 of this scan for what §2 actually contains) |
| `references/label_rules.py` | operational tier definitions (§2 of the spec defers to v0.9, which is not in the repo) |
| `src.cohort.rna_model_ids()` / `crispr_model_ids()` | **row-membership probes only**; no gene values read |

Prohibitions honoured: `CRISPRGeneEffect.csv` gene effects never read (only the
first field of each row, via the sanctioned `crispr_model_ids()` membership probe);
`src.dependency_stats` never imported; the expression matrix touched only for
`ModelID` membership (`rna_model_ids()`), never for gene values. No existing file
modified. This is the only file written.

---

## 1. Cohort arithmetic

### 1.1 Stepwise filter trace

Universe = the 35 `cohort == "DepMap_26Q1"` manifest rows. (Cross-check:
`src.cohort.osteosarcoma_models()` returns exactly 35 `Model.csv` rows and the
two sets are identical — zero symmetric difference. The manifest's DepMap stratum
*is* the DepMap OS panel.)

| # | Step | N |
|---|---|---|
| 0 | `master_manifest.csv` rows, all cohorts | 115 |
| 1 | `cohort == "DepMap_26Q1"` | **35** |
| 2 | matched 1:1 to `Model.csv` osteosarcoma rows | 35 |
| 3 | `alt_call_direction ∈ {positive, negative}` (drops 4 `unknown`) | **31** |
| 4 | … **and** present in the expression matrix (RNA requirement) | **31** |
| 5a | → ALT+ arm (RNA-only cohort) | **18** |
| 5b | → ALT− arm (RNA-only cohort) | **13** |
| 6 | step 4 **and** present in `CRISPRGeneEffect.csv` | **25** |
| 7a | → ALT+ arm (RNA+CRISPR cohort) | **16** |
| 7b | → ALT− arm (RNA+CRISPR cohort) | **9** |

### 1.2 Membership counts over the 35 rows

| Property | True | False |
|---|---|---|
| `alt_call_direction ∈ {positive, negative}` | 31 | 4 (all `alt_unknown`) |
| in expression matrix (1,719 ModelIDs) | 34 | 1 |
| in CRISPR matrix (1,208 ModelIDs) | 26 | 9 |

RNA × CRISPR within the 31 labelled rows:

| | in CRISPR | not in CRISPR | All |
|---|---|---|---|
| **in RNA** | 25 | 6 | **31** |
| **not in RNA** | 0 | 0 | 0 |
| All | 25 | 6 | 31 |

**Every labelled line has RNA.** The single line without RNA (`ACH-001094`, HTK)
is `alt_unknown` and is excluded by the label filter anyway, so the RNA requirement
costs nothing. The CRISPR requirement is the only binding constraint.

### 1.3 Resulting arm sizes

| Requirement | ALT+ | ALT− | Total | Unknown excluded |
|---|---|---|---|---|
| **(a) RNA only** | **18** | **13** | **31** | 4 |
| **(b) RNA + CRISPR** | 16 | 9 | 25 | 1 screened + 9 unscreened |

Row (b) reproduces the Stage 1 figures exactly (16 / 9 / 1 unknown among the 26
CRISPR-screened OS lines), which validates the trace. Dropping the CRISPR
requirement is a **+6 lines (+24%)** gain: ALT+ 16→18 (+2), ALT− 9→13 (+4). The
gain is concentrated in the smaller arm, improving balance from 16:9 (1.78:1) to
18:13 (1.38:1).

### 1.4 The 6 ModelIDs gained by dropping the CRISPR requirement

| ModelID | Cell line | Arm | `alt_status` | `source_confidence` | Source institution | Sex | Age | Primary/Met |
|---|---|---|---|---|---|---|---|---|
| ACH-000613 | HOS | negative | alt_negative_high_conf | peer_reviewed | ATCC | Female | 13 | Primary |
| ACH-001813 | C242 | negative | alt_silver | computational | Montefiore Einstein (Gorlick) | Unknown | — | Primary |
| ACH-001817 | C393 | negative | alt_silver | computational | Montefiore Einstein (Gorlick) | Unknown | — | Metastatic |
| ACH-002779 | NCC-ESOS1-C1 | negative | alt_silver | computational | NCC Research Institute (Kondo) | Unknown | — | — |
| ACH-002178 | NY | positive | alt_gold_positive | peer_reviewed | JCRB | Male | 15 | Primary |
| ACH-002689 | OS525 | positive | alt_gold_positive | peer_reviewed | UCSF (Sweet-Cordero) | Female | 20 | Metastatic |

Composition of the gain: **2 gold positives, 1 high-confidence negative, 3 silver
negatives.** Consequences to carry forward:

- The negative arm's silver fraction rises from 2/9 (22%) to 5/13 (38%).
- 3 of the 6 gained lines have `Sex = Unknown` and missing `Age`; covariate
  completeness therefore degrades (see §3.4).
- `ACH-000613` (HOS) shares `PatientID = PT-epgMPc` with `ACH-001001` (143B),
  already in the ALT− arm. Expanding the cohort puts **two lines from the same
  patient in the same arm** — a non-independence that did not exist at N=9. See §5.4.

### 1.5 The 4 excluded `unknown` rows

| ModelID | Cell line | `alt_status` | `alt_source` | in RNA | in CRISPR |
|---|---|---|---|---|---|
| ACH-001094 | HTK | alt_unknown | pending_literature_merge | No | No |
| ACH-001712 | Hs 860.T | alt_unknown | pending_literature_merge | Yes | No |
| ACH-001818 | C396 | alt_unknown | pending_literature_merge | Yes | No |
| ACH-002778 | NCC-OS1-X2-C1 | alt_unknown | pending_literature_merge | Yes | Yes |

All four are unlabelled for the same reason (`pending_literature_merge`), not for
assay failure. Three have RNA: **if a label were ever obtained for them, the
RNA-only ceiling would be 34.**

---

## 2. Label structure

### 2.1 `alt_status` × `alt_call_direction`, all 35 DepMap_26Q1 rows

| `alt_status` | negative | positive | unknown | All |
|---|---|---|---|---|
| alt_gold_positive | 0 | **12** | 0 | 12 |
| alt_negative_high_conf | **8** | 0 | 0 | 8 |
| alt_silver | **5** | **6** | 0 | 11 |
| alt_unknown | 0 | 0 | 4 | 4 |
| **All** | **13** | **18** | 4 | 35 |

The table is clean: gold and high-confidence tiers are direction-deterministic,
and `alt_silver` is the only tier that splits across both directions — which is
why `src.cohort` resolves arms from `alt_call_direction` and never from `alt_status`.

### 2.2 Tier composition of each arm

| Arm | Total | Gold / high-confidence | Silver | Silver % |
|---|---|---|---|---|
| ALT+ (RNA-only, N=18) | 18 | 12 `alt_gold_positive` | 6 `alt_silver` | 33% |
| ALT− (RNA-only, N=13) | 13 | 8 `alt_negative_high_conf` | 5 `alt_silver` | 38% |
| ALT+ (RNA+CRISPR, N=16) | 16 | 10 | 6 | 38% |
| ALT− (RNA+CRISPR, N=9) | 9 | 7 | 2 | 22% |

**Direct answers.** Of the 18 positives: **12 gold, 6 silver.** Of the 13
negatives: **8 high-confidence, 5 silver.** A gold-only / high-confidence-only
sensitivity cohort (`src.cohort.gold_only`) is therefore **12 vs 8 = 20** under the
RNA-only requirement, versus 10 vs 7 = 17 under RNA+CRISPR.

### 2.3 Tier ↔ provenance identity in this cohort

In the DepMap_26Q1 stratum the tier is a perfect function of the evidence source:

| `alt_status` | `source_confidence` | `alt_source` |
|---|---|---|
| alt_gold_positive (12) | peer_reviewed (12/12) | `ALTitude_2024` (12) |
| alt_negative_high_conf (8) | peer_reviewed (8/8) | `ALTitude_2024` (7) + `ALTitude_2024; PMID:33344901; ALT_assay=C-circle; ...Negative (low C-circle signal)` (1, = HOS) |
| alt_silver (11) | computational (11/11) | `ALTitude_2024_predicted` (11) |
| alt_unknown (4) | unknown (4/4) | `pending_literature_merge` (4) |

So "silver" ≡ "ALTitude prediction, no assay" and "gold/high-conf" ≡ "ALTitude
curated assay call" in this cohort. `source_confidence` and `alt_status` are
redundant covariates here, not independent ones.

### 2.4 Quoted tier definitions

**Important provenance caveat.** `references/label_spec_v1.0.md` §2 ("ALT Status")
is a **placeholder**: its entire body is the line

> `[Copy from v0.9 Section 2 — no changes]`

and `label_spec_v0.9` is **not present anywhere in the repository**. The only
normative text §2 adds is the Aim-2 implementation note:

> "**Aim 2 implementation in v1.0:** With N=48 positive / N=11 negative, train a
> regularized signature score (elastic net or LASSO logistic regression with strong
> regularization) rather than a binary classifier with confidence intervals. …
> Do NOT report calibrated probabilities at small N."

The tier semantics must therefore be taken from the frozen builder code, which is
the executable definition of record. One line each, from
`references/label_rules.py` (`ALTITUDE_ALT_MAP`, L253-256; `derive_alt_call_direction`,
L265-285) and `references/manifest_builder.py`:

| Tier | Operational definition (as coded) |
|---|---|
| `alt_gold_positive` | ALTitude `Curated.Marker == "ALT+ (high confidence)"` — i.e. an assay-backed positive ALT call; direction is deterministic (`positive`) from the tier alone. Other cohorts reach the same tier via a direct assay (telomere-FISH "Positive", C-circle, TRAP/ALT binary = 1). |
| `alt_silver` | ALTitude `Curated.Marker == "ALT+"`, or a **prediction-only** ALTitude row (`altitude_evidence_tier == "prediction_only"` with `ALT_pred_class ∈ {ALT+, ALT−}`); the tier does **not** encode direction, which must be resolved from `altitude_pred_class`. In DepMap_26Q1 all 11 silver rows are prediction-only (`alt_source == "ALTitude_2024_predicted"`, `source_confidence == "computational"`). |
| `alt_negative_high_conf` | ALTitude `Curated.Marker == "ALT−"`, or an explicit assay-negative result in a literature cohort (telomere-FISH "Negative", `ALT == "no"`, `alt_binary == 0`); direction is deterministic (`negative`) from the tier alone. |

`src.cohort` admits exactly these three tiers (`ADMITTED_STATUSES`) and splits on
`alt_call_direction`; `ALT_POSITIVE_STATUSES` and `ALT_NEGATIVE_STATUSES` both
contain `alt_silver` precisely because the tier is direction-agnostic.

Binding limitation from `label_spec_v1.0.md` §6 relevant to this iteration:

> "1. **DepMap 26Q1 lacks structural variant calls.** Architecture calls for DepMap
> cell lines are derived from copy-number segments and treated as `unknown` for
> primary architecture analysis. Three CN-derived `disrupted_high_conf` calls
> (Saos-2, NRH-OS1, NRH-OS4) are documented with caveats …"
>
> "2. **CCLE 2019 SvABA data has minimal overlap with modern DepMap OS cohort.**
> … only 1 of 35 DepMap OS lines overlaps, yielding zero locus-window SVs."

---

## 3. Covariate inventory

Criteria: computed over the **35 DepMap_26Q1 OS lines**; a column is retained if
it is (i) non-constant among non-missing values and (ii) < 30% missing. Empty
strings and `nan` count as missing; the literal strings `Unknown` / `unknown` are
**not** counted as missing (they are DepMap sentinel levels) but are flagged
separately because they are analytically useless.

### 3.1 `Model.csv` — retained columns (17 of 49)

| Column | Distinct | % missing | Values (count) | Confounder risk for a WRAP53-expression → ALT-status association |
|---|---|---|---|---|
| `ModelID` | 35 | 0 | identifier | — (key) |
| `PatientID` | 33 | 0 | 33 patients for 35 lines; `PT-epgMPc` ×3 (HOS, 143B, HTK) | **Non-independence**, not confounding. See §5.4. |
| `CellLineName` / `StrippedCellLineName` | 35 | 0 | identifiers | — (key) |
| `CCLEName` | 32 | 8.6 | present 32, absent 3 | **Era proxy** — legacy CCLE naming. Weak. |
| `DepmapModelType` | 2 | 0 | OS 34, OSOS 1 | Low (n=1 minority). |
| `OncotreeSubtype` | 2 | 0 | Osteosarcoma 34, Osteoblastic Osteosarcoma 1 | Low (n=1 minority). |
| `OncotreeCode` | 2 | 0 | OS 34, OSOS 1 | Duplicate of above. |
| `Age` | 14 | 22.9 | n=27; min 9, median 14, mean 16.3, max 38 | **Plausible** — ALT prevalence varies with age in OS; and WRAP53/telomere biology is age-linked. |
| `AgeCategory` | 3 | 0 | Pediatric 22, Adult 8, Unknown 5 | **Plausible** — same channel as `Age`, coarser but complete. |
| `Sex` | 3 | 0 | Male 17, Female 14, Unknown 4 | **Plausible** — standard demographic confounder; WRAP53 is autosomal (17p13.1) so a direct mechanism is weak, but report it. |
| `PrimaryOrMetastasis` | 3 | 8.6 | Primary 25, Metastatic 5, Unknown 2, missing 3 | **Plausible** — disease-stage axis. |
| `SampleCollectionSite` | 4 | 5.7 | bone 25, soft_tissue 4, Unknown 3, rib 1, missing 2 | **Plausible, and entangled with source** — all 4 `soft_tissue` lines are Oslo NRH lines. See §5.2. |
| `SourceType` | 14 | 0 | ATCC 7, UCSF 6, Montefiore Einstein 4, Oslo/Radium 4, JCRB 3, RIKEN 2, NCC Japan 2, ECACC 1, Broad 1, Dana-Farber 1, DSMZ 1, CCLF 1, MD Anderson 1, Clínica Univ. de Navarra 1 | **HIGH — the main batch axis.** See §5. |
| `GrowthPattern` | 2 | 0 | Adherent 31, Mixed 4 | **Plausible** — culture phenotype; all 4 `Mixed` are the Oslo NRH lines, so it is a source alias here. |
| `OnboardedMedia` | 14 | 2.9 | MF-002-001 9, MF-001-001 5, MF-015-008 4, MF-015-001 3, MF-023-001 2, MF-007-006 2, MF-001-003 2, then 7 singletons + 1 missing | **Plausible — culture/media batch**, and strongly source-aligned. |
| `FormulationID` | 14 | 2.9 | DMEM+10%FBS 9, RPMI+10%FBS 5, EMEM+20%FBS+pyruvate+NEAA 4, EMEM+10%FBS 3, McCoy's 5A+10%FBS 2, DMEM:F12+10%FBS+GlutaMAX 2, RPMI+10%FBS+2mM Gln 2, + 7 singletons | **Plausible** — human-readable alias of `OnboardedMedia`; 1:1 with it. |

### 3.2 `Model.csv` — dropped columns (32 of 49), with reason

Constant (no stratification possible): `OncotreeLineage` (all `Bone`),
`OncotreePrimaryDisease` (all `Osteosarcoma`), `ModelType` (all `Cell Line`),
`PediatricModelType` (all `True`), `SerumFreeMedia` (all `False`).

**100% missing — no covariate available at all:** `TissueOrigin`, `ModelTreatment`,
`PatientTreatmentType`, `PatientTreatmentDetails`, `Stage`, `StagingSystem`,
`PatientTumorGrade`, `PlateCoating`, `EngineeredModel`, `EngineeredModelDetails`,
`CulturedResistanceDrug`, `PublicComments`, `ModelSubtypeFeatures`, `ModelIDAlias`.

≥30% missing (or near-singleton): `PatientSubtypeFeatures` (97.1%),
`PatientRace` (74.3%; only caucasian / black_or_african_american),
`PatientTreatmentStatus` (85.7%; 5 of 35 populated: Active treatment / Pre-treatment),
`PatientTreatmentResponse` (97.1%; one `Partial Response`),
`RRID` (40.0%), `SourceDetail` (42.9%), `CatalogNumber` (65.7%),
`ModelDerivationMaterial` (80.0%; the populated ones are `PDX`), `HCMIID` (97.1%),
`ModelAvailableInDbgap` (42.9%), `WTSIMasterCellID` (71.4%),
`SangerModelID` (62.9%), `COSMICID` (71.4%).

**Prior-treatment and clinical-stage covariates are effectively unavailable.**
`PatientTreatmentStatus` is 85.7% missing, `PatientTreatmentType`/`Details`/
`Stage`/`StagingSystem`/`PatientTumorGrade` are 100% missing, and
`PatientTreatmentResponse` has one populated cell. Treatment exposure, chemotherapy
history, grade and stage therefore **cannot be adjusted for or even described** in
this cohort. Note that `SourceDetail` (42.9% missing) is dropped by the stated rule
even though it carries the depositing-investigator name, which is the sharpest
batch handle available — it is used descriptively in §5 but should not be a model term.

### 3.3 `master_manifest.csv` — retained columns (35 of 46)

| Column | Distinct | % missing | Values (count) | Note / confounder risk |
|---|---|---|---|---|
| `sample_id` | 35 | 0 | identifier | key (= `ModelID`) |
| `has_rna` | 2 | 0 | True 34, False 1 | matches the live RNA probe exactly |
| `has_wes` | 2 | 0 | True 34, False 1 | — |
| `has_wgs` | 2 | 0 | True 33, False 2 | — |
| `tp53_status` | 3 | 0 | not_altered 26, altered 8, unknown 1 | **HIGH — the canonical OS confounder.** 17p13.1 shares the locus with WRAP53. |
| `tp53_evidence` | 4 | 0 | `OmicsSomaticMutationsMatrixDamaging value=0; no CN check yet` 26, `value=2` 7, `value=1` 1, `ACH-001094 not in damaging matrix` 1 | Shows `tp53_status` is mutation-only; **no CN check performed**. |
| `tcab1_rna_quantile` | 34 | 0 | 0.0 ×2, else distinct; one `unknown` | **This is the exposure, not a covariate** — WRAP53 expression, pan-cohort quantile. |
| `tcab1_rna_value` | 34 | 0 | log2(TPM+1); 0.0 ×2 (`ACH-002834`/OS384 gold+, `ACH-003185`/NRH-OS4 silver+), one `unknown` (`ACH-001094`) | **Exposure.** Frozen WRAP53 values are already in the manifest. |
| `tcab1_rna_quantile_within_cohort` | 34 | 0 | 0.0–1.0 in 33 steps, one `unknown` | **Exposure** (within-DepMap rank). |
| `tcab1_genomic_status` | 2 | 0 | unknown 32, disrupted_high_conf 3 | Near-constant; 3 CN-derived calls flagged in spec §6.1. |
| `tcab1_genomic_evidence` | 2 | 0 | pending_sv_review 32, Manual review notebook 02 3 | — |
| `tcab1_functional_status` | 3 | 0 | functional_intact_supported 25, functional_loss_possible 9, unknown 1 | **Derived from the exposure** (bottom-quartile WRAP53 RNA). Not an independent covariate. |
| `tcab1_functional_evidence` | 34 | 0 | free text restating the RNA quantile | Derived from exposure. |
| `tcab1_overall_status` | 2 | 0 | unknown 26, candidate_disrupted 9 | Derived from exposure. |
| `tcab1_confidence` | 2 | 0 | 0.0 ×26, 0.5 ×9 | Derived from exposure. |
| `alt_status` | 4 | 0 | see §2.1 | **This is the outcome/label.** |
| `alt_source` | 4 | 0 | ALTitude_2024 19, ALTitude_2024_predicted 11, pending_literature_merge 4, ALTitude+PMID:33344901 1 | **Label-provenance axis — see §5.3.** |
| `source_confidence` | 3 | 0 | peer_reviewed 20, computational 11, unknown 4 | Redundant with `alt_status` (§2.3). |
| `telomerase_status` | 3 | 0 | indeterminate 26, positive 8, unknown 1 | **Derived from TERT RNA** — a correlate of the ALT label, not an exogenous covariate. See §5.5. |
| `tert_expression` | 21 | 0 | log2(TPM+1); 0.0 ×15, max 3.40 | Same channel; RNA-derived. |
| `telomerase_evidence` | 20 | 0 | free text restating TERT value | Same channel. |
| `atrx_status` | 3 | 0 | not_altered 33, altered 1, unknown 1 | **Mechanistically central but uninformative here** (1 altered line, `ACH-003171`/NRH-OS1). |
| `daxx_status` | 2 | 0 | not_altered 34, unknown 1 | Near-constant; unusable. |
| `smarcal1_status` | 2 | 0 | not_altered 34, unknown 1 | Near-constant; unusable. |
| `top3a_status` | 2 | 0 | not_altered 34, unknown 1 | Near-constant; unusable. |
| `include_aim1_primary_strict` | 2 | 0 | False 33, True 2 | Stale Stage-0 flag; `src.cohort` deliberately does not trust these. |
| `include_aim1_primary_pragmatic` | 2 | 0 | False 33, True 2 | idem |
| `include_aim1_primary_broad` | 2 | 0 | True 19, False 16 | idem |
| `include_aim2_strict` | 2 | 0 | True 20, False 15 | = gold + high-conf tiers |
| `include_aim2_pragmatic` | 2 | 0 | True 31, False 4 | = the 31 labelled lines |
| `include_aim2_broad` | 2 | 0 | True 31, False 4 | idem |
| `include_tp53_sensitivity` | 2 | 0 | False 34, True 1 | unusable |
| `sv_architecture_auto` | 2 | 0 | Both-single-SV 33, Neither 2 | **Artefactual** — derived from CN segments, not SVs (spec §6.1); the 33 "Both-single-SV" calls come from a single CN record spanning both WRAP53 and TP53. Do not use. |
| `sv_auto_reason` | 34 | 0 | `Records hitting both WRAP53_coding and TP53_coding: [<id>]` ×33, `No SVs in locus window` ×2 | Confirms the above. |
| `alt_call_direction` | 3 | 0 | positive 18, negative 13, unknown 4 | **This is the arm variable.** |

Dropped from the manifest (constant or all-missing): `cohort` (all DepMap_26Q1),
`sample_type` (all cell_line), `has_sc` (all False), `has_clinical` (all True),
`sv_architecture_17p13_1` (all `unknown`), `sv_architecture_evidence` (all
`pending_sv_review`), `alt_assay_evidence` (all `unknown`), `purity_estimate` and
`purity_quartile` (all `not_applicable`), `include_telomerase_negative_secondary`
(all False), `notes` (100% empty — note that spec §6.1 claims the three CN-derived
caveats are "documented … in the notes column"; in this frozen copy the column is
empty for all 35 rows).

### 3.4 Covariate completeness cost of the expansion

Count of missing-or-`Unknown` cells, 31-line RNA-only cohort vs 25-line RNA+CRISPR:

| Column | bad / 31 (RNA-only) | bad / 25 (RNA+CRISPR) |
|---|---|---|
| `Sex` | 4 (12.9%) | 1 (4.0%) |
| `Age` | 7 (22.6%) | 4 (16.0%) |
| `AgeCategory` | 4 (12.9%) | 2 (8.0%) |
| `PrimaryOrMetastasis` | 3 (9.7%) | 2 (8.0%) |
| `SampleCollectionSite` | 3 (9.7%) | 2 (8.0%) |
| `GrowthPattern` | 0 | 0 |
| `OnboardedMedia` | 0 | 0 |
| `SourceType` | 0 | 0 |
| `CCLEName` | 3 | 2 |
| `RRID` | 14 (45.2%) | 10 (40.0%) |

The +6 lines bring disproportionately incomplete demographics: 3 of 6 are
`Sex = Unknown` with missing `Age`. The absolute demographic coverage still rises
(Sex known 24→27, Age known 21→24), but any Sex- or Age-adjusted model loses 4
and 7 lines respectively from N=31.

---

## 4. Usable binary splits

Group sizes on the **31-line RNA-only labelled cohort**. `n missing/Unknown` is
excluded from the two informative levels. Flag = both informative arms ≥ 4.

| Candidate binary split | Group sizes | n missing/Unknown | Both arms ≥ 4 |
|---|---|---|---|
| `Sex`: Male vs Female | 14 / 13 | 4 | **YES** |
| `AgeCategory`: Pediatric vs Adult | 20 / 7 | 4 | **YES** |
| `Age` ≥ 14 (median) vs < 14 | 13 / 11 | 7 | **YES** |
| `Age` ≥ 18 vs < 18 | 7 / 17 | 7 | **YES** |
| `PrimaryOrMetastasis`: Primary vs Metastatic | 24 / 4 | 3 | **YES** (at the floor) |
| `SampleCollectionSite`: bone+rib vs soft_tissue | 24 / 4 | 3 | **YES** (at the floor) |
| `GrowthPattern`: Adherent vs Mixed | 27 / 4 | 0 | **YES** (at the floor) |
| `SourceType`: commercial repository (ATCC/JCRB/RIKEN/DSMZ/ECACC) vs academic/hospital | 13 / 18 | 0 | **YES** |
| Media: DMEM-family `FormulationID` vs other | 11 / 20 | 0 | **YES** |
| Media: `OnboardedMedia == MF-002-001` vs other | 8 / 23 | 0 | **YES** |
| ModelID block: `ACH-000/001xxx` vs `ACH-002/003xxx` (onboarding-era proxy) | 13 / 18 | 0 | **YES** |
| `RRID` present vs absent | 17 / 14 | 0 | **YES** |
| `SangerModelID` present vs absent | 12 / 19 | 0 | **YES** |
| `ModelAvailableInDbgap`: CCLE-approved vs HMB-MDS | 7 / 8 (+2 In Progress) | 14 | **YES** (but 45% unusable) |
| CRISPR screened vs not | 25 / 6 | 0 | **YES** |
| `tp53_status`: altered vs not_altered | 7 / 24 | 0 | **YES** |
| `source_confidence` / tier: peer_reviewed(gold+high-conf) vs computational(silver) | 20 / 11 | 0 | **YES** |
| `tcab1_functional_status`: loss_possible vs intact_supported | 7 / 24 | 0 | **YES** — but this is a function of the exposure |
| `telomerase_status`: positive vs indeterminate | 8 / 23 | 0 | **YES** — but RNA-derived and label-entangled (§5.5) |
| `CCLEName` present vs absent | 28 / 3 | 0 | no (3) |
| `OncotreeSubtype`: Osteosarcoma vs Osteoblastic OS | 30 / 1 | 0 | no (1) |
| `DepmapModelType`: OS vs OSOS | 30 / 1 | 0 | no (1) |
| `atrx_status`: altered vs not_altered | 1 / 30 | 0 | no (1) |
| `daxx_status` / `smarcal1_status` / `top3a_status` | 31 / 0 | 0 | no (constant) |
| `has_wgs`: True vs False | 30 / 1 | 0 | no (1) |
| `tcab1_genomic_status`: disrupted_high_conf vs unknown | 3 / 28 | 0 | no (3) |
| `sv_architecture_auto`: Both-single-SV vs Neither | 30 / 1 | 0 | no (1) |
| `include_tp53_sensitivity` | 1 / 30 | 0 | no (1) |

**Summary.** 19 of 28 candidate splits clear the n ≥ 4 floor on both arms, but
four of them (`PrimaryOrMetastasis`, `SampleCollectionSite`, `GrowthPattern`, and
`ModelAvailableInDbgap`) sit exactly at or barely above it, and three of them
(`tcab1_functional_status`, `telomerase_status`, and `source_confidence`/tier) are
not exogenous to the exposure or the label. ATRX — the single most mechanistically
relevant stratifier for ALT — is unusable at n=1.

---

## 5. Batch / provenance risk

The brief asks whether ALT+ and ALT− lines came from systematically different
sources, eras or institutions. Descriptive cross-tabulations only; **no test of
association was performed.**

### 5.1 Institution × arm (31-line RNA-only cohort)

| `SourceType` | ALT− | ALT+ | All |
|---|---|---|---|
| UCSF Helen Diller (Sweet-Cordero) | 1 | **5** | 6 |
| ATCC | 4 | 3 | 7 |
| Oslo Univ. Hosp. / Norwegian Radium (Myklebost) | 1 | **3** | 4 |
| Montefiore Einstein (Gorlick) | **3** | 0 | 3 |
| JCRB | 1 | 2 | 3 |
| RIKEN | 1 | 1 | 2 |
| NCC Research Institute Japan (Kondo) | 1 | 0 | 1 |
| The Broad Institute (Kadoch) | 1 | 0 | 1 |
| CCLF | 0 | 1 | 1 |
| Clínica Universitaria de Navarra | 0 | 1 | 1 |
| DSMZ | 0 | 1 | 1 |
| MD Anderson (Gorlick) | 0 | 1 | 1 |
| **All** | **13** | **18** | **31** |

**There is clear institutional clustering by arm.** The three multi-line academic
depositors each lean one way: UCSF/Sweet-Cordero is 5:1 positive (and all 5
positives are `alt_gold_positive`), Oslo/NRH is 3:1 positive (all 4 silver),
Montefiore/Gorlick is 0:3 negative. Collapsing to commercial repository vs
academic/hospital removes the signal (13 vs 18 overall; 6:7 and 7:11 by arm), so
the clustering is at the level of *individual depositor*, not of repository class.

Of the 6 lines gained by dropping the CRISPR requirement, **3 of the 4 gained
negatives come from two investigator blocks** (Montefiore ×2, NCC Japan ×1). The
expansion therefore *increases* the institutional imbalance of the negative arm,
not just its size.

### 5.2 Site, growth pattern and media are source aliases

| | ALT− | ALT+ |
|---|---|---|
| `SampleCollectionSite = soft_tissue` | 1 | 3 |
| `GrowthPattern = Mixed` | 1 | 3 |

Both cells are the **same four lines** — the Oslo NRH block (NRH-GCT1 negative;
NRH-OS1, NRH-OS2, NRH-OS4 positive). `SampleCollectionSite`, `GrowthPattern` and
`SourceType` are therefore not independent axes in this cohort; adjusting for one
partly adjusts for the others.

Media shows the same pattern. `OnboardedMedia` × arm:

| `OnboardedMedia` | ALT− | ALT+ | Dominant depositor |
|---|---|---|---|
| MF-002-001 (DMEM+10%FBS) | 2 | **6** | 5 of 8 are UCSF |
| MF-015-008 (EMEM+20%FBS+pyruvate+NEAA) | **3** | 0 | all 3 Montefiore |
| MF-001-001 (RPMI+10%FBS) | 2 | 3 | mixed |
| MF-023-001/002 (McCoy's 5A) | 0 | 3 | all ATCC (U2OS, Saos-2, G-292) |
| MF-015-001 (EMEM+10%FBS) | 2 | 1 | ATCC/JCRB |
| remaining 9 formulations | 4 | 5 | singletons |

Media assignment tracks the depositing institution, so "media effect" and
"institution effect" cannot be separated at this N.

### 5.3 Label-provenance imbalance

| `alt_source` | ALT− | ALT+ |
|---|---|---|
| `ALTitude_2024` (curated assay) | 7 | 12 |
| `ALTitude_2024_predicted` (prediction only) | 5 | 6 |
| `ALTitude_2024; PMID:33344901; C-circle Negative` | 1 | 0 |

Label provenance is *not* badly skewed between arms (silver 38% of negatives vs
33% of positives), which is reassuring. But both arms derive from a **single
upstream label source, ALTitude 2024** — 31 of 31 lines. There is no independent
label replication within this cohort, so any systematic error in ALTitude's calls
(or in its ALT+/ALT− predictor for the 11 silver lines) propagates undetected
into both arms. The one externally corroborated call (HOS, PMID:33344901 C-circle)
is explicitly a *low-signal* negative.

### 5.4 Patient-level non-independence

`PatientID` has 33 distinct values for 35 lines. One patient contributes three
models:

| ModelID | Line | `PatientID` | Arm | Source |
|---|---|---|---|---|
| ACH-000613 | HOS | PT-epgMPc | **negative** | ATCC |
| ACH-001001 | 143B | PT-epgMPc | **negative** | ATCC |
| ACH-001094 | HTK | PT-epgMPc | unknown (excluded, no RNA) | ECACC |

143B is a known HOS derivative. Under the RNA+CRISPR requirement only 143B was in
the cohort; under the RNA-only requirement **both HOS and 143B enter the same
(negative) arm**, so 2 of 13 negatives are the same patient. This is a
pseudo-replication issue introduced *by* the expansion and should be handled
either by collapsing the pair or by a leave-pair-out sensitivity run.

### 5.5 Variables that are not exogenous covariates

Three retained manifest columns are derived from the same RNA matrix as the
proposed exposure, or from the label itself, and would be colliders/mediators
rather than confounders:

- `telomerase_status` / `tert_expression` / `telomerase_evidence` — thresholded
  TERT expression. Cross-tab: `positive` 8/13 of ALT−, 0/18 of ALT+;
  `indeterminate` 5/13 of ALT−, 18/18 of ALT+. This near-separation is expected
  (ALT and telomerase are alternative TMMs) and makes the variable a proxy for
  the label, not a control for it.
- `tcab1_rna_*`, `tcab1_functional_status`, `tcab1_overall_status`,
  `tcab1_confidence` — all functions of WRAP53 expression, i.e. of the exposure.
- `source_confidence` — deterministic in `alt_status` (§2.3).

### 5.6 What I could **not** determine

Stated plainly:

- **No era/date information exists in this cohort.** `Model.csv` has **no**
  onboarding date, collection date, derivation year, or DepMap-release-of-first-
  appearance column (the only "date-like" column name is `OnboardedMedia`, which
  is a media formulation ID). I used the `ACH-0000xx…ACH-0035xx` ModelID block as
  an *unvalidated* ordinal proxy (ACH-000xxx 3−/3+; ACH-001xxx 5−/2+;
  ACH-002xxx 4−/8+; ACH-003xxx 1−/5+ — i.e. positives skew toward the
  higher-numbered, later-onboarded blocks). I cannot confirm that ACH numbering
  is chronological, so this is suggestive only.
- **No sequencing / profiling batch information.** `data/raw/depmap/` contains only
  `Model.csv`, `CRISPRGeneEffect.csv`, `OmicsCNGeneWGS.csv`,
  `OmicsCNSegmentsWGS.csv` and the expression matrix. There is no
  `OmicsProfiles.csv`, `OmicsDefaultModelProfiles.csv` or sequencing-metadata
  table, so library prep, sequencing platform, flowcell/lane, read length,
  alignment pipeline version and profiling date are **unknown and unknowable from
  this repository**. An RNA-expression contrast between two label groups cannot be
  checked for sequencing-batch confounding with the data on hand.
- **No culture-passage, mycoplasma, STR-authentication or doubling-time fields.**
- **Prior treatment is unrecoverable** (§3.2): 85.7% missing for
  `PatientTreatmentStatus`, 100% missing for treatment type/details, stage, staging
  system and grade.
- **No independent ALT assay.** All 31 labels trace to ALTitude 2024; I cannot
  cross-validate them against a second assay within this cohort.
- **I did not look at expression values**, so I cannot say whether the UCSF /
  Montefiore / Oslo blocks differ in global expression profile (library depth,
  gene-detection rate, 3′ bias). That check is the natural next step and is
  outside this scan's permissions.

---

## LIMITATIONS OF THIS SCAN

1. **Metadata only; no outcome data touched.** No CRISPR gene effects were read
   (CRISPR entered only as a row-membership set), no expression values were read
   (the expression matrix entered only as a `ModelID` set), and `src.dependency_stats`
   was never imported. Consequently nothing here speaks to effect size, power, or
   whether any covariate actually predicts anything.
2. **No statistical testing was performed, by instruction.** Every imbalance
   reported in §5 is a raw count. Cells of size 3–6 are routine in these tables;
   an apparently striking 5:1 or 0:3 split is entirely compatible with chance at
   this N. Treat §5 as hypothesis-generating about *confounding structure*, not as
   evidence of it.
3. **The tier definitions are reconstructed from code, not from the spec.**
   `label_spec_v1.0.md` §2 is a placeholder pointing at v0.9, which is absent from
   the repository. §2.4 quotes `label_rules.py` / `manifest_builder.py` instead.
   If the prose v0.9 definitions differ from the implemented ones, this scan
   follows the implementation.
4. **The `notes` column is empty for all 35 rows**, although spec §6.1 says the
   three CN-derived `disrupted_high_conf` caveats are documented there. Provenance
   caveats recorded in a different copy of the manifest are not visible here.
5. **Missingness convention is a judgement call.** `Unknown` / `unknown` sentinel
   strings were counted as present-but-uninformative rather than missing. Under the
   stricter convention (sentinels = missing), `Sex` becomes 11.4% missing (still
   retained), `PrimaryOrMetastasis` 14.3% (retained), `SampleCollectionSite` 14.3%
   (retained), and `AgeCategory` 14.3% (retained) — no column's retain/drop verdict
   flips, but several manifest `*_status` columns (`atrx`, `daxx`, `smarcal1`,
   `top3a`, `tcab1_genomic_status`) become effectively constant rather than
   two-level. The 30%-missing rule was applied to the **35-line** universe, not to
   the 31-line analysis cohort; §3.4 gives the cohort-specific completeness.
6. **The 30%/non-constant rule mechanically discards informative-but-sparse
   fields.** `SourceDetail` (42.9% missing) is the sharpest depositing-investigator
   handle in the file and `RRID` (40.0%) is the cleanest era/authentication proxy;
   both are dropped by the rule. §5 uses them descriptively.
7. **Arm-level N is small and the covariate space is large.** With ALT− = 13 and
   roughly 19 usable binary splits plus a continuous `Age`, there is ample room for
   a spurious adjustment to change any downstream estimate. The number of
   *available* covariates reported here is not the number that can be *used*
   simultaneously.
8. **Single-release snapshot.** Everything is DepMap 26Q1 as frozen in this repo.
   Label and metadata drift across releases was not assessed, and the manifest is
   a frozen copy whose own provenance I did not verify beyond reading
   `references/MANIFEST.sha256` paths.
9. **Alias collapsing not audited.** `CCLEName` values such as `PSS089_BONE`
   (= OS186) and `SYO1_SOFT_TISSUE` indicate legacy naming from other panels;
   I did not check for the same biological line entering twice under two ModelIDs
   beyond the `PatientID` check in §5.4.
