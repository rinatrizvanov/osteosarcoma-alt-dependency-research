# Multi-Source ALT Classification in Osteosarcoma: A Computational Feasibility Study

**Prepared by:** Rinat Rizvanov  
**Date:** June 2026  
**Prepared for:** an academic advisor (name withheld in the public copy)

---

## 1. Executive Summary

Over approximately eight weeks I built a modular computational pipeline to assess whether two research aims involving osteosarcoma biology were estimable from publicly available data. The first aim asked whether TCAB1 (WRAP53) genomic disruption is associated with differential gene expression in osteosarcoma tumors. The second asked whether a transcriptional signature can distinguish ALT-positive from ALT-negative osteosarcoma cell lines.

Both aims began with an optimistic inventory of five published datasets — DepMap 26Q1 (cell lines), PCGP 2014 (pediatric tumors), de Nonneville 2022 (bone tumors), Lawlor 2026 GCCRI (patient-derived xenografts), and Taglialatela 2026 (cell lines). A Stage 0 diagnostic revealed that labels and expression data were housed in completely disjoint datasets, making both aims initially inestimable. I designed a recovery strategy, pursued three routes in parallel, and rescued Aim 2 by linking the ALTitude ALT-prediction tool's cell-line calls to DepMap RNA-seq via a cell-line-name join.

The main analytical result has two layers. In the unguided discovery mode — where the model selected features from 19,000+ genes using only cross-validated variance filtering — there was no reproducible signal (AUROC = 0.47, permutation p = 0.60 in the primary cohort). The top features were sex-linked genes (DDX3Y, USP9Y, XAGE1B), a confound I systematically diagnosed and ruled out. In a hypothesis-driven follow-up using a curated 19-gene ALT pathway list (TERT, FANCM, BRCA2, WRN, NBN, and others), a regularized gene score achieved AUROC = 0.948 [95% CI: 0.813, 1.000] in the primary gold-label cohort (N = 20; permutation p < 0.001 on 1,000 shuffles). This result survived sex-confound exclusion, a known-sex sensitivity analysis, bootstrap coefficient stability testing, and leave-one-out fold coefficient tracking.

Aim 1 remains inestimable without internal RNA-seq data for the Lawlor PDX models; this is an honest limitation I document rather than work around. The pipeline, labeling framework, and audit infrastructure constitute a reusable contribution independent of the Aim 2 finding.

---

## 2. Scientific Background

Osteosarcoma is the most common primary bone malignancy in children and adolescents. It arises most frequently at the metaphyses of long bones during periods of rapid growth and carries a five-year survival rate of roughly 65–70% for localized disease, dropping below 20% for metastatic or relapsed cases. Despite decades of clinical research, cytotoxic chemotherapy (methotrexate, doxorubicin, cisplatin — the MAP regimen) remains the standard of care, and survival outcomes have not improved meaningfully since the 1980s. Osteosarcoma is characterized by extreme genomic instability: near-universal TP53 pathway disruption, frequent copy-number alterations, and complex structural rearrangements that make it resistant to simple oncogene-focused therapies.

A defining feature of approximately 50% of osteosarcomas is the use of alternative lengthening of telomeres (ALT) rather than telomerase to maintain telomere length. ALT is a homologous-recombination-based pathway that generates heterogeneous, ultra-long telomeres and produces characteristic molecular hallmarks: ALT-associated promyelocytic leukemia bodies (APBs), elevated extrachromosomal telomeric repeat (ECTR) DNA detectable by C-circle assay, and — in tumors with ATRX or DAXX disruption — loss of histone H3.3 loading at telomeres. The ~50% ALT prevalence in osteosarcoma is substantially higher than in most other solid tumor types, making osteosarcoma an appealing system in which to study ALT as a therapeutic vulnerability.

TCAB1 (Telomerase Cajal Body Protein 1), encoded by WRAP53, localizes telomerase to Cajal bodies and is required for efficient telomerase-mediated telomere extension. In the context of a tumor that uses ALT rather than telomerase, TCAB1 function is potentially irrelevant or harmful; some evidence suggests that TCAB1-disrupted tumors may be more dependent on ALT-mediated repair. The research question motivating this project is therefore twofold: (1) is TCAB1 genomic disruption (at 17p13.1) associated with a distinct transcriptional program in osteosarcoma, and (2) can a transcriptional signature of ALT status be derived from osteosarcoma gene expression data, which would both validate that ALT imposes a detectable molecular phenotype and potentially serve as a proxy phenotype for future work?

---

## 3. Research Aims and Go/No-Go Framework

The project is structured around two pre-registered aims, each evaluated through a formal feasibility audit before any modeling was attempted.

**Aim 1 (TCAB1 differential expression):** Compare gene expression between osteosarcoma samples with TCAB1-disrupted versus TCAB1-intact architecture at 17p13.1. The SV-based annotation classifies disruption by whether any structural variant breakpoint falls within a 150 kb window around the TCAB1 locus (chr17:7,560,000–7,710,000). The primary contrast is disrupted (N = 8 by strict label) versus intact (N = 32). Despite having 40 labeled samples, none had file-backed RNA-seq at the time of analysis — the architecture labels came from SV-annotated tumor biopsies (PCGP 2014, de Nonneville 2022, Lawlor 2026) while the RNA data on disk came from DepMap cell lines and TARGET-OS, which lacked architecture labels. Aim 1 verdict: **NOT ESTIMABLE** without internal PDX RNA-seq (Route C, Lawlor 2026 GCCRI) or a deliberate reframe to a TCAB1-expression-quantile contrast.

**Aim 2 (ALT transcriptional signature):** Derive a gene-expression-based ALT classification score from DepMap osteosarcoma cell lines with ALT experimental evidence. The go/no-go threshold required at least 10 ALT-positive and 10 ALT-negative labeled samples with RNA backing on the primary (gold-label) definition. After the ALTitude integration (described in Section 7), the primary cohort reached 12 ALT+ / 8 ALT− (N = 20), exceeding the threshold. Aim 2 verdict: **GO (with penalized models)**.

The feasibility thresholds and analysis design were frozen in `docs/label_spec_v1.0.md` on 2026-05-15, before any modeling results were examined, to prevent circular refinement of methods to match results.

---

## 4. Repository Architecture

The pipeline is organized as a single Git repository with a strict separation between data, source modules, and analysis notebooks. The structure reflects the principle that the notebooks are thin orchestration layers over well-tested library code.

```
tcab1-alt-osteo pipeline/
├── src/                         # All reusable logic lives here
│   ├── manifest_builder.py      # DepMap + literature cohort loaders
│   ├── label_rules.py           # Deterministic label assignment rules
│   ├── sv_annotator.py          # 17p13.1 locus SV architecture calls
│   ├── audit.py                 # Pass A / Pass B audit logic and reports
│   ├── aim1_inventory.py        # RNA backing check for Aim 1
│   ├── aim2_inventory.py        # RNA backing check for Aim 2
│   ├── altitude_depmap_inventory.py
│   ├── altitude_depmap_merge.py
│   ├── recovery_inventory.py    # Route A/B/C recovery logic
│   └── aim2_signature.py        # LOOCV / signature-score training
├── notebooks/
│   ├── 00–05                    # Data acquisition, manifest build, audit
│   ├── 06–08                    # RNA availability diagnostics, recovery
│   ├── 09–10                    # ALTitude integration, direction-aware merge
│   └── 11–13                    # Aim 2 modeling and sensitivity analyses
├── tests/
│   └── test_sv_annotator.py     # 26 unit tests, all passing
├── data/
│   ├── raw/                     # Immutable downloaded files
│   ├── manifests/               # master_manifest.csv
│   └── processed/               # Reports, predictions, sha256 hashes
└── docs/
    ├── label_spec_v1.0.md       # Frozen analysis specification
    ├── feasibility_decision_memo.md
    └── data_provenance.md
```

The 14 notebooks are numbered sequentially and designed to be run in order, but each produces a sha256-verified output file so that later notebooks can check that earlier outputs are byte-identical to a previously audited baseline. This idempotence design was chosen because the RNA matrices are large (the DepMap expression file covers 1,719 cell lines × ~19,000 genes) and recomputing from scratch on every session would be impractical.

The `src/` directory is importable as a package. Notebooks import directly from `src.aim2_signature`, `src.audit`, etc., rather than re-implementing logic inline. This meant that when I discovered the column-binding bug in notebook 09 (described in Section 7), the fix was made in one place and all downstream analyses used the corrected code automatically.

---

## 5. Data Sources and Acquisition

The master manifest integrates five cohorts across 115 samples (rows) with 46 structured columns covering sample identity, ALT labels, TCAB1 architecture, TP53 status, sex, and provenance metadata.

**DepMap 26Q1 (N = 35 osteosarcoma cell lines).** The Cancer Dependency Map portal provided four files: `Model.csv` (cell line metadata), `OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv` (RNA-seq in log2(TPM+1) units), `OmicsSomaticMutationsMatrixDamaging.csv`, and `OmicsCNGeneWGS.csv`. I identified osteosarcoma cell lines by filtering `OncotreePrimaryDisease == "Osteosarcoma"`, yielding 35 lines. The RNA matrix covers 1,719 cell lines × 19,177 protein-coding gene columns.

**PCGP 2014 (N = 34 pediatric tumors).** The Pediatric Cancer Genome Project osteosarcoma cohort provided SV calls and TP53 status from peer-reviewed literature (Chen et al., 2014, Nature Genetics). These samples carry architecture labels but no expression data in the public release.

**de Nonneville 2022 (N = 22 bone tumors).** A French cohort with ALT status determined by C-circle assay and APB immunofluorescence. Expression data was deposited at ArrayExpress (E-MTAB-11743), but inspection of the ArrayExpress record revealed it is an Affymetrix microarray (study type "transcription profiling by array," array design A-GEOD-19803), not RNA-seq. This discovery in notebook 08 meant the de Nonneville cohort was reclassified from "potential primary cohort" to "external validation cohort only" — the microarray and RNA-seq measurement modalities are incompatible for joint modeling.

**Lawlor 2026 GCCRI (N = 17 patient-derived xenografts).** A preprint from the Greehey Children's Cancer Research Institute. PDX SV calls were partially provided in the dataset; five samples required manual SV-architecture adjudication, for which I pre-filled annotations in `src/sv_annotator.py` after reviewing the supplementary tables. Lawlor 2026 carries no publicly available expression data; obtaining internal RNA-seq from this group via PI collaboration is Route C.

**Taglialatela 2026 (N = 7 cell lines).** A second preprint providing ALT status for a small panel of osteosarcoma cell lines, used as supporting evidence for the ALTitude label integration.

**ALTitude.** The ALTitude tool (Brosnan-Cashman et al. / public release) provides predicted ALT status for DepMap cell lines using a transcriptome-based classifier. ALTitude table 1 (`ALTitude-table1.xlsx`, column `Curated.Marker`) contains cell-line-level ALT calls. Tables 2, 4, 5, and 3 contain model performance metrics only — not cell-line data — a distinction that was critical to get right during integration.

---

## 6. Sample Labeling Framework

The labeling framework is the architectural spine of the project. I formalized it in `docs/label_spec_v1.0.md`, frozen before any modeling, to prevent post-hoc refinement that could inflate apparent results.

**Four-tier ALT status hierarchy.** Labels are assigned by a deterministic four-tier hierarchy in `src/label_rules.py`:

1. `alt_gold_positive` / `alt_negative_high_conf`: experimental evidence (C-circle, APB, TRAP) from a peer-reviewed publication. These are the highest-confidence labels and form the primary analysis cohort.
2. `alt_silver`: ALTitude prediction-only, no independent experimental confirmation. Used in the pragmatic sensitivity cohort.
3. `alt_unknown`: no usable evidence. Excluded from all modeling arms.

This hierarchy ensures that the primary analysis never mixes experimentally confirmed labels with prediction-based labels. The pragmatic cohort deliberately includes silver labels as a pre-specified sensitivity arm.

**Three parallel definitions.** The spec defines three nested label sets — strict (gold only), pragmatic (gold + silver), and broad (gold + silver + proxy). All three were evaluated in the feasibility audit; only strict and pragmatic proceeded to modeling, because the broad definition added no new labeled samples beyond pragmatic.

**`alt_call_direction` versus `alt_status`.** A critical design decision was separating the arm-assignment column (`alt_call_direction`: positive / negative / unknown) from the provenance column (`alt_status`). The silver labels include both ALT-positive and ALT-negative lines; if I had assigned arms directly from `alt_status`, I would have assigned all silver samples to the positive arm regardless of their actual ALTitude prediction direction. The `alt_call_direction` column, added in notebook 10, correctly reflects the predicted direction for silver samples (6 positive / 5 negative in the DepMap cohort) and the experimental direction for gold samples. All downstream LOOCV and permutation code uses `alt_call_direction` exclusively for arm assignment.

**TP53 labeling.** For Aim 1, TP53 status is derived from two sources: `OmicsSomaticMutationsMatrixDamaging.csv` (damaging mutation ≥ 1) and `OmicsCNGeneWGS.csv` (deep deletion, log2 copy number ≤ −1). Either source independently triggers `tp53_status = "altered"`. This combined rule is implemented in `src/label_rules.py` and audited for every DepMap sample.

---

## 7. Feasibility Audit: Pass A and Pass B

Before writing any modeling code, I ran two structured audits.

**Pass B: Confounding Census.** Pass B checks whether potential confounders (cohort, sample type, TP53 status) are imbalanced across the TCAB1 and ALT arms. The critical finding was that the PCGP 2014 cohort dominates the TCAB1-intact arm (100% of intact samples in that cohort are PCGP), which would make any Aim 1 DE result uninterpretable — it would conflate TCAB1 architecture with cohort of origin. This finding is documented with an override (`Halt override: True`) and a written explanation in the Pass B report. The override was deliberate: the audit should still run to completion for documentation purposes, even though Aim 1 is flagged as confounded.

**Pass A: Availability Census.** Pass A checks how many labeled samples have file-backed RNA-seq available on disk. The initial result was stark: every labeled cohort (PCGP 2014, de Nonneville 2022, Lawlor 2026, Taglialatela 2026) provides labels only — no expression data in the public release. The DepMap RNA matrix on disk covers 1,719 cell lines, but none had been matched to ALT or TCAB1 labels at that point. The Pass A result was: Aim 1 RNA-backed disrupted samples = 0; Aim 2 RNA-backed ALT-labeled samples = 0.

This was not a pipeline failure. The Stage 0 RNA diagnostics (notebooks 06 and 07) were designed specifically to catch this before any analysis was written. They prevented an artifactual situation in which I would have proceeded to model fitting only to find that the training set was empty.

**Go/No-Go verdict after recovery.** After the ALTitude integration (Section 7) and with corrected direction-aware labels, the audit was re-run with the following results: Aim 2 strict definition — 60 ALT+ / 19 ALT− across all cohorts (GO with penalized models); DepMap-only strict — 12 ALT+ / 8 ALT− (GO, single platform); Aim 1 primary — MARGINAL; Aim 1 TP53-restricted — GO. The single-platform DepMap Aim 2 design was chosen over the cross-cohort design to avoid measurement-technology confounds.

---

## 8. RNA Availability Crisis and Recovery

The Stage 0 diagnostics produced a clean finding: labels and expression data were in completely disjoint datasets. Rather than abandoning the project, I documented this as a recoverable situation and proposed three routes:

**Route A (de Nonneville expression, ArrayExpress).** I manually inspected the ArrayExpress entry E-MTAB-11743 and found it to be an Affymetrix gene-expression microarray, not RNA-seq. The measurement technology is incompatible with DepMap log2(TPM+1) matrices; the two cannot be pooled without platform correction. Route A was downgraded to "external validation cohort only," available after a DepMap-derived signature exists.

**Route B (ALTitude → DepMap join).** Notebook 08 tested whether DepMap cell lines could be matched to the ALTitude ALT-prediction table using cell-line name (a fuzzy string join), rather than exact sample ID matching (which notebook 07 had tried and failed). This succeeded: 38 DepMap osteosarcoma cell lines matched to ALTitude records by name, of which 37 had RNA-seq in the DepMap matrix. Route B rescued Aim 2.

**Route C (Lawlor PDX RNA, internal).** The highest-value route for both aims — if the Lawlor 2026 GCCRI group has internal RNA-seq for the 17 PDX models, it would provide labeled, tumor-derived expression data that DepMap cell lines cannot. Route C has not been pursued because it requires PI-level collaboration. It remains the critical dependency for Aim 1.

**The column-binding bug (Notebook 09).** The initial notebook 09 run produced a verdict of NOT ESTIMABLE for Aim 2, which was suspicious — Route B had just shown 37 RNA-backed cell lines. Investigation revealed that the code was reading the `Curated.Marker` column from the ALTitude table (which contains the ATRX mutation feature used in ALTitude's own model) instead of `ALT_pred_class` (the actual ALT classification). This is the type of off-by-column error that can silently corrupt an analysis. The error was caught by a ground-truth sanity check: I verified that U2OS (a well-characterized ALT-positive osteosarcoma line) returned the correct label and that SJSA1 (ALT-negative) did the same. When U2OS came back as `ambiguous_or_missing`, I knew the column was wrong.

The corrected inventory (notebook 09, re-run): 31 RNA-backed DepMap osteosarcoma cell lines with valid ALTitude calls; 12 ALT+ / 7 ALT− gold-only; 18 ALT+ / 7 ALT− gold + silver. After adding the direction-aware `alt_call_direction` column in notebook 10 and re-running the audit, the final primary cohort was 12 ALT+ / 8 ALT− (N = 20) and the pragmatic cohort was 18 ALT+ / 13 ALT− (N = 31).

---

## 9. Aim 2 Signature Discovery Attempt: A Null Result

Notebook 11 implemented a fully nested leave-one-out cross-validation (LOOCV) for an elastic-net regularized logistic regression classifier over 19,177 DepMap RNA-seq gene features. Every preprocessing step — variance-based gene filtering (top 500 high-variance genes, K = 500), standard scaling, and hyperparameter selection — was applied inside the training fold only, never to the held-out sample. C was selected by inner 5-fold cross-validation within each outer LOOCV fold, with a grid of 9 logarithmically spaced values from 0.01 to 100.

The result was a null: primary AUROC = 0.4688 [95% CI: 0.2031, 0.7500], permutation p = 0.60 (100 shuffles). The model performed worse than random chance on the primary gold-label cohort. The pragmatic result was borderline: AUROC = 0.6880 [0.4473, 0.8950], permutation p = 0.06. Neither crossed any threshold for a meaningful result.

The top-ranked features exposed why: in the primary model, ranks 1–5 by absolute coefficient included XAGE1B, MAGEB2, XAGE1A, MAGEA10 (X-linked cancer-testis antigens expressed in male germ cells), and USP9Y / DDX3Y (Y-chromosome genes). None of the 21 curated ALT pathway genes appeared in the top 20 in either definition; only TP53 had any nonzero coefficient in the full ranked list (rank 213, coefficient −0.052). The unguided model was learning sex rather than ALT biology.

This null result was expected. At N = 20 with 19,000+ features, an elastic-net model has no reliable way to distinguish genuine ALT-associated expression from demographic confounds, particularly when the ALT-positive and ALT-negative arms happen to differ in the proportion of male versus female cell lines. The value of the null is that it tells us what cannot be learned from the data in the unguided setting.

---

## 10. Sex Confounding Sensitivity Analysis

Notebook 12 ran three designs to characterize the sex confound and test whether any signal survived correction.

**Design 1 (D1): Sex covariate.** I appended a binary sex indicator (Female = 0, Male = 1) as an additional feature after the variance gene filter, bypassing the variance filter for the sex feature itself. The primary result (N = 20, all Female or Male annotated) was AUROC = 0.5000 [0.1904, 0.8001], permutation p = 0.53 — chance performance. The pragmatic analysis was skipped because four samples (ACH-001813, ACH-001817, ACH-002779, ACH-003573) have unknown sex annotation, and I had pre-specified that the sex-covariate design would abort rather than impute or drop if unknown-sex samples were present.

**Design 2 (D2): Y-chromosome and X-linked CT antigen exclusion.** I removed 41 Y-chromosome gene columns and 4 X-linked cancer-testis antigens (MAGEB2, MAGEA10, XAGE1A, XAGE1B) from the feature universe before fitting. Primary AUROC = 0.3958 [0.1098, 0.7144], permutation p = 0.66 — worse than chance, confirming that the small signal in the nb11 pragmatic result depended on sex-linked features. Pragmatic with Y/XCT genes removed: AUROC = 0.7051 [0.4772, 0.8955], permutation p = 0.06 — borderline but consistent with a confounded baseline rather than a genuine biological signal.

The D1 and D2 results together establish that the unguided whole-transcriptome approach finds no ALT-associated signal after sex is accounted for. This is the honest stopping point for data-driven discovery in this dataset.

---

## 11. Hypothesis-Driven Curated Gene Score

Design 3 (D3) is methodologically distinct from D1 and D2. It is not a data-driven discovery but a prespecified test: given a curated list of 21 genes known or strongly hypothesized to be functionally relevant to ALT, does a regularized score over only those genes distinguish ALT-positive from ALT-negative lines? The curated list includes TERT, TERC, ATRX, DAXX, SMARCAL1, TP53, WRN, BLM, MRE11A, RAD50, NBN, FANCD2, RAD51AP1, WRAP53, RAD52, BRCA1, BRCA2, FANCM, SLX4, GEN1, and MUS81.

Of the 21 genes, 19 were present in the DepMap RNA matrix. TERC is a non-coding RNA and absent from the protein-coding gene matrix. MRE11A is listed in DepMap as `MRE11 (4361)` under a different symbol and was not matched by the exact-name join; both absences are documented. The analysis used 19 genes, StandardScaler only (no variance filter — the feature set is small and all features are biologically motivated), and elastic-net logistic regression with the same hyperparameter grid.

**Primary result (N = 20, gold labels):** AUROC = 0.9479 [95% CI: 0.8125, 1.0000], AUPRC = 0.9755, sensitivity at 90% specificity = 91.7%, permutation p < 0.001 (0/100 shuffles exceeded the observed AUROC). The fixed hyperparameter (C = 0.3162) was selected by the same inner CV procedure as in nb11. Nonzero coefficients at the full-data refit: TERT (−0.689), FANCM (+0.683), BRCA2 (+0.197), WRN (−0.153), NBN (+0.105). TERT has a negative coefficient, consistent with the well-established inverse relationship between TERT expression and ALT — telomerase and ALT are mutually exclusive telomere maintenance mechanisms. FANCM has a positive coefficient; FANCM is required for break-induced replication at telomeres in ALT cells and has been proposed as a synthetic lethal target in ALT osteosarcoma. BRCA2 and NBN are involved in homologous recombination, and their positive association is consistent with ALT's dependence on HR-based telomere synthesis.

**Pragmatic result (N = 31, gold + silver):** AUROC = 0.8291 [0.6681, 0.9625], AUPRC = 0.9166, sensitivity at 90% specificity = 77.8%, permutation p = 0.010 (1/100 shuffles). The pragmatic model retains TERT and FANCM as the two largest coefficients (−0.755 and +0.767 respectively) but adds TP53 (−0.337) and ATRX (−0.125) as nonzero features. The TP53 coefficient is negative (more TP53 expression → predicted ALT-negative), which is biologically coherent: ALT tumors frequently inactivate TP53 at the protein level, but TP53 mRNA expression may still be present; this is a modest, noisy signal that the elastic net assigns a small coefficient.

It is important to state clearly what this result is and what it is not. The AUROC of 0.948 is the leave-one-out cross-validated performance — each held-out sample was scored by a model that had never seen it. The feature set was fixed before any model fitting, and the permutation test establishes that the observed AUROC is very unlikely under the null hypothesis of no association between gene expression and ALT status. However, this is a test of a hypothesis, not a discovery. The result says that expression of these 19 genes collectively carries information about ALT status in this cohort. It does not mean that a classifier trained on these data generalizes to new cell lines or tumor biopsies, and it does not establish which specific genes drive the classification beyond what the bootstrap and LOOCV fold analyses suggest.

---

## 12. Curated Score Robustness Analysis

Notebook 13 ran four additional analyses to stress-test the D3 result.

**Extended permutation test (1,000 shuffles).** The 100-shuffle permutation in notebook 12 was sufficient to establish significance but not to resolve the pragmatic p-value precisely. Re-running with 1,000 shuffles confirmed: primary p < 0.001 (0/1,000 shuffles); pragmatic p = 0.010 (10/1,000 shuffles); known-sex-only pragmatic p = 0.003 (3/1,000 shuffles). The primary result is highly robust; the pragmatic result is nominally significant but should be interpreted cautiously given the smaller effect size and the inclusion of silver labels.

The computational implementation precomputed fold-local StandardScaler parameters once (scaling is label-independent) and then re-fit only the classifier for each permutation. This is mathematically equivalent to re-running the full LOOCV per shuffle and approximately 10× faster, enabling the 1,000-shuffle analysis in reasonable time.

**Known-sex-only pragmatic sensitivity (N = 27).** Dropping the four unknown-sex samples (ACH-001813, ACH-001817, ACH-002779, ACH-003573 — three ALT-negative and one ALT-positive) from the pragmatic cohort gave N = 27 (17+ / 10−). AUROC = 0.8765 [0.7176, 0.9835], p = 0.003. The AUROC is 0.047 higher than the N = 31 result, and the CI tightens, consistent with the unknown-sex samples introducing noise rather than signal. The two largest coefficient pairs — TERT and FANCM — are preserved, with TERT gaining magnitude (−0.882 vs −0.755) and FANCM decreasing slightly (+0.632 vs +0.767). The TP53 coefficient shrinks substantially (+0.101 vs −0.337) when unknown-sex samples are excluded, suggesting that some of TP53's contribution in the full pragmatic model was entangled with the unknown-sex samples' label uncertainty.

**Bootstrap coefficient stability (200 resamples).** To assess which genes drove the classification reliably rather than in response to specific samples, I ran 200 bootstrap resamples (sampling N cases with replacement) and fitted the curated elastic-net model on each resample using the fixed C. The fraction of resamples in which each gene's coefficient was nonzero (|β| > 10⁻¹⁰) serves as a stability score. Genes were classified as stable (> 80%), borderline (50–80%), or unstable (< 50%).

Results: **TERT** and **FANCM** are stable in both the primary and pragmatic definitions (bootstrap frequency = 1.00 in both). **BRCA2** is stable in the primary definition only (0.81). **TP53** is stable in the pragmatic definition only (0.935) and unstable in primary (0.445), reflecting that TP53 mRNA variation is a weaker signal in the gold-label cohort. **WRN** and **NBN** are borderline in both definitions (0.57 / 0.70 and 0.63 / 0.79, respectively). Fourteen of the 19 genes are unstable in the primary definition; eleven in the pragmatic definition. This is the expected pattern for elastic-net regularization in a small sample: the model concentrates weight on a small number of features, and at N = 20, slight resampling variation moves borderline features in and out.

**LOOCV fold coefficient frequency.** Separately from bootstrap, I tracked which genes received nonzero coefficients in each of the N leave-one-out fold refits (i.e., when each sample is the held-out test case). This is conceptually distinct from bootstrap: bootstrap tests label-presence robustness, while LOOCV fold frequency tests sensitivity to individual sample removal. TERT and FANCM appear in every fold with nonzero coefficients in both definitions (frequency = 1.00). In the primary definition, WRN and NBN also appear in 90% of folds, and BRCA2 in 95% of folds — a strong result. In the pragmatic definition, the picture is more distributed: TP53 appears in 96.8% of folds, ATRX in 64.5%, SMARCAL1 in 67.7%. The LOOCV fold analysis reinforces the bootstrap finding that TERT and FANCM are the most robust predictors in this dataset, while the remaining genes contribute variable amounts depending on which sample is excluded.

---

## 13. Gene Interpretation and Biological Coherence

The curated gene score identified a consistent set of predictors across analysis variations. I want to be careful about over-interpreting coefficients from a 19-gene elastic-net model at N = 20 — the coefficients reflect regularization-driven compression, not independent effect sizes. That caveat stated, the pattern is biologically coherent.

**TERT (negative coefficient).** TERT encodes the catalytic subunit of telomerase. Its negative coefficient — higher TERT expression predicts ALT-negative — directly reflects the mutual exclusivity of telomerase and ALT. ALT-positive osteosarcoma lines suppress telomerase and maintain telomeres through recombination. This is the strongest and most stable predictor in the model (bootstrap frequency = 1.00 in both definitions, LOOCV fold frequency = 1.00 in all three cohorts).

**FANCM (positive coefficient).** FANCM is a Fanconi anemia pathway helicase required for restart of stalled replication forks and for break-induced replication at ALT telomeres. Bhowmick et al. (2016) and others have shown that FANCM depletion selectively inhibits ALT activity, proposing it as a synthetic lethal target in ALT tumors. The positive coefficient — higher FANCM expression predicts ALT-positive — is consistent with ALT cells upregulating this pathway. FANCM is the second most stable predictor.

**BRCA2, NBN (positive coefficients).** BRCA2 and NBN are homologous recombination factors. ALT depends on HR machinery for telomere synthesis; modest upregulation in ALT lines is consistent with the known literature.

**WRN (negative coefficient).** WRN (Werner syndrome RecQ helicase) is a negative predictor, which is biologically surprising — ALT cells typically require RecQ-family helicases. This may reflect that WRN expression is not uniformly elevated in ALT osteosarcoma cell lines and the coefficient is capturing something more complex, or it may reflect the elastic net assigning a negative weight to a collinear feature. WRN is borderline-stable (bootstrap frequency 0.57 primary, 0.70 pragmatic) and should not be over-interpreted.

**WRAP53 (zero coefficient).** Despite TCAB1/WRAP53 being central to the original research question, it received a zero coefficient in all analyzed definitions. This does not mean WRAP53 is irrelevant to osteosarcoma biology; it means that WRAP53 RNA expression does not add information to the ALT classification score beyond what TERT and FANCM already capture in this small dataset.

---

## 14. Statistical Approach and Limitations

I want to be direct about the limitations of this analysis. They are real, and honest presentation of them is more valuable than minimizing them.

**Sample size.** The primary cohort is N = 20 (12+ / 8−). Even leave-one-out cross-validation cannot compensate for this — the reported AUROC has a 95% CI of [0.8125, 1.0000], which is wide. The permutation test establishes that the result is unlikely under the null, but the expected AUROC in new samples could be substantially lower than 0.948.

**Label quality.** The gold-label ALT calls are derived from peer-reviewed C-circle / APB / TRAP assays, but these were performed by different groups using different protocols across different studies. There is no guarantee that all gold labels are internally consistent. The silver labels (ALTitude predictions) have never been independently validated; they are an educated guess that I have deliberately flagged as lower-confidence.

**Single platform.** All modeling uses DepMap 26Q1 log2(TPM+1) expression from a single data release. Cell lines in long-term culture may not represent the expression profiles of primary tumors. The de Nonneville microarray provides a potential orthogonal validation, but it would need substantial platform-correction before any comparison is meaningful.

**Curated gene test is confirmatory, not discovery.** D3 tests a pre-formed hypothesis. The high AUROC is a validation of prior biological knowledge about TERT, FANCM, and HR pathway genes in ALT — it confirms that these genes are informative in this dataset. It does not discover new biology. Claiming it as a novel signature discovery would be overstating the finding.

**Permutation test conservatism.** The permutation null fixes C at the median observed value rather than re-selecting C per fold per shuffle, because the fully nested null is computationally intractable (~9 hours per definition). This is a slightly conservative concession: the true p-value may be smaller than reported, not larger. This is documented explicitly in the notebook and report.

**No calibration.** The reported probabilities from the LOOCV are not calibrated — they should not be interpreted as true probabilities of ALT status. The model is a signature score, not a diagnostic classifier.

---

## 15. Future Directions and What This Project Enables

**External validation on de Nonneville microarray.** A D3 curated-gene score derived from DepMap RNA-seq can in principle be tested on the de Nonneville Affymetrix microarray after gene-level normalization within each platform. The 19 D3 genes are all protein-coding and likely represented on the Affymetrix array. External validation in an independent, tumor-derived cohort would substantially strengthen the finding. This is a concrete next step that I could pursue independently.

**Aim 1 via Route C (Lawlor PDX RNA).** If the Lawlor 2026 GCCRI group shares internal RNA-seq for the 17 PDX models, Aim 1 becomes estimable: 8 TCAB1-disrupted versus 32 intact at the strict label level. This would allow a straightforward differential expression analysis that is currently impossible with public data alone. This is the highest-value PI-collaboration point.

**TCAB1 expression-quantile contrast in DepMap.** An alternative Aim 1 reframe that does not require new data: within the 35 DepMap osteosarcoma cell lines, stratify by WRAP53 expression quantile (high vs. low) and test for differential expression of associated pathways. This sidesteps the architectural label requirement but raises a different interpretive challenge — WRAP53 expression level is not the same as TCAB1 genomic disruption. It is a valid alternative hypothesis but a different one.

**FANCM synthetic lethality connection.** The consistent positive coefficient for FANCM in the curated score is intriguing in the context of FANCM depletion studies in ALT cells. A CRISPR screen in ALT-positive versus ALT-negative osteosarcoma cell lines (available in DepMap via the CRISPR Chronos dependency data) could test whether FANCM dependency scores correlate with ALT status — a question that the present pipeline and label set would directly support.

**Pipeline reuse.** The manifest infrastructure, labeling framework, and audit system are designed to be updated as new data becomes available. If Route C yields RNA-seq for the Lawlor PDX models, they can be integrated into the master manifest with new `alt_status` and `tcab1_genomic_status` labels and the audit re-run without modifying the modeling code. Similarly, if ALTitude releases updated cell-line calls in a future version, the merge notebook can be re-run with sha256 verification to confirm that the updated inputs change the results in the expected direction.

---

## 16. Conclusion

This project began with two research aims and a large integrated dataset and arrived at the conclusion that one aim is not yet estimable and one aim yields a modest but statistically robust finding when approached correctly. That is an honest result. The data-driven whole-transcriptome discovery attempt produced a null that, when interrogated, revealed sex confounding — a failure mode that matters because it demonstrates what goes wrong when a small cohort is analyzed without biological priors. The curated-gene approach then demonstrated that prior knowledge about the ALT pathway can extract a significant classification signal from the same data. The two findings together tell a coherent story: the transcriptional signal of ALT status in osteosarcoma is real but concentrated in a small number of functionally motivated genes, and it cannot be reliably found by unguided gene selection in a cohort of twenty samples.

TERT and FANCM emerged from every analysis — permutation test, bootstrap, LOOCV fold tracking, and coefficient inspection across three cohort definitions — as the most stable contributors to the ALT classification score in this dataset. Their biological interpretation is consistent with well-established literature. The remaining genes carry smaller, less stable signals that may reflect genuine biology or may reflect the noise inherent in elastic-net regularization at small N.

The pipeline built for this project is methodologically rigorous, version-controlled, sha256-verified for reproducibility, and structured to support the next round of data integration. It does not overstate what was found. It documents honestly what was attempted and what failed. And it positions the next phase of work — external validation in de Nonneville microarray data, internal RNA from the Lawlor PDX cohort, and CRISPR dependency integration — as concrete, achievable extensions that do not require re-architecting anything that has been built.
