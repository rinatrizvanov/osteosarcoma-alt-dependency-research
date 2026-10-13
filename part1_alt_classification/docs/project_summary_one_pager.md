# Project Summary: ALT Classification in Osteosarcoma

**Rinat Rizvanov — June 2026**

---

## What I Built

I built a modular computational pipeline to assess whether a transcriptional signature can distinguish ALT-positive from ALT-negative osteosarcoma cell lines using publicly available data. Starting from five published datasets — DepMap 26Q1, PCGP 2014, de Nonneville 2022, Lawlor 2026 GCCRI, and Taglialatela 2026 — I built a 115-sample master manifest with a four-tier ALT label hierarchy, a two-stage feasibility audit (confounding census, then RNA-availability check), and a recovery strategy when labels and expression data proved to be in completely disjoint datasets. The recovery linked ALTitude's cell-line ALT predictions to DepMap RNA-seq via a cell-line-name join, yielding a usable primary cohort of 20 DepMap osteosarcoma lines (12 ALT+, 8 ALT−) with gold-quality experimental labels.

---

## What I Found

**The data-driven discovery attempt failed honestly.** A fully nested leave-one-out cross-validation of elastic-net logistic regression over 19,177 RNA-seq features returned AUROC = 0.47 (permutation p = 0.60) in the primary cohort. The top features were Y-chromosome and X-linked cancer-testis antigen genes — the model was learning sex, not ALT biology. Systematic sensitivity analyses (sex covariate, Y-chromosome gene exclusion) confirmed that no ALT signal survived sex-confound correction in the unguided setting.

**A hypothesis-driven curated gene score showed strong signal.** Using a pre-specified 19-gene ALT pathway list (TERT, FANCM, BRCA2, WRN, NBN, and 14 others), the same LOOCV procedure achieved AUROC = 0.948 [95% CI: 0.813, 1.000] in the primary cohort and 0.829 [0.668, 0.963] in the pragmatic (N = 31) cohort. The permutation p-value on 1,000 shuffles was < 0.001 for the primary cohort. Bootstrap across 200 resamples showed TERT (−) and FANCM (+) appearing with nonzero coefficients in 100% of resamples in both definitions; they also appeared in every leave-one-out fold. The result survived dropping the four unknown-sex samples (AUROC = 0.877) and was consistent in coefficient direction with established ALT biology: TERT negative (telomerase and ALT are mutually exclusive), FANCM positive (required for break-induced replication at ALT telomeres).

This is a confirmatory test of prior biological knowledge, not a novel signature discovery. The finding is that the prior knowledge is recoverable from this dataset.

---

## Honest Limitations

The primary cohort is N = 20 — small enough that the wide CI matters and the permutation test, while clearing the bar, should not be over-interpreted. ALTitude silver labels (used in the pragmatic cohort) are prediction-based and not independently validated. All modeling uses a single DepMap RNA-seq release; generalization to primary tumors is unknown. Aim 1 (TCAB1 differential expression) is not estimable from public data alone because architecture labels and RNA-seq are in disjoint datasets — this limitation is documented rather than worked around.

---

## What I Need and What Comes Next

The most valuable next step is external validation on the de Nonneville microarray (E-MTAB-11743, 22 osteosarcoma tumors with C-circle ALT calls) using the 19-gene curated score after within-platform normalization. This I can pursue independently. Aim 1 requires internal RNA-seq from the Lawlor 2026 GCCRI PDX models — a PI-level collaboration request. The pipeline infrastructure (manifest, audit, LOOCV framework, sha256-verified outputs) is designed to incorporate new data without re-architecting any existing analysis. If the Lawlor RNA data becomes available, running Aim 1 is a matter of adding samples to the manifest and re-running the audit.
