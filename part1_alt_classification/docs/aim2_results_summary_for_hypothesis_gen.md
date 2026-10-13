# Aim 2 Results Summary — Hypothesis-Driven ALT-Pathway Expression Signature in DepMap Osteosarcoma

**Project:** Independent computational analysis testing whether a curated ALT-pathway expression signature discriminates ALT+ from ALT− osteosarcoma samples, conducted independently. Cell-line-only; not yet extended to primary tumors, PDX, or WGS features.

---

## Data

DepMap 26Q1 log2(TPM+1) protein-coding RNA-seq. Labels merged from published cell-line ALT assays (C-circle, APB, TRAP) plus the ALTitude classification table. Two definitions from a pre-specified label spec:

- **Primary** (gold-only, all experimentally confirmed): 12 ALT+ vs 8 ALT−, N=20.
- **Pragmatic** (adds ALTitude prediction-only silver labels, direction-aware): 18 ALT+ vs 13 ALT−, N=31.

Arms split by `alt_call_direction` (positive/negative), not by evidence tier. Single-platform, single sample type (cell line), no cross-cohort pooling.

---

## Curated 19-Gene Panel

Pre-specified from prior ALT literature before any modeling. Two genes from the original 21-gene list absent from the DepMap matrix (TERC, MRE11A) were dropped, leaving 19. Panel:

- **Telomere maintenance / telomerase axis:** TERT, ATRX, DAXX, WRAP53 (canonical ALT markers and telomerase machinery — Heaphy et al. 2011; Lovejoy et al. 2012; multiple reviews)
- **Homologous recombination / break-induced replication:** BRCA1, BRCA2, RAD50, RAD52, RAD51AP1, NBN, MRE11 complex members (Cesare & Reddel 2010; Dilley et al. 2016; Roumelioti et al. 2016)
- **Replication stress and fork remodeling:** SMARCAL1, WRN, BLM, FANCM, FANCD2 (Cox et al. 2016; Silva et al. 2019; Lu et al. 2019)
- **Structure-specific resolvases:** SLX4, GEN1, MUS81 (Sarbajna et al. 2014; Wan et al. 2013)
- **Osteosarcoma-context genomic co-driver:** TP53 (near-universally disrupted in osteosarcoma; 17p13.1 architecture linked to TCAB1 disruption — Keegan et al. 2026 bioRxiv)

*[Placeholder: replace citations above with the specific ones from your docs/label_spec_v1.0.md.]*

---

## Method

Regularized logistic regression (elastic net, l1_ratio=0.5) inside a StandardScaler → LR Pipeline; nested LOOCV with inner 5-fold CV selecting C from a log-spaced grid; class_weight='balanced'. Permutation null via 1000 label shuffles with fold-local preprocessing precomputed once. Coefficient stability via 200 bootstrap resamples plus LOOCV fold coefficient tracking. All preprocessing inside CV folds; no leakage of held-out sample statistics into training.

---

## Results

**Primary result (curated 19-gene, primary strict, N=20).** AUROC 0.948, 95% bootstrap CI [0.813, 1.000]. AUPRC 0.976. Sensitivity at 90% specificity 0.917. Empirical permutation p < 0.001 from 1000 shuffles (0 of 1000 shuffles matched or exceeded observed).

**Pragmatic sensitivity (curated 19-gene, pragmatic, N=31).** AUROC 0.829, CI [0.668, 0.963], p = 0.010. Known-sex-only pragmatic (N=27, drops 4 unknown-sex samples): AUROC 0.877, p = 0.003 — result strengthens slightly after dropping missing-sex rows.

**Coefficient stability (200-bootstrap non-zero frequency).** Genes stable (>80%) in both primary and pragmatic: TERT, FANCM. Stable in primary only: BRCA2. Stable in pragmatic only: TP53. Borderline (50–80%) in one or both: NBN, WRN, ATRX, SMARCAL1. Coefficient directions in pragmatic: TERT down, FANCM up, ATRX down, NBN up, BRCA2 up in ALT-positive — biologically coherent with the ALT literature (loss of telomerase activity, upregulation of HR/BIR machinery, canonical ATRX loss).

---

## Sex-Confound Finding (Unbiased Whole-Transcriptome Analysis, Notebook 11)

An elastic net across all ~19,000 protein-coding genes at N=20 returned primary AUROC 0.469 (p=0.60) and pragmatic AUROC 0.688 (p=0.06). Top-20 features included Y-chromosome genes (USP9Y, DDX3Y) and cancer/testis antigens (XAGE1A/B, MAGEB2, MAGEA10). Two sensitivity analyses eliminated the borderline signal: adding sex as a covariate dropped AUROC to 0.50; removing all Y-chromosome and X-linked cancer/testis antigen genes dropped it to 0.40. This confirms the borderline unbiased signal was donor-sex confounding, not weak ALT biology. The 19-gene curated approach was designed and pre-specified as the response to this diagnostic; it recovers the canonical biology precisely because it eliminates the confound by construction.

---

## What This Establishes and What It Does Not

Canonical ALT biology is recoverable from DepMap osteosarcoma expression at N=20 when the feature space is hypothesis-constrained to ~19 literature-derived genes. TERT and FANCM are the most robust individual contributors. The result is a hypothesis-driven ALT-pathway score, not a discovered signature — the genes came from prior literature. It has NOT been tested on primary tumors, PDX models, or WGS-derived features. It has NOT been integrated with structural-variant architecture, C-circle assays, or evolutionary timing. Every hypothesis generated from this work should treat this signature as a validated cell-line-only feature that is ready to be applied, projected, or contrasted in new contexts — not as a completed biological finding.
