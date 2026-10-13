# Notebook 02 — Strand-awareness check on the TCAB1/WRAP53 stratifier

Evaluates the branch rule committed in `docs/analytical_prespec.md` §4c. Reads no
gene-effect values; the gate check asserts `src.dependency_stats` is not loaded.

## Risk being tested

WRAP53 exon 1 carries a ~1,419 bp head-to-head antisense overlap with TP53. On an
unstranded library, reads there cannot be assigned to the correct gene, which would
contaminate WRAP53 TPM with TP53 and invalidate the stratifier.

## Documentary check: INCONCLUSIVE

| Source | Finding |
|---|---|
| Broad `DepMap_processing_pipeline.md` | GTEx pipeline, STAR + RSEM, GENCODE v38 / GRCh38. Strandedness not stated. |
| DepMap 26Q1 release announcement | No expression-pipeline or strandedness detail. |
| `depmap_omics` per-quarter configs | Not publicly reachable. |
| Third-party data dictionaries (non-authoritative) | RSEM "unstranded mode"; batch-corrected variant "stranded or unstranded depending on each sample's sequencing protocol". |
| DepMap 25Q2 release notes | STAR 2.7.11b; quantification moved to Salmon v1.10.0. |

No authoritative source states strandedness. Best available indication is unstranded,
possibly mixed per sample. Recorded as an unresolved assumption, not a finding.

## Empirical diagnostics

- rho(WRAP53, TP53) across 1719 lines = +0.3539 (p = 6.89e-52)
- rho(WRAP53, TP53) within 34 osteosarcoma lines = +0.2443 (p = 0.164)
- Expression-matched null over 1000 genes: median +0.1344, p95 +0.3663, p99 +0.5007
- Observed rho sits at the 94.40th percentile of that null
- rho(WRAP53 expr, WRAP53 CN) in osteosarcoma lines = +0.5537 (p = 0.00083)

| Quantity | Value |
|---|---|
| Lines (all) | 1719 |
| Lines (osteosarcoma) | 34 |
| rho(WRAP53, TP53), all lines | +0.3539 (p=6.89e-52) |
| rho(WRAP53, TP53), osteosarcoma | +0.2443 (p=0.164) |
| Expression-matched null genes | 1000 |
| Null median / p95 / p99 | +0.1344 / +0.3663 / +0.5007 |
| Observed percentile in null | 94.40 |
| rho(WRAP53 expr, WRAP53 CN), osteosarcoma | +0.5537 (p=0.00083) |
| Exceeds null 99th percentile | False |
| CN tracks better than TP53 | True |

## Branch selected

**`direct`**

Use WRAP53 log2(TPM+1) directly as the Section 4 stratifier.

DIRECT: rho(WRAP53, TP53) is within what expression-matched gene pairs show (below the 99th percentile). No evidence of antisense cross-assignment. Use WRAP53 log2(TPM+1) directly as the Section 4 stratifier.

## Stratifier candidates

Direct and TP53-residualised WRAP53 computed for all 34 osteosarcoma
lines with RNA, written to `wrap53_stratifier_candidates.csv`. Rank agreement between
the two (Spearman): 0.9287.

## Figure

`figures/strand_check_null.png` — observed rho against the expression-matched null.
