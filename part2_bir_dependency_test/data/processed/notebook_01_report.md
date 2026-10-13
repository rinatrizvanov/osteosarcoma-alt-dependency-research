# Notebook 01 — Data acquisition and verification

DepMap release **26Q1**. All inputs already on disk; nothing downloaded.

## File integrity

| File | MB | sha256 match | Read-only |
|---|---|---|---|
| CRISPRGeneEffect.csv | 440.6 | yes | yes |
| OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv | 305.0 | yes | yes |
| Model.csv | 0.7 | yes | yes |
| OmicsCNGeneWGS.csv | 409.0 | yes | yes |
| OmicsCNSegmentsWGS.csv | 63.7 | yes | yes |

## Cohort availability

- osteosarcoma lines in `Model.csv`: 35
- with a CRISPR screen: 26
- in the expression matrix: 34
- with both: 26

## Gene availability (outcome genes + stratifier)

| Gene | Missing across the 26 screened osteosarcoma lines |
|---|---|
| POLD3 | 0 |
| RAD52 | 0 |
| SMARCAL1 | 0 |
| BRCA1 | 0 |
| BRCA2 | 0 |
| RAD51 | 0 |
| WRAP53 | 0 |

No imputation is needed or permitted: every primary gene is fully observed.

## Note

This notebook reads gene-effect *availability* only. No gene-effect value is
retained or inspected; the outcome is first loaded in notebook 05, after the
pre-specification freeze.
