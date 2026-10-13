# Notebook 09 — ALTitude × DepMap ALT Inventory
_Generated: 2026-05-19T23:13:48Z_
Read-only diagnostic quantifying Aim 2 class balance from ALTitude curated ALT calls for RNA-backed DepMap OS cell lines.
## Data loading
- Manifest rows: **115**
- DepMap RNA ModelIDs (index only): **1719**
- ALTitude panel rows: **959**  unique ModelIDs: **958**
- Route B RNA-backed unique ModelIDs: **32**
## Non-OS lines in RNA-backed set
- **ACH-000592** (TM31): glioblastoma (TM31) — final_alt=ambiguous_or_missing. **Excluded from OS blocks.**
- **ACH-000309** (SKLU1): lung_carcinoma — maps via Route B but **absent from DepMap RNA matrix**; excluded from all blocks.
## Ambiguous or missing ALT calls
- **ACH-000592** (TM31): manifest=alt_unknown, altitude=ambiguous_or_missing
- **ACH-000613** (HOS): manifest=alt_unknown, altitude=ambiguous_or_missing
- **ACH-001275** (SYO1): manifest=alt_unknown, altitude=ambiguous_or_missing
- **ACH-001813** (C242): manifest=alt_unknown, altitude=ambiguous_or_missing
- **ACH-001817** (C393): manifest=alt_unknown, altitude=ambiguous_or_missing
- **ACH-002779** (NCCESOS1C1): manifest=alt_unknown, altitude=ambiguous_or_missing
- **ACH-003176** (NRHGCT1): manifest=alt_unknown, altitude=ambiguous_or_missing
## Disagreements (gold vs ALTitude)
None.
## Block counts summary
| Block | n_positive | n_negative | n_ambiguous | min_arm |
|---|---|---|---|---|
| block1_os_gold | 12 | 7 | 6 | 7 |
| block2_os_sensitivity | 18 | 7 | 6 | 7 |
| block3_all_gold | 12 | 7 | 7 | 7 |
_Within-cohort minimum per arm = 3_
## Provenance breakdown (OS lines only)
| final_provenance         |   alt_gold_positive |   alt_negative_high_conf |   alt_silver |   ambiguous_or_missing |
|:-------------------------|--------------------:|-------------------------:|-------------:|-----------------------:|
| altitude_experimental    |                  12 |                        7 |            0 |                      0 |
| altitude_prediction_only |                   0 |                        0 |            6 |                      0 |
| no_call                  |                   0 |                        0 |            0 |                      6 |
## AIM 2 VIA ALTITUDE-DEPMAP VERDICT
**WITHIN-COHORT VIABLE**
> OS-only gold labels: 12 ALT-positive and 7 ALT-negative RNA-backed DepMap lines. Both arms ≥ 3; all from the same cohort/platform (DepMap 26Q1 log2TPM). Sufficient for within-cohort Aim 2 analysis.
Positive arm: **12**  Negative arm: **7**  Block: block1_os_gold
## Next steps
1. Review and accept this class-balance report.
2. If verdict is WITHIN-COHORT VIABLE, proceed to Aim 2 ALT classifier design.
3. If verdict is POOLED-EXPLORATORY or NOT ESTIMABLE, request Route A/C data.
4. No DE or modeling code has been written.
