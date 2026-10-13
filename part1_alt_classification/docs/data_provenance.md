# Data Provenance

Generated and maintained by notebook 00_download_data.ipynb.
Last updated: see notebook output for timestamp.

---

## User-downloaded files

### DepMap 26Q1

| File | Source | Access Date | SHA256 | Notes |
|---|---|---|---|---|
| `data/raw/depmap/Model.csv` | manually downloaded by user (DepMap 26Q1) | see file mtime | filled by notebook 00 | |
| `data/raw/depmap/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv` | manually downloaded by user (DepMap 26Q1) | see file mtime | filled by notebook 00 | replaces OmicsExpressionProteinCodingGenesTPMLogp1.csv |
| `data/raw/depmap/OmicsSomaticMutations.csv` | manually downloaded by user (DepMap 26Q1) | see file mtime | filled by notebook 00 | |
| `data/raw/depmap/OmicsSomaticMutationsMatrixDamaging.csv` | manually downloaded by user (DepMap 26Q1) | see file mtime | filled by notebook 00 | |
| `data/raw/depmap/OmicsCNGeneWGS.csv` | manually downloaded by user (DepMap 26Q1) | see file mtime | filled by notebook 00 | replaces OmicsCNGene.csv |
| `data/raw/depmap/OmicsCNSegmentsWGS.csv` | manually downloaded by user (DepMap 26Q1) | see file mtime | filled by notebook 00 | |
| `data/raw/depmap/OmicsFusionFiltered.csv` | manually downloaded by user (DepMap 26Q1) | see file mtime | filled by notebook 00 | |
| `data/raw/depmap/CRISPRGeneEffect.csv` | manually downloaded by user (DepMap 26Q1) | see file mtime | filled by notebook 00 | |

**DepMap release:** 26Q1
**Note on SV data:** OmicsStructuralVariants.csv is not available in 26Q1 public release.
Cell-line SV-equivalent data reconstructed from OmicsFusionFiltered.csv (WRAP53/TP53 fusions)
and OmicsCNSegmentsWGS.csv (focal CN segments at locus chr17:7,560,000-7,710,000).

### Literature tables

| File | Source | Access Date | SHA256 | Confidence |
|---|---|---|---|---|
| `data/raw/literature_tables/de_Nonneville_2022_table_EV3.csv` | de Nonneville et al., EMBO Mol Med 2022, Table EV3 | see file mtime | filled by notebook 00 | peer_reviewed |
| `data/raw/literature_tables/tcab1_pdx_labels.csv` | Lawlor et al. 2026, bioRxiv 2026.01.02.691864, Fig 1 | see file mtime | filled by notebook 00 | preprint |
| `data/raw/literature_tables/tcab1_pdx_sv_calls.csv` | Lawlor et al. 2026, bioRxiv 2026.01.02.691864, Fig 2 | see file mtime | filled by notebook 00 | preprint |
| `data/raw/literature_tables/taglialatela_2026_celllines.csv` | Taglialatela et al. 2026, bioRxiv 10.64898/2026.02.15.704070 | see file mtime | filled by notebook 00 | preprint |
| `data/raw/literature_tables/ALTitude-table1.xlsx` | ALTitude paper (peer-reviewed) | see file mtime | filled by notebook 00 | peer_reviewed |
| `data/raw/literature_tables/ALTitude-table2,4,5.xlsx` | ALTitude paper (peer-reviewed) | see file mtime | filled by notebook 00 | peer_reviewed |
| `data/raw/literature_tables/ALTitude-table3.xlsx` | ALTitude paper (peer-reviewed) | see file mtime | filled by notebook 00 | peer_reviewed |
| `data/raw/literature_tables/Top3A_EV3.xlsx` | de Nonneville et al. 2022, fallback xlsx | see file mtime | filled by notebook 00 | peer_reviewed |
| `data/raw/literature_tables/chen_2014_tableS3_SVs.xlsx` | Chen et al. 2014, *Cell Reports*, Tables S3a-e | 2026-05-01 | `b2b4bba0dded…` (filled by notebook 00) | peer_reviewed |
| `data/raw/literature_tables/chen_2014_tableS4_TP53_SVs.xlsx` | Chen et al. 2014, *Cell Reports*, Table S4a-b | 2026-05-01 | `d06157011b5e…` (filled by notebook 00) | peer_reviewed |
| `data/raw/literature_tables/chen_2014_tableS6_ATRX_ALT.xlsx` | Chen et al. 2014, *Cell Reports*, Table S6 | 2026-05-01 | `2930b7b81d9c…` (filled by notebook 00) | peer_reviewed |

### CCLE 2019 SvABA structural variants (Ghandi et al., *Nature*)

Substitutes for the missing ``OmicsStructuralVariants.csv`` in the DepMap 26Q1 public release.
Reference genome **GRCh37 / hg19** — coordinates are lifted to GRCh38 in-pipeline by
``src.sv_annotator.annotate_ccle_svs`` (using ``pyliftover`` against the UCSC hg19→hg38
chain file). Liftover failures are logged explicitly and excluded from annotation; failed
lifts are never silently treated as zero-coordinate.

| File | Source | Caller | Reference genome | Access Date | SHA256 | Confidence |
|---|---|---|---|---|---|---|
| `data/raw/ccle/CCLE_translocations_SvABA_20181221.xlsx` | CCLE 2019 (Ghandi et al., *Nature* 2019) | SvABA | GRCh37 | 2026-05-01 | `0fca02bf2e9c…` (filled by notebook 00) | peer_reviewed |

### Chen 2014 PCGP osteosarcoma cohort

Adds 34 patient tumor samples (20 discovery + 14 validation) with peer-reviewed TP53 SV
breakpoints, ATRX status, and telomere-FISH ALT calls. Reference genome **GRCh37**;
coordinates are classified IN-FRAME by ``src.sv_annotator.classify_chen_sv_architecture``
and only the categorical architecture call is stored in the manifest (no liftover is
performed for Chen 2014, by design — see deviation note below).

Sample IDs in S3a/S3b use a ``_M``/``_D`` suffix (sample type) while S6 uses unsuffixed
patient stems; the loader resolves this via stem fallback. One patient (SJOS010) has
two S6 entries (``SJOS010_1``/``SJOS010_2``) which do not unambiguously map to
``SJOS010_M``/``SJOS010_D``; both end up with ``atrx_status='unknown'`` and this is
documented as a known limitation.

### ALTitude table content notes (inspected 2026-04-28)

- **ALTitude-table1.xlsx (Sheet1):** Primary data table, 959 cell lines × 120 features.
  Contains `Curated.Marker` (peer-reviewed gold-standard ALT call used for model training),
  `ALT_pred_class` (ensemble predicted label), `ALT_pred_prob`, `ALT source`, `ALT_assay type`,
  `ALT assay strength`, ATRX mutation details, DAXX mutation details, SMARCAL1 expression
  and dependency, telomere features. **This is the source for peer-reviewed silver-tier ALT labels
  for DepMap cell lines.**
- **ALTitude-table2,4,5.xlsx:** Three sheets — model performance metrics (LOOCV/test AUC,
  sensitivity, specificity for training set n=82), comparison vs Wu Y et al. 2025, and DeLong
  statistical comparisons. No cell-line ALT calls.
- **ALTitude-table3.xlsx (ablation_results_full):** Ablation study results — feature-set × model
  ROC/PR AUC across 181 feature combinations. No cell-line ALT calls.
- **ALTitude CSV exports** (sh1, sh2, sh3): Correspond to the three sheets in table2,4,5.xlsx
  (performance metrics, comparison, statistics). Not cell-line data.

---

## Downloaded by pipeline (notebook 00)
| 2026-04-29 | CRISPRGeneEffect.csv | First column unnamed; ModelIDs in index (standard DepMap convention for gene × model matrices) | Schema check updated to verify file by index content (`ACH-` prefix) rather than column header. Data is correct; file passes verification. |

### TARGET-OS (UCSC Xena)

| File | Source URL | Download Date | SHA256 | Notes |
|---|---|---|---|---|
| `data/raw/target_os/` | UCSC Xena TARGET hub | filled by notebook 00 | filled by notebook 00 | Open-access |

### GDC API

| File | Source | GDC Release | Download Date | SHA256 |
|---|---|---|---|---|
| `data/raw/target_os/` | api.gdc.cancer.gov (open-access tier) | filled by notebook 00 | filled by notebook 00 | filled by notebook 00 |

### Reference files

| File | Source | Download Date | SHA256 |
|---|---|---|---|
| `data/raw/reference/gencode.v45.annotation.gtf.gz` | GENCODE v45, GRCh38 | filled by notebook 00 | filled by notebook 00 |
| `data/raw/reference/hg38.chrom.sizes` | UCSC | filled by notebook 00 | filled by notebook 00 |

---

## Deviations

*Document any schema changes, missing files, or download failures here.*

| Date | File | Issue | Action taken |
|---|---|---|---|
| 2026-05-01 | `data/raw/ccle/CCLE_translocations_SvABA_20181221.xlsx` | Coordinates are GRCh37/hg19 (per file's ``descriptions`` sheet) but the project working frame is GRCh38. | Lifted in-pipeline by ``src.utils.liftover_breakpoint`` (UCSC hg19→hg38 chain via pyliftover ≥ 0.4.1). ``src.sv_annotator.annotate_ccle_svs`` filters chr17 endpoints to ``GRCH37_LOCUS_WINDOW`` *before* lifting, then re-evaluates locus + features in GRCh38. Liftover failures are logged and the failing endpoint is dropped (the partner endpoint is still emitted if it lifts). |
| 2026-05-01 | `data/raw/literature_tables/chen_2014_tableS{3,4,6}_*.xlsx` | Coordinates are GRCh37/hg19. | Classified IN-FRAME (no liftover) by ``src.sv_annotator.classify_chen_sv_architecture`` against ``GRCH37_TP53_CODING``, ``GRCH37_TP53_INTRON1``, ``GRCH37_WRAP53_CODING`` (all verified 2026-05-01 against Ensembl GRCh37 ENST00000269305 / ENSG00000141499 + NCBI Gene). Only the categorical architecture call is stored in the manifest; coordinates never leave GRCh37. |
| 2026-05-01 | CCLE × DepMap OS overlap | Of 32 osteosarcoma ``CCLEName``s in DepMap 26Q1, only ``U2OS_BONE`` (ACH-000364) appears in the CCLE 2019 SvABA cohort, and U2OS has zero chr17 SVs in the TP53/WRAP53 locus window after liftover. | The CCLE integration runs and produces zero new locus calls for the OS cohort. This is recorded as a real, auditable result rather than a bug; cell lines without CCLE coverage retain their CN-segment-based auto-call. |
| 2026-05-01 | `data/raw/literature_tables/chen_2014_tableS6_ATRX_ALT.xlsx` (sample SJOS010) | Patient SJOS010 has two distinct primary tumors recorded as ``SJOS010_1`` and ``SJOS010_2`` in S6, but the discovery cohort lists ``SJOS010_M`` and ``SJOS010_D``. The mapping is not unambiguous. | Both ``SJOS010_M`` and ``SJOS010_D`` are emitted with ``atrx_status='unknown'`` and ``alt_status='alt_unknown'``. This is a known limitation; the actual count of PCGP samples with ATRX/ALT data populated is 18, not the 19 implied by the original spec. |
