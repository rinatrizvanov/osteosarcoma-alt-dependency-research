# Essentiality scan — iteration 03, BIR module coherence

Feasibility and **marginal description only**. This file exists so that
`prespec.md` can be written with realistic expectations about missingness,
variance and dynamic range in Chronos gene effect for the candidate gene sets.
It contains **no cross-gene statistic of any kind** — see
[§6 "What I did not compute, by instruction"](#6-what-i-did-not-compute-by-instruction).

- Date of scan: 2026-10-12 (UTC timestamp of the run: 2026-10-12T06:49:56Z)
- Release: DepMap 26Q1
- git HEAD at scan time: `cf73508afd9cf191ddd416ced7ee23ff9097bd01`
- Interpreter: `uv run python` in the project venv — Python 3.12.13,
  pandas 2.3.3, numpy 2.5.4, scipy 1.18.1
- Loaders reused verbatim, no reimplementation: `src.cohort.build_cohort`,
  `src.cohort.assert_all`, `src.dependency_stats.crispr_gene_columns`,
  `src.dependency_stats.load_all_gene_effects_for_line`,
  `src.utils.verify_expected`
- Sources read: `data/raw/depmap/Model.csv`,
  `data/raw/depmap/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv`
  (membership probe only, via `build_cohort`), `references/master_manifest.csv`,
  `data/raw/depmap/CRISPRGeneEffect.csv`
- `CRISPRGeneEffect.csv` sha256 verified against `src.utils.EXPECTED_SHA256`:
  `e610a4cefb13a82b5b256b47eb08b63ff14843f8dbd0fb164bc0a32688e5b89e` — **match = True**
  (440,646,050 bytes)

Sign convention, binding (prespec, unchanged): Chronos gene effect — **more
negative = more essential**; 0 ≈ non-essential, −1 ≈ median common essential.

---

## 1. What I ran

1. `cohort = src.cohort.build_cohort()` then `src.cohort.assert_all(cohort)` —
   both ground-truth anchors (ACH-000364 U2OS positive, ACH-000748 SJSA1
   negative), the expected set sizes and the expected tier × direction counts
   all passed. No exception raised.
2. `src.utils.verify_expected(src.utils.CRISPR_PATH)` — digest match (above).
3. `src.dependency_stats.crispr_gene_columns([g])` once **per gene**, for each
   of the 35 requested symbols, to resolve the exact `SYMBOL (entrez)` header
   name and to detect absence. Called one symbol at a time so that a single
   missing gene would report as absent rather than aborting the whole scan
   (the helper raises `KeyError` listing every missing symbol).
4. `src.dependency_stats.load_all_gene_effects_for_line(pos_ids + neg_ids)` —
   the 25 screened osteosarcoma rows, all 18,531 gene columns. This is the
   permitted loader used because `load_gene_effects` **raises** on
   any NaN ("Stage 1 does not impute the outcome") and therefore cannot be used
   to *characterise* missingness, which is the question here.
5. Per-gene marginals computed from **one column at a time**
   (`eff.loc[ids, col]`, a 1-D Series), separately for the 16 ALT+ and the
   9 ALT− lines. Two independent runs of the scan script gave byte-identical
   numbers.

### Cohort trace (membership only, reproduced from `cohort.trace`)

| Step | N |
|---|---|
| Model.csv rows matching 'steosarcoma' | 35 |
| present in CRISPRGeneEffect.csv (screened) | 26 |
| present in expression matrix | 26 |
| joined to frozen manifest (cohort == 'DepMap_26Q1') | 26 |
| alt_status in ('alt_gold_positive', 'alt_silver', 'alt_negative_high_conf') | 25 |
| alt_call_direction in ('positive', 'negative') | 25 |
| ALT+ (primary contrast set) | 16 |
| ALT- (negative-control set) | 9 |

Line sets used below: **ALT+ n = 16**, **ALT− n = 9**,
`excluded_unknown` n = 1 (not used).

---

## 2. Gene presence in `CRISPRGeneEffect.csv`

**All 35 requested genes are present.** No symbol required an alias, and no
symbol was absent: `crispr_gene_columns` resolved every one of them (35/35).
The exact column names are carried in the tables in §3 and §4, so the prespec
can quote them verbatim.

Genes absent from the matrix: **none** (`[]`).

---

## 3. ALT+ lines (n = 16) — per-gene marginals

Each row is computed from that gene's column alone. `n < −0.5` and `n < −1.0`
are counts of lines at or beyond those Chronos thresholds. SD is `ddof = 1`.

| Gene | Column in `CRISPRGeneEffect.csv` | n | miss | mean | SD | median | q25 | q75 | IQR | min | max | n < −0.5 | n < −1.0 |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| **BIR module (Stage 1)** | | | | | | | | | | | | | |
| POLD3 | `POLD3 (10714)` | 16 | 0 | -1.7238 | 0.3191 | -1.6664 | -1.9935 | -1.5094 | 0.4841 | -2.4133 | -1.3121 | 16 | 16 |
| RAD52 | `RAD52 (5893)` | 16 | 0 | +0.0171 | 0.1727 | -0.0100 | -0.0235 | +0.0991 | 0.1227 | -0.4162 | +0.3782 | 0 | 0 |
| SMARCAL1 | `SMARCAL1 (50485)` | 16 | 0 | -0.4907 | 0.3950 | -0.4624 | -0.7767 | -0.2218 | 0.5549 | -1.1035 | +0.2238 | 8 | 3 |
| **HR comparator (Stage 1)** | | | | | | | | | | | | | |
| BRCA1 | `BRCA1 (672)` | 16 | 0 | -0.5731 | 0.4395 | -0.5548 | -0.6868 | -0.3009 | 0.3859 | -1.8413 | +0.0152 | 9 | 2 |
| BRCA2 | `BRCA2 (675)` | 16 | 0 | -0.5372 | 0.3749 | -0.4064 | -0.8579 | -0.2595 | 0.5983 | -1.2213 | -0.0170 | 7 | 1 |
| RAD51 | `RAD51 (5888)` | 16 | 0 | -1.1850 | 0.3168 | -1.1274 | -1.3069 | -1.0124 | 0.2945 | -1.8214 | -0.6105 | 16 | 13 |
| **Extended BIR / fork restart** | | | | | | | | | | | | | |
| MUS81 | `MUS81 (80198)` | 16 | 0 | -0.3532 | 0.1528 | -0.3778 | -0.4582 | -0.1993 | 0.2589 | -0.6123 | -0.1329 | 3 | 0 |
| PIF1 | `PIF1 (80119)` | 16 | 0 | -0.1087 | 0.1182 | -0.0930 | -0.1681 | -0.0441 | 0.1240 | -0.3393 | +0.0699 | 0 | 0 |
| POLA1 | `POLA1 (5422)` | 16 | 0 | -0.9127 | 0.3493 | -0.8941 | -1.1941 | -0.7383 | 0.4557 | -1.5419 | -0.3345 | 14 | 5 |
| RAD18 | `RAD18 (56852)` | 16 | 0 | -0.1095 | 0.1775 | -0.0971 | -0.1808 | -0.0045 | 0.1763 | -0.4522 | +0.2158 | 0 | 0 |
| SMARCAD1 | `SMARCAD1 (56916)` | 16 | 0 | +0.0610 | 0.0995 | +0.0587 | +0.0109 | +0.1092 | 0.0983 | -0.1174 | +0.3099 | 0 | 0 |
| ERCC4 | `ERCC4 (2072)` | 16 | 0 | -0.2868 | 0.1883 | -0.2913 | -0.3777 | -0.1538 | 0.2238 | -0.6606 | -0.0093 | 2 | 0 |
| EXD2 | `EXD2 (55218)` | 16 | 0 | +0.1112 | 0.1469 | +0.1429 | +0.0071 | +0.1875 | 0.1804 | -0.1506 | +0.4721 | 0 | 0 |
| BLM | `BLM (641)` | 16 | 0 | -0.0183 | 0.1649 | -0.0150 | -0.1945 | +0.1392 | 0.3337 | -0.2340 | +0.2460 | 0 | 0 |
| FANCM | `FANCM (57697)` | 16 | 0 | -1.0045 | 0.3875 | -1.0202 | -1.1539 | -0.7137 | 0.4403 | -1.9730 | -0.5007 | 16 | 8 |
| POLQ | `POLQ (10721)` | 16 | 0 | -0.3825 | 0.1780 | -0.4089 | -0.5298 | -0.2151 | 0.3147 | -0.6501 | -0.1140 | 5 | 0 |
| PRIMPOL | `PRIMPOL (201973)` | 16 | 0 | +0.1233 | 0.1488 | +0.1323 | +0.0134 | +0.2257 | 0.2123 | -0.2362 | +0.3558 | 0 | 0 |
| REV1 | `REV1 (51455)` | 16 | 0 | +0.0159 | 0.2861 | -0.0341 | -0.1391 | +0.2901 | 0.4292 | -0.5202 | +0.4896 | 1 | 0 |
| REV3L | `REV3L (5980)` | 16 | 0 | -0.6782 | 0.3674 | -0.5777 | -0.9797 | -0.3720 | 0.6077 | -1.2516 | -0.2385 | 9 | 4 |
| POLD1 | `POLD1 (5424)` | 16 | 0 | -2.1912 | 0.4051 | -2.2499 | -2.5383 | -1.8408 | 0.6974 | -2.7758 | -1.3815 | 16 | 16 |
| POLD2 | `POLD2 (5425)` | 16 | 0 | -1.0516 | 0.1983 | -1.0583 | -1.1855 | -0.9326 | 0.2529 | -1.3940 | -0.7306 | 16 | 9 |
| LIG3 | `LIG3 (3980)` | 16 | 0 | -0.1551 | 0.1358 | -0.1432 | -0.3010 | -0.0392 | 0.2618 | -0.3357 | +0.0447 | 0 | 0 |
| SLX4 | `SLX4 (84464)` | 16 | 0 | -0.3279 | 0.2318 | -0.3721 | -0.4830 | -0.1392 | 0.3438 | -0.7545 | +0.1016 | 4 | 0 |
| GEN1 | `GEN1 (348654)` | 16 | 0 | -0.1311 | 0.1882 | -0.0864 | -0.1877 | -0.0284 | 0.1593 | -0.5459 | +0.0973 | 2 | 0 |
| WRN | `WRN (7486)` | 16 | 0 | -0.2693 | 0.2493 | -0.2968 | -0.3878 | -0.1244 | 0.2634 | -0.6308 | +0.2830 | 4 | 0 |
| RAD51AP1 | `RAD51AP1 (10635)` | 16 | 0 | +0.0748 | 0.1275 | +0.0454 | -0.0063 | +0.1711 | 0.1774 | -0.1176 | +0.2811 | 0 | 0 |
| NBN | `NBN (4683)` | 16 | 0 | -0.0571 | 0.1632 | -0.0661 | -0.1724 | +0.0094 | 0.1819 | -0.3283 | +0.2717 | 0 | 0 |
| MRE11 | `MRE11 (4361)` | 16 | 0 | -0.4346 | 0.2584 | -0.4724 | -0.5834 | -0.2257 | 0.3577 | -0.9213 | -0.0184 | 7 | 0 |
| RAD50 | `RAD50 (10111)` | 16 | 0 | -0.1018 | 0.2418 | -0.0338 | -0.2407 | +0.0723 | 0.3130 | -0.7889 | +0.1994 | 1 | 0 |
| ATRX | `ATRX (546)` | 16 | 0 | -0.0594 | 0.3690 | -0.0002 | -0.3731 | +0.1966 | 0.5697 | -0.7342 | +0.4976 | 2 | 0 |
| DAXX | `DAXX (1616)` | 16 | 0 | -0.4343 | 0.2968 | -0.3877 | -0.5309 | -0.2388 | 0.2921 | -1.1625 | -0.0527 | 4 | 1 |
| TERT | `TERT (7015)` | 16 | 0 | +0.0071 | 0.1138 | +0.0019 | -0.0460 | +0.0451 | 0.0910 | -0.1794 | +0.2584 | 0 | 0 |
| SETX | `SETX (23064)` | 16 | 0 | -0.2751 | 0.0910 | -0.2704 | -0.3300 | -0.2238 | 0.1063 | -0.4080 | -0.0969 | 0 | 0 |
| PARP1 | `PARP1 (142)` | 16 | 0 | -0.2970 | 0.1816 | -0.2964 | -0.4050 | -0.1886 | 0.2164 | -0.7186 | -0.0192 | 2 | 0 |
| PARP2 | `PARP2 (10038)` | 16 | 0 | -0.0097 | 0.1183 | -0.0298 | -0.0482 | +0.0546 | 0.1028 | -0.2178 | +0.2847 | 0 | 0 |

---

## 4. ALT− lines (n = 9) — per-gene marginals

Same construction, computed independently on the 9 ALT− lines.

| Gene | Column in `CRISPRGeneEffect.csv` | n | miss | mean | SD | median | q25 | q75 | IQR | min | max | n < −0.5 | n < −1.0 |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| **BIR module (Stage 1)** | | | | | | | | | | | | | |
| POLD3 | `POLD3 (10714)` | 9 | 0 | -1.6913 | 0.3563 | -1.7968 | -1.8589 | -1.4861 | 0.3729 | -2.1272 | -0.9629 | 9 | 8 |
| RAD52 | `RAD52 (5893)` | 9 | 0 | -0.0364 | 0.1745 | -0.1042 | -0.1644 | +0.0544 | 0.2188 | -0.2109 | +0.3217 | 0 | 0 |
| SMARCAL1 | `SMARCAL1 (50485)` | 9 | 0 | +0.0557 | 0.1370 | +0.0660 | +0.0073 | +0.1278 | 0.1206 | -0.2106 | +0.2791 | 0 | 0 |
| **HR comparator (Stage 1)** | | | | | | | | | | | | | |
| BRCA1 | `BRCA1 (672)` | 9 | 0 | -0.5380 | 0.2023 | -0.5334 | -0.7318 | -0.4058 | 0.3260 | -0.8169 | -0.2666 | 5 | 0 |
| BRCA2 | `BRCA2 (675)` | 9 | 0 | -0.5092 | 0.2496 | -0.5498 | -0.7600 | -0.2782 | 0.4818 | -0.7909 | -0.1449 | 5 | 0 |
| RAD51 | `RAD51 (5888)` | 9 | 0 | -1.1340 | 0.4406 | -1.1689 | -1.4747 | -0.9067 | 0.5679 | -1.7429 | -0.3152 | 8 | 6 |
| **Extended BIR / fork restart** | | | | | | | | | | | | | |
| MUS81 | `MUS81 (80198)` | 9 | 0 | -0.2782 | 0.1964 | -0.2707 | -0.3851 | -0.2346 | 0.1505 | -0.5528 | +0.1288 | 1 | 0 |
| PIF1 | `PIF1 (80119)` | 9 | 0 | -0.0892 | 0.1239 | -0.1435 | -0.1463 | -0.0246 | 0.1217 | -0.3042 | +0.0952 | 0 | 0 |
| POLA1 | `POLA1 (5422)` | 9 | 0 | -0.6986 | 0.3241 | -0.8404 | -0.9169 | -0.6039 | 0.3130 | -1.0728 | -0.1037 | 7 | 1 |
| RAD18 | `RAD18 (56852)` | 9 | 0 | -0.0497 | 0.1543 | -0.0470 | -0.1408 | +0.0114 | 0.1522 | -0.3305 | +0.1795 | 0 | 0 |
| SMARCAD1 | `SMARCAD1 (56916)` | 9 | 0 | +0.0356 | 0.1188 | +0.0692 | -0.0239 | +0.0968 | 0.1207 | -0.1514 | +0.2135 | 0 | 0 |
| ERCC4 | `ERCC4 (2072)` | 9 | 0 | -0.2702 | 0.1560 | -0.2421 | -0.3508 | -0.2178 | 0.1330 | -0.5567 | +0.0188 | 1 | 0 |
| EXD2 | `EXD2 (55218)` | 9 | 0 | +0.0614 | 0.1130 | +0.0296 | -0.0481 | +0.1566 | 0.2048 | -0.0636 | +0.2205 | 0 | 0 |
| BLM | `BLM (641)` | 9 | 0 | -0.2396 | 0.2858 | -0.1723 | -0.3467 | -0.0680 | 0.2787 | -0.8643 | +0.1238 | 1 | 0 |
| FANCM | `FANCM (57697)` | 9 | 0 | -0.2498 | 0.2387 | -0.3420 | -0.3982 | -0.1336 | 0.2646 | -0.5132 | +0.2640 | 1 | 0 |
| POLQ | `POLQ (10721)` | 9 | 0 | -0.1888 | 0.1341 | -0.2329 | -0.2655 | -0.0966 | 0.1689 | -0.3583 | -0.0077 | 0 | 0 |
| PRIMPOL | `PRIMPOL (201973)` | 9 | 0 | +0.0415 | 0.1034 | +0.0607 | -0.0141 | +0.1301 | 0.1442 | -0.1229 | +0.1691 | 0 | 0 |
| REV1 | `REV1 (51455)` | 9 | 0 | -0.1751 | 0.2622 | -0.1296 | -0.4263 | +0.0259 | 0.4523 | -0.6119 | +0.1161 | 1 | 0 |
| REV3L | `REV3L (5980)` | 9 | 0 | -0.6020 | 0.1266 | -0.6318 | -0.6971 | -0.5243 | 0.1728 | -0.7441 | -0.3848 | 7 | 0 |
| POLD1 | `POLD1 (5424)` | 9 | 0 | -2.2923 | 0.3928 | -2.2603 | -2.3027 | -2.1403 | 0.1624 | -3.1433 | -1.7650 | 9 | 9 |
| POLD2 | `POLD2 (5425)` | 9 | 0 | -1.1768 | 0.2674 | -1.1830 | -1.3762 | -1.0682 | 0.3080 | -1.4975 | -0.7287 | 9 | 7 |
| LIG3 | `LIG3 (3980)` | 9 | 0 | -0.3615 | 0.2256 | -0.3001 | -0.5754 | -0.2558 | 0.3196 | -0.7227 | -0.0621 | 3 | 0 |
| SLX4 | `SLX4 (84464)` | 9 | 0 | -0.2774 | 0.2058 | -0.2301 | -0.3808 | -0.1746 | 0.2062 | -0.6960 | +0.0014 | 1 | 0 |
| GEN1 | `GEN1 (348654)` | 9 | 0 | -0.1515 | 0.0934 | -0.1135 | -0.2106 | -0.0744 | 0.1362 | -0.2884 | -0.0511 | 0 | 0 |
| WRN | `WRN (7486)` | 9 | 0 | -0.1303 | 0.1454 | -0.0945 | -0.1865 | -0.0401 | 0.1463 | -0.3666 | +0.0844 | 0 | 0 |
| RAD51AP1 | `RAD51AP1 (10635)` | 9 | 0 | +0.0596 | 0.1019 | +0.0574 | +0.0219 | +0.0828 | 0.0609 | -0.0897 | +0.2759 | 0 | 0 |
| NBN | `NBN (4683)` | 9 | 0 | -0.0256 | 0.2007 | -0.0236 | -0.0408 | +0.1124 | 0.1532 | -0.4756 | +0.1999 | 0 | 0 |
| MRE11 | `MRE11 (4361)` | 9 | 0 | -0.7183 | 0.5174 | -0.6060 | -1.0623 | -0.5138 | 0.5485 | -1.4692 | +0.0070 | 7 | 3 |
| RAD50 | `RAD50 (10111)` | 9 | 0 | -0.2643 | 0.4723 | -0.1524 | -0.4313 | +0.0653 | 0.4965 | -1.0170 | +0.3509 | 2 | 1 |
| ATRX | `ATRX (546)` | 9 | 0 | -0.1603 | 0.2252 | -0.1278 | -0.3971 | +0.0574 | 0.4545 | -0.4838 | +0.0881 | 0 | 0 |
| DAXX | `DAXX (1616)` | 9 | 0 | -0.3096 | 0.2334 | -0.3201 | -0.3534 | -0.1734 | 0.1801 | -0.7623 | +0.0426 | 2 | 0 |
| TERT | `TERT (7015)` | 9 | 0 | -0.0306 | 0.1747 | -0.0705 | -0.1631 | -0.0025 | 0.1606 | -0.2494 | +0.2536 | 0 | 0 |
| SETX | `SETX (23064)` | 9 | 0 | -0.2378 | 0.1678 | -0.2372 | -0.2615 | -0.1746 | 0.0869 | -0.5098 | +0.0521 | 1 | 0 |
| PARP1 | `PARP1 (142)` | 9 | 0 | -0.2506 | 0.1751 | -0.2747 | -0.3019 | -0.1624 | 0.1395 | -0.5913 | +0.0513 | 1 | 0 |
| PARP2 | `PARP2 (10038)` | 9 | 0 | -0.0453 | 0.0681 | -0.0553 | -0.0686 | -0.0021 | 0.0665 | -0.1882 | +0.0448 | 0 | 0 |

---

## 5. Missingness and the genome-wide reference pool

### 5.1 Missingness in the 35 candidate genes

**Zero.** Every one of the 35 genes has a non-missing Chronos value in all 16
ALT+ lines and all 9 ALT− lines (`miss = 0` in every row of §3 and §4). There is
no candidate gene whose availability constrains the iteration, and no
imputation question arises for this gene set.

### 5.2 Reference-pool size (counts of *columns*, nothing more granular)

Total gene columns in `CRISPRGeneEffect.csv`: **18,531**.
Columns with **no missing value across the whole line set**:

| Line set | n lines | columns fully non-missing | columns with ≥ 1 missing value |
|---|---:|---:|---:|
| ALT+ | 16 | **17,525** | 1,006 |
| ALT− | 9 | **17,526** | 1,005 |
| ALT+ and ALT− combined | 25 | **17,525** | 1,006 |

So a genome-wide pool of candidate genes for a matched-triplet calibration is
**≈17.5k genes**, i.e. ~94.6% of the matrix, and it is essentially the same
pool in both arms. The three counts pin that down exactly: the union of the two
arms' incomplete columns is 1,006, which equals the ALT+ incomplete count, so
the 1,005 columns incomplete in ALT− are a **subset** of the 1,006 incomplete in
ALT+ — the two arms differ by a single column. Missingness here is a property of
the gene, not of the arm. Pool size is not a constraint: 17,525 complete
columns support ~9 x 10^11 unordered triplets, far more than any calibration
would draw.

### 5.3 Supplementary — pool size after a per-gene variance floor

Added beyond what was asked, because a matched-triplet calibration will likely
want to exclude genes with no usable variance. Each number below is **a count
of columns** for one line set; per-gene SD across lines is a permitted per-gene
statistic, no gene is named, no ModelID appears, and nothing joint across genes
was computed. Drop this subsection if the prespec does not need it.

| Line set | complete | SD > 0.05 | SD > 0.10 | SD > 0.15 | SD > 0.25 | median per-gene SD | median per-gene range |
|---|---:|---:|---:|---:|---:|---:|---:|
| ALT+ (16) | 17,525 | 17,518 | 15,371 | 7,402 | 1,769 | 0.1404 | 0.5211 |
| ALT− (9) | 17,526 | 17,305 | 12,110 | 5,441 | 1,844 | 0.1213 | 0.3815 |
| combined (25) | 17,525 | 17,525 | 15,476 | 6,647 | 1,885 | 0.1361 | 0.5734 |

Read-across for the prespec: the **typical** gene in this cohort has a
per-gene SD of ~0.14 across the 16 ALT+ lines (median range ~0.52). A candidate
gene whose own SD is at or below that median is not distinguishable from the
genome-wide background in terms of spread, so a matched reference pool chosen on
SD will be large at any threshold the prespec is likely to pick (7,402 genes
even at SD > 0.15 in ALT+).

---

## 6. What I did not compute, by instruction

The outcome variable of iteration 03 is the **cross-gene correlation structure**
of Chronos gene effect. Nothing in this file touches it. Specifically, I did
**not**:

1. compute any **correlation, covariance, or any other statistic involving two
   or more genes jointly** — every number in §3 and §4 comes from a single
   1-column Series; the only quantity that ranges over genes at all is a *count
   of columns* in §5.2/§5.3;
2. emit any **per-line table** — no row of this file corresponds to a cell line,
   no ModelID appears anywhere beside the two ground-truth anchor IDs named in
   §1 as assertion inputs (no value is attached to either), and there is **no
   wide matrix and no multi-gene per-line listing**;
3. run **any hypothesis test, permutation or p-value**, of any kind, for any
   gene or pair of genes;
4. compute the **Stage 1 composite** `mean(POLD3,RAD52,SMARCAL1) −
   mean(BRCA1,BRCA2,RAD51)`, or **any** composite, difference, sum, ratio or
   z-score across genes — `src.dependency_stats.composite_score` and
   `standardize_within_line` were never called;
5. **rank or order cell lines** by any gene's value, and I do not name which
   lines are the most or least dependent on anything. (Within-gene order
   statistics — median, q25, q75, min, max — are reported as permitted, but they
   are unlabelled: which line sits at a quantile is not stated and was not
   looked at.)
6. compare ALT+ to ALT− **inferentially**. §3 and §4 are two independent
   descriptions printed side by side; any apparent arm difference in them is
   uncalibrated at n = 16 vs 9 and is **not** a finding of this scan.

Left out because I was unsure it was permitted, per instruction:

- any measure of how a gene's value in one line relates to another gene's value
  in the same line, including scatter summaries, residuals after regressing one
  gene on another, and any pairwise distance between genes;
- any clustering, PCA, factor or network summary over the gene set (all of these
  are functions of the cross-gene covariance, i.e. the outcome);
- any "module score" or within-set average, even descriptive;
- `load_gene_effects(...)` on multi-gene symbol lists for the purpose of
  producing a joint frame to describe (it was never called at all);
- any count of lines satisfying a condition on **two** genes at once.

If the prespec needs any of the above, it must specify it before the analysis
notebook runs — it is not in this file and I did not look at it.

---

## 7. Feasibility verdict

Every gene is **present and fully observed in both arms**, so nothing is
excluded on availability or missingness. The only axis on which a gene can be
unsuitable here is its own **marginal shape**: too little spread to carry a
signal at n = 16, or so uniformly essential that its variation is a
common-essential artefact. The four classes below are assigned from the ALT+
marginals in §3 only, are **non-overlapping**, and are assigned by the priority
saturated → flat → non-essential → usable. Thresholds are stated so the prespec
can override them: "flat" = ALT+ SD < 0.12 (the genome-wide median
per-gene SD in this arm is 0.1404, §5.3);
"saturated" = ALT+ median < −1.0; "non-essential" = 0 of 16 ALT+ lines below −0.5.

### 7.1 Usable, no reservation — 17 genes

Real spread, at least one ALT+ line reaching −0.5, and not pinned to the
essential floor:

`SMARCAL1`, `BRCA1`, `BRCA2`, `MUS81`, `POLA1`, `ERCC4`, `POLQ`, `REV1`, `REV3L`, `SLX4`, `GEN1`, `WRN`, `MRE11`, `RAD50`, `ATRX`, `DAXX`, `PARP1`

Of the six Stage 1 primary genes, three fall here: `SMARCAL1`, `BRCA1`, `BRCA2`.

### 7.2 Usable but **flat** — ALT+ SD below 0.12 (5 genes)

Measured, but line-to-line spread at or below the genome-wide median, so
anything estimated from them at n = 16 will be noise-dominated:

- `SETX` — ALT+ SD 0.0910, median -0.2704, range 0.3111
- `SMARCAD1` — ALT+ SD 0.0995, median +0.0587, range 0.4273; also flat in ALT− (SD 0.1188)
- `TERT` — ALT+ SD 0.1138, median +0.0019, range 0.4379
- `PIF1` — ALT+ SD 0.1182, median -0.0930, range 0.4092
- `PARP2` — ALT+ SD 0.1183, median -0.0298, range 0.5026; also flat in ALT− (SD 0.0681)

`TERT` is in this class and is worth noting for interpretation: in ALT+ lines it
is ALT+ median +0.0019 with SD 0.1138 and 0/16 lines below −0.5 — exactly what a telomerase-independent cohort should look
like, and a useful internal sanity check, but not a usable analysis channel.

### 7.3 Usable but **non-essential throughout** — 0 of 16 ALT+ lines below −0.5 (8 genes)

These have ordinary spread but the spread sits entirely inside the
non-essential band around 0; no ALT+ line shows a dependency on them at all:

- `RAD52` — ALT+ median -0.0100, SD 0.1727, min -0.4162, 0/16 below −0.5
- `RAD18` — ALT+ median -0.0971, SD 0.1775, min -0.4522, 0/16 below −0.5
- `EXD2` — ALT+ median +0.1429, SD 0.1469, min -0.1506, 0/16 below −0.5
- `BLM` — ALT+ median -0.0150, SD 0.1649, min -0.2340, 0/16 below −0.5
- `PRIMPOL` — ALT+ median +0.1323, SD 0.1488, min -0.2362, 0/16 below −0.5
- `LIG3` — ALT+ median -0.1432, SD 0.1358, min -0.3357, 0/16 below −0.5
- `RAD51AP1` — ALT+ median +0.0454, SD 0.1275, min -0.1176, 0/16 below −0.5
- `NBN` — ALT+ median -0.0661, SD 0.1632, min -0.3283, 0/16 below −0.5

**`RAD52` is in this class, and it is a Stage 1 BIR-module member.** ALT+ median -0.0100, SD 0.1727, most negative value -0.4162, **0 of 16** ALT+ lines below −0.5 and 0 below
−1.0. In the ALT− arm it is the same picture (median -0.1042, 0/9 below −0.5). Whatever variance `RAD52` contributes to a module endpoint is
variance within the non-essential band; it is not dependency signal. The prespec
should state this expectation in advance rather than meet it afterwards.

### 7.4 Usable but **near-saturated essential** — ALT+ median < −1.0 (5 genes)

- `POLD1` — ALT+ median -2.2499, SD 0.4051, 16/16 ALT+ lines below −1.0, max -1.3815
- `POLD3` — ALT+ median -1.6664, SD 0.3191, 16/16 ALT+ lines below −1.0, max -1.3121
- `RAD51` — ALT+ median -1.1274, SD 0.3168, 13/16 ALT+ lines below −1.0, max -0.6105
- `POLD2` — ALT+ median -1.0583, SD 0.1983, 9/16 ALT+ lines below −1.0, max -0.7306
- `FANCM` — ALT+ median -1.0202, SD 0.3875, 8/16 ALT+ lines below −1.0, max -0.5007

These are not degenerate — SD 0.20–0.41, wider than the genome-wide median — but
their variation is variation *among already-essential* lines. At −1.0 ≈ the
median common essential, a gene with 16/16 lines past that point is behaving as
a housekeeping dependency in this cohort, and its line-to-line spread is as
plausibly screen quality or proliferation rate as pathway biology.

**`POLD3` is in this class, and it is a Stage 1 BIR-module member**: ALT+ median
-1.6664, 16/16 ALT+ lines below −1.0, least-essential line still -1.3121. `RAD51`, the Stage 1 HR comparator, is in the
same class (median -1.1274, 13/16 below −1.0), as are `POLD1`, `POLD2` and `FANCM` from the extended set.

### 7.5 Genes dropped

**None.** 35/35 resolved to a column, 0/25 missing values, and no gene is
constant: the flattest gene in the ALT+ arm, `SETX`, still has SD 0.0910 over a range
of 0.3111. Classes 7.2–7.4
are **flagged, not dropped** — which of them a prespec admits is a
pre-specification decision, and this scan deliberately does not make it.

### 7.6 Bottom line for the prespec

- Availability imposes **no** constraint on iteration 03: 35/35 genes present,
  zero missing values in either arm, ~17.5k complete genome-wide reference
  columns (17,525 in ALT+), so a matched-triplet calibration has an ample pool.
- The binding constraints are **dynamic range and n**, not missingness. The
  three Stage 1 BIR-module members land in **three different classes**:
  `POLD3` near-saturated essential (§7.4), `RAD52` non-essential throughout
  (§7.3), `SMARCAL1` the only one with an informative mid-range spread
  (§7.1 — ALT+ median -0.4624, SD 0.3950, 8/16 below −0.5). The HR comparator is likewise split:
  `BRCA1` and `BRCA2` mid-range, `RAD51` near-saturated. A prespec that assumes
  these six channels are comparable in scale should say why.
- The ALT− arm is n = 9. Its §4 marginals are descriptive only. An SD computed
  from 9 values has a sampling SE of roughly `sd / sqrt(2(n−1))` ≈ 25% of itself
  (≈18% at n = 16), so an arm-to-arm difference in spread of less than about a
  third is not interpretable from this file, and no calibration of it is
  provided here.
- Whether these genes behave as a coherent **module** is the iteration's outcome
  variable. It is **not answered here and cannot be inferred from this file** —
  §6 lists exactly what was withheld to keep it untouched.

---

## 8. Appendix — reproduction

The scan ran from a throwaway script in the session scratchpad, not from a file
in the repo (the instruction was to write exactly one file). It is short enough
to reproduce verbatim; this is the whole of it, minus printing:

```python
from src import cohort as C, dependency_stats as D, utils as U

coh = C.build_cohort(); C.assert_all(coh)
pos = coh.alt_positive["ModelID"].astype(str).tolist()   # 16
neg = coh.alt_negative["ModelID"].astype(str).tolist()   #  9
U.verify_expected(U.CRISPR_PATH)

colmap, absent = {}, []
for g in ALL_35_SYMBOLS:                 # one symbol at a time, so that an
    try:                                 # absent gene is reported, not fatal
        colmap.update(D.crispr_gene_columns([g]))
    except KeyError:
        absent.append(g)                 # -> [] , all 35 resolved

eff = D.load_all_gene_effects_for_line(pos + neg)        # 25 x 18531

# reference pool: counts of COLUMNS only
notna = eff.notna()
n_pool_pos = int((notna.loc[pos].sum(axis=0) == len(pos)).sum())   # 17525
n_pool_neg = int((notna.loc[neg].sum(axis=0) == len(neg)).sum())   # 17526

# per-gene marginals: ONE column at a time, never a frame of two genes
for g, col in colmap.items():
    for ids in (pos, neg):
        x = eff.loc[ids, col].to_numpy(float)            # 1-D, single gene
        ...  # n, n_missing, mean, sd(ddof=1), median, q25, q75, min, max,
             # (x < -0.5).sum(), (x < -1.0).sum()
```

Determinism: the script was run twice; every number in §3, §4 and §5 was
identical between runs (only the recorded timestamp and git HEAD differed).
`src/` was not modified — no file in the repo was changed by this scan except
this one.

---

## 9. Follow-up — BTR complex members TOP3A, RMI1, RMI2

Added after the literature scan nominated the **BTR complex (BLM, TOP3A, RMI1,
RMI2)** as the candidate assay-sensitivity positive control. Run under the same
prohibitions as the rest of this file: per-gene marginals only, one column at a
time, nothing joint across genes. `BLM` was already covered in §3/§4 and is
restated below for completeness; the Pol δ comparison genes are likewise
restated from §3, not recomputed.

- Run: 2026-10-12, git HEAD `1cb2d870cd62d425453203e3c26e74d5491f3c6a`, same
  loaders and the same cohort objects as §1; `crispr_gene_columns` resolved all
  three symbols (absent: `[]`). Script run twice, numbers identical.

### 9.1 ALT+ lines (n = 16)

| Gene | Column in `CRISPRGeneEffect.csv` | n | miss | mean | SD | median | q25 | q75 | IQR | min | max | n < −0.5 | n < −1.0 |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| TOP3A | `TOP3A (7156)` | 16 | 0 | -1.0505 | 0.2196 | -1.0780 | -1.1881 | -0.8778 | 0.3103 | -1.4713 | -0.6322 | 16 | 10 |
| RMI1 | `RMI1 (80010)` | 16 | 0 | -0.4215 | 0.3703 | -0.4045 | -0.5433 | -0.1392 | 0.4041 | -1.4009 | +0.0906 | 4 | 1 |
| RMI2 | `RMI2 (116028)` | 16 | 0 | -0.1681 | 0.1707 | -0.2145 | -0.2389 | -0.0979 | 0.1410 | -0.4593 | +0.1657 | 0 | 0 |

Restated from §3 for comparison, not recomputed:

| Gene | Column in `CRISPRGeneEffect.csv` | n | miss | mean | SD | median | q25 | q75 | IQR | min | max | n < −0.5 | n < −1.0 |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| BLM | `BLM (641)` | 16 | 0 | -0.0183 | 0.1649 | -0.0150 | -0.1945 | +0.1392 | 0.3337 | -0.2340 | +0.2460 | 0 | 0 |
| POLD1 | `POLD1 (5424)` | 16 | 0 | -2.1912 | 0.4051 | -2.2499 | -2.5383 | -1.8408 | 0.6974 | -2.7758 | -1.3815 | 16 | 16 |
| POLD2 | `POLD2 (5425)` | 16 | 0 | -1.0516 | 0.1983 | -1.0583 | -1.1855 | -0.9326 | 0.2529 | -1.3940 | -0.7306 | 16 | 9 |
| POLD3 | `POLD3 (10714)` | 16 | 0 | -1.7238 | 0.3191 | -1.6664 | -1.9935 | -1.5094 | 0.4841 | -2.4133 | -1.3121 | 16 | 16 |

### 9.2 ALT− lines (n = 9)

| Gene | Column in `CRISPRGeneEffect.csv` | n | miss | mean | SD | median | q25 | q75 | IQR | min | max | n < −0.5 | n < −1.0 |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| TOP3A | `TOP3A (7156)` | 9 | 0 | -1.2887 | 0.3001 | -1.3248 | -1.5605 | -1.1833 | 0.3772 | -1.6404 | -0.7143 | 9 | 8 |
| RMI1 | `RMI1 (80010)` | 9 | 0 | -0.7700 | 0.3019 | -0.8435 | -0.9113 | -0.6087 | 0.3026 | -1.1422 | -0.1422 | 8 | 2 |
| RMI2 | `RMI2 (116028)` | 9 | 0 | -0.2111 | 0.1251 | -0.1850 | -0.2910 | -0.1270 | 0.1640 | -0.3852 | -0.0090 | 0 | 0 |

Restated from §4 for comparison, not recomputed:

| Gene | Column in `CRISPRGeneEffect.csv` | n | miss | mean | SD | median | q25 | q75 | IQR | min | max | n < −0.5 | n < −1.0 |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| BLM | `BLM (641)` | 9 | 0 | -0.2396 | 0.2858 | -0.1723 | -0.3467 | -0.0680 | 0.2787 | -0.8643 | +0.1238 | 1 | 0 |
| POLD1 | `POLD1 (5424)` | 9 | 0 | -2.2923 | 0.3928 | -2.2603 | -2.3027 | -2.1403 | 0.1624 | -3.1433 | -1.7650 | 9 | 9 |
| POLD2 | `POLD2 (5425)` | 9 | 0 | -1.1768 | 0.2674 | -1.1830 | -1.3762 | -1.0682 | 0.3080 | -1.4975 | -0.7287 | 9 | 7 |
| POLD3 | `POLD3 (10714)` | 9 | 0 | -1.6913 | 0.3563 | -1.7968 | -1.8589 | -1.4861 | 0.3729 | -2.1272 | -0.9629 | 9 | 8 |

### 9.3 Class assignment (§7 definitions, ALT+ marginals, same priority)

| Gene | ALT+ median | ALT+ SD | n < −0.5 | Class |
|---|---:|---:|---:|---|
| TOP3A | -1.0780 | 0.2196 | 16/16 | **near-saturated essential (§7.4)** |
| RMI1 | -0.4045 | 0.3703 | 4/16 | **unreserved (§7.1)** |
| RMI2 | -0.2145 | 0.1707 | 0/16 | **non-essential throughout (§7.3)** |
| BLM (from §7) | -0.0150 | 0.1649 | 0/16 | **non-essential throughout (§7.3)** |

The four BTR members therefore fall into **three different classes**: one
near-saturated essential (`TOP3A`), one unreserved (`RMI1`), two non-essential
throughout (`RMI2`, `BLM`).

### 9.4 Marginal adequacy of BTR as a sensitivity control — variance only

**On marginal grounds alone, all four BTR members are non-degenerate, so BTR is
not disqualified — but two of its four members sit in the band where spread is
hardest to distinguish from screen noise.** Every member has an ALT+ SD above
the genome-wide median per-gene SD of 0.1404 (§5.3): `RMI1` 0.3703, `TOP3A`
0.2196, `RMI2` 0.1707, `BLM` 0.1649; none is flat by the §7.2 threshold and
none is constant, so across-line variance exists in all four and a correlation
among them is not ruled out *a priori* by degeneracy. The caution is that the
genome-wide median SD is itself set largely by non-essential genes whose spread
is mostly measurement noise, and `RMI2` (0.1707) and `BLM` (0.1649) sit only
~0.02–0.03 above it with **0 of 16** ALT+ lines reaching −0.5 — i.e. their
variance is variance inside the non-essential band, which is exactly where
noise-dominated spread is expected. `BLM`'s own location (ALT+ median -0.0150)
means it is not a dependency in this cohort at all, whatever its physical
partners do. The Pol δ alternative has the **larger** spread in every member
(`POLD1` 0.4051, `POLD3` 0.3191, `POLD2` 0.1983) and no member near the noise
floor, but all three are near-saturated (medians -2.2499, -1.6664, -1.0583;
§7.4), so their spread is spread *among already-essential* lines and is as
plausibly a global screen-quality or proliferation axis as a complex-specific
one. Neither candidate is clean: **BTR risks a floor of noise in two members,
Pol δ risks a shared screen-quality confound in all three.** Judged only on
whether each member has non-degenerate across-line variance — the only question
asked here — Pol δ is the stronger of the two, and BTR passes with `RMI2` and
`BLM` as its weak links. Which to
pre-specify, and whether a control whose members' variance may be confounded
can serve as a gate at all, is a prespec decision this scan does not make.

**Not computed, by instruction:** no correlation or mean pairwise rho among BTR
members, among Pol δ members, or between the two sets; no composite; no test.
Choosing a positive control by its own correlation result is exactly what this
section is written to avoid.
