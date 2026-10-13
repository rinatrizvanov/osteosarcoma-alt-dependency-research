# Osteosarcoma ALT: classification, and a pre-registered dependency test

Two linked computational projects on **alternative lengthening of telomeres (ALT)** in
osteosarcoma, built on the public DepMap 26Q1 cell-line resource. The first builds a
transcriptomic classifier for ALT status. The second uses its frozen output to test a
specific mechanistic hypothesis about DNA-repair dependency — and does not confirm it.

Both were run under pre-registration: the analysis plan, the success criteria **and the
failure criteria** were committed to version control before any outcome value was read. The
results are a strong classifier with unstable internals, and a well-characterised null.

---

## Part 1 — ALT classification signature

`part1_alt_classification/`

A 19-gene curated ALT signature (19 of 21 intended genes present in DepMap 26Q1; `TERC`
and `MRE11A` absent) trained on osteosarcoma expression data against manually curated ALT
labels.

| Cohort | N | AUROC | 95% CI | Permutation p |
|---|---|---|---|---|
| primary label definition | 20 | **0.948** | 0.813–1.000 | < 0.001 |
| pragmatic definition | 31 | 0.829 | 0.668–0.963 | 0.010 |

Discrimination is strong. **The coefficients are not.** Bootstrap resampling (200 draws)
found only **3 of 19 genes stable** at the primary N = 20, with 14 below 50% selection
frequency; just **`FANCM` and `TERT` were stable in both cohorts**. `TERT`'s coefficient is
negative, the expected direction for a telomerase-negative phenotype. So the model separates
classes well while the genes it leans on at this sample size are mostly not reproducible —
reported in full rather than omitted.

## Part 2 — Pre-registered BIR dependency test

`part2_bir_dependency_test/`

**Hypothesis.** Within ALT+ osteosarcoma lines, those with low TCAB1/WRAP53 expression
should depend more on the break-induced-replication module (POLD3, RAD52, SMARCAL1) than on
canonical homologous recombination (BRCA1, BRCA2, RAD51). The mechanism: 17p13.1 structural
variants disrupt TCAB1 and WRAP53β together, so a cell losing both the telomerase escape
route and the HR scaffold should be pushed onto BIR.

**Result.** Cliff's δ = **+0.281** (95% CI −0.313 to +0.844), p = 0.382, median split 8 v 8.
The prediction was δ < 0. The estimate runs the opposite way, and the sign is stable across
all 16 leave-one-out refits, both tertile variants and the continuous version. Minimum
BH-adjusted q across the six pre-specified genes: 0.959.

**Follow-up campaign.** Four further pre-registered iterations tested alternative models
from the recent literature under Benjamini–Hochberg correction across every iteration's
primary endpoint — the bar was q < 0.1 at the campaign level, not p < 0.05 within an
iteration. Final q-values: **0.420, 0.420, 1.000**. No positive. It closed when every
remaining candidate required a variable the cohort does not contain.

**The most useful finding is a construct post-mortem.** Two of the three BIR genes carry no
dependency variation in these lines — POLD3 is essential in 16/16, RAD52 in 0/16 — and the
HR comparator arm supplied **68.2%** of the composite endpoint's variance while RAD52
contributed **−2.7%**. The endpoint was mostly measuring its own comparator, and the arm it
was meant to contrast against is not a coherent module here. Separately, a positive control
(the BTR complex, T = +0.526, p = 0.0001) confirms the cohort *can* recover module structure
at n = 16, so these nulls are informative within their stated detection floors rather than
simply underpowered.

## How the two parts connect

Part 1 produced `master_manifest.csv`, the frozen ALT label set. Part 2 consumed it
read-only and **never re-derived a label** — it is vendored into
`part2_bir_dependency_test/references/` with sha256 pins, which is why several Part 1 files
appear twice here. The duplication is deliberate: it is the provenance record, and
deduplicating it would break the checksums the design depends on.

---

## What this is NOT

- **This is not a refutation of the biological hypothesis.** No claim is made that
  TCAB1-stratified BIR dependency does not exist. The claim is narrower and entirely about
  measurability: the prediction, as operationalised on this cohort, could not have been
  detected — the exposure had no natural low group, and the endpoint averaged genes in three
  incompatible dependency regimes. A non-confirmation, not a falsification.
- **It is not a claim that the hypothesis space is exhausted.** One candidate (cross-arm
  antagonism between HR and alternative restart routes) remains genuinely **untested** on
  data already here. An earlier draft of the campaign's conclusion wrongly closed it; the
  correction is recorded as a visible erratum in
  `part2_bir_dependency_test/autoresearch/FINAL.md` rather than quietly rewritten.
- **It is not a validated biomarker.** Part 1's classifier is a within-cohort result at
  N = 20–31 with mostly unstable coefficients and no external validation.
- **It is not a clinical or diagnostic tool.**
- **It does not redistribute DepMap or any publisher material.** Code, pre-registrations and
  derived summaries only.

## How to reproduce

Requires Python 3.12 and [`uv`](https://docs.astral.sh/uv/). Raw data is **not** included
and must be fetched from [depmap.org](https://depmap.org/portal/download/) — release
**26Q1 specifically**, because the ALT labels were derived against it and a different
release silently changes which lines are screened.

Place the five files in `part2_bir_dependency_test/data/raw/depmap/` and verify:

```
e610a4cefb13a82b5b256b47eb08b63ff14843f8dbd0fb164bc0a32688e5b89e  CRISPRGeneEffect.csv                                 440646050 bytes
0377be80c525fde98cbd2c6e8b06bdf2a4014a9683eb70182c1f8649d711021a  OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv   305007605 bytes
ea4e0b2a3bc806f81df62689a5ae75f1a100135727a3d7b8a4c7ccc8815183f8  Model.csv                                               697455 bytes
d2fe65b516802fb2b355b953fcc43dcb908840f275e9166c0dbd7defadca8a1d  OmicsCNGeneWGS.csv                                   409005582 bytes
690e861b3a35db89b6a6b1122edd842c1adb090ec52ec1b13b0ee91d369a142f  OmicsCNSegmentsWGS.csv                                63718976 bytes
```

Check with `shasum -a 256 <file>` (macOS) or `sha256sum <file>` (Linux). The frozen label
input is pinned the same way:
`d51c5f180391ffbce9fd4d5d66ed449532a796859e2bd186a73049c8cbc68b31  master_manifest.csv`.

Then `uv sync && uv run pytest` in either part. Part 2's suite is 185 tests. Notebooks are
idempotent: re-running writes a byte-identical `notebook_NN_sha256.json`.

## Repository map

| Path | What it is |
|---|---|
| `docs/project_narrative.md` | The whole arc in plain English — start here |
| `docs/methodology_summary.md` | How the project was run: pre-registration, FDR, campaign structure |
| `docs/PREREGISTRATION_INDEX.md` | Every frozen pre-specification, outcome memo and ledger, indexed |
| `docs/project_history.md` | How this repository was consolidated from two parent projects |
| `part1_alt_classification/` | ALT signature: notebooks, label spec, audits |
| `part2_bir_dependency_test/` | Stage 1 test, the four-iteration campaign, literature index |

## Citation and licence

MIT licence (`LICENSE`). Citation metadata in `CITATION.cff`.
