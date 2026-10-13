# Part 1 — ALT classification signature

> **Part of a two-project repository.** See the [top-level README](../README.md) for how this
> connects to Part 2, and [`docs/project_narrative.md`](../docs/project_narrative.md) for the
> whole arc in plain English. This part produced the frozen ALT labels that Part 2 consumes
> read-only.

A curated transcriptomic classifier for alternative lengthening of telomeres (ALT) status in
DepMap osteosarcoma cell lines, with a pre-registered label specification and an honest
coefficient-stability analysis.

## Result

19 of 21 intended curated genes were present in DepMap 26Q1 (`TERC` and `MRE11A` absent).

| Cohort | N | AUROC | 95% CI | AUPRC | Sens@90Spec | Permutation p (1000) |
|---|---|---|---|---|---|---|
| primary | 20 | **0.948** | 0.813–1.000 | 0.976 | 0.917 | < 0.001 |
| pragmatic | 31 | 0.829 | 0.668–0.963 | 0.917 | 0.778 | 0.010 |
| pragmatic, known-sex only | 27 | 0.877 | 0.718–0.984 | 0.937 | 0.824 | 0.003 |

**The coefficients are mostly unstable, and that is the more important finding.** Bootstrap
resampling (200 draws, non-zero threshold 1e-10):

| Cohort | N | stable (>80%) | borderline (50–80%) | unstable (<50%) |
|---|---|---|---|---|
| primary | 20 | 3 genes | 2 | 14 |
| pragmatic | 31 | 3 genes | 5 | 11 |

Only **`FANCM` and `TERT`** are stable in both, both at 100% selection frequency. `TERT`
carries a negative coefficient — lower telomerase, more likely ALT — which is the expected
direction. `SMARCAL1`, which becomes the protagonist of Part 2, is itself unstable here
(16% in primary, 58% in pragmatic).

Read the discrimination and the stability together: at N = 20–31 the model separates the
classes well, while the particular gene weights it uses are largely a function of the sample.
This is a within-cohort result with no external validation and is **not** a validated
biomarker.

## Pre-registration

The label specification was frozen before modelling, with an explicit versioning rule:
thresholds may not be adjusted after differential-expression results are seen, and any change
after the v1.0 freeze requires a new version number and a written rationale.

- [`docs/label_spec_v1.0.md`](docs/label_spec_v1.0.md) — the frozen contract
- [`docs/label_spec_v0.9.md`](docs/label_spec_v0.9.md) — provisional predecessor, kept for provenance
- [`docs/feasibility_decision_memo.md`](docs/feasibility_decision_memo.md) — the sample-size and confounding review that gated the freeze

## Layout

| Path | Contents |
|---|---|
| `notebooks/00`–`03` | Data download, manifest construction, SV review, ALT label merge |
| `notebooks/04`–`05` | Confounding and availability audits |
| `notebooks/06`–`10` | RNA inventories, recovery routes, ALTitude/DepMap merge |
| `notebooks/11`–`13` | The signature, sex-adjusted sensitivity, curated robustness |
| `src/` | Label rules, manifest builder, signature fitting, audit helpers |
| `data/manifests/master_manifest.csv` | **The frozen label output.** Consumed read-only by Part 2. |
| `data/processed/` | Notebook reports, predictions, fitted models, sha256 manifests |
| `docs/project_slides.*`, `docs/project_summary_one_pager.md` | Presentation material |

## Running it

Raw DepMap data is not included. See [How to reproduce](../README.md#how-to-reproduce) in the
top-level README for the required release and exact SHA256 checksums, then:

```sh
uv sync
uv run pytest
uv run jupyter lab
```

Note that one gene pair in the primary signature (`USP9Y`, `DDX3Y`) is Y-linked, which is why
`notebooks/12` repeats the analysis with sex handled explicitly — as a covariate, and with
Y-chromosome genes dropped.
