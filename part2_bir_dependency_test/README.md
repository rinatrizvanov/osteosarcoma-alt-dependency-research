# Pre-registered test of a TCAB1-stratified BIR dependency in ALT+ osteosarcoma

> **Part of a two-project repository.** See the [top-level README](../README.md) for how this
> connects to Part 1, [`docs/project_narrative.md`](../docs/project_narrative.md) for the whole
> arc in plain English, and [`docs/PREREGISTRATION_INDEX.md`](../docs/PREREGISTRATION_INDEX.md)
> for every frozen pre-specification indexed in one place. The ALT labels this part consumes
> were produced by Part 1 and are never re-derived here.
>
> **The campaign that followed Stage 1 is closed.** Its conclusion of record, including a
> published erratum, is [`autoresearch/FINAL.md`](autoresearch/FINAL.md).

**A null result, reported in full.** The hypothesis was pre-registered before any outcome
value was read, the analysis ran as frozen, and it did not confirm.

[![result](https://img.shields.io/badge/result-non--confirmation-B2182B)](docs/stage1_outcome_memo.md)
[![verdict](https://img.shields.io/badge/verdict-indeterminate-666666)](data/processed/stage1_verdict.json)
[![prespec](https://img.shields.io/badge/prespec-v1.2%20frozen-2166AC)](docs/analytical_prespec.md)

---

## Headline result

> **Cliff's δ = +0.281** (95% CI −0.313 to +0.844), **p = 0.382**, n = 8 v 8.
> The prediction was δ **< 0**. The estimate runs the **opposite** way.
> Minimum BH-adjusted q across the six pre-specified genes: **0.959**.

Eight of ten pre-specified checks agree, including all 16 leave-one-out refits. The two
that favour the hypothesis are the two smallest comparisons, and one of them is the
negative control.

| | |
|---|---|
| **Verdict** | `indeterminate` — [`stage1_verdict.json`](data/processed/stage1_verdict.json) |
| **Two-page summary** | [`docs/stage1_outcome_memo.md`](docs/stage1_outcome_memo.md) |
| **Frozen pre-specification** | [`docs/analytical_prespec.md`](docs/analytical_prespec.md) (v1.2) |
| **Figures** | [`figures/`](figures/) |

## Hypothesis

Within ALT+ osteosarcoma cell lines, those with low TCAB1/WRAP53 expression should show
selectively elevated CRISPR essentiality for the **BIR module** (POLD3, RAD52, SMARCAL1)
relative to **canonical HR** (BRCA1, BRCA2, RAD51), compared with ALT+ lines with intact
TCAB1.

**Mechanism.** 17p13.1 structural variants disrupt TCAB1 and WRAP53β together (Keegan et
al. 2026). TCAB1 loss blocks telomerase; WRAP53β loss cripples classical HR (Henriksson &
Farnebo 2015). An ALT+ cell losing both the telomerase escape route and the classical HR
scaffold should be forced onto non-classical BIR machinery. Bennett et al. 2026 (ALTitude)
established SMARCAL1 as a broad ALT dependency *without* stratifying by TCAB1; this adds
that stratification.

**Endpoint.** One composite per line, Chronos gene effect (more negative = more essential):

```
composite = mean(POLD3, RAD52, SMARCAL1) − mean(BRCA1, BRCA2, RAD51)
```

## What was done

DepMap 26Q1. 35 osteosarcoma lines → 26 CRISPR-screened → **16 ALT+ / 9 ALT− / 1 excluded**.
ALT labels taken frozen from prior work ([`references/`](references/)), never re-derived.
TCAB1 stratum from WRAP53 log2(TPM+1); two-sided exact Mann-Whitney; Cliff's δ with
bootstrap CIs and Hodges-Lehmann shifts; BH-FDR across the six genes.

## How the pre-registration was enforced

Not by good intentions — by construction:

- `src/dependency_stats.py` is the only module that loads gene-effect values, and
  `load_gene_effects` is first called in **notebook 05**.
- Notebooks 02–04 carry a **gate check** asserting that module is absent from
  `sys.modules`. All three passed. (Notebook 01 imports it for column names and
  missingness counts only, loading no values — stated precisely in prespec §0.)
- `src/criteria.py` evaluates the success / inconclusive / failure rules **mechanically**,
  so the verdict comes from the rules as written, not from reading the numbers.
- The freeze commit hash is recorded in `data/processed/notebook_04_report.md`.

Two things were fixed *before* the freeze because exact-null calculations exposed them:
the §8.5 selectivity criterion (the original would have blocked a true positive 20–30% of
the time at the real ALT− group sizes), and §5a, which records in advance that at 8v8 the
smallest |δ| reaching p<0.05 is **0.594** — so this study was only ever powered for a very
large effect.

## Reproducing

```bash
uv sync                                      # Python 3.12, pinned
uv run pytest                                # 79 tests
uv run jupyter nbconvert --to notebook --execute --inplace notebooks/0*.ipynb
uv run python scripts/make_publication_figures.py
```

Notebooks are **idempotent**: re-running an unchanged notebook reproduces its
`notebook_NN_sha256.json` byte-for-byte. Run timestamps go to the gitignored
`data/processed/provenance.jsonl`, deliberately not into report bodies, so hash equality
is a real test rather than a formality. Verified end to end on a fresh clone.

Note that git tracks the executable bit but **not** the write bit, so the `chmod -w` on
`references/` does not survive `git clone`. Notebook 00 restores mode 444 when it finds a
reference file writable; the sha256 check against `references/MANIFEST.sha256` is the
hard integrity guarantee.

**Data.** DepMap 26Q1 is not redistributed here (`data/raw/` is gitignored). Fetch the
five files named in [`docs/data_provenance.md`](docs/data_provenance.md) from the DepMap
**26Q1 archive specifically** — not "latest", which would change which lines are screened
and silently break the label join — and place them in `data/raw/depmap/`. Every load is
sha256-verified against `src.utils.EXPECTED_SHA256`.

## Layout

```
docs/analytical_prespec.md     the frozen contract
docs/stage1_outcome_memo.md    two-page citable summary
docs/preprint_outline.md       short-preprint outline
references/                    read-only prior artefacts + frozen ALT labels
src/                           cohort · stratify · strand_check · dependency_stats · criteria · plots
notebooks/                     00–07, idempotent, sha256-verified
scripts/                       publication figure generation
data/processed/                run reports, result JSON, sha256 manifests (committed)
figures/                       publication figures (PNG 300dpi + PDF)
```

## Incidental observation

**OS384** (ACH-002834) shows zero expression of **both WRAP53 and TP53** with WRAP53 copy
number ≈ 0 — the 17p13.1 co-deletion architecture, visible in public data. Recorded in
prespec §4d *before* the outcome was read, so it is not post-hoc. It stands independently
of the dependency result.

## Citation

Rizvanov R. *Pre-registered test of a TCAB1-stratified break-induced replication
dependency in ALT+ osteosarcoma cell lines: a non-confirmation.* 2026.
Release `stage1-v1.0-null`.
