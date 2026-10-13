# PROGRAM.md — standing rules for the autoresearch loop

> **STATUS: campaign closed 2026-10-12 on stop condition (c), out of literature-grounded candidates. Four iterations ran. Read `FINAL.md` first.** The rules below are unchanged and are kept as the record of what the loop ran under; they are no longer an instruction to continue. Reopening the loop would need a new mandate and, on the evidence in `FINAL.md`, a different cohort.

**This file is the loop's constitution.** A resumed session reads
`PROGRAM.md`, `LEDGER.md`, `NULL_SPACE.md` and `SUMMARY.md` and continues
correctly from those four alone. Do not change the rules below mid-campaign; if
a rule is wrong, stop and report.

## Mission

Given what the literature now says about why Stage 1 returned null, find the
best computational model consistent with this cohort — or characterise the
space as un-modelable from this data. Iterate
(hypothesise → prespec → test → evaluate → revise) until a stop condition fires.

## Context priming — read only these

1. `Literature/INDEX.md` — the literature interface (85 papers, Priority
   Reading, and the **Prespec Conflicts** table: the primary hypothesis source).
2. `docs/analytical_prespec.md` — frozen Stage 1 contract v1.2.
3. `docs/stage1_outcome_memo.md` — what Stage 1 found.
4. `Literature/openalex_scan/next_hypotheses.md` — candidates H1–H5.
5. `README.md` — orientation.

No PDFs at startup. Pull one paper only when the iteration under evaluation
needs its methods or figures; use `pdftotext`, not direct PDF reads; stop
referencing it when the iteration closes.

## Data — this is all there is

`data/raw/depmap/` (DepMap 26Q1, gitignored, sha256-verified on load):
`CRISPRGeneEffect.csv`, `OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv`,
`Model.csv`, `OmicsCNGeneWGS.csv`, `OmicsCNSegmentsWGS.csv`.

No somatic mutations, proteomics, RNAi, PRISM drug response, methylation or SV
calls. **A hypothesis needing data outside that list is skipped with a note in
`NULL_SPACE.md`** — no new release, no substitute proxy.

ALT labels come frozen from `references/master_manifest.csv` (mode 444, 115
rows, 35 `DepMap_26Q1`) and are **never** re-derived.

## Iteration procedure

One hypothesis per iteration, in `autoresearch/iteration_NN_<slug>/` containing
`prespec.md`, `notebooks/`, `data/`, `outcome_memo.md`. Sequential, never
parallel — each informs the next.

1. **State the hypothesis** with its mechanistic basis cited to a specific
   paper indexed in `Literature/INDEX.md`. No hypothesis without literature
   grounding.
2. **Scan** with three parallel subagents (below) to establish feasibility.
3. **Write and commit `prespec.md`** — directional prediction, exposure
   definition, endpoint, test, effect size, success/failure/inconclusive
   criteria via `src.verdict.evaluate_iteration`, a minimum-cell guard, and the
   minimum detectable |delta| from `src.verdict.power_note(n, m)`.
   **Commit before any outcome variable is touched.**
4. **Run** the notebooks; re-run and confirm each `notebook_NN_sha256.json` is
   byte-identical.
5. **Evaluate** mechanically. Never round a result toward a nearby verdict.
6. **Write `outcome_memo.md`**, append to `LEDGER.md`, add to `NULL_SPACE.md`
   if unsupported, update `SUMMARY.md`. Commit.

No post-hoc rule changes inside an iteration. A revised idea is a **new
iteration with its own prespec**, not an amendment.

## The verdict partition

`src/verdict.py` (not `src/criteria.py`, which stays frozen at
`stage1-v1.0-null`). Six named outcomes tile the (direction, p, |delta|) space
exhaustively:

| Outcome | Cell |
|---|---|
| `null` | \|aligned\| < 0.20, either direction, any p |
| `supported` | aligned > 0, p < 0.05, all pre-specified conditions met |
| `partial_support` | aligned > 0, p < 0.05, some condition failed |
| `inconclusive_underpowered` | aligned > 0, p >= 0.05 |
| `refuted` | aligned < 0, p < 0.05 |
| `counter_signal` | aligned < 0, p >= 0.05, \|aligned\| >= 0.20 |

`aligned = delta * predicted_sign`, so `aligned > 0` always means "as
predicted". `indeterminate` is reserved for missing or non-finite inputs and is
never the answer for an in-range result. Stage 1's own result
(delta = +0.281, p = 0.382, predicted negative) is `counter_signal`.

## Three anti-bias rules — the point of the exercise

**(a) Campaign-wide multiplicity.** `LEDGER.md` counts **every** statistical
test run across all iterations. At completion apply BH across all iterations'
primary endpoints (`src.verdict.campaign_fdr`). `SUMMARY.md` reports raw p and
experiment-level q for each. **The formal success criterion is q < 0.1 at the
experiment level, not p < 0.05 at the iteration level.** Any positive is
reported with the number of preceding tests and labelled exploratory.

**(b0) Scoping output must not pair IDs with outcome labels.** Learned in
iteration 01. A feasibility check run for iteration *n+1* can contaminate
iteration *n*: an ATRX scoping table printed for iteration 02 carried the full
ALT label vector for 25 lines, one join away from iteration 01's exposure table.
Write scoping output that pairs ModelIDs with outcome labels to a file instead
of displaying it, or defer it until the iteration that needs it has frozen.
Record any such exposure in the iteration's `DISCLOSURES.md` at the moment it
is noticed.

**(b) Exposure-derivation isolation.** A **new** exposure variable (re-derived
ATRX call, ALT subtype label, branch assignment) is derived in its own numbered
notebook that **does not import `dependency_stats`** or otherwise load
gene-effect values; its threshold is set on the widest non-outcome superset
available (normally the 34 RNA-backed OS lines, not the analysis subset); it is
frozen with its own sha256; and **its commit precedes the analysis notebook's**
so `git log` shows the order. A label fitted on the data that supplies the
outcome is circular and inadmissible.

**(c) Mandatory subtype secondary.** Whenever a composite endpoint is used,
also report the per-subtype split (type I-like vs type II-like, by whatever
marker the iteration can defend) as a required secondary. If no
literature-grounded subtype marker exists for that hypothesis, say so
explicitly as a limitation in the memo.

## Subagents

Parallel, **scan step only**, launched in one message, three per iteration:

- **Essentiality scanner** — gene-effect distributions for candidate genes over
  the relevant lines. Distributions and group sizes only, no hypothesis tests.
- **Metadata scanner** — candidate stratification axes in `Model.csv`, the
  frozen manifest, CN and expression; reports which yield arms >= 4.
- **Literature checker** — is this specific hypothesis supported, contradicted
  or already published? `Literature/INDEX.md` plus at most one or two PDFs.

Written files only, under the iteration folder. No shared state. The main agent
synthesises; subagents do not decide.

## Literature constraints — hard

- **SMARCAL1 is likely an ALT suppressor in osteosarcoma** (Liu 2023
  `10.1093/neuonc/noad022`; Brosnan-Cashman 2021 `10.1038/s41379-021-00841-7`).
  Frame it as a suppressor, not a dependency, absent specific justification in
  that iteration's prespec.
- **WRAP53beta loss impairs HR AND NHEJ non-selectively** (Henriksson 2014
  `10.1101/gad.246546.114`) — upstream at RNF8/ubiquitin signalling, reducing
  53BP1, BRCA1 and RAD51 recruitment together. "Selective HR disabling channels
  cells into BIR" is **not** an available premise.
- **BRCA2 loss drives ALT/BIR** (Lee 2024 `10.1093/nar/gkae251`). **Exclude
  BRCA2 from any HR-neutral control set.**
- **ATRX-gated stratification is the field's working axis** but needs proper
  re-derivation here, subject to the iteration_02 feasibility gate.
- **17p13.1-architecture hypotheses require SV calls DepMap 26Q1 lacks.** The
  expression proxy was already shown insufficient. Skip with a note; invent no
  new proxy.
- **Type I vs Type II ALT may average to null in any single-mode composite**
  (Broderick 2023 `10.1038/s41467-023-38029-z`; Zhang 2026
  `10.1002/bies.70187`; Folini 2026). A composite must select one subtype or
  stratify — see rule (c).

## Stop conditions — stop and report; do not continue

(a) A clean positive: success criterion fires, directionally consistent,
    robust, **and q < 0.1 experiment-wide**.
(b) **Five consecutive** iterations inconclusive or null with no new
    mechanistic angle. Firm: write `FINAL.md` and stop. **Do not fabricate
    further hypotheses to keep the loop alive.**
(c) Out of literature-grounded candidates.
(d) A hypothesis requires collaborator-held data.
(e) A hypothesis requires experimental validation to interpret.
(f) Any single iteration exceeds 3 hours wall-clock.
(g) Two weeks wall-clock without (a).
(h) 10 iterations total.

Whichever comes first. **If a stop condition fires, stop. No "just one more".**

## Hard rules and out of scope

Do not re-run Stage 1 under a different name. Do not fabricate or cherry-pick
essentiality patterns. Do not duplicate published work verbatim without a
specific reason — osteosarcoma-specific replication is justified; a naive
re-test of Bennett's SMARCAL1-across-all-ALT finding is not.

The loop does **not**: write the preprint (manual, separate, sacred); touch the
Stage 1 analysis, its outputs in `data/processed/`, its notebooks or its tag;
propose collaborator data asks; modify the frozen manifest, `references/` or any Aim 2
output. Never `git add Literature/new/` wholesale (503 MB of publisher
supplementary material).

## On completion

`autoresearch/FINAL.md` answers exactly three questions:
**(a)** Did any iteration produce a result that survives experiment-wide FDR?
**(b)** What does the pattern of nulls tell us about the hypothesis space —
which mechanistic framings are ruled out, which remain untestable from this
data?
**(c)** What is the honest next step: preprint with loop results as
supplementary characterisation, pursue a specific hypothesis in a different
cohort, or stop the project entirely?
