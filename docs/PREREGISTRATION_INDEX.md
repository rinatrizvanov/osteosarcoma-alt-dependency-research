# Index of pre-registrations, outcome memos and ledgers

These documents are the substance of the project. Everything else is code that executes
them. Each pre-specification was committed **before** the outcome variable it governs was
read; the freeze commit hash is cited inside each file, and the ordering evidence lives in
the parent repository's `git log` (see `project_history.md`).

## Part 2 — the dependency test and campaign

### Stage 1

| Document | What it is |
|---|---|
| [`analytical_prespec.md`](../part2_bir_dependency_test/docs/analytical_prespec.md) | The frozen Stage 1 contract, v1.2. Prediction, exposure rule, endpoint, test, success **and** failure criteria, minimum-cell guard, exact minimum detectable effect. |
| [`stage1_outcome_memo.md`](../part2_bir_dependency_test/docs/stage1_outcome_memo.md) | What Stage 1 found. Two pages. |
| [`stage1_verdict.json`](../part2_bir_dependency_test/data/processed/stage1_verdict.json) | The verdict as computed by code, not chosen. |

### Campaign governance

| Document | What it is |
|---|---|
| [`PROGRAM.md`](../part2_bir_dependency_test/autoresearch/PROGRAM.md) | The campaign's constitution: mission, permitted data, iteration procedure, anti-bias rules, the eight stop conditions. Carries a closed-status banner. |
| [`LEDGER.md`](../part2_bir_dependency_test/autoresearch/LEDGER.md) | Append-only, one row per iteration, with the cumulative test count. Committed rows are never edited. |
| [`SUMMARY.md`](../part2_bir_dependency_test/autoresearch/SUMMARY.md) | Standing state: ranking, campaign FDR, the formal success criterion quoted verbatim. |
| [`NULL_SPACE.md`](../part2_bir_dependency_test/autoresearch/NULL_SPACE.md) | Every hypothesis tried and why it is closed, split into `tested-and-null`, `not-estimable` and `skipped-missing-data`. Includes the binding constraints later iterations had to respect — and a closing erratum withdrawing one entry. |
| [`FINAL.md`](../part2_bir_dependency_test/autoresearch/FINAL.md) | **The conclusion of record.** Answers three questions: did anything survive experiment-wide FDR, what the pattern of nulls means, and the honest next step. Read its erratum section too; it is part of the conclusion. |
| [`candidate_screen_iteration05.md`](../part2_bir_dependency_test/autoresearch/candidate_screen_iteration05.md) | Why there is no fifth iteration. All thirteen candidates, the filter each failed, and the numbers. |
| [`hypothesis_vs_test_postmortem.md`](../part2_bir_dependency_test/autoresearch/hypothesis_vs_test_postmortem.md) | Method post-mortem on the distinction between a wrong hypothesis, a blind test and a broken construct. The most portable document here. |

### Per-iteration

Each directory holds `prespec.md` (frozen), `outcome_memo.md`, `notebooks/` and `data/`.
A `prespec.DRAFT.md` alongside the frozen version is the pre-scan draft, kept to establish
ordering. `DISCLOSURES.md`, where present, declares contamination found before the iteration
was designed.

| # | Directory | Hypothesis | Verdict |
|---|---|---|---|
| 01 | [`iteration_01_tcab1_predicts_alt_status`](../part2_bir_dependency_test/autoresearch/iteration_01_tcab1_predicts_alt_status/) | WRAP53 mRNA predicts ALT *status* rather than dependency within ALT | `counter_signal` |
| 02 | [`iteration_02_atrx_gated_smarcal1`](../part2_bir_dependency_test/autoresearch/iteration_02_atrx_gated_smarcal1/) | SMARCAL1 dependency is ATRX-gated | `not-estimable` — killed at the feasibility gate, no outcome read |
| 03 | [`iteration_03_bir_module_coherence`](../part2_bir_dependency_test/autoresearch/iteration_03_bir_module_coherence/) | POLD3/RAD52/SMARCAL1 do not behave as a module | `not-estimable` by frozen gate / mechanical `null` |
| 04 | [`iteration_04_smarcal1_dependency_model`](../part2_bir_dependency_test/autoresearch/iteration_04_smarcal1_dependency_model/) | PRIMPOL expression predicts SMARCAL1 dependency | `null` |

Iteration 01 and 04 carry `DISCLOSURES.md`. Iteration 02 has an outcome memo but no prespec,
because it never reached a freeze — the feasibility gate stopped it first, which is the
intended behaviour and not an omission.

## Part 1 — label specification and audits

| Document | What it is |
|---|---|
| [`label_spec_v1.0.md`](../part1_alt_classification/docs/label_spec_v1.0.md) | **The frozen label contract.** Category definitions and thresholds for TCAB1 and ALT status. Thresholds may not be adjusted after differential-expression results are seen; changes after the v1.0 freeze require a new version and a written rationale. |
| [`label_spec_v0.9.md`](../part1_alt_classification/docs/label_spec_v0.9.md) | The provisional predecessor, retained to show how the spec evolved before freezing. A `.docx` original sits alongside it. |
| [`feasibility_decision_memo.md`](../part1_alt_classification/docs/feasibility_decision_memo.md) | The sample-size and confounding review that gated the v1.0 freeze, and the decision to adopt the pragmatic definition as primary. |
| [`audit_passA_report.md`](../part1_alt_classification/data/processed/audit_passA_report.md) | Availability audit — what data exists per line. |
| [`audit_passB_report.md`](../part1_alt_classification/data/processed/audit_passB_report.md) | Confounding audit. |
| [`audit_worklog.md`](../part1_alt_classification/docs/audit_worklog.md) | Running log of audit decisions. |
| [`notebook_13_report.md`](../part1_alt_classification/data/processed/notebook_13_report.md) | The robustness analysis: 1000-permutation test, bootstrap coefficient stability, LOOCV fold frequency. The source of the "only FANCM and TERT are stable" finding. |

## Reading order

1. `project_narrative.md` — the arc, no jargon
2. `part2_bir_dependency_test/autoresearch/FINAL.md` — the conclusion, including its erratum
3. `part2_bir_dependency_test/docs/analytical_prespec.md` — what was actually promised in advance
4. `hypothesis_vs_test_postmortem.md` — what went wrong and what transfers
