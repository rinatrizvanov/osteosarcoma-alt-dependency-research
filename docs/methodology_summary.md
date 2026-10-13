# How this project was run

A technical summary of the pre-registration discipline, the multiplicity correction and the
campaign structure. For the science itself see `project_narrative.md`.

## 1. Pre-registration, and what made it binding

Every analysis in Part 2 was frozen before any outcome value was read. A frozen
pre-specification names: the directional prediction and its sign; the exposure definition
and the exact threshold rule; the endpoint; the test statistic; the minimum detectable
effect at the available n; and — the part usually omitted — explicit **failure** and
**inconclusive** criteria, so that a non-result has a predetermined label rather than a
negotiated one.

Three mechanisms made this more than a promise:

- **Commit ordering.** The prespec is committed before the analysis code, which is committed
  before execution. `git log` is the evidence, and a draft-before-scan commit establishes
  ordering independently of anyone's account of it.
- **Structural isolation.** Exactly one module may load outcome data. Analysis notebooks
  assert it is absent from `sys.modules`. Exposure variables are derived in a separate
  notebook that cannot import it, with thresholds set on the widest *non-outcome* superset
  available (the 34 RNA-backed osteosarcoma lines) and committed before the analysis runs.
- **Mechanical verdicts.** `src/verdict.py` implements an exhaustive partition — `supported`,
  `partial_support`, `inconclusive_underpowered`, `counter_signal`, `refuted`, `null`, plus
  `indeterminate` only for missing inputs. The verdict is computed from the numbers by code.
  There is no step at which a human reads a p-value and chooses a label.

Reproducibility is checked rather than asserted: re-running a notebook must write a
byte-identical `notebook_NN_sha256.json`. Part 2's test suite is 185 tests, including
validation of the bootstrap against the closed form `delta = 2U/(nm) - 1` and of the
Benjamini-Hochberg implementation against `scipy.stats.false_discovery_control`.

## 2. The multiplicity problem, and the correction

The campaign tested five hypotheses in sequence against **one fixed cohort**, and each was
chosen partly *because* the previous one failed. Per-iteration p-values understate the
error rate badly under those conditions. With ten iterations at roughly seven tests each,
three to four nominally significant results would be expected under a complete null.

So the success criterion was set at the campaign level, and stated verbatim in the
repository before any iteration ran:

> A result counts as positive only if its pre-specified success criterion fires, it is
> directionally consistent and robust, **and its primary endpoint survives Benjamini-Hochberg
> FDR at q < 0.1 computed across every iteration's primary endpoint in the campaign** — not
> merely p < 0.05 within its own iteration.

Final values across the three iterations that ran a primary endpoint: **q = 0.420, 0.420,
1.000**. Cumulative test count: **101**.

Worth noting as a property rather than a bug: adding the fourth iteration's primary *raised*
q for the first and third from 0.280 to 0.420. Each additional pre-registered primary makes
the family-wide bar stricter for every member. An analysis that reported only the best
iteration would have shown 0.280 and called it suggestive.

## 3. Three anti-bias rules

- **Campaign-wide FDR** across all primary endpoints (above).
- **Exposure-derivation isolation.** A new exposure variable is derived in its own notebook
  that does not import the outcome loader, thresholded on the widest non-outcome superset,
  and committed before the analysis notebook exists.
- **Contamination disclosure.** When an earlier step had already revealed a value needed by
  a later design, that was declared in a `DISCLOSURES.md` written *before* the next
  iteration was designed, and the affected candidate was withdrawn. This happened twice; in
  one case it cost the campaign its own leading candidate (an HR/BRCA2-centred endpoint),
  because a feasibility scan had already printed the relevant group means.

## 4. Stop conditions, declared in advance

The campaign could end eight ways, all specified before it started — including a clean
positive at q < 0.1, five consecutive non-positive iterations with no new mechanistic angle,
running out of literature-grounded candidates, or needing data the cohort does not contain.
It ended on **running out of candidates**: all thirteen remaining options failed at least one
of four admission filters, and the governing rule was explicit that fabricating a further
hypothesis to extend the loop was not permitted.

## 5. Two frozen rules that fired against the evidence

Recorded because pre-registration is only meaningful if it binds when inconvenient.

- One iteration's assay-sensitivity gate named a **single** positive control, which failed —
  while two other control complexes passed decisively, demonstrating the design was in fact
  sensitive. The frozen rule was applied anyway and the iteration's verdict of record is
  `not-estimable`. The defect was documented as found-by-data, and the identical test was
  deliberately **not** re-run under a corrected label, because the statistic and p-value
  would be unchanged and only the interpretation would differ.
- Another iteration's support conditions were written *relative to the primary's magnitude*.
  The primary returned exactly zero, so three of four failed arithmetically while carrying
  no information — the negative control's own p was 0.540, behaving exactly as intended. The
  outcome memo says so rather than banking "three of four checks failed" as evidence.

Both are written up in `part2_bir_dependency_test/autoresearch/hypothesis_vs_test_postmortem.md`,
alongside a third error of the opposite kind: a hypothesis wrongly declared *not estimable*
on an argument the project's own positive control refutes. That one is corrected by a visible
erratum in `FINAL.md`, not a rewrite. The general lesson — a not-estimable verdict owes a
positive demonstration, not an inference from a classification — is the most portable thing
the campaign produced.

## 6. What this discipline cost and bought

It did not produce a positive result, and it was not supposed to guarantee one. What it
bought is that the nulls are interpretable: each has a pre-registered detection floor, a
calibration control showing the cohort can recover structure of the relevant kind
(BTR complex, T = +0.526, p = 0.0001 at n = 16), and a documented reason why the design could
or could not have seen the effect. The distinction between *the hypothesis is wrong*, *the
test could not see it* and *the construct does not measure what it is named after* is kept
explicit throughout, and the project made errors in both directions on exactly that
distinction before getting it right.
