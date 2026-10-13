# Iteration 03 pre-specification — DRAFT

**Hypothesis H4.** POLD3, RAD52 and SMARCAL1 do not behave as a coherent module in
ALT+ osteosarcoma cells.

**Status: DRAFT.** Written and committed **before any of the three scan reports for this
iteration existed on disk.** No gene-effect value has been read. The draft is committed
separately from the frozen `prespec.md` so that `git log` shows which design decisions
predate the scans and which were made in light of them — the same discipline iteration 01
used (`07220fd` draft, `8e6ed51` freeze).

Three things in this draft are deliberately left open, because they are exactly what the
scans are for: the final control gene sets (§4), the subtype secondary (§7), and the
missingness-driven line count (§3). Everything else below is fixed as written.

---

## 1. Hypothesis and directional prediction

**The premise under test is the Stage 1 endpoint's own premise**, not its negation. Stage 1
scored each line with

```
composite = mean(POLD3, RAD52, SMARCAL1) − mean(BRCA1, BRCA2, RAD51)
```

and averaging three genes into one score presumes they measure one thing. If their
essentiality profiles across lines are mutually uncorrelated, the average is not a module
score, it is three unrelated measurements added together, and the Stage 1 endpoint was not
a legitimate endpoint regardless of the biology.

So the prediction registered here is the **module-premise** direction:

> Across ALT+ osteosarcoma lines, the Chronos gene-effect profiles of POLD3, RAD52 and
> SMARCAL1 are **positively** co-dependent: mean pairwise Spearman rho **T > 0**.

`predicted_sign = +1`. A `null` or `counter_signal` verdict is then evidence **against the
Stage 1 endpoint's construct validity**, which is the finding this iteration exists to
produce. Registering it this way round matters: if the prediction were "T ≈ 0", a null
would be scored as a success, and a confirmed null is the one result a loop like this is
most at risk of manufacturing.

**Mechanistic basis.** Batenburg et al. (2023) *Int J Mol Sci*,
https://doi.org/10.3390/ijms241512419, gives the canonical BIR module as
**MUS81–RAD52–POLD3**; Richards et al. (2023) *Nat Commun*,
https://doi.org/10.1038/s41467-023-40071-w, shows RAD52 and PolD3 are mechanistically
coupled in break-induced fork recovery, with site-specific ADP-ribosylation of PolD3
required for BIR. Those two papers are the strongest available support for treating
POLD3 and RAD52 as one functional unit, and they are the reason the premise is testable
rather than rhetorical.

**Counter-evidence, recorded before the test.** Three strands, all from
`Literature/INDEX.md`:

1. Broderick et al. (2023) *Nat Commun*, https://doi.org/10.1038/s41467-023-38029-z — ALT
   runs **two distinct BIR processes**, with EXD2 directing which. A composite averaging
   across both mixes pathways that may carry opposite dependencies.
2. Carson & Flynn (2023) *Curr Opin Pharmacol*,
   https://doi.org/10.1016/j.coph.2023.102380 — RAD51 and RAD52 regulate **distinct
   branches** of ALT; the BRCA axis is absent from the ALT-vulnerability frame entirely.
3. Wu et al. (2023) *Nat Commun*, https://doi.org/10.1038/s41467-023-35992-5 — POLD3
   essentiality reports a general oncogenic replication-stress requirement, not an
   ALT/BIR-specific one. SMARCAL1, meanwhile, has contested polarity in ALT
   (Liu 2023 https://doi.org/10.1093/neuonc/noad022) and is not constitutively essential —
   germline biallelic loss is viable in humans (Sobieszek & Idzik 2026).

Prespec Conflicts **#1** and **#4** in `Literature/INDEX.md` are the consolidated form of
this. The honest prior is that the premise fails; the test is pre-registered in its
favour anyway.

**Not Stage 1 renamed, and not iteration 01 renamed.** Stage 1 split 16 ALT+ lines on
WRAP53 expression and compared a composite across the two arms. This iteration performs
**no split, uses no exposure variable at all**, and tests the covariance structure
*among* the endpoint's constituent genes. There is no stratifier to tune and no arm to
assign.

## 2. Primary endpoint and statistic

**Statistic.** `T` = mean Spearman rho over the three gene pairs
(POLD3–RAD52, POLD3–SMARCAL1, RAD52–SMARCAL1), computed on Chronos gene effect across the
ALT+ CRISPR-screened osteosarcoma lines. Implemented as
`src.coherence.mean_pairwise_rho`, committed with this draft.

`T` is bounded on [−1, +1], the same scale as Cliff's delta, so it feeds
`src.verdict.evaluate_iteration` unmodified and inherits the campaign's
`NULL_CEILING = 0.200`.

**Test.** Two-sided permutation test, `src.coherence.N_PERMUTATIONS = 10_000`,
`SEED = 20261012`, both fixed in the module rather than the notebook so neither can be
raised after seeing a p-value. Each gene's values are permuted **independently** across
lines, which destroys cross-gene association while leaving every gene's marginal
distribution exactly intact. p is `(#{|T_null| >= |T_obs|} + 1) / (N + 1)`; the +1 keeps a
finite permutation set from ever reporting p = 0.

**One primary test. No multiplicity correction within the iteration**, and the primary
p enters the campaign-wide BH family in `LEDGER.md` as iteration 03's single entry.

## 3. Cohort

The **16 ALT+ CRISPR-screened osteosarcoma lines** from `src.cohort.build_cohort()`, with
`assert_all` run first. ALT labels are taken frozen from `references/master_manifest.csv`
and never re-derived.

Lines with a missing value for **any** of the three genes are dropped before the
statistic is computed, and both the count and the identity of any dropped line are
recorded. `src.coherence.mean_pairwise_rho` raises on non-finite input rather than
imputing. The realised line count goes in the frozen prespec once the essentiality scan
reports missingness; the **minimum-cell guard is `n >= 10` lines**, below which the
primary is declared not estimable and the iteration is logged to `NULL_SPACE.md` as a
feasibility failure in the manner of iteration 02.

## 4. Detection floor — registered before the outcome is read

Because `T` is a function of within-gene **ranks**, and the null permutes those ranks
independently, the null distribution of `T` is **distribution-free given (n_lines,
n_genes)** in the absence of ties. It therefore does not depend on a single gene-effect
value, and the detection floor can be pre-registered rather than reported afterwards.
Computed by `src.coherence.detection_floor` at `SEED = 20261012`:

| n lines | k genes | smallest \|T\| reaching two-sided p < 0.05 |
|---|---|---|
| **16** | **3** | **0.289** |
| 16 | 2 | 0.494 |
| 16 | 5 | 0.156 |
| 9 | 3 | 0.384 |
| 9 | 2 | 0.683 |

Null |T| at 16×3 has median 0.103 and SD 0.149. Indicative p-values: T = 0.20 → 0.181,
T = 0.30 → 0.042, T = 0.40 → 0.005.

**Consequence, stated in advance.** This design is markedly better powered than Stage 1's,
whose 8v8 median split needed |delta| >= 0.594 (prespec §5a). The reason is structural, not
clever: coherence uses all 16 lines at once instead of splitting them into two arms of
eight. A moderate true coherence of T ≈ 0.35–0.40 is detectable here, where the equivalent
effect was invisible to Stage 1.

**The band 0.200 <= |T| < 0.289 is real and is not a null.** It sits above the campaign's
clean-null ceiling but below what this n can resolve, and `src.verdict` will return
`inconclusive_underpowered` there. That is the correct answer in that band and it will be
reported as such, not rounded to either neighbour.

**Tie fallback.** If any gene shows ties among the analysed lines — `Coherence.has_ties`
reports this mechanically — the distribution-free guarantee lapses, and the null is
recomputed by permuting the **observed** tied rank vectors instead of `1..n`. The
pre-registered floor above is then superseded by the recomputed one and both are reported.
Chronos values are continuous floats, so this is expected to be moot; the rule exists so
that the decision is not made by judgement after the fact.

## 5. Assay-sensitivity gate — pre-registered, and it can void the primary

A coherence test that cannot detect coherence where coherence certainly exists tells us
nothing about the BIR triplet. So, in the **same 16 lines, with the same statistic and the
same null**, a positive control is run:

**Positive control: POLD1, POLD2, POLD3** — the catalytic, structural and accessory
subunits of DNA polymerase delta, an obligate physical complex. If co-essentiality across
cell lines detects anything, it detects shared membership of an obligate complex; this is
the standard logic by which DepMap co-dependency analysis recovers protein complexes.

> **Gate.** If the positive control returns `T_pos < 0.40`, the primary result is reported
> as **"design insensitive at this n"** and the iteration's verdict is recorded as a
> feasibility failure, *not* as evidence against the BIR module. The threshold 0.40 is
> fixed here, before either value is computed, and is chosen as comfortably above the
> 0.289 floor so that passing the gate means the design has real headroom rather than
> bare significance.

Note what this gate costs: it is a second test on the same 16 lines, and it is counted in
the campaign test total. It is worth that because without it a null primary is
uninterpretable, which is the failure mode stop condition (e) exists to prevent.

The negative-control set is left open in this draft pending the literature scan, under one
hard constraint already binding: **BRCA2 is excluded from any HR-neutral control set**,
because BRCA2 loss itself drives ALT/BIR (Lee et al. 2024 *NAR*,
https://doi.org/10.1093/nar/gkae251; Prespec Conflict #3).

## 6. Support conditions

Passed to `src.verdict.evaluate_iteration` as `conditions`. They separate `supported` from
`partial_support` and are reported in full either way; by construction they cannot leave a
result unassigned.

| Key | Condition |
|---|---|
| `all_pairs_positive` | all three pairwise rho > 0 |
| `min_pair_magnitude` | the weakest pair has rho >= 0.20 |
| `exceeds_random_triplets` | `T` exceeds the 95th percentile of `T` over random gene triplets drawn from the genome-wide reference pool in the same 16 lines (empirical calibration; pool and draw count fixed in the frozen prespec) |
| `assay_sensitive` | positive control `T_pos >= 0.40` (§5) |
| `loo_no_sign_flip` | no leave-one-line-out refit changes the sign of `T` |

## 7. Secondary analyses — pre-specified

Each is reported with its own effect and p where applicable, and **none can change the
primary verdict**. All are counted in the campaign test total.

- **S1 — per-pair detail.** The three pairwise rho values individually, each with a
  2000-resample bootstrap 95% CI.
- **S2 — the Stage 1 endpoint audit.** Coherence of the HR comparator arm, and the
  variance contribution of each of the six Stage 1 genes to the Stage 1 composite. If one
  gene dominates the composite's variance, the Stage 1 result was effectively a single-gene
  result, which bears directly on its R8 check.
- **S3 — is incoherence ALT-specific?** The same `T` in the 9 ALT− lines. Floor there is
  0.384 (§4), so this is descriptive.
- **S4 — extended BIR set.** Coherence over a literature-nominated extended
  fork-restart set, gene list to be fixed in the frozen prespec from the scans. At k = 5
  the floor falls to 0.156, so this is the better-powered question, and it is registered as
  secondary precisely so it cannot be promoted to primary after the fact.
- **S5 — subtype secondary, required by rule (c).** Standing rule (c) of `PROGRAM.md`
  requires a type I-like vs type II-like split whenever a composite endpoint is used, or an
  explicit statement that no defensible marker exists. Which of those two this iteration
  does is **the one question the frozen prespec will answer from the scans**, and if the
  answer is "no defensible marker", that is stated as a limitation in `outcome_memo.md`
  rather than worked around with an invented proxy.

## 8. Verdict mapping

`src.verdict.evaluate_iteration(delta=T, p_value=p_perm, predicted_sign=+1, conditions=...)`.
No new thresholds, no iteration-local evaluator. In plain terms:

| Verdict | Reading for the Stage 1 endpoint |
|---|---|
| `supported` / `partial_support` | the three genes **do** form a module; the Stage 1 composite was a legitimate construct and its null stands as a biological null |
| `inconclusive_underpowered` | positive but unresolvable at n = 16 |
| `null` (\|T\| < 0.200) | **the composite was not a module score.** Stage 1's endpoint lacked construct validity, and its null is uninterpretable as a biological result |
| `counter_signal` / `refuted` | the genes are **anti**-correlated — worse than incoherent, since averaging them actively cancels signal |

The `null` row is the outcome the counter-evidence in §1 predicts, and it is a substantive
methodological finding about this project's own prior work, not an absence of one.

## 9. What this iteration cannot do

It cannot rescue Stage 1. A coherent module would not make the Stage 1 null a positive, and
an incoherent one does not make the underlying biology wrong — only the measurement of it.
It also cannot speak to whether POLD3, RAD52 and SMARCAL1 are coupled *mechanistically*;
co-essentiality across 16 cell lines is a statement about cross-line variance, and a
pathway whose members are uniformly essential, or uniformly inert, in all 16 lines will
show no coherence while being perfectly coupled in biochemistry. That asymmetry is recorded
here so it cannot be forgotten at interpretation time.
