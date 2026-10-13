# Iteration 03 pre-specification — FROZEN

**Hypothesis H4.** POLD3, RAD52 and SMARCAL1 do not behave as a coherent module in
ALT+ osteosarcoma cells, so no composite over them is a legitimate endpoint.

**Status:** frozen before any cross-gene quantity was computed.
**Supersedes:** `prespec.DRAFT.md` (committed `cf73508`, written before any scan report
existed on disk — the `ls` in that commit's own tool call shows the folder held only
`data/`, `notebooks/` and `src/`).
**Scans informing it:** `scan_essentiality.md`, `scan_metadata.md`, `scan_literature.md`.
**Statistic and detection floor:** `src/coherence.py`, committed with the draft and
validated by `tests/test_coherence_iteration03.py` (58 tests; suite 185).

---

## 1. Hypothesis and directional prediction

**The premise under test is the Stage 1 endpoint's own premise, not its negation.**
Stage 1 scored each line with

```
composite = mean(POLD3, RAD52, SMARCAL1) − mean(BRCA1, BRCA2, RAD51)
```

Averaging three genes into one score presumes they measure one thing. If their
essentiality profiles across lines are mutually uncorrelated, that average is not a
module score — it is three unrelated measurements added together, and Stage 1's endpoint
lacked construct validity regardless of the biology.

> **Prediction.** Across ALT+ osteosarcoma lines, the Chronos gene-effect profiles of
> POLD3, RAD52 and SMARCAL1 are **positively** co-dependent: mean pairwise Spearman rho
> **T > 0**. `predicted_sign = +1`.

Registered in this direction deliberately. If the prediction were "T ≈ 0", a null would
score as a success, and a manufactured null is the single result a loop like this is most
at risk of producing. Here a `null` verdict is a **finding against this project's own
prior endpoint**, and it is the outcome the counter-evidence below predicts.

**Mechanistic basis for the premise.** Batenburg et al. (2023) *Int J Mol Sci*,
https://doi.org/10.3390/ijms241512419, gives the canonical BIR module as
**MUS81–RAD52–POLD3**; Richards et al. (2023) *Nat Commun*,
https://doi.org/10.1038/s41467-023-40071-w, couples RAD52 assembly to PolD3
ADP-ribosylation in break-induced fork recovery.

**Counter-evidence, recorded before the test.** `scan_literature.md` establishes that the
premise has **no co-dependency support anywhere** — no paper in the 87-PDF corpus names
this triplet as one module, and Bennett et al. 2026 (ALTitude,
https://doi.org/10.64898/2026.03.04.709438), the source of our silver ALT labels and of
the SMARCAL1 dependency claim, mentions POLD3, RAD52 and BIR **zero times**. It is
contradicted with a sign:

1. **Broderick et al. (2023)** *Nat Commun*, https://doi.org/10.1038/s41467-023-38029-z —
   ALT runs two distinct BIR processes, and the RAD52-**independent** branch "depends
   specifically on the fork remodeler SMARCAL1 and the MUS81 nuclease". RAD52 and
   SMARCAL1 therefore mark a **switch**, predicting rho(RAD52, SMARCAL1) <= 0.
2. **Wu et al. (2023)** *Nat Commun*, https://doi.org/10.1038/s41467-023-35992-5 — POLD3
   reports general oncogenic replication stress, and both BIR branches share
   PCNA–RFC–Pol delta, so rho(POLD3, ·) ≈ 0.
3. **SMARCAL1 polarity is contested** (Liu 2023 https://doi.org/10.1093/neuonc/noad022;
   Brosnan-Cashman 2021 https://doi.org/10.1038/s41379-021-00841-7; Cox, Maréchal & Flynn
   2016 https://doi.org/10.1016/j.celrep.2016.01.011 runs the opposite way), and germline
   biallelic SMARCAL1 loss is viable in humans, so SMARCAL1 is not constitutively
   essential.

Consolidated as Prespec Conflicts **#1** and **#4** in `Literature/INDEX.md`. The honest
prior is that the premise fails. It is registered in its favour anyway.

**The crux, stated precisely.** Mechanistic coupling within a pathway and correlated
gene-effect across cell lines are **different claims**. Co-essentiality correlation
recovers obligate complexes (Pan et al. 2018, https://doi.org/10.1016/j.cels.2018.04.011;
Wainberg et al. 2021, https://doi.org/10.1038/s41588-021-00840-z; Chen et al. 2023,
https://doi.org/10.1038/s41467-023-39990-5); pathway coupling at a branch point predicts
epistasis or redundancy, not correlation. This iteration tests the *correlation* claim,
because that is the claim a composite endpoint relies on. It does not test, and cannot
refute, the biochemistry.

**Already published?** The method is established (the three citations above). The result is
not: nobody has published the co-dependency structure of this triplet, pan-cancer or in
ALT+ lines. The nearest work, Taglialatela et al. 2026 *Genes Dev*
(https://doi.org/10.1101/gad.353377.125), clusters *cell lines* over top-20 differential
dependencies in DepMap 23Q4 — not *genes* by gene-effect correlation — and POLD3 and RAD52
are not in its gene set. This iteration is osteosarcoma-specific and is a post-mortem on
the loop's own endpoint, so it falls on the permitted side of the no-duplication rule.

**Not Stage 1 renamed, and not iteration 01 renamed.** Stage 1 split 16 ALT+ lines on
WRAP53 expression and compared a composite across arms. This iteration performs **no
split and uses no exposure variable at all** — there is no stratifier to tune and no arm to
assign. Iteration 01 used the ALT label as its outcome; this one uses Chronos gene effect.

## 2. What the scans changed before freezing

Recorded so the result cannot be re-framed afterwards.

**2.1 The three genes occupy three different dependency regimes** (`scan_essentiality.md`
§7, marginals only, no cross-gene quantity computed). Over the 16 ALT+ lines:

| Gene | median | SD | lines < −0.5 | lines < −1.0 | regime |
|---|---|---|---|---|---|
| POLD3 | −1.666 | 0.319 | **16/16** | **16/16** | pan-essential in this cohort |
| RAD52 | −0.010 | 0.173 | **0/16** | 0/16 | essential in no line |
| SMARCAL1 | −0.462 | 0.395 | 8/16 | 3/16 | genuinely variable |

The HR comparator splits the same way: RAD51 median −1.127 (16/16 below −0.5, 13/16 below
−1.0), against BRCA1 −0.555 and BRCA2 −0.406 mid-range.

**This is a construct-validity problem that holds deterministically, before any
correlation or p-value**: the Stage 1 BIR arm averaged a gene essential in every line with
a gene essential in none. §7's variance decomposition quantifies it and does not depend on
the primary test's outcome. It is recorded here, pre-freeze, so that it stands as a
pre-registered observation rather than a post-hoc rationalisation of whatever T turns out
to be.

**2.2 Missingness is zero**, for all 35 scanned genes in all 16 ALT+ and all 9 ALT− lines.
So **n = 16** lines, the §4 minimum-cell guard is satisfied with margin, and the
pre-registered floor of 0.289 applies unchanged. The complete-column reference pool is
**17,525** genes in the ALT+ arm.

**2.3 No gene is degenerate.** The flattest gene scanned (SETX) still has SD 0.091 over a
range of 0.311; the genome-wide median per-gene SD in this arm is 0.1404. RAD52's SD
(0.173) and POLD3's (0.319) both exceed that median, so neither is so flat that
correlation is unestimable. Range restriction will attenuate T, and §5's matched null is
how that attenuation is accounted for rather than argued about.

**2.4 No continuous ALT-activity measure exists in this cohort** (`scan_metadata.md`):
`alt_assay_evidence` is `unknown` for all 35 rows, there is no C-circle, APB, TRAP or
telomere-length column, and `telomerase_status` is constant (`indeterminate`) inside the
ALT+ set. So ALT activity cannot be used as a covariate, and the §7 S5 statement is forced.

**2.5 All 16 ALT+ lines are distinct patients**, so the permutation test faces no
pseudo-replication. (A duplicated-patient pair exists only if the cohort is widened to
the 34 RNA-backed lines, which this iteration does not do.)

**2.6 A free internal check on the frozen labels.** TERT in the ALT+ arm has median
+0.0019, SD 0.114, and **0/16** lines below −0.5 — exactly what a telomerase-independent
cohort should look like. This corroborates the frozen ALT labels from the dependency side,
independently of how they were derived. It is an observation, not a test, and it enters no
criterion.

## 3. Cohort

The **16 ALT+ CRISPR-screened osteosarcoma lines** from `src.cohort.build_cohort()`, with
`assert_all()` run first. ALT labels are taken frozen from `references/master_manifest.csv`
and never re-derived.

Lines missing any analysed gene are dropped before the statistic is computed, with the
count and identity recorded; `mean_pairwise_rho` raises on non-finite input rather than
imputing. **Minimum-cell guard: n >= 10 lines.** Below that the primary is not estimable
and the iteration is logged to `NULL_SPACE.md` as a feasibility failure in the manner of
iteration 02. Per §2.2 the realised n is 16.

## 4. Primary endpoint, test, and the detection floor

**Statistic.** `T` = mean Spearman rho over the three gene pairs (POLD3–RAD52,
POLD3–SMARCAL1, RAD52–SMARCAL1) on Chronos gene effect across the 16 lines
(`src.coherence.mean_pairwise_rho`). Bounded on [−1, +1], the same scale as Cliff's delta,
so it feeds `src.verdict.evaluate_iteration` unmodified and inherits the campaign's
`NULL_CEILING = 0.200`.

**Test.** Two-sided permutation, `N_PERMUTATIONS = 10_000`, `SEED = 20261012`, both fixed
in the module rather than the notebook so neither can be raised after a p-value is seen.
Each gene's values are permuted **independently** across lines, destroying cross-gene
association while leaving every marginal exactly intact. p = `(#{|T_null| >= |T_obs|} + 1) /
(N + 1)`; the +1 stops a finite permutation set ever reporting p = 0.

**One primary test.** No within-iteration correction; the primary p is iteration 03's
single entry in the campaign-wide BH family (§8).

**Detection floor, pre-registered.** `T` is a function of within-gene **ranks**, and the
null permutes those ranks independently, so the null distribution is **distribution-free
given (n_lines, n_genes)** absent ties. It therefore depends on no gene-effect value, and
the floor is genuinely pre-registerable — a property asserted exactly, not assumed, in
`tests/test_coherence_iteration03.py::TestNullIsCentredAndDistributionFree`:

| n lines | k genes | smallest \|T\| reaching two-sided p < 0.05 |
|---|---|---|
| **16** | **3** | **0.289** |
| 16 | 4 | 0.203 |
| 16 | 2 | 0.494 |
| 9 | 3 | 0.384 |

Null \|T\| at 16x3: median 0.103, SD 0.149. Indicative: T = 0.20 → p = 0.181;
T = 0.30 → 0.042; T = 0.40 → 0.005.

**This design is markedly better powered than Stage 1's**, which needed |delta| >= 0.594 at
its 8v8 median split (Stage 1 prespec §5a). The reason is structural, not clever: coherence
uses all 16 lines at once instead of splitting them into two arms of eight.

**The band 0.200 <= |T| < 0.289 is real and is not a null.** It lies above the campaign's
clean-null ceiling but below what n = 16 can resolve, and `src.verdict` returns
`inconclusive_underpowered` there. That is the correct answer in that band and will be
reported as such, not rounded to either neighbour.

**Tie fallback.** If `Coherence.has_ties` is true the distribution-free guarantee lapses,
and the null is recomputed by permuting the **observed** tied rank vectors instead of
`1..n`; the floor above is then superseded by the recomputed one and both are reported.
Chronos values are continuous floats so this is expected to be moot; the rule exists so
the decision is mechanical rather than a judgement made after the fact.

## 5. Calibration against matched random triplets

A raw T near zero could mean "not a module" or merely "these marginals cannot produce
correlation at n = 16". The discriminator, following Wainberg et al. 2021, is a
**mean- and variance-matched** null — an unmatched random null is an artefact generator:

For each of POLD3, RAD52 and SMARCAL1, candidate replacements are drawn from the 17,525
complete-column reference pool with ALT+ mean within **±0.15** and ALT+ SD within **±0.05**
of the target gene. **2,000** matched triplets are sampled (`SEED = 20261012`), T computed
for each, and the observed T placed as a percentile of that distribution. Matching
tolerances and draw count are fixed here, before any triplet is drawn.

If fewer than 50 candidates exist for any gene, the tolerance for that gene is widened in
steps of ±0.05 (mean) and ±0.02 (SD) until 50 are available, and the realised tolerance is
reported. This is a stated fallback, not a free parameter.

## 6. Assay-sensitivity gate — pre-registered, and it can void the primary

A coherence test that cannot detect coherence where coherence certainly exists says nothing
about the BIR triplet. So positive controls run in the **same 16 lines, same statistic,
same null**.

**Primary gate: Pol delta — POLD1, POLD2, POLD3.** An obligate physical complex, and
critically it **contains POLD3 itself**. If POLD3's profile cannot correlate even with its
own obligate complex partners in these 16 lines, then POLD3's non-correlation with RAD52
and SMARCAL1 is uninformative about modularity. That makes this the sharpest available
sensitivity test for this particular primary. Marginals (`scan_essentiality.md` §7.4):
POLD1 median −2.250 SD 0.405, POLD2 −1.058 SD 0.198, POLD3 −1.666 SD 0.319 — all three
near-saturated, all with SD above the genome-wide median.

> **Gate.** The design counts as sensitive only if Pol delta reaches **T_pos >= 0.40 on
> both raw and within-line-standardised gene effect** (§7 S2b). If it fails either, the
> primary is reported as **"design insensitive at this n"**, the verdict is recorded as a
> feasibility failure, and the result is **not** read as evidence against the BIR module.

**Why the gate requires both scales.** Pan-essential genes can correlate across lines
through line-level screen quality rather than shared function — every gene looks more
essential in a line with a stronger screen. A gate passing only on raw values would be
passing through that artefact, which is not the kind of sensitivity the primary needs.
Within-line standardisation (Stage 1 prespec §6.3's transform, z-scoring each line's
effects across all ~18.5k genes) removes the line-level location and scale; genuine
co-essentiality of an obligate complex survives it. The threshold 0.40 is fixed here,
before either value is computed, and sits comfortably above the 0.289 floor so that
passing means real headroom rather than bare significance.

**Secondary control: BTR — BLM, TOP3A, RMI1, RMI2.** Reported for calibration, **with no
gate threshold**, because the follow-up scan shows BTR spans the same regime heterogeneity
as the BIR triplet (TOP3A near-saturated 16/16 below −0.5; RMI1 mid-range; RMI2 and BLM
essential in no line). That makes it a fairer structural analogue of the triplet but a
weaker sensitivity instrument, so it informs interpretation and gates nothing.

**Third module: HR — BRCA1, BRCA2, RAD51,** tested as a *module*, **never as a neutral
comparator**. BRCA2 loss itself drives ALT/BIR (Lee et al. 2024 *NAR*,
https://doi.org/10.1093/nar/gkae251; Prespec Conflict #3), RAD51 is branch-specific, and
BRCA1 suppresses alternative restart routes. The Stage 1 composite's subtrahend is
therefore audited, not trusted.

**Pre-committed interpretive rule** (flagged by `scan_literature.md` as load-bearing): if
**all** positive controls fail, the verdict is `inconclusive_underpowered` — design
insensitive — **not** `refuted` or `null`. At n = 16, "the premise fails" and "the design
is underpowered" are otherwise indistinguishable, and this rule is what keeps them apart.

## 7. Support conditions and secondary analyses

**Support conditions**, passed to `evaluate_iteration` as `conditions`. They separate
`supported` from `partial_support`, are reported in full either way, and by construction
cannot leave a result unassigned:

| Key | Condition |
|---|---|
| `all_pairs_positive` | all three pairwise rho > 0 |
| `min_pair_magnitude` | the weakest pair has rho >= 0.20 |
| `exceeds_matched_triplets` | observed T above the 95th percentile of the §5 matched-triplet distribution |
| `assay_sensitive` | §6 gate: Pol delta T_pos >= 0.40 on raw **and** standardised |
| `loo_no_sign_flip` | no leave-one-line-out refit changes the sign of T |

**Secondaries.** Each reported with its own effect and p where applicable; **none can
change the primary verdict**; all counted in the campaign test total (§8).

- **S1 — per-pair detail.** The three pairwise rho individually, each with a 2000-resample
  bootstrap 95% CI and its own permutation p.
- **S1b — Broderick's directional prediction.** Broderick 2023 places SMARCAL1 in the
  RAD52-independent branch, predicting **rho(RAD52, SMARCAL1) <= 0**. Registered as a
  directional secondary so this iteration can corroborate a specific published claim even
  when its own primary premise fails.
- **S2a — variance decomposition of the Stage 1 composite.** Each of the six Stage 1 genes'
  contribution to the variance of `mean(POLD3,RAD52,SMARCAL1) − mean(BRCA1,BRCA2,RAD51)`
  across the 16 lines. If one gene dominates, Stage 1's result was effectively a
  single-gene result, which bears directly on its R8 check (dropping SMARCAL1 moved delta
  from +0.281 to +0.156). Deterministic, no p-value.
- **S2b — within-line-standardised repeat.** The primary recomputed on gene effects
  z-scored within each line across all genes, per Stage 1 §6.3. Also supplies the gate's
  second scale (§6).
- **S3 — is incoherence ALT-specific?** The same T in the 9 ALT− lines. Floor there is
  0.384, so descriptive only.
- **S4 — the literature's actual BIR module.** Coherence of **MUS81, RAD52, POLD3**
  (Batenburg's module, which names MUS81 and not SMARCAL1). MUS81 marginals are mid-range
  (ALT+ median −0.378, SD 0.153). Descriptive, and registered as secondary precisely so it
  cannot be promoted to primary after the fact.
- **S4b — extended fork-restart set.** Coherence over MUS81, RAD52, POLD3, SMARCAL1 and
  PIF1 (k = 5, floor 0.156). Secondary for the same reason.
- **S5 — subtype secondary: no defensible marker exists, and this is the explicit
  statement rule (c) requires.** Both scans converge. Type I/II ALT is defined in **mouse**
  ESCs via the ~7 kb mTALT template, and Lee et al. 2026's own supplementary material
  (https://doi.org/10.1016/j.mocell.2026.100405) reports that **no compact human genomic
  segment resembling that template was identified in GRCh38**. The three human two-mode
  models (Broderick; Zhang 2026 https://doi.org/10.1002/bies.70187; Folini 2026
  https://doi.org/10.21203/rs.3.rs-11116403/v1) all require activity, imaging or C-circle
  readouts, none of which exist in DepMap 26Q1 (§2.4). Over the 34 RNA-backed lines, no
  candidate marker gene is bimodal — the largest guard-respecting interior gap is below 1
  SD for all 16 tested (gap/SD 0.22–0.90) — so any cut would be a researcher choice, and
  with no ALT assay in the cohort the resulting label would be **unfalsifiable**. The only
  expression/CN-obtainable partition is ATRX/DAXX-driven vs SMARCAL1-driven, which is not
  type I/II and which iteration 02 already found non-estimable. **No invented proxy will be
  substituted**; this limitation is restated in `outcome_memo.md`.
- **S6 — batch robustness.** T refit excluding the 6 commercially-deposited lines
  (6/10 split, the one exogenous batch axis `scan_metadata.md` found usable), to check that
  coherence or its absence is not a depositor artefact.

## 8. Test accounting and verdict mapping

**Tests this iteration contributes to the campaign count: 29.** Enumerated so the total is
auditable rather than asserted: primary (1); Pol delta gate raw and standardised (2); BTR
(1); HR module (1); S1 three pairwise permutation p (3); S2b standardised primary (1);
S3 ALT− (1); S4 (1); S4b (1); S6 (1); leave-one-out refits (16). The matched-triplet
calibration (§5) and the S2a variance decomposition are not hypothesis tests and are not
counted. Bootstrap CIs are interval estimates, not tests, and are not counted.

**Only the primary p enters the campaign BH family.** `LEDGER.md` records all 29 in the
running total, per anti-bias rule (a).

**Verdict.** `evaluate_iteration(delta=T, p_value=p_perm, predicted_sign=+1,
conditions=...)`. No new thresholds and no iteration-local evaluator:

| Verdict | Reading for the Stage 1 endpoint |
|---|---|
| `supported` / `partial_support` | the three genes **do** form a module; the Stage 1 composite was a legitimate construct and its null stands as a biological null |
| `inconclusive_underpowered` | positive but unresolvable at n = 16 — or, per §6, the design was insensitive |
| `null` (\|T\| < 0.200) | **the composite was not a module score.** Stage 1's endpoint lacked construct validity and its null is uninterpretable as a biological result |
| `counter_signal` / `refuted` | the genes are **anti**-correlated — worse than incoherent, since averaging them actively cancels signal. This is what Broderick 2023 predicts |

## 9. What this iteration cannot do

It cannot rescue Stage 1: a coherent module would not make the Stage 1 null a positive, and
an incoherent one does not make the underlying biology wrong — only the measurement of it.

It cannot speak to whether these genes are coupled **mechanistically**. Co-essentiality
across 16 lines is a statement about cross-line variance; a pathway whose members are
uniformly essential, or uniformly inert, in all 16 lines will show no coherence while being
perfectly coupled in biochemistry. §2.1 shows this cohort is exactly that situation for two
of the three genes, which is why §5's matched null and §6's gate exist and why the
mechanistic claim is out of scope. Recorded here so it cannot be forgotten at
interpretation time.

It cannot test ALT subtype (§7 S5), ATRX gating (iteration 02, `NULL_SPACE.md`), or
anything requiring SV calls, proteomics or an ALT-activity assay.
