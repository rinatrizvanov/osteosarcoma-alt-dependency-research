# A pre-registered TCAB1-stratified BIR dependency in ALT+ osteosarcoma: non-confirmation

**Rinat Rizvanov** · 2026-10-11 · Release `stage1-v1.0-null`
Pre-specification: [`analytical_prespec.md`](analytical_prespec.md) v1.2, frozen `81745ed`, amended `a63d56c`
Verdict: `indeterminate` — computed by `src.criteria.evaluate`, not asserted

---

## Summary

A pre-registered test of whether TCAB1/WRAP53 status modulates break-induced replication
(BIR) dependency in ALT+ osteosarcoma cell lines. **It did not confirm.** The primary
estimate ran opposite to the prediction (Cliff's δ = +0.281, p = 0.382) and no gene
approached significance (minimum BH q = 0.959). Eight of ten pre-specified checks agree.
The analysis plan, including the success, inconclusive and failure criteria, was committed
to version control before any gene-effect value was read.

## Background

Alternative Lengthening of Telomeres (ALT) is the telomerase-independent telomere
maintenance mechanism used by a large fraction of osteosarcomas. Keegan et al. (2026)
describe 17p13.1 structural variants disrupting TCAB1 and WRAP53β simultaneously. TCAB1
loss blocks telomerase; WRAP53β loss cripples classical homologous recombination
(Henriksson & Farnebo 2015, *Front Genet* 6:91). An ALT+ cell losing both the telomerase
escape route and the classical HR scaffold would, on this reasoning, be forced onto
non-classical BIR machinery. Bennett et al. (2026, ALTitude) established SMARCAL1 as a
broad ALT dependency but did not stratify by TCAB1 status or 17p13.1 architecture; this
study adds that stratification using public data.

## Hypothesis

Within ALT+ osteosarcoma cell lines, those with **low** TCAB1/WRAP53 expression show
selectively elevated CRISPR essentiality for BIR-module genes (POLD3, RAD52, SMARCAL1)
relative to canonical HR genes (BRCA1, BRCA2, RAD51), compared with ALT+ lines with intact
TCAB1.

Chronos gene effect is more negative when a gene is more essential, so with

```
composite = mean(POLD3, RAD52, SMARCAL1) − mean(BRCA1, BRCA2, RAD51)
```

the prediction is a **more negative** composite in TCAB1-low lines: **Cliff's δ < 0**.

## Methods

**Data.** DepMap 26Q1 (Chronos `CRISPRGeneEffect`, `OmicsExpressionTPMLogp1…`, `Model`,
`OmicsCNGeneWGS`), sha256-verified at every load. 26Q1 specifically, because the ALT
labels were derived against it.

**Cohort.** 35 osteosarcoma lines in `Model.csv` → 26 CRISPR-screened and expressed →
**16 ALT+**, **9 ALT−** (negative control), 1 excluded as `alt_unknown`. ALT labels were
taken frozen from prior work and never re-derived; ground-truth anchors (U2OS positive,
SJSA1 negative) were asserted before analysis.

**Stratification.** TCAB1-low vs TCAB1-intact by WRAP53 log2(TPM+1). The pre-specified
bottom-tertile rule, with its cutpoint estimated on all 34 osteosarcoma lines with RNA so
it could not be tuned to the outcome, placed only 3 of 16 ALT+ lines in the low arm. The
pre-specified minimum-cell guard (<4 per arm) therefore fired and the primary fell back to
a **median split (8 v 8)**, with the superseded tertile (3 v 13) retained as a parallel
sensitivity analysis.

**Measurement check.** WRAP53 exon 1 overlaps TP53 antisense by ~1,419 bp, so an
unstranded library would contaminate the stratifier, and DepMap does not document
strandedness. A pre-committed rule was evaluated empirically: ρ(WRAP53, TP53) = +0.354 sat
at the 94.4th percentile of an expression-matched null, below the committed
99th-percentile flag; within the 34 osteosarcoma lines ρ was only +0.244 (p = 0.164); and
WRAP53 tracked its own **copy number** more strongly (ρ = +0.554, p = 0.0008) than it
tracked TP53 — the opposite of what read cross-assignment produces. Direct quantification
was used.

**Statistics.** Two-sided exact Mann-Whitney on the composite (one pre-specified test, no
correction). Cliff's δ with 2000-resample bootstrap CIs, validated against the identity
δ = 2U/(nm) − 1, plus Hodges-Lehmann shifts. Six per-gene tests with BH-FDR across exactly
those six. Robustness checks R1–R8 and a within-stratum permutation test for selectivity.

**Power, stated in advance.** §5a records that at 8 v 8 the smallest |δ| able to reach
two-sided p < 0.05 is **0.594**, and ≈0.90 at 3 v 13: the study was only ever powered for
a very large effect. Committed before any result was seen, not offered afterwards.

## Results

| Analysis | n | Cliff's δ | 95% CI | p | Predicted direction |
|---|---|---|---|---|---|
| **Primary** (median 8v8) | 8 v 8 | **+0.281** | −0.313 to +0.844 | 0.382 | no |
| Tertile sensitivity | 3 v 13 | +0.538 | −0.231 to +1.000 | 0.189 | no |
| R3 tertile within ALT+ | 5 v 11 | +0.564 | +0.055 to +0.964 | 0.090 | no |
| R4 continuous Spearman | 16 | ρ = −0.442 | — | 0.087 | no |
| R5 ATRX-wildtype | 7 v 8 | +0.179 | −0.500 to +0.786 | 0.613 | no |
| R8 BIR without SMARCAL1 | 8 v 8 | +0.156 | −0.438 to +0.750 | 0.645 | no |
| §6.3 within-line z-scored | 8 v 8 | +0.281 | −0.313 to +0.844 | 0.382 | no |
| R1 gold-only labels | 3 v 7 | −0.619 | −1.000 to +0.143 | 0.183 | yes |
| R6 ALT− control | 4 v 5 | −0.400 | −1.000 to +0.500 | 0.413 | yes |

**Per-gene (primary contrast).** All six BH q = 0.959. POLD3 δ = −0.031, RAD52 −0.344,
SMARCAL1 +0.125, BRCA1 −0.031, BRCA2 −0.250, RAD51 −0.188. Every interval spans zero
(Figure 2).

**Leave-one-out (R7).** None of 16 refits reaches the predicted direction; δ range
[+0.179, +0.464]. The wrong-direction estimate is not driven by any single line.

**Selectivity (§6.2).** The ALT− stratum showed no significant shift (p = 0.413);
interaction permutation p = 0.176.

**The two favourable checks are the two smallest comparisons.** At 3v7 and 4v5 the exact
null median |δ| is 0.33 and 0.30, with P(|δ| ≥ 0.474) of 27% and 29% — both consistent
with noise. R6 is moreover the *negative control*, where the hypothesis predicts no
effect, so a predicted-direction result there argues against a TCAB1-specific account
rather than for one.

## Verdict

`indeterminate`. Three of five §8 support conditions failed. The result is not formally
`falsified` only because |δ| = 0.281 exceeds the 0.20 "clean null" ceiling frozen in §10,
while the sign is wrong and p is far from 0.05. **§11 (negative-result fallback) was
therefore not authorized and no §11 analysis was run.** In plain language: the hypothesis
was not supported.

## What this does and does not establish

**No support.** Nothing in the predicted direction came near significance, and the
better-powered analyses consistently point the other way.

**Not a decisive refutation.** The primary 95% CI (−0.313 to +0.844) contains δ values as
negative as −0.31, so a small-to-moderate effect in the predicted direction is not
excluded.

**But low power does not explain a wrong-direction estimate.** Underpowering explains a
failure to detect. It does not explain why the primary, the tertile sensitivity, the
within-ALT+ tertile, the continuous correlation and all 16 leave-one-out refits agree on
the opposite sign.

## Incidental observation

Two ALT+ lines carry homozygous WRAP53 loss, recorded in prespec §4d *before* the outcome
was read:

| Line | ALT tier | WRAP53 log2(TPM+1) | TP53 log2(TPM+1) | WRAP53 CN |
|---|---|---|---|---|
| ACH-002834 / **OS384** | gold+ | 0.000 | 0.000 | 0.00009 |
| ACH-003185 / **NRHOS4** | silver+ | 0.000 | 0.106 | 0.00001 |

OS384 has zero expression of **both** WRAP53 and TP53 with WRAP53 CN ≈ 0 — the 17p13.1
co-deletion architecture described by Keegan et al., visible in a public cell line. This
stands independently of the dependency result.

Separately, **SMARCAL1 ran backwards** (δ +0.125 at 8v8, +0.436 at 3v13): the WRAP53-null
lines were *less* SMARCAL1-dependent. Dropping SMARCAL1 (R8) moved the composite δ from
+0.281 to +0.156. Against Bennett et al.'s broad ALT SMARCAL1 dependency this is a
counter-observation, but it is post-hoc, underpowered and non-significant.

## Limitations

1. **N = 16 ALT+ lines.** Minimum detectable |δ| 0.594 at 8v8 (§5a).
2. **WRAP53 expression is a proxy for 17p13.1 structural disruption, not a measurement of
   it.** DepMap 26Q1 ships no SV calls. This is the central limitation.
3. **The primary split is the median fallback**, so the low arm contains mid-range-WRAP53
   lines that may not be TCAB1-disrupted, biasing toward the null. The pre-registered §6.4
   check found no concentration in the extreme-loss lines (+0.257, below the 0.300
   threshold) — moot given the sign.
4. **DepMap strandedness is undocumented**; the direct-quantification branch rests on the
   empirical diagnostics above.
5. **ATRX confounding is uncontrolled in practice** — only 1 of 16 ALT+ lines is
   ATRX-altered, so R5 is uninformative. Anticipated in §7, not found afterwards.
6. **6 of 16 ALT+ lines carry prediction-only (silver) ALT labels**; R1 is the guard, and
   it is the one favourable check — at n = 3 in the low arm.
7. **Cell lines are not tumours.**

## Methods lesson, recorded not retrofitted

§10's clean-null ceiling (|δ| < 0.20) was calibrated for a null *in the predicted
direction* and did not anticipate a non-significant estimate of moderate magnitude
pointing the **wrong** way — so the result fell into `indeterminate` rather than a named
outcome. A v2 pre-specification should define the falsification region over the *signed*
effect. **The rule was left exactly as frozen**, not amended after seeing the result it
applies to; the `indeterminate` bucket existed precisely so an unanticipated pattern would
be reported rather than rounded into the nearest verdict.

## Availability

Code, frozen pre-specification, run reports and figures are in the repository at release
`stage1-v1.0-null`; the primary result, per-gene results, robustness tables and the
machine-readable verdict are under `data/processed/`, and both figures under `figures/`
as 300 dpi PNG and PDF. DepMap 26Q1 is not redistributed — `docs/data_provenance.md` names
the five required files with their checksums and the 26Q1-archive fetch policy.
