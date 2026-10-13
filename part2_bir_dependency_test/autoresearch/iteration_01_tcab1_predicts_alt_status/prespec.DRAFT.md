# Iteration 01 pre-specification — DRAFT, NOT YET FROZEN

> **Status: DRAFT.** Sections marked `[[SCAN]]` are filled from the three scan
> reports before freezing. This draft is committed as a draft so the two-pass
> construction is auditable: everything below that is *not* marked `[[SCAN]]`
> was written before any scan report was read, and the scans report on the
> **exposure and metadata only**, never on the exposure-outcome relationship.
> The frozen version is `prespec.md`; this file is retained as provenance.

## 1. Hypothesis

**H1.** TCAB1/WRAP53 loss predicts ALT **status** — whether an osteosarcoma
line uses Alternative Lengthening of Telomeres at all — rather than predicting
which recombination machinery an already-ALT+ line depends on.

**Directional prediction.** ALT-positive osteosarcoma lines have **lower**
WRAP53 log2(TPM+1) than ALT-negative lines.

With Cliff's delta computed as δ(ALT+ vs ALT−) on WRAP53 expression, the
prediction is **δ < 0**, so `predicted_sign = -1`.

## 2. Mechanistic basis

Keegan et al. (2026) *BMC Cancer*, https://doi.org/10.1186/s12885-026-16250-5
— whole-genome sequencing of osteosarcoma identified structural variants in
TP53 intron 1 that simultaneously inactivate TP53 and TCAB1, exploiting the
head-to-head 1,419-nt TP53/WRAP53 overlap at 17p13.1. These TCAB1 SVs were
present in **~40% of ALT-positive osteosarcoma tumours**, and the paper frames
functional inactivation of the telomerase holoenzyme as "an early and
previously unrecognized event contributing to the **activation** of the ALT
pathway."

`Literature/INDEX.md` Prespec Conflict **#8** records the consequence for
Stage 1: if TCAB1 loss *initiates* ALT, then conditioning on ALT+ status
conditions on the mediator, and a within-ALT+ contrast on TCAB1 is a
collider-adjacent design. H1 is the direct test of that diagnosis — it moves
the contrast to the stratum where the literature actually places the effect.

**This is not Stage 1 under another name.** Stage 1 tested, within ALT+ lines,
whether WRAP53 expression modulates a BIR-minus-HR Chronos dependency
composite. H1 tests, across ALT+ and ALT− lines, whether WRAP53 expression
discriminates ALT status. Different cohort, different contrast, different
outcome variable. **No gene-effect value is read in this iteration at all.**

## 3. Known information — declared, not concealed

This is the one place where H1's pre-registration is imperfect, and the
imperfection is declared rather than hidden.

Stage 1's prespec §4d recorded, **before any outcome was read**, that two lines
have WRAP53 log2(TPM+1) = 0: ACH-002834 (OS384, ALT gold-positive) and
ACH-003185 (NRHOS4, ALT silver-positive). Both are ALT-positive. That is
2 of the ~31 lines in this iteration's cohort, and it is consistent with the
direction H1 predicts.

Consequences, accepted in advance:

1. The test is **partially informed**. 2 of ~31 label–exposure pairs were known
   before this prespec was written.
2. A pre-specified sensitivity analysis **excludes both lines** (§7 R3). If the
   result depends on them, that is reported as the primary interpretation.
3. Marginal distributions of the outcome (how many lines are ALT+ vs ALT−) were
   known and must be, in order to pre-specify group sizes and power — Stage 1's
   prespec stated its 16/9 split the same way. What was **not** observed is the
   *joint* distribution of WRAP53 expression against the ALT label. The exposure
   scan for this iteration was explicitly prohibited from cross-tabulating the
   two, and states which files it did not open.

## 4. Cohort

`[[SCAN]]` — filter trace and arm sizes from `scan_metadata.md`.

Definition, fixed before the scan: all `cohort == "DepMap_26Q1"` rows of the
frozen `references/master_manifest.csv` with `alt_call_direction ∈
{positive, negative}` **and** present in the expression matrix. CRISPR
screening is **not** required, because this iteration reads no gene-effect
data — which expands the cohort beyond Stage 1's 16/9.

ALT labels are taken frozen from the manifest and never re-derived.
Ground-truth anchors asserted before analysis, aborting on failure:
ACH-000364 (U2OS) = positive, ACH-000748 (SJSA1) = negative.

## 5. Exposure

WRAP53 log2(TPM+1) from `OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv`,
direct quantification. The strandedness branch decision ("direct, use WRAP53
TPM") is inherited from Stage 1 notebook 02, together with its documented
limitation that Broad never states library strandedness. No new exposure
variable is derived in this iteration, so anti-bias rule (b) imposes no extra
notebook here.

`[[SCAN]]` — distributional shape and whether a WRAP53-low group is naturally
separable, from `scan_exposure.md`.

## 6. Primary analysis

- **Endpoint:** WRAP53 log2(TPM+1), continuous, per line.
- **Contrast:** ALT-positive vs ALT-negative.
- **Test:** two-sided Mann-Whitney U, exact where no ties, asymptotic with
  `mwu_method` recorded otherwise. One test; no correction within the iteration.
- **Effect size:** Cliff's delta with a 2000-resample bootstrap 95% CI, plus the
  Hodges-Lehmann median shift in log2(TPM+1) units.
- **Discrimination:** AUROC, reported as `(δ + 1) / 2`, because "predicts
  status" is a discrimination claim and AUROC is the quantity a reader will
  want. It is a monotone transform of δ, not a second test.
- **Multiplicity:** the campaign-level BH across all iterations' primaries
  (anti-bias rule (a)) is the binding correction. Experiment-level q < 0.1 is
  the success criterion, not this iteration's raw p.

## 7. Pre-specified secondary and robustness analyses

| # | Check | Pre-specified expectation if H1 is true |
|---|---|---|
| R1 | Gold-only labels (drop `alt_silver`) | same direction |
| R2 | Restrict to the 25 Stage 1 CRISPR-screened lines | same direction; shows the result is not an artefact of cohort expansion |
| R3 | **Exclude ACH-002834 and ACH-003185** (the two known WRAP53-zero lines) | same direction; this is the §3 disclosure guard |
| R4 | WRAP53 gene-level copy number as an orthogonal, strand-agnostic exposure | same direction |
| R5 | TP53 log2(TPM+1) as the exposure | a 17p13.1 co-deletion architecture predicts a same-direction TP53 effect; a WRAP53-specific effect predicts a weaker one. **Descriptive, non-gating** |
| R6 | Leave-one-line-out on the primary | no sign flip across refits |

R5 is deliberately not a gate. The two genes overlap head-to-head, so they are
not independent, and Stage 1 measured ρ(WRAP53, TP53) = +0.354 across all lines.
Reporting it as a gate would penalise H1 for a correlation that the locus
architecture guarantees.

## 8. Criteria

Evaluated mechanically by `src.verdict.evaluate_iteration` with
`predicted_sign = -1`. The six outcomes partition the (direction, p, |δ|) space
exhaustively, so every result reaches exactly one named verdict; no result can
fall into a gap as Stage 1's did.

Support conditions, all of which must pass for `supported` rather than
`partial_support`:

| Key | Condition |
|---|---|
| `primary` | two-sided p < 0.05 **and** δ ≤ −0.330 |
| `r1_gold_only` | R1 same direction (δ < 0) |
| `r3_drop_known` | R3 same direction (δ < 0) — the known-information guard |
| `r6_loo` | no sign flip on any leave-one-line-out refit |
| `discrimination` | AUROC ≥ 0.70 |

Outcome mapping, for the avoidance of doubt:

- `supported` — δ < 0, p < 0.05, every condition above met.
- `partial_support` — δ < 0, p < 0.05, some condition failed.
- `inconclusive_underpowered` — δ < 0 (as predicted), p ≥ 0.05.
- `null` — |δ| < 0.200 either way.
- `counter_signal` — δ > 0 (ALT+ lines have *higher* WRAP53), p ≥ 0.05, δ ≥ 0.200.
- `refuted` — δ > 0, p < 0.05.

A `supported` verdict here is **not** a campaign positive until it also clears
experiment-wide BH at q < 0.1.

## 9. Minimum detectable effect, stated in advance

`[[SCAN]]` — the table is filled with the realised group sizes. Computed by
exhaustive enumeration of the exact Mann-Whitney null
(`src.verdict.min_detectable_delta`), the same routine that independently
reproduces all six values in the frozen Stage 1 table.

Known already, independent of the scan: at 18v13 the smallest |δ| able to reach
two-sided exact p < 0.05 is **0.427**; at Stage 1's 16v9 it is **0.486**; at a
gold-only 12v8 it is **0.542**. So even the most favourable realisable cohort
here can only detect a large effect, and a non-significant result will be
genuinely uninformative about small-to-moderate ones.

## 10. Minimum-cell guard

If either arm has n < 4 after the §4 filters, H1 is declared not estimable, the
reason and the numbers are logged to `NULL_SPACE.md` as `not-estimable`, and the
iteration closes without a dependency read. Inherited from Stage 1
`stratify.MIN_CELL_SIZE = 4`.

## 11. Rule (c) — mandatory subtype secondary

Anti-bias rule (c) binds only when a **composite** endpoint is used. This
iteration's endpoint is a single continuous variable, not a composite, so rule
(c) does not bind.

Recorded separately, because it will matter for later iterations: no defensible
type I vs type II ALT marker exists in DepMap 26Q1. The one paper that directly
profiles the two subtypes (Lee et al. 2026 *Mol Cells*,
https://doi.org/10.1016/j.mocell.2026.100405) could not be retrieved —
`Literature/INDEX.md` lists it as the highest-value gap, supplementary files
only.

## 12. What would make this iteration uninformative

Stated in advance so it cannot be argued afterwards:

1. **Power.** See §9. A null is uninformative about effects below ~0.43.
2. **Range restriction.** Mahmoudi et al. (2011), 10.1038/cddis.2010.90, report
   WRAP53 is typically *over*-expressed in cancer cell lines and that knockdown
   triggers massive apoptosis, with tumour cells more sensitive than normal
   ones. Viable screened lines may therefore be depleted of genuinely
   WRAP53-low states, compressing the exposure range in both arms.
3. **Cell lines are not tumours.** Keegan's ~40% prevalence is a **tumour** WGS
   figure. Nothing guarantees the same prevalence survives into an immortalised
   panel; selection during line establishment runs the other way (point 2).
4. **Expression is not structural disruption.** The central Stage 1 limitation
   carries over unchanged. DepMap 26Q1 ships no SV calls.

## 13. Changelog

| Version | Date | Notes |
|---|---|---|
| draft | 2026-10-12 | Written before the scan reports were read; `[[SCAN]]` sections pending. |
