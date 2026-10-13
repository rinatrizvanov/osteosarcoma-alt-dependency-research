# Iteration 02 outcome — H2: ATRX-gated SMARCAL1 dependency in ALT+ osteosarcoma

**Verdict: `not-estimable`.** Killed at the pre-specified feasibility gate in
~15 minutes, before any dependency value was read. No prespec was frozen,
because there is no estimable contrast to pre-specify.

**Hypothesis.** Within ALT+ osteosarcoma lines, SMARCAL1 essentiality is
selectively elevated in ATRX-deficient lines, and TCAB1/WRAP53 status adds
nothing once ATRX is accounted for. An osteosarcoma-specific replication of
Brown et al. 2026 *Neuro-Oncology*,
https://doi.org/10.1093/neuonc/noaf300 — justified because replicating a
published effect makes a null interpretable.

**Cohort it would have used.** The 16 ALT+ CRISPR-screened lines (H2's outcome
is Chronos SMARCAL1 gene effect, so CRISPR screening is required — iteration
01's RNA-only expansion is unavailable here).

---

## The gate: no defensible ATRX rule yields both arms ≥ 4

The inherited `atrx_status` is point-mutation-only, so the brief called for
re-derivation including copy-number loss and expression silencing. Four
candidate rules, evaluated on the 16 ALT+ lines:

| Rule | ATRX-altered / 16 | Clears `MIN_CELL_SIZE = 4`? |
|---|---|---|
| Manifest `atrx_status` (point mutations) | **1** (NRHOS1) | no |
| Expression, `ATRX log2(TPM+1) < 3.0` | **2** | no |
| Gene-level CN, `ATRX CN < 0.01` | **4** | yes — **but see below** |
| Segment-level focal homozygous loss (`RelCopyRatio < 0.10`) | **1** (CAL72) | no |

The only rule that clears the guard is the gene-level copy-number one, and it
is wrong.

## Why the gene-level CN rule is unusable

Four lines carry gene-level `ATRX CN ≈ 0` in `OmicsCNGeneWGS.csv`: U2OS, CAL72,
NOS1, OS457. Two checks were run.

**It is not a chromosome-X scaling artefact.** In those same four lines the
flanking chrX genes read normally — KDM6A 0.49–1.21, STAG2 0.37–0.50, RPS4X
0.48–0.70 — so the ATRX calls are focal, not chrX-wide.

**But three of the four are contradicted by the segment file.** Taking ATRX at
hg38 chrX:77,504,878–77,786,233 and reading `OmicsCNSegmentsWGS.csv`:

| Line | gene-level ATRX CN | ATRX-spanning segment `RelCopyRatio` | focal event? |
|---|---|---|---|
| CAL72 | 0.000 | **1e-06** over chrX:77,555,000–77,736,000, flanked by 0.445 / 0.437 | **yes — genuine focal intragenic deletion** |
| NOS1 | 0.001 | 0.406 over a single 92-Mb segment | no |
| OS457 | 0.000 | 0.424 over a single 47-Mb segment | no |
| U2OS | 0.005 | 0.671 over a single 30-Mb segment | no |

For NOS1, OS457 and U2OS the ATRX locus sits inside one broad segment at
0.41–0.67 — roughly single-copy on X, unremarkable in a male line — with no
focal deletion anywhere near the gene. The gene-level "homozygous deletion"
calls are not supported by the data they are derived from.

Across all 16 ALT+ lines, **exactly one (CAL72) has a segment-validated focal
homozygous ATRX loss.**

**The two n=1 methods do not even agree on which line.** The manifest's
point-mutation call is NRHOS1; the segment-validated deletion is CAL72. So four
methods give n ∈ {1, 2, 4}, and the two that give 1 disagree on the identity.

The discordance also explains itself: the four gene-level-CN-zero lines are the
four *highest* ATRX expressers (4.93–6.88, U2OS at 6.88, the cohort maximum).
That pattern is incoherent for a true homozygous deletion and is what a
mis-called locus looks like.

## Power, had the split existed

Even granting the unreliable 4v12 split, the minimum detectable |Cliff's δ| at
4v12 is **0.708** (exact Mann-Whitney null). It would need near-perfect
separation of SMARCAL1 dependency between arms to register at all.

## Three independent reasons H2 is closed on this cohort

1. **Group sizes.** No defensible rule yields both arms ≥ 4. The two credible
   rules give n = 1.
2. **Data quality.** The one rule that reaches n = 4 is refuted for 3 of its 4
   calls by the segment-level file.
3. **Power.** MDE 0.708 even at the unreliable 4v12.

## What this closes, and what it does not

**Closed:** ATRX-gated stratification of SMARCAL1 dependency in DepMap 26Q1
osteosarcoma. The field's working axis (Brown 2026, Liu 2023, Diplas 2018) is
real but **not testable in this lineage on this release** — osteosarcoma ALT is
not ATRX-driven, which is the substance of Stage 1's limitation #5 and is now
quantified rather than asserted.

**Not closed:** Brown et al.'s glioma result stands untouched. A pan-cancer
version using Bennett/ALTitude's expanded ALT+ set across lineages would have
the group sizes — but that needs ALT labels for non-osteosarcoma lines, which
the loop's hard rule forbids deriving ("ALT labels never reclassified — the
frozen manifest is authoritative", and the manifest covers osteosarcoma only).
**That is a scope decision for the user, not for the loop.**

## Data-quality finding worth carrying

**`OmicsCNGeneWGS.csv` gene-level ATRX values disagree with
`OmicsCNSegmentsWGS.csv` in 3 of 4 apparent homozygous deletions.** The naive
"CN loss OR expression silencing" re-derivation the brief specified would have
produced a 4v12 ATRX-gated analysis with three of four exposed lines
misclassified, and nothing downstream would have caught it. Any future
iteration using gene-level copy number from this release must validate against
the segment file first.

This is the third defect of its kind in this project: Stage 1's
`load_wrap53_cn` layout bug, the manifest's artefactual `sv_architecture_auto`
(33/35 rows from one CN segment spanning both WRAP53 and TP53), and now
gene-level ATRX.

## Campaign accounting

Statistical tests on an outcome variable: **0.** No Chronos value was read. The
gate used exposure-side data only (ATRX/KDM6A/STAG2/RPS4X copy number, ATRX
expression, CN segments). Nothing enters campaign BH.
