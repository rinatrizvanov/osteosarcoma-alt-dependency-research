# Iteration 01 — information disclosures before the prespec freeze

Written **before** `prespec.md` was frozen and before the per-line exposure
table in `scan_exposure.md` was read. The point of this file is that leakage is
recorded at the moment it is noticed, not reconstructed afterwards.

## D1 — Stage 1 §4d: two WRAP53-zero lines, both ALT+

Stage 1's prespec §4d recorded, before any Stage 1 outcome was read, that
ACH-002834 (OS384, gold+) and ACH-003185 (NRH-OS4, silver+) have WRAP53
log2(TPM+1) = 0. Both are ALT-positive. **2 known label–exposure pairs.**

Consistent with H1's predicted direction. Guard: R3 drops both lines.

## D2 — Ground-truth anchors are label-bearing, and the scan table is exposure-bearing

The exposure scanner disclosed that `src/cohort.py` and `src/stratify.py`
module docstrings state aggregate cohort constants (ALT+ n = 16, ALT− n = 9)
and the two ground-truth anchors with their directions: ACH-000364 (U2OS) =
positive, ACH-000748 (SJSA1) = negative.

**Ruling on the scanner's conduct: not material.** Marginal outcome counts must
be known in order to pre-specify group sizes and power — Stage 1's prespec
states its own 16/9 split for exactly that reason, and §3 of this prespec
declares it. The anchors are literature-standard reference facts used as
abort-on-failure integrity checks, not data-derived findings. Neither entered
any statistic, table or cutpoint in the deliverable, and neither reveals any
line's WRAP53 value.

**But the combination does.** `scan_exposure.md` contains a per-line table of
WRAP53 values for all 34 lines, keyed by ModelID. Read alongside the anchor
definitions, it yields the (WRAP53, ALT) pair for U2OS and for SJSA1.
**2 further known pairs, so 4 of ~31 in total.** Guard: R3b drops all four.

## D3 — The main agent holds the full ALT label vector for the 25 Stage 1 lines

This is the significant one, and it is the main agent's exposure, not a
subagent's.

Earlier in this session, while scoping the feasibility of a *different*
hypothesis (ATRX-gated SMARCAL1, iteration 02), the main agent printed a
per-line table of ATRX / DAXX / SMARCAL1 expression and copy number for the 25
labelled CRISPR-screened lines, **grouped and sorted by ALT arm**. That output
contains the complete ALT label vector for those 25 lines, by ModelID.

It contains **no WRAP53 column**. The only WRAP53 figure the main agent saw at
that point was an aggregate `describe()` over all 34 lines (min 0.00, median
3.60, max 5.89) — no per-line values.

So the joint distribution of WRAP53 × ALT label is **not** known to the main
agent. But it is now one join away: the label vector is in this session's
history and the per-line exposure table is in `scan_exposure.md`.

**Action taken.** The main agent has deliberately **not read the per-line table
in `scan_exposure.md`**, and will not, until `prespec.md` is committed. Only
these aggregate outputs from that scan were used to write the prespec:

- membership arithmetic (35 → 34 with RNA → 33 with CN, and the two ModelIDs lost)
- distribution summary (min/quartiles/median/max/mean/SD, 2 lines at exactly 0)
- shape: one continuous mode over 2.041–5.888 plus a discrete 2-line spike at 0,
  no low-but-nonzero values, largest interior gap 2.041 units
- candidate cutpoints and their blind group sizes
- rho(WRAP53, WRAP53 CN) = +0.554, p = 8.30e-4; rho(WRAP53, TP53) = +0.244,
  p = 0.164 — both independently reproducing Stage 1's reported values
- the inequality caveat: the tertile cut 3.335580 is an observed value, so
  `<` gives 11/23 and `<=` gives 12/22

None of those reveals a label–exposure pairing.

**Why this is tolerable rather than disqualifying.** What pre-registration
protects is the ordering: the analysis plan must be fixed before the
exposure–outcome relationship is seen. That ordering holds here and is
verifiable in `git log` — `prespec.DRAFT.md` was committed before any scan
report existed, and `prespec.md` is committed before the per-line table is
read. The four known pairs (D1, D2) are a real, bounded imperfection, declared
and guarded by R3/R3b.

**Why it is worth recording anyway.** This is precisely the failure mode that
anti-bias rule (b) exists to prevent, and it arose in the loop's own first
iteration from an incidental scoping check run for a *different* hypothesis.
The operational lesson, carried into `PROGRAM.md` for later iterations: a
feasibility scan for iteration *n+1* can contaminate iteration *n*. Scoping
output that pairs ModelIDs with outcome labels should be written to a file and
not displayed, or deferred until the iteration that needs it has frozen.

## Consequence for the prespec

§3 is amended from "2 of ~31" to "4 of ~31" known label–exposure pairs, and a
second guard is added:

- **R3** — exclude ACH-002834 and ACH-003185 (the D1 pair)
- **R3b** — exclude all four known pairs: ACH-002834, ACH-003185, ACH-000364,
  ACH-000748

If the primary verdict and either guard disagree in direction, the guard is
reported as the primary interpretation.
