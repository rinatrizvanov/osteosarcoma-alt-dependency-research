# Iteration 04 disclosures — information already seen, before this prespec is written

Anti-bias rule (b0) requires that any outcome information already visible to this session
be declared before an iteration is designed, so the design cannot quietly exploit it.

## 1. ALT+ vs ALT− marginals for 38 genes are already known

`iteration_03_bir_module_coherence/scan_essentiality.md` §4 reports, **separately for the
16 ALT+ and the 9 ALT− lines**, the mean, SD, median, quartiles, min, max and
threshold counts of Chronos gene effect for 38 genes: the six Stage 1 genes, the extended
fork-restart set, and the BTR members.

**Consequence, binding on this iteration.** Any primary endpoint that is a simple
**ALT+ vs ALT− group difference** on one of those 38 genes is **contaminated** and cannot
be pre-registered as a blind test. I have, for example, already seen that BRCA2 means are
nearly identical between the arms (−0.537 ALT+ vs −0.509 ALT−), which disposes of the
BRCA2-centred contrast that iteration 03's memo listed as a candidate. That candidate is
withdrawn here rather than run as though blind.

The 38 genes: POLD3, RAD52, SMARCAL1, BRCA1, BRCA2, RAD51, MUS81, PIF1, POLA1, RAD18,
SMARCAD1, ERCC4, EXD2, BLM, FANCM, POLQ, PRIMPOL, REV1, REV3L, POLD1, POLD2, LIG3, SLX4,
GEN1, WRN, RAD51AP1, NBN, MRE11, RAD50, ATRX, DAXX, TERT, SETX, PARP1, PARP2, TOP3A,
RMI1, RMI2.

## 2. What is *not* contaminated

The scan was explicitly forbidden from emitting per-line values, from ranking lines, and
from computing any cross-gene quantity, and its §6 records those prohibitions. So:

- **No line-level gene-effect value has been seen** for any gene. Which lines are
  SMARCAL1-dependent is unknown.
- **No exposure–outcome relationship has been seen.** No expression or copy-number
  variable has been related to any gene effect.
- Cross-gene correlations are known only for the sets iteration 03 computed and reported:
  the BIR triplet, Pol delta, BTR, HR, MUS81/RAD52/POLD3, and
  MUS81/RAD52/POLD3/SMARCAL1/PIF1.

An iteration whose primary is a **correlation between an exposure variable and a
line-level gene effect** is therefore uncontaminated, because neither side's per-line
values nor their relationship has been observed.

## 3. Iteration 03's own coherence results are known

T values and pairwise rho for the six module sets above are known (`h4_modules.csv`). Any
iteration reusing one of those exact sets as its primary would be re-reporting a computed
result. Iteration 04 must not do that.
