## 2026-05-01 — Pass B halt override

Pass B fired 7 cohort-dominance warnings, all driven by structural data-availability
mismatch rather than biological confounding:

- Strict/Pragmatic disrupted_high_conf: N=2, both DepMap (insufficient sample size,
  not a confounder).
- Broad disrupted_high_conf: N=3 disrupted + 16 unknown, all DepMap (Aim 1 inclusion
  flags require has_rna=True, which currently collapses the universe to cell lines).
- alt_unknown dominance warnings reflect DepMap's lack of ALT labels (relevant for
  Aim 2 scope, not for confounding within the analyzed set).

Decision: override HALT_DECISION = False to allow Pass A to quantify feasibility.
The cohort dominance is a data-type mismatch, not a within-cohort stratification
issue. Pass A sample-size counts are needed to choose between three project shapes:
  (A) Aim 2 primary, Aim 1 exploratory/supplementary
  (B) Redefine Aim 1 around TCAB1 expression proxy + integrate TARGET-OS
  (C) Pursue de novo SV calling on raw BAMs to fill the architecture gap

Pending Pass A output before commitment.