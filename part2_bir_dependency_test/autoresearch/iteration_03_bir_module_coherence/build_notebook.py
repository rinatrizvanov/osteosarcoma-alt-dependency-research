"""Builds notebooks/01_h4_module_coherence.ipynb from the cell sources below."""
import json, pathlib

CELLS = []
def md(s): CELLS.append(("markdown", s.strip("\n")))
def co(s): CELLS.append(("code", s.strip("\n")))

md("""
# 01 — H4: do POLD3, RAD52 and SMARCAL1 behave as a module?

Implements `prespec.md`, frozen at commit `b24a8ed`, before any cross-gene quantity was
computed. The statistic and its detection floor were committed earlier still, with the
draft (`cf73508`), in a module that contains no path to the outcome.

**The premise under test is the Stage 1 endpoint's own premise**, so the prediction is
`T > 0`, `predicted_sign = +1`. A `null` verdict here is a finding **against** this
project's prior endpoint, not a success.

Pre-registered detection floor at 16 lines x 3 genes: **|T| = 0.289**, against Stage 1's
0.594 at its 8v8 split. The band 0.200 <= |T| < 0.289 resolves to
`inconclusive_underpowered` and will not be rounded to either neighbour.

Assay-sensitivity gate (prespec §6): Pol delta must reach **T >= 0.40 on both raw and
within-line-standardised** effect, or the primary is reported as design-insensitive rather
than as evidence against the module.
""")

md("## 1 — Header, cohort and gene-set resolution")
co("""
import sys, random, warnings
from pathlib import Path

import numpy as np

random.seed(42)
np.random.seed(42)

_cwd = Path('.').resolve()
ROOT = next((p for p in [_cwd, *_cwd.parents] if (p / 'src').is_dir() and (p / 'references').is_dir()), None)
assert ROOT is not None, f'Cannot locate repo root from {_cwd}'
ITER = ROOT / 'autoresearch' / 'iteration_03_bir_module_coherence'
for p in (str(ROOT), str(ITER)):
    if p not in sys.path:
        sys.path.insert(0, p)
print(f'ROOT: {ROOT}')

warnings.filterwarnings('ignore', category=FutureWarning)

import pandas as pd
from src.utils import (log_versions, run_timestamp, record_provenance, write_report,
                       write_sha256_manifest, git_head)

RUN_TIMESTAMP = run_timestamp()
VERSIONS = log_versions(['pandas', 'numpy', 'scipy'])
assert sys.version_info[:2] == (3, 12), f'Expected Python 3.12, got {sys.version_info[:2]}'
print(f'Python: {sys.version.split()[0]}')

PRESPEC_FREEZE = 'b24a8ed'   # prespec.md frozen
STATISTIC_COMMIT = 'cf73508' # src/coherence.py + draft prespec
print(f'Prespec freeze: {PRESPEC_FREEZE} | statistic: {STATISTIC_COMMIT} | runtime HEAD: {git_head()}')

OUT = ITER / 'data'
OUT.mkdir(parents=True, exist_ok=True)
""")

co("""
import h4_analysis as h4
from src.cohort import build_cohort, assert_all
from src.dependency_stats import (load_gene_effects, load_all_gene_effects_for_line,
                                  standardize_within_line, crispr_gene_columns)
from src.verdict import evaluate_iteration, Condition, power_note

cohort = build_cohort()
assert_all(cohort)
ALT_POS = [str(m) for m in cohort.alt_positive['ModelID']]
ALT_NEG = [str(m) for m in cohort.alt_negative['ModelID']]
print(f'ALT+ {len(ALT_POS)} | ALT- {len(ALT_NEG)} | sizes {cohort.sizes}')

# Prespec §3 minimum-cell guard, checked before anything is computed.
assert len(ALT_POS) >= h4.MIN_LINES, f'guard: {len(ALT_POS)} < {h4.MIN_LINES}'

ALL_GENES = sorted({*h4.BIR_TRIPLET, *h4.POL_DELTA, *h4.BTR, *h4.HR_MODULE,
                    *h4.BATENBURG_BIR, *h4.EXTENDED_RESTART})
print(f'{len(ALL_GENES)} genes: {ALL_GENES}')
crispr_gene_columns(ALL_GENES)  # raises naming any absent symbol
""")

md("""
## 2 — The outcome is read here

Everything above this point is membership and column resolution. The next cell is the
first to load a Chronos **value**. Both the primary and every control come out of the same
single read, so no gene set can be added after a result is seen without re-running the
whole notebook.
""")
co("""
effects = load_gene_effects(ALT_POS, ALL_GENES)
effects_neg = load_gene_effects(ALT_NEG, ALL_GENES)
print(f'ALT+ {effects.shape} | ALT- {effects_neg.shape}')
assert not effects.isna().any().any()

# Prespec §2.1 recorded these marginals pre-freeze from the scan; re-derived here so the
# notebook is self-contained and the pre-registered numbers are checkable against it.
marginals = pd.DataFrame({
    'median': effects.median(), 'sd': effects.std(ddof=1),
    'n_lt_-0.5': (effects < -0.5).sum(), 'n_lt_-1.0': (effects < -1.0).sum(),
}).loc[list(h4.BIR_TRIPLET) + list(h4.STAGE1_HR)]
marginals.to_csv(OUT / 'h4_marginals.csv')
marginals
""")

md("## 3 — Primary: coherence of the Stage 1 BIR triplet")
co("""
primary = h4.module_coherence(effects, h4.BIR_TRIPLET, 'BIR triplet (Stage 1)')
print(f'T = {primary.statistic:+.4f}   p_perm = {primary.p_value:.4f}')
print(f'floor at n={primary.n_lines}, k={primary.n_genes}: {primary.floor:.4f}  '
      f'(reached: {abs(primary.statistic) >= primary.floor})')
print(f'ties present: {primary.has_ties}')
for (a, b), r in primary.pairwise.items():
    print(f'  rho({a}, {b}) = {r:+.4f}')
""")

md("""
### 3b — Per-pair detail with bootstrap CIs (S1) and Broderick's prediction (S1b)

Broderick 2023 places SMARCAL1 in the RAD52-**independent** BIR branch, predicting
`rho(RAD52, SMARCAL1) <= 0`. Registered as a directional secondary so this iteration can
corroborate a specific published claim even where its own premise fails.
""")
co("""
pairs_ci = h4.pairwise_bootstrap_ci(effects, h4.BIR_TRIPLET)
pairs_p = h4.pairwise_permutation_p(effects, h4.BIR_TRIPLET)
pairs = pairs_ci.merge(pairs_p[['pair', 'p_perm', 'floor_p05_k2']], on='pair')
pairs.to_csv(OUT / 'h4_pairwise.csv', index=False)

rho_rad52_smarcal1 = primary.pairwise[('RAD52', 'SMARCAL1')]
broderick_supported = rho_rad52_smarcal1 <= 0
print(f'S1b Broderick: rho(RAD52, SMARCAL1) = {rho_rad52_smarcal1:+.4f} '
      f'-> prediction (<= 0) {"SUPPORTED" if broderick_supported else "NOT supported"}')
pairs
""")

md("""
## 4 — Assay-sensitivity gate (§6)

Pol delta contains POLD3 itself. If POLD3's profile cannot correlate with its own obligate
complex partners in these 16 lines, its non-correlation with RAD52 and SMARCAL1 is
uninformative about modularity.

The gate needs **both** scales. Pan-essential genes can correlate through line-level screen
quality rather than shared function, so a gate passing only on raw values would be passing
through that artefact. Within-line standardisation removes each line's location and scale.
""")
co("""
pol_raw = h4.module_coherence(effects, h4.POL_DELTA, 'Pol delta (gate, raw)')
print(f'Pol delta raw:          T = {pol_raw.statistic:+.4f}  p = {pol_raw.p_value:.4f}')

# Within-line standardisation needs each line's genome-wide mean and SD (Stage 1 §6.3).
all_effects = load_all_gene_effects_for_line(ALT_POS)
z_all = standardize_within_line(all_effects)
z_cols = {g: c for g, c in crispr_gene_columns(ALL_GENES).items()}
z_effects = pd.DataFrame({g: z_all[c] for g, c in z_cols.items()}, index=z_all.index)
z_effects = z_effects.loc[ALT_POS]

pol_z = h4.module_coherence(z_effects, h4.POL_DELTA, 'Pol delta (gate, standardised)')
print(f'Pol delta standardised: T = {pol_z.statistic:+.4f}  p = {pol_z.p_value:.4f}')

GATE_PASSED = (pol_raw.statistic >= h4.GATE_MIN_T) and (pol_z.statistic >= h4.GATE_MIN_T)
print(f'\\nGate threshold T >= {h4.GATE_MIN_T} on BOTH scales -> '
      f'{"PASSED" if GATE_PASSED else "FAILED (design insensitive)"}')
""")

md("## 5 — Other control modules: BTR, HR, and the literature's own BIR module")
co("""
controls = [
    pol_raw, pol_z,
    h4.module_coherence(effects, h4.BTR, 'BTR (BLM/TOP3A/RMI1/RMI2)'),
    h4.module_coherence(effects, h4.HR_MODULE, 'HR module (Stage 1 subtrahend)'),
    h4.module_coherence(effects, h4.BATENBURG_BIR, 'Batenburg BIR (MUS81/RAD52/POLD3)'),
    h4.module_coherence(effects, h4.EXTENDED_RESTART, 'Extended fork restart (k=5)'),
]
table = pd.DataFrame([primary.row()] + [c.row() for c in controls])
table.to_csv(OUT / 'h4_modules.csv', index=False)
table
""")

md("""
## 6 — Matched-triplet calibration (§5)

A raw T near zero could mean "not a module" or merely "these marginals cannot produce
correlation at n = 16". The discriminator is a **mean- and variance-matched** null: random
triplets whose marginals mimic POLD3's, RAD52's and SMARCAL1's. An unmatched null would be
an artefact generator, because the three genes sit in three different dependency regimes.
""")
co("""
pool = all_effects.dropna(axis=1, how='any')
pool.columns = [c.split(' (')[0] for c in pool.columns]
pool = pool.loc[:, ~pool.columns.duplicated()]
print(f'complete-column reference pool: {pool.shape[1]} genes x {pool.shape[0]} lines')

matched_draws, matched_diag = h4.matched_triplet_null(pool, h4.BIR_TRIPLET)
matched_pct = float((matched_draws < primary.statistic).mean() * 100)
matched_q95 = float(np.quantile(matched_draws, h4.MATCHED_PERCENTILE / 100))
exceeds_matched = primary.statistic > matched_q95

matched_diag.to_csv(OUT / 'h4_matched_diagnostics.csv', index=False)
pd.DataFrame({'T': matched_draws}).to_csv(OUT / 'h4_matched_draws.csv', index=False)

print(f'matched null: mean {matched_draws.mean():+.4f}, sd {matched_draws.std():.4f}, '
      f'q95 {matched_q95:+.4f}')
print(f'observed T = {primary.statistic:+.4f} -> {matched_pct:.1f}th percentile; '
      f'exceeds q95: {exceeds_matched}')
matched_diag
""")

md("""
## 7 — S2a: variance decomposition of the Stage 1 composite

Deterministic given these 16 lines, and independent of the primary's outcome. If one gene
dominates the composite's variance, Stage 1's result was effectively a single-gene result —
which bears directly on its R8 check, where dropping SMARCAL1 moved delta from +0.281 to
+0.156.
""")
co("""
vardec = h4.variance_decomposition(effects)
vardec.to_csv(OUT / 'h4_variance_decomposition.csv', index=False)
print(f"composite variance = {vardec.attrs['composite_variance']:.5f}")
vardec.sort_values('share_of_total', ascending=False)
""")

md("## 8 — S2b, S3, S6: standardised primary, the ALT− arm, and batch robustness")
co("""
primary_z = h4.module_coherence(z_effects, h4.BIR_TRIPLET, 'BIR triplet (standardised)')
primary_neg = h4.module_coherence(effects_neg, h4.BIR_TRIPLET, 'BIR triplet (ALT-)')
print(f'S2b standardised: T = {primary_z.statistic:+.4f}  p = {primary_z.p_value:.4f}')
print(f'S3 ALT- (n={primary_neg.n_lines}, floor {primary_neg.floor:.3f}): '
      f'T = {primary_neg.statistic:+.4f}  p = {primary_neg.p_value:.4f}')

# S6 batch robustness: drop commercially-deposited lines, the one exogenous batch axis
# scan_metadata.md found usable.
from src.cohort import osteosarcoma_models
meta = osteosarcoma_models().set_index('ModelID')
source = meta.loc[ALT_POS, 'SourceDetail'].astype(str) if 'SourceDetail' in meta.columns else pd.Series('NA', index=ALT_POS)
commercial = source.str.contains('ATCC|DSMZ|JCRB|RIKEN|ECACC|KCLB', case=False, na=False)
keep = [m for m in ALT_POS if not commercial.get(m, False)]
print(f'S6: dropping {int(commercial.sum())} commercially-deposited lines -> n={len(keep)}')
primary_s6 = (h4.module_coherence(effects.loc[keep], h4.BIR_TRIPLET, 'BIR triplet (non-commercial)')
              if len(keep) >= h4.MIN_LINES else None)
if primary_s6 is not None:
    print(f'   T = {primary_s6.statistic:+.4f}  p = {primary_s6.p_value:.4f}  '
          f'floor {primary_s6.floor:.3f}')
else:
    print('   below the §3 minimum-cell guard, not computed')
""")

md("## 9 — Leave-one-out")
co("""
loo = h4.leave_one_out(effects, h4.BIR_TRIPLET)
loo.to_csv(OUT / 'h4_loo.csv', index=False)
loo_no_flip = bool((np.sign(loo['T']) == np.sign(primary.statistic)).all())
print(f'T range over {len(loo)} refits: [{loo["T"].min():+.4f}, {loo["T"].max():+.4f}]')
print(f'no sign flip: {loo_no_flip}')
loo
""")

md("""
## 10 — Mechanical verdict

`src.verdict.evaluate_iteration` is the campaign's exhaustive partition. The support
conditions can only move a result between `supported` and `partial_support`; they cannot
leave it unassigned.

Prespec §6's pre-committed interpretive rule applies first: if **all** positive controls
fail, the verdict is design-insensitivity, not evidence against the module.
""")
co("""
conditions = [
    Condition(key='all_pairs_positive',
              description='all three pairwise rho > 0',
              passed=bool(all(r > 0 for r in primary.pairwise.values())),
              detail='; '.join(f'{a}-{b}={r:+.3f}' for (a, b), r in primary.pairwise.items())),
    Condition(key='min_pair_magnitude',
              description=f'weakest pair rho >= {h4.MIN_PAIR_RHO}',
              passed=bool(min(primary.pairwise.values()) >= h4.MIN_PAIR_RHO),
              detail=f'weakest {min(primary.pairwise.values()):+.4f}'),
    Condition(key='exceeds_matched_triplets',
              description=f'T above the q{h4.MATCHED_PERCENTILE} of mean/SD-matched random triplets',
              passed=bool(exceeds_matched),
              detail=f'observed {primary.statistic:+.4f} vs q95 {matched_q95:+.4f} '
                     f'({matched_pct:.1f}th percentile)'),
    Condition(key='assay_sensitive',
              description=f'Pol delta T >= {h4.GATE_MIN_T} on raw and standardised',
              passed=bool(GATE_PASSED),
              detail=f'raw {pol_raw.statistic:+.4f}, standardised {pol_z.statistic:+.4f}'),
    Condition(key='loo_no_sign_flip',
              description='no leave-one-out refit flips the sign',
              passed=bool(loo_no_flip),
              detail=f'range [{loo["T"].min():+.4f}, {loo["T"].max():+.4f}]'),
]

v = evaluate_iteration(delta=primary.statistic, p_value=primary.p_value,
                       predicted_sign=+1, conditions=conditions,
                       attributes={
                           'statistic': 'mean pairwise Spearman rho',
                           'n_lines': primary.n_lines,
                           'detection_floor': round(primary.floor, 4),
                           'gate_passed': bool(GATE_PASSED),
                           'matched_percentile': round(matched_pct, 1),
                       })

print(f'OUTCOME: {v.outcome}   rule = {v.rule}')
print(v.headline)
for c in conditions:
    print(f'  {"PASS" if c.passed else "FAIL"}  {c.key}: {c.description}')
if not GATE_PASSED:
    print('\\nNOTE: §6 rule -- positive control failed, so this is design insensitivity, '
          'not evidence against the module.')
""")

md("## 11 — Outputs, report and sha256 manifest")
co("""
import json

result = {
    'iteration': '03_bir_module_coherence',
    'prespec_freeze': PRESPEC_FREEZE,
    'statistic': 'mean pairwise Spearman rho',
    'genes': list(h4.BIR_TRIPLET),
    'n_lines': primary.n_lines,
    'T': round(primary.statistic, 6),
    'p_perm': round(primary.p_value, 6),
    'predicted_sign': 1,
    'detection_floor': round(primary.floor, 6),
    'pairwise': {f'{a}|{b}': round(r, 6) for (a, b), r in primary.pairwise.items()},
    'gate_pol_delta_raw': round(pol_raw.statistic, 6),
    'gate_pol_delta_standardised': round(pol_z.statistic, 6),
    'gate_passed': bool(GATE_PASSED),
    'matched_percentile': round(matched_pct, 2),
    'matched_q95': round(matched_q95, 6),
    'exceeds_matched': bool(exceeds_matched),
    'standardised_T': round(primary_z.statistic, 6),
    'alt_negative_T': round(primary_neg.statistic, 6),
    'broderick_s1b_rho_rad52_smarcal1': round(rho_rad52_smarcal1, 6),
    'broderick_s1b_supported': bool(broderick_supported),
    'loo_T_min': round(float(loo['T'].min()), 6),
    'loo_T_max': round(float(loo['T'].max()), 6),
    'verdict': v.outcome,
}
(OUT / 'h4_primary_result.json').write_text(json.dumps(result, indent=2) + '\\n')
(OUT / 'h4_verdict.json').write_text(json.dumps(v.as_dict(), indent=2) + '\\n')

lines = [
    f'Prespec freeze: {PRESPEC_FREEZE} | statistic commit: {STATISTIC_COMMIT}',
    '',
    f'PRIMARY  T = {primary.statistic:+.4f}  p_perm = {primary.p_value:.4f}  '
    f'(n = {primary.n_lines}, floor {primary.floor:.4f})',
    f'VERDICT  {v.outcome}  (rule: {v.rule})',
    f'         {v.headline}',
    '',
    'Pairwise rho:',
    *[f'  {a}-{b}: {r:+.4f}' for (a, b), r in primary.pairwise.items()],
    '',
    f'Gate (Pol delta): raw {pol_raw.statistic:+.4f}, standardised {pol_z.statistic:+.4f}, '
    f'threshold {h4.GATE_MIN_T} -> {"PASSED" if GATE_PASSED else "FAILED"}',
    f'Matched triplets: observed at the {matched_pct:.1f}th percentile, q95 {matched_q95:+.4f}',
    f'S2b standardised T: {primary_z.statistic:+.4f}',
    f'S3 ALT- T: {primary_neg.statistic:+.4f} (n = {primary_neg.n_lines}, floor {primary_neg.floor:.3f})',
    f'S1b Broderick rho(RAD52,SMARCAL1) = {rho_rad52_smarcal1:+.4f} '
    f'({"supported" if broderick_supported else "not supported"})',
    f'LOO T range: [{loo["T"].min():+.4f}, {loo["T"].max():+.4f}], no sign flip: {loo_no_flip}',
    '',
    'Support conditions:',
    *[f'  {"PASS" if c.passed else "FAIL"} {c.key}' for c in conditions],
    '',
    'Module table:',
    table.to_string(index=False),
    '',
    'Variance decomposition of the Stage 1 composite:',
    vardec.to_string(index=False),
]
report = write_report('01', lines, outdir=OUT)

outputs = sorted(str(p) for p in OUT.glob('h4_*.csv'))
outputs += [str(OUT / 'h4_primary_result.json'), str(OUT / 'h4_verdict.json'), str(report)]
manifest = write_sha256_manifest('01', outputs, outdir=OUT)
record_provenance([{'notebook': '01_h4_module_coherence', 'run_timestamp': RUN_TIMESTAMP,
                    'git_head': git_head(), 'versions': VERSIONS, 'outcome': v.outcome}],
                  OUT / 'provenance.jsonl')
print(f'Report:   {report}')
print(f'Manifest: {manifest}')
print(f'OUTCOME:  {v.outcome}')
""")

nb = {
    "cells": [
        {"cell_type": t, "metadata": {}, "source": s.splitlines(keepends=True),
         **({"outputs": [], "execution_count": None} if t == "code" else {})}
        for t, s in CELLS
    ],
    "metadata": {
        "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
        "language_info": {"name": "python", "version": "3.12"},
    },
    "nbformat": 4, "nbformat_minor": 5,
}
out = pathlib.Path("notebooks/01_h4_module_coherence.ipynb")
out.parent.mkdir(parents=True, exist_ok=True)
out.write_text(json.dumps(nb, indent=1) + "\n")
print(f"wrote {out} with {len(CELLS)} cells")
