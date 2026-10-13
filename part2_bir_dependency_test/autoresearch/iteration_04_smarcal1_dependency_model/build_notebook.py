"""Builds notebooks/01_h6_primpol_smarcal1.ipynb from the cell sources below."""
import json, pathlib

CELLS = []
def md(s): CELLS.append(("markdown", s.strip("\n")))
def co(s): CELLS.append(("code", s.strip("\n")))

md("""
# 01 — H6: does PRIMPOL expression predict SMARCAL1 dependency in ALT+ osteosarcoma?

Implements `prespec.md`, frozen at `5f2f039`, before any exposure–outcome quantity was
computed.

Prediction: **rho < 0** — more PRIMPOL, more negative SMARCAL1 Chronos, more dependent.
`predicted_sign = -1`. Mechanism: Taglialatela 2026 (`10.1101/gad.353377.125`), where the
telomeric ssDNA that accumulates on SMARCAL1 loss is PRIMPOL-repriming-dependent.

**Power, pre-registered.** At n = 16 the smallest |rho| reaching two-sided p < 0.05 is
**0.500**, and the campaign bar (p < 0.033 at three primaries) needs **0.535**. A
non-significant result in the predicted direction is *uninformative about the mechanism*
and must not be reported as evidence against it (prespec §4).

**The exposure is the mechanism's weaker arm.** BLM — full suppression of SMARCAL1-dependent
growth defects in Taglialatela — failed the §3 gate on its distribution in this cohort
(prespec §2) and is reported here as explicitly uninterpretable.
""")

md("## 1 — Header, cohort, exposures")
co("""
import sys, random, warnings
from pathlib import Path

import numpy as np

random.seed(42)
np.random.seed(42)

_cwd = Path('.').resolve()
ROOT = next((p for p in [_cwd, *_cwd.parents] if (p / 'src').is_dir() and (p / 'references').is_dir()), None)
assert ROOT is not None, f'Cannot locate repo root from {_cwd}'
ITER = ROOT / 'autoresearch' / 'iteration_04_smarcal1_dependency_model'
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

PRESPEC_FREEZE = '5f2f039'
print(f'Prespec freeze: {PRESPEC_FREEZE} | runtime HEAD: {git_head()}')

OUT = ITER / 'data'
OUT.mkdir(parents=True, exist_ok=True)
""")

co("""
import h6_analysis as h6
from src.cohort import build_cohort, assert_all
from src.stratify import load_expression, osteosarcoma_rna_models
from src.verdict import evaluate_iteration, Condition

cohort = build_cohort()
assert_all(cohort)
ALT_POS = [str(m) for m in cohort.alt_positive['ModelID']]
ALT_NEG = [str(m) for m in cohort.alt_negative['ModelID']]
print(f'ALT+ {len(ALT_POS)} | ALT- {len(ALT_NEG)}')

EXPOSURES = [h6.PRIMARY_EXPOSURE, h6.NEGATIVE_CONTROL, h6.DESCRIPTIVE_ONLY]
expr_all = load_expression(None, EXPOSURES)
RNA34 = [m for m in osteosarcoma_rna_models() if m in expr_all.index]
expr16 = expr_all.loc[[m for m in ALT_POS if m in expr_all.index]]
expr9 = expr_all.loc[[m for m in ALT_NEG if m in expr_all.index]]
print(f'reference set n={len(RNA34)} | ALT+ with RNA n={len(expr16)} | ALT- with RNA n={len(expr9)}')
assert len(expr16) >= h6.MIN_LINES, f'guard: {len(expr16)} < {h6.MIN_LINES}'
""")

md("""
## 2 — The §3 exposure-variance gate, evaluated before the outcome is read

No Chronos value has been loaded at this point. The gate is the rule that disqualified BLM
in the prespec; it is re-evaluated here so the disqualification is reproducible from code
rather than asserted in prose.
""")
co("""
gates = [h6.evaluate_gate(g, expr16[g].values) for g in EXPOSURES]
gate_table = pd.DataFrame([g.row() for g in gates])
gate_table.to_csv(OUT / 'h6_exposure_gate.csv', index=False)

PRIMARY_GATE = next(g for g in gates if g.gene == h6.PRIMARY_EXPOSURE)
assert PRIMARY_GATE.passed, 'primary exposure fails its own pre-registered gate'
BLM_GATE = next(g for g in gates if g.gene == h6.DESCRIPTIVE_ONLY)
print(f'{h6.PRIMARY_EXPOSURE} gate passed: {PRIMARY_GATE.passed}')
print(f'{h6.DESCRIPTIVE_ONLY} gate passed: {BLM_GATE.passed}  <- prespec §2, expected to fail')
gate_table
""")

md("""
## 3 — The outcome is read here

First load of a Chronos value in this iteration. The primary, every control and the
screen-amplitude diagnostic all come from this one read.
""")
co("""
from src.dependency_stats import load_gene_effects, load_all_gene_effects_for_line

outcome16 = load_gene_effects(list(expr16.index), [h6.OUTCOME_GENE])[h6.OUTCOME_GENE]
outcome9 = load_gene_effects(list(expr9.index), [h6.OUTCOME_GENE])[h6.OUTCOME_GENE]
print(f'SMARCAL1 ALT+ n={len(outcome16)} median={outcome16.median():+.4f} sd={outcome16.std(ddof=1):.4f}')
""")

md("## 4 — Primary")
co("""
primary = h6.correlate(expr16[h6.PRIMARY_EXPOSURE], outcome16,
                       label='primary', exposure_name=h6.PRIMARY_EXPOSURE,
                       outcome_name=h6.OUTCOME_GENE)
print(f'rho = {primary.rho:+.4f}   p_perm = {primary.p_value:.4f}   n = {primary.n}')
print(f'floor for p<0.05 at this n: {primary.floor:.4f}  (reached: {abs(primary.rho) >= primary.floor})')
print(f'aligned (rho * predicted_sign) = {primary.rho * h6.PREDICTED_SIGN:+.4f}')
print(f'campaign bar needs |rho| >= 0.535 (p < {h6.CAMPAIGN_BAR_P})')
""")

md("## 5 — Controls and secondaries (S1, S2, S2b, S3, S4)")
co("""
from src.dependency_stats import standardize_within_line

# S2: screen amplitude -- each line's genome-wide mean gene effect.
all_effects16 = load_all_gene_effects_for_line(list(expr16.index))
screen_amplitude = all_effects16.mean(axis=1)
screen_amplitude.name = 'screen_mean'

z16 = standardize_within_line(all_effects16)
smarcal1_col = [c for c in z16.columns if c.split(' (')[0] == h6.OUTCOME_GENE][0]
outcome16_z = z16[smarcal1_col]
outcome16_z.name = h6.OUTCOME_GENE

analyses = [
    primary,
    h6.correlate(expr16[h6.NEGATIVE_CONTROL], outcome16, label='S1 negative control',
                 exposure_name=h6.NEGATIVE_CONTROL, outcome_name=h6.OUTCOME_GENE),
    h6.correlate(expr16[h6.PRIMARY_EXPOSURE], screen_amplitude, label='S2 screen amplitude',
                 exposure_name=h6.PRIMARY_EXPOSURE, outcome_name='genome-wide mean effect'),
    h6.correlate(expr16[h6.PRIMARY_EXPOSURE], outcome16_z, label='S2b standardised',
                 exposure_name=h6.PRIMARY_EXPOSURE, outcome_name=f'{h6.OUTCOME_GENE} (z within line)'),
    h6.correlate(expr9[h6.PRIMARY_EXPOSURE], outcome9, label='S3 ALT- specificity',
                 exposure_name=h6.PRIMARY_EXPOSURE, outcome_name=h6.OUTCOME_GENE),
    h6.correlate(expr16[h6.DESCRIPTIVE_ONLY], outcome16, label='S4 BLM (UNINTERPRETABLE)',
                 exposure_name=h6.DESCRIPTIVE_ONLY, outcome_name=h6.OUTCOME_GENE),
]
table = pd.DataFrame([a.row() for a in analyses])
table.to_csv(OUT / 'h6_analyses.csv', index=False)

S1, S2, S2b, S3, S4 = analyses[1], analyses[2], analyses[3], analyses[4], analyses[5]
print(f'S4 BLM rho = {S4.rho:+.4f} (p={S4.p_value:.4f}) -- FAILS the §3 gate '
      f'(Shapiro p={BLM_GATE.shapiro_p:.4f}); 2 outliers + 14 lines in a 1.5-unit band. '
      f'NOT evidence either way.')
table
""")

md("## 6 — Leave-one-out")
co("""
loo = h6.leave_one_out(expr16[h6.PRIMARY_EXPOSURE], outcome16)
loo.to_csv(OUT / 'h6_loo.csv', index=False)
loo_no_flip = bool((np.sign(loo['rho']) == np.sign(primary.rho)).all())
print(f'rho range over {len(loo)} refits: [{loo["rho"].min():+.4f}, {loo["rho"].max():+.4f}]')
print(f'no sign flip: {loo_no_flip}')
loo
""")

md("## 7 — Mechanical verdict")
co("""
conditions = [
    Condition(key='exposure_variance_gate',
              description='primary exposure passes the §3 gate',
              passed=bool(PRIMARY_GATE.passed),
              detail=f'sd={PRIMARY_GATE.sd:.3f}, shapiro p={PRIMARY_GATE.shapiro_p:.3f}, '
                     f'guarded gap={PRIMARY_GATE.guarded_gap:.3f}'),
    Condition(key='negative_control_null',
              description='ZRANB3 shows no relationship and a smaller |rho| than the primary',
              passed=bool(S1.p_value >= 0.05 and abs(S1.rho) < abs(primary.rho)),
              detail=f'ZRANB3 rho={S1.rho:+.4f}, p={S1.p_value:.4f}'),
    Condition(key='not_screen_confounded',
              description='|rho(exposure, screen amplitude)| < |primary rho|',
              passed=bool(abs(S2.rho) < abs(primary.rho)),
              detail=f'screen-amplitude rho={S2.rho:+.4f} vs primary {primary.rho:+.4f}'),
    Condition(key='loo_no_sign_flip',
              description='no leave-one-out refit flips the sign',
              passed=bool(loo_no_flip),
              detail=f'range [{loo["rho"].min():+.4f}, {loo["rho"].max():+.4f}]'),
]

v = evaluate_iteration(delta=primary.rho, p_value=primary.p_value,
                       predicted_sign=h6.PREDICTED_SIGN, conditions=conditions,
                       attributes={
                           'statistic': 'Spearman rho (PRIMPOL expression vs SMARCAL1 Chronos)',
                           'n_lines': primary.n,
                           'detection_floor': round(primary.floor, 4),
                           'campaign_bar_rho': 0.535,
                           'exposure_swapped_from': 'BLM (failed §3 gate)',
                       })
print(f'OUTCOME: {v.outcome}   rule = {v.rule}')
print(v.headline)
for c in conditions:
    print(f'  {"PASS" if c.passed else "FAIL"}  {c.key}: {c.detail}')
""")

md("## 8 — Outputs, report and sha256 manifest")
co("""
import json

result = {
    'iteration': '04_smarcal1_dependency_model',
    'prespec_freeze': PRESPEC_FREEZE,
    'statistic': 'Spearman rho',
    'exposure': h6.PRIMARY_EXPOSURE,
    'exposure_swapped_from': 'BLM',
    'outcome': h6.OUTCOME_GENE,
    'n_lines': primary.n,
    'rho': round(primary.rho, 6),
    'p_perm': round(primary.p_value, 6),
    'predicted_sign': h6.PREDICTED_SIGN,
    'aligned_delta': round(primary.rho * h6.PREDICTED_SIGN, 6),
    'detection_floor': round(primary.floor, 6),
    'S1_zranb3_rho': round(S1.rho, 6), 'S1_p': round(S1.p_value, 6),
    'S2_screen_amplitude_rho': round(S2.rho, 6),
    'S2b_standardised_rho': round(S2b.rho, 6), 'S2b_p': round(S2b.p_value, 6),
    'S3_alt_negative_rho': round(S3.rho, 6), 'S3_p': round(S3.p_value, 6),
    'S4_blm_rho_UNINTERPRETABLE': round(S4.rho, 6), 'S4_p': round(S4.p_value, 6),
    'S4_blm_gate_passed': bool(BLM_GATE.passed),
    'loo_rho_min': round(float(loo['rho'].min()), 6),
    'loo_rho_max': round(float(loo['rho'].max()), 6),
    'verdict': v.outcome,
}
(OUT / 'h6_primary_result.json').write_text(json.dumps(result, indent=2) + '\\n')
(OUT / 'h6_verdict.json').write_text(json.dumps(v.as_dict(), indent=2) + '\\n')

lines = [
    f'Prespec freeze: {PRESPEC_FREEZE}',
    '',
    f'PRIMARY  rho({h6.PRIMARY_EXPOSURE}, {h6.OUTCOME_GENE}) = {primary.rho:+.4f}  '
    f'p_perm = {primary.p_value:.4f}  n = {primary.n}  floor = {primary.floor:.4f}',
    f'         aligned = {primary.rho * h6.PREDICTED_SIGN:+.4f}  (predicted_sign = {h6.PREDICTED_SIGN})',
    f'VERDICT  {v.outcome}  (rule: {v.rule})',
    f'         {v.headline}',
    '',
    'Exposure gate (evaluated before the outcome was read):',
    gate_table.to_string(index=False),
    '',
    'Analyses:',
    table.to_string(index=False),
    '',
    f'LOO rho range: [{loo["rho"].min():+.4f}, {loo["rho"].max():+.4f}], no sign flip: {loo_no_flip}',
    '',
    'Support conditions:',
    *[f'  {"PASS" if c.passed else "FAIL"} {c.key}: {c.detail}' for c in conditions],
    '',
    'S4 BLM is UNINTERPRETABLE by prespec §2/§5: it fails the exposure gate',
    f'(Shapiro p = {BLM_GATE.shapiro_p:.4f}); its distribution in these lines is 2 low',
    'outliers plus 14 lines inside a 1.496-unit band. A large value is not evidence.',
]
report = write_report('01', lines, outdir=OUT)
outputs = sorted(str(p) for p in OUT.glob('h6_*.csv'))
outputs += [str(OUT / 'h6_primary_result.json'), str(OUT / 'h6_verdict.json'), str(report)]
manifest = write_sha256_manifest('01', outputs, outdir=OUT)
record_provenance([{'notebook': '01_h6_primpol_smarcal1', 'run_timestamp': RUN_TIMESTAMP,
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
out = pathlib.Path("notebooks/01_h6_primpol_smarcal1.ipynb")
out.parent.mkdir(parents=True, exist_ok=True)
out.write_text(json.dumps(nb, indent=1) + "\n")
print(f"wrote {out} with {len(CELLS)} cells")
