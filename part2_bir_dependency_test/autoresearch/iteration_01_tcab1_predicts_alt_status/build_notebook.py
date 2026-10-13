"""Builds notebooks/01_h1_alt_status.ipynb from the cell sources below."""
import json, pathlib

CELLS = []
def md(s): CELLS.append(("markdown", s.strip("\n")))
def co(s): CELLS.append(("code", s.strip("\n")))

md("""
# 01 — H1: does WRAP53 expression discriminate ALT status?

Implements `prespec.md` **v1.0.1**, frozen at commit `2d68239` (erratum) on top of
`8e6ed51` (v1.0 freeze), both strictly before any exposure–outcome join.

**Reads no gene-effect value.** The outcome is the frozen ALT label.
`dependency_stats.load_gene_effects` is monkeypatched to raise, so any attempt to
reach Chronos data fails loudly rather than silently.

Prediction: ALT+ lines have **lower** WRAP53, i.e. Cliff's δ < 0, `predicted_sign = -1`.
Power, pre-registered: at 18v12 the minimum detectable |δ| is **0.435**; at Keegan's own
35% prevalence with perfect transcript separation this design has **27% power**. A
non-significant result in the predicted direction is therefore **uninformative** and must
not be read as evidence against H1 (prespec §8).
""")

md("## 1 — Header and gate")
co("""
import sys, random, warnings
from pathlib import Path

import numpy as np

random.seed(42)
np.random.seed(42)

_cwd = Path('.').resolve()
ROOT = next((p for p in [_cwd, *_cwd.parents] if (p / 'src').is_dir()), None)
assert ROOT is not None, f'Cannot locate src/ from {_cwd}'
ITER = ROOT / 'autoresearch' / 'iteration_01_tcab1_predicts_alt_status'
for p in (str(ROOT), str(ITER)):
    if p not in sys.path:
        sys.path.insert(0, p)
print(f'ROOT: {ROOT}')

warnings.filterwarnings('ignore', category=FutureWarning)

import pandas as pd
from src.utils import log_versions, run_timestamp, record_provenance, write_report, write_sha256_manifest, git_head

RUN_TIMESTAMP = run_timestamp()
VERSIONS = log_versions(['pandas', 'numpy', 'scipy'])
assert sys.version_info[:2] == (3, 12), f'Expected Python 3.12, got {sys.version_info[:2]}'
print(f'Python: {sys.version.split()[0]}')

PRESPEC_FREEZE = '8e6ed51'      # prespec v1.0 frozen
PRESPEC_ERRATUM = '2d68239'     # v1.0.1 AUROC transform erratum, pre-analysis
print(f'Prespec freeze: {PRESPEC_FREEZE} | erratum: {PRESPEC_ERRATUM} | runtime HEAD: {git_head()}')

OUT = ITER / 'data'
OUT.mkdir(parents=True, exist_ok=True)
""")

md("""
### Outcome gate

This iteration must not touch Chronos gene effect. Rather than merely asserting the
module is unimported — it is imported, for its pure statistical helpers — the data
loader itself is replaced with a raiser. That is the stronger guarantee.
""")
co("""
import src.dependency_stats as _ds

_CALLS = []
def _forbidden(*a, **k):
    _CALLS.append(1)
    raise AssertionError('OUTCOME GATE: load_gene_effects() is forbidden in iteration 01')
_ds.load_gene_effects = _forbidden
print('Outcome gate armed: dependency_stats.load_gene_effects now raises.')

import h1_analysis as H
print(f'Analysis module loaded. predicted_sign = {H.PREDICTED_SIGN}')
""")

md("## 2 — Cohort (prespec §4)")
co("""
cohort_full, trace = H.build_h1_cohort()
print('\\n'.join(trace.as_markdown()))

H.assert_anchors(cohort_full)
H.assert_no_forbidden_columns(cohort_full)
print('\\nAnchors asserted (U2OS positive, SJSA1 negative); no forbidden columns present.')

cohort = H.collapse_patients(cohort_full)
n_pos = int((cohort.arm == 'ALT+').sum()); n_neg = int((cohort.arm == 'ALT-').sum())
print(f'Primary (patient-collapsed): {n_pos} ALT+ / {n_neg} ALT-  N={len(cohort)}')
print(f'Dropped for patient collapse: {H.PATIENT_COLLAPSE_DROP} (shares PatientID with {H.PATIENT_COLLAPSE_KEEP})')

assert n_pos >= 4 and n_neg >= 4, 'prespec §10 minimum-cell guard fired'
print('§10 minimum-cell guard: not triggered.')

cohort_full.to_csv(OUT / 'h1_cohort_uncollapsed.csv')
cohort.to_csv(OUT / 'h1_cohort.csv')
""")

md("## 3 — Exposure")
co("""
from src.stratify import load_expression

expr = load_expression(list(cohort_full.index), symbols=[H.TCAB1_GENE, H.TP53_GENE])
wrap53 = expr[H.TCAB1_GENE]
tp53 = expr[H.TP53_GENE]
print(f'WRAP53 loaded for {wrap53.notna().sum()} lines; {(wrap53 == 0).sum()} at exactly 0.')

cn = H.pd.Series(dtype=float)
try:
    from src.stratify import load_wrap53_cn
    cn = load_wrap53_cn(list(cohort_full.index))
    print(f'WRAP53 CN loaded for {cn.notna().sum()} of {len(cohort_full)} lines.')
except Exception as exc:
    print(f'WRAP53 CN unavailable: {exc}')
""")

md("## 4 — Primary analysis (prespec §6)")
co("""
primary = H.contrast(wrap53, cohort, 'primary_18v12')
for k, v in primary.items():
    print(f'  {k:24s} {v}')
""")

md("## 5 — Robustness R1–R9 (prespec §7)")
co("""
rows = [primary]

gold = cohort[cohort.alt_status.isin(H.GOLD_TIERS)]
rows.append(H.contrast(wrap53, gold, 'R1_gold_only'))

from src.cohort import crispr_model_ids
_crispr = crispr_model_ids()
rows.append(H.contrast(wrap53, cohort_full[cohort_full.index.isin(_crispr)], 'R2_stage1_lines'))

rows.append(H.contrast(wrap53, cohort.drop(index=[i for i in H.KNOWN_D1 if i in cohort.index]), 'R3_drop_D1'))
rows.append(H.contrast(wrap53, cohort.drop(index=[i for i in H.KNOWN_ALL if i in cohort.index]), 'R3b_drop_all_known'))
if len(cn):
    rows.append(H.contrast(cn, cohort, 'R4_wrap53_copy_number'))
rows.append(H.contrast(tp53, cohort, 'R5_tp53_exposure_descriptive'))
rows.append(H.contrast(wrap53, cohort_full, 'R7_uncollapsed_18v13'))

robust = H.pd.DataFrame(rows)
print(robust[['name', 'n_pos', 'n_neg', 'cliffs_delta', 'ci_low', 'ci_high', 'p_value',
              'auroc', 'min_detectable_delta', 'mwu_method']].to_string(index=False))
robust.to_csv(OUT / 'h1_robustness.csv', index=False)
""")

co("""
loo = H.leave_one_out_deltas(wrap53, cohort)
loo_ok = all(d < 0 for d in loo) or all(d > 0 for d in loo)
print(f'R6 leave-one-out: {len(loo)} refits, delta range [{min(loo):+.3f}, {max(loo):+.3f}], '
      f'no sign flip = {loo_ok}')
H.pd.Series(loo, index=cohort.index, name='loo_cliffs_delta').to_csv(OUT / 'h1_loo_deltas.csv')

lodo = H.leave_one_depositor_out(wrap53, cohort)
_fin = lodo.cliffs_delta.dropna()
lodo_ok = bool(len(_fin)) and (all(_fin < 0) or all(_fin > 0))
print(f'\\nR8 leave-one-depositor-out: {len(lodo)} depositors, no sign flip = {lodo_ok}')
print(lodo.to_string(index=False))
lodo.to_csv(OUT / 'h1_leave_one_depositor_out.csv', index=False)
""")

co("""
TERTILE_CUT = 3.335580   # prespec §5, strict '<', from scan_exposure.md
r9 = H.tertile_enrichment(wrap53, cohort, TERTILE_CUT)
print('R9 bottom-tertile enrichment (DESCRIPTIVE ONLY -- cannot reach p<0.05):')
for k, v in r9.items():
    print(f'  {k:22s} {v}')
""")

md("## 6 — Verdict (prespec §8), computed mechanically")
co("""
from src.verdict import Condition, evaluate_iteration

def _d(name):
    m = robust[robust.name == name]
    return float(m.cliffs_delta.iloc[0]) if len(m) else float('nan')

conds = [
    Condition('primary', 'p<0.05 and delta <= -0.330',
              primary['p_value'] < 0.05 and primary['cliffs_delta'] <= -0.330,
              f"p={primary['p_value']:.4g}, delta={primary['cliffs_delta']:+.3f}"),
    Condition('r1_gold_only', 'R1 same direction', _d('R1_gold_only') < 0,
              f"delta={_d('R1_gold_only'):+.3f}"),
    Condition('r3_drop_known', 'R3 same direction', _d('R3_drop_D1') < 0,
              f"delta={_d('R3_drop_D1'):+.3f}"),
    Condition('r3b_drop_all_known', 'R3b same direction', _d('R3b_drop_all_known') < 0,
              f"delta={_d('R3b_drop_all_known'):+.3f}"),
    Condition('r6_loo', 'no sign flip on leave-one-out', loo_ok,
              f"range [{min(loo):+.3f}, {max(loo):+.3f}]"),
    Condition('r8_loDo', 'no sign flip on leave-one-depositor-out', lodo_ok,
              f"{len(_fin)} evaluable depositors"),
    Condition('discrimination', 'AUROC >= 0.70 (delta <= -0.400)',
              primary['auroc'] >= 0.70, f"AUROC={primary['auroc']:.3f}"),
]

v = evaluate_iteration(
    delta=primary['cliffs_delta'], p_value=primary['p_value'],
    predicted_sign=H.PREDICTED_SIGN, conditions=conds,
    attributes={'auroc': primary['auroc'],
                'min_detectable_delta': primary['min_detectable_delta'],
                'mwu_method': primary['mwu_method'],
                'n_pos': primary['n_pos'], 'n_neg': primary['n_neg']},
)
print(f'OUTCOME: {v.outcome}  (rule: {v.rule})')
print(v.headline)
print()
for c in conds:
    print(f"  [{'x' if c.passed else ' '}] {c.key:20s} {c.detail}")

import json as _json
(OUT / 'h1_verdict.json').write_text(_json.dumps(v.as_dict(), indent=2, sort_keys=True) + '\\n')
(OUT / 'h1_primary_result.json').write_text(_json.dumps(primary, indent=2, sort_keys=True) + '\\n')
""")

md("## 7 — Report and sha256 manifest")
co("""
L = []
L.append('# Iteration 01 — H1: WRAP53 expression vs ALT status')
L.append('')
L.append(f'Pre-specification: `prespec.md` v1.0.1, frozen `{PRESPEC_FREEZE}`, erratum `{PRESPEC_ERRATUM}`.')
L.append('Outcome variable: the frozen ALT label. **No gene-effect value was read** — '
         '`load_gene_effects` was replaced with a raiser and never fired '
         f'({len(_CALLS)} attempts).')
L.append('')
L.append('## Cohort (prespec §4)')
L.append('')
L += trace.as_markdown()
L.append('')
L.append(f'Primary is patient-collapsed: **{n_pos} ALT+ / {n_neg} ALT−**, N={len(cohort)}. '
         f'`{H.PATIENT_COLLAPSE_DROP}` dropped as a same-patient derivative of '
         f'`{H.PATIENT_COLLAPSE_KEEP}`, by the rule fixed in §4.')
L.append('')
L.append('## Primary result (prespec §6)')
L.append('')
L.append('| Quantity | Value |')
L.append('|---|---|')
for k in ['n_pos', 'n_neg', 'median_pos', 'median_neg', 'cliffs_delta', 'ci_low', 'ci_high',
          'hodges_lehmann', 'auroc', 'u_statistic', 'p_value', 'mwu_method', 'min_detectable_delta']:
    val = primary[k]
    L.append(f'| `{k}` | {val:.4g} |' if isinstance(val, float) else f'| `{k}` | {val} |')
L.append('')
L.append('## Robustness (prespec §7)')
L.append('')
_cols = ['name', 'n_pos', 'n_neg', 'cliffs_delta', 'ci_low', 'ci_high', 'p_value', 'auroc',
         'min_detectable_delta', 'mwu_method']
L.append('| ' + ' | '.join(_cols) + ' |')
L.append('|' + '---|' * len(_cols))
for _, r in robust.iterrows():
    cells = []
    for c in _cols:
        x = r[c]
        cells.append(f'{x:.4g}' if isinstance(x, float) else str(x))
    L.append('| ' + ' | '.join(cells) + ' |')
L.append('')
L.append(f'**R6 leave-one-out:** {len(loo)} refits, δ range '
         f'[{min(loo):+.3f}, {max(loo):+.3f}], no sign flip = {loo_ok}.')
L.append('')
L.append(f'**R8 leave-one-depositor-out:** {len(lodo)} depositors, no sign flip = {lodo_ok}.')
L.append('')
L.append('| depositor | n dropped | δ | p | note |')
L.append('|---|---|---|---|---|')
for _, r in lodo.iterrows():
    dd = f"{r['cliffs_delta']:+.3f}" if H.pd.notna(r['cliffs_delta']) else 'n/a'
    pp = f"{r['p_value']:.4g}" if H.pd.notna(r['p_value']) else 'n/a'
    L.append(f"| {r['depositor']} | {r['n_dropped']} | {dd} | {pp} | {r['note']} |")
L.append('')
L.append(f"**R9 bottom-tertile enrichment** (cut {r9['cut']:.6f}, strict `<`, descriptive only): "
         f"low arm {r9['table_low_pos_neg']} (ALT+/ALT−), high arm {r9['table_high_pos_neg']}, "
         f"OR {r9['odds_ratio']:.3f}, Fisher two-sided p = {r9['fisher_p_two_sided']:.4g}.")
L.append('')
L.append('## Verdict (prespec §8)')
L.append('')
L.append(f'**`{v.outcome}`** — rule `{v.rule}`.')
L.append('')
L.append(v.headline)
L.append('')
L.append('| Condition | Met | Detail |')
L.append('|---|---|---|')
for c in conds:
    L.append(f"| {c.description} | {'yes' if c.passed else 'NO'} | {c.detail} |")
L.append('')
L.append('## Power, as pre-registered (prespec §9)')
L.append('')
L.append(f"Minimum detectable |δ| at {n_pos}v{n_neg} is "
         f"**{primary['min_detectable_delta']:.3f}**. At Keegan's own 35% lesion prevalence "
         f"with perfect transcript separation this design has 27% power; under 13% at the "
         f"realistic attenuated magnitude (|δ| < 0.2). Per §8, a non-significant result in "
         f"the predicted direction is uninformative and is not evidence against H1.")

report = write_report('01', L, outdir=OUT)
outputs = sorted(str(p) for p in OUT.glob('h1_*')) + [str(report)]
manifest = write_sha256_manifest('01', outputs, outdir=OUT)
record_provenance([{'notebook': '01_h1_alt_status', 'run_timestamp': RUN_TIMESTAMP,
                    'git_head': git_head(), 'versions': VERSIONS,
                    'outcome': v.outcome}], OUT / 'provenance.jsonl')
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
out = pathlib.Path("notebooks/01_h1_alt_status.ipynb")
out.parent.mkdir(parents=True, exist_ok=True)
out.write_text(json.dumps(nb, indent=1) + "\n")
print(f"wrote {out} with {len(CELLS)} cells")
