"""Iteration 01 — H1: does WRAP53 expression discriminate ALT status?

Implements `prespec.md` v1.0.1 exactly. Reads no gene-effect value: the outcome
here is the frozen ALT label, not Chronos essentiality.

Sign convention, stated once. `dependency_stats.cliffs_delta(a, b)` returns
P(a > b) − P(a < b). With ``a`` = ALT-positive WRAP53 and ``b`` = ALT-negative
WRAP53, H1 predicts ALT+ lines are *lower*, hence **δ < 0** and
``predicted_sign = -1``. The discrimination metric for the classifier
"lower WRAP53 ⇒ ALT+" is therefore **AUROC = (1 − δ)/2** (prespec §14 erratum).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from src.dependency_stats import bootstrap_delta_ci, cliffs_delta, hodges_lehmann
from src.utils import RAW_DEPMAP_DIR, ROOT
from src import verdict as V

MANIFEST = ROOT / "references" / "master_manifest.csv"
MODEL_CSV = RAW_DEPMAP_DIR / "Model.csv"

TCAB1_GENE = "WRAP53"
TP53_GENE = "TP53"
PREDICTED_SIGN = -1

#: prespec §4 — patient collapse. 143B (ACH-001001) is a HOS (ACH-000613)
#: derivative; both are ALT-negative. Retain the Stage 1 cohort member.
PATIENT_COLLAPSE_DROP = "ACH-000613"
PATIENT_COLLAPSE_KEEP = "ACH-001001"

#: prespec §3 — label–exposure pairs known before the freeze.
KNOWN_D1 = ("ACH-002834", "ACH-003185")              # WRAP53 == 0, both ALT+
KNOWN_D2 = ("ACH-000364", "ACH-000748")              # U2OS +, SJSA1 −
KNOWN_ALL = KNOWN_D1 + KNOWN_D2

#: prespec §4 — abort-on-failure ground truth.
ANCHORS = {"ACH-000364": "positive", "ACH-000748": "negative"}

GOLD_TIERS = ("alt_gold_positive", "alt_negative_high_conf")

#: prespec §11 item 6 — forbidden variables.
FORBIDDEN_COLUMNS = ("sv_architecture_auto", "telomerase_status")


@dataclass
class FilterTrace:
    steps: list[tuple[str, int]] = field(default_factory=list)

    def add(self, label: str, n: int) -> None:
        self.steps.append((label, n))

    def as_markdown(self) -> list[str]:
        out = ["| Step | N |", "|---|---|"]
        out += [f"| {label} | {n} |" for label, n in self.steps]
        return out


def build_h1_cohort() -> tuple[pd.DataFrame, FilterTrace]:
    """prespec §4. Labelled DepMap_26Q1 osteosarcoma lines with RNA.

    CRISPR screening is deliberately *not* required — this iteration reads no
    gene-effect data, so requiring it would discard 6 informative lines.
    """
    from src.cohort import rna_model_ids

    trace = FilterTrace()
    man = pd.read_csv(MANIFEST, low_memory=False)
    trace.add("master_manifest rows", len(man))

    dep = man[man["cohort"] == "DepMap_26Q1"].copy()
    trace.add("cohort == DepMap_26Q1", len(dep))

    informative = dep[dep["alt_call_direction"].isin(["positive", "negative"])].copy()
    trace.add("alt_call_direction in {positive, negative}", len(informative))

    rna = rna_model_ids()
    cohort = informative[informative["sample_id"].isin(rna)].copy()
    trace.add("present in the expression matrix", len(cohort))

    cohort = cohort.rename(columns={"sample_id": "ModelID"})
    cohort["arm"] = np.where(
        cohort["alt_call_direction"] == "positive", "ALT+", "ALT-"
    )
    cohort = cohort[["ModelID", "arm", "alt_status", "alt_call_direction"]]

    meta = pd.read_csv(MODEL_CSV, low_memory=False)
    keep = ["ModelID", "StrippedCellLineName", "PatientID", "SourceType"]
    cohort = cohort.merge(meta[keep], on="ModelID", how="left")
    cohort = cohort.set_index("ModelID").sort_index()
    trace.add(
        f"arms: {(cohort.arm == 'ALT+').sum()} ALT+ / {(cohort.arm == 'ALT-').sum()} ALT-",
        len(cohort),
    )
    return cohort, trace


def assert_anchors(cohort: pd.DataFrame) -> None:
    """prespec §4 — abort rather than analyse a mis-joined cohort."""
    for model_id, direction in ANCHORS.items():
        if model_id not in cohort.index:
            raise AssertionError(f"anchor {model_id} absent from the cohort")
        got = cohort.loc[model_id, "alt_call_direction"]
        if got != direction:
            raise AssertionError(
                f"anchor {model_id}: expected {direction}, manifest says {got}"
            )


def assert_no_forbidden_columns(frame: pd.DataFrame) -> None:
    """prespec §11 item 6."""
    present = [c for c in FORBIDDEN_COLUMNS if c in frame.columns]
    if present:
        raise AssertionError(f"forbidden columns reached the analysis frame: {present}")


def collapse_patients(cohort: pd.DataFrame) -> pd.DataFrame:
    """prespec §4 — one line per patient, by the rule fixed in advance."""
    if PATIENT_COLLAPSE_DROP in cohort.index:
        return cohort.drop(index=PATIENT_COLLAPSE_DROP)
    return cohort


def auroc_from_delta(delta: float) -> float:
    """prespec §6 + §14 erratum: classifier is 'lower WRAP53 => ALT+'."""
    return (1.0 - float(delta)) / 2.0


def contrast(
    values: pd.Series, cohort: pd.DataFrame, name: str, *, seed: int = 42
) -> dict:
    """Two-sided Mann-Whitney + Cliff's δ + HL shift + AUROC for one contrast."""
    from scipy import stats

    ids = [i for i in cohort.index if i in values.index and np.isfinite(values.get(i, np.nan))]
    sub = cohort.loc[ids]
    pos = values.loc[[i for i in ids if sub.loc[i, "arm"] == "ALT+"]].to_numpy(float)
    neg = values.loc[[i for i in ids if sub.loc[i, "arm"] == "ALT-"]].to_numpy(float)

    if pos.size < 1 or neg.size < 1:
        return {"name": name, "n_pos": int(pos.size), "n_neg": int(neg.size),
                "cliffs_delta": float("nan"), "p_value": float("nan")}

    pooled = np.concatenate([pos, neg])
    has_ties = len(np.unique(pooled)) < pooled.size
    method = "asymptotic" if has_ties else "exact"
    res = stats.mannwhitneyu(pos, neg, alternative="two-sided", method=method)
    delta = cliffs_delta(pos, neg)
    ci_low, ci_high = bootstrap_delta_ci(pos, neg, n_boot=2000, seed=seed)
    return {
        "name": name,
        "n_pos": int(pos.size),
        "n_neg": int(neg.size),
        "median_pos": float(np.median(pos)),
        "median_neg": float(np.median(neg)),
        "u_statistic": float(res.statistic),
        "p_value": float(res.pvalue),
        "mwu_method": method,
        "cliffs_delta": float(delta),
        "ci_low": float(ci_low),
        "ci_high": float(ci_high),
        "hodges_lehmann": float(hodges_lehmann(pos, neg)),
        "auroc": auroc_from_delta(delta),
        "min_detectable_delta": V.min_detectable_delta(int(pos.size), int(neg.size)),
    }


def leave_one_out_deltas(values: pd.Series, cohort: pd.DataFrame) -> list[float]:
    """prespec R6."""
    out = []
    for drop in cohort.index:
        r = contrast(values, cohort.drop(index=drop), f"loo_{drop}")
        out.append(r["cliffs_delta"])
    return out


def leave_one_depositor_out(values: pd.Series, cohort: pd.DataFrame) -> pd.DataFrame:
    """prespec R8 — SourceType is the depositor handle (complete, 14 levels)."""
    rows = []
    for dep in sorted(cohort["SourceType"].dropna().unique()):
        sub = cohort[cohort["SourceType"] != dep]
        if (sub.arm == "ALT+").sum() < 2 or (sub.arm == "ALT-").sum() < 2:
            rows.append({"depositor": dep, "n_dropped": int((cohort.SourceType == dep).sum()),
                         "cliffs_delta": float("nan"), "p_value": float("nan"),
                         "note": "skipped: an arm fell below n=2"})
            continue
        r = contrast(values, sub, f"loDo_{dep}")
        rows.append({"depositor": dep, "n_dropped": int((cohort.SourceType == dep).sum()),
                     "cliffs_delta": r["cliffs_delta"], "p_value": r["p_value"],
                     "n_pos": r["n_pos"], "n_neg": r["n_neg"], "note": ""})
    return pd.DataFrame(rows)


def tertile_enrichment(values: pd.Series, cohort: pd.DataFrame, cut: float) -> dict:
    """prespec R9 — descriptive only; cannot reach significance (§7)."""
    from scipy import stats

    ids = [i for i in cohort.index if i in values.index]
    sub = cohort.loc[ids].copy()
    sub["low"] = values.loc[ids] < cut          # strict '<', pinned per §5
    tab = pd.crosstab(sub["low"], sub["arm"])
    for col in ("ALT+", "ALT-"):
        if col not in tab.columns:
            tab[col] = 0
    for row in (False, True):
        if row not in tab.index:
            tab.loc[row] = 0
    table = [[int(tab.loc[True, "ALT+"]), int(tab.loc[True, "ALT-"])],
             [int(tab.loc[False, "ALT+"]), int(tab.loc[False, "ALT-"])]]
    odds, p = stats.fisher_exact(table, alternative="two-sided")
    return {"cut": float(cut), "table_low_pos_neg": table[0],
            "table_high_pos_neg": table[1], "odds_ratio": float(odds),
            "fisher_p_two_sided": float(p)}
