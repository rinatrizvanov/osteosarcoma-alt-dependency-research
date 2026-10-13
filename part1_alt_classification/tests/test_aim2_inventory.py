"""
Tests for ``src/aim2_inventory.py``.

Mirrors the structure of ``test_aim1_inventory.py`` — synthetic inventories
exercise each branch of the decision rule, plus a few glue tests for
``select_aim2_contrast_rows`` and ``build_inventory_aim2``. The cheap
file-backed probes themselves live in ``aim1_inventory`` and are already
covered by ``test_aim1_inventory.py``.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent))

from src.aim2_inventory import (  # noqa: E402
    AIM2_ARMS,
    ALT_NEGATIVE_LABEL,
    ALT_POSITIVE_LABEL,
    WITHIN_COHORT_PER_ARM_MIN,
    assess_confounding_aim2,
    build_inventory_aim2,
    decide_aim2_design,
    select_aim2_contrast_rows,
)


# ---------------------------------------------------------------------------
# select_aim2_contrast_rows
# ---------------------------------------------------------------------------

def test_select_aim2_contrast_keeps_only_gold_pos_and_neg_high_conf():
    manifest = pd.DataFrame([
        {"sample_id": "a", "alt_status": "alt_gold_positive",     "cohort": "X"},
        {"sample_id": "b", "alt_status": "alt_negative_high_conf","cohort": "Y"},
        {"sample_id": "c", "alt_status": "alt_unknown",           "cohort": "X"},
        {"sample_id": "d", "alt_status": "alt_silver",            "cohort": "Z"},
        {"sample_id": "e", "alt_status": "alt_proxy",             "cohort": "Z"},
    ])
    sub = select_aim2_contrast_rows(manifest)
    assert set(sub["sample_id"]) == {"a", "b"}
    assert set(sub["arm"]) == set(AIM2_ARMS)


def test_select_aim2_contrast_preserves_alt_status_column():
    manifest = pd.DataFrame([
        {"sample_id": "a", "alt_status": "alt_gold_positive", "cohort": "X"},
    ])
    sub = select_aim2_contrast_rows(manifest)
    assert sub.loc[0, "arm"] == "alt_gold_positive"
    assert sub.loc[0, "alt_status"] == "alt_gold_positive"


# ---------------------------------------------------------------------------
# build_inventory_aim2
# ---------------------------------------------------------------------------

def test_build_inventory_aim2_routes_depmap_membership():
    contrast = pd.DataFrame([
        {"sample_id": "ACH-1",  "cohort": "DepMap_26Q1",       "sample_type": "cell_line",
         "has_rna": True,  "arm": "alt_gold_positive"},
        {"sample_id": "ACH-99", "cohort": "DepMap_26Q1",       "sample_type": "cell_line",
         "has_rna": True,  "arm": "alt_negative_high_conf"},
        {"sample_id": "P-1",    "cohort": "PCGP_2014",         "sample_type": "tumor",
         "has_rna": False, "arm": "alt_negative_high_conf"},
    ])
    inv = build_inventory_aim2(
        contrast,
        depmap_models={"ACH-1"},
        target_os_samples=set(),
    )
    assert list(inv.columns) == [
        "sample_id", "cohort", "alt_status", "sample_type",
        "manifest_has_rna", "file_backed_rna", "expression_source", "data_type",
    ]
    idx = inv.set_index("sample_id")
    assert bool(idx.loc["ACH-1",  "file_backed_rna"]) is True
    assert bool(idx.loc["ACH-99", "file_backed_rna"]) is False
    assert bool(idx.loc["P-1",    "file_backed_rna"]) is False
    assert idx.loc["ACH-1",  "alt_status"] == "alt_gold_positive"
    assert idx.loc["P-1",    "alt_status"] == "alt_negative_high_conf"


# ---------------------------------------------------------------------------
# Decision rule
# ---------------------------------------------------------------------------

def _inv(rows: list[dict]) -> pd.DataFrame:
    base = {
        "sample_id": "",
        "cohort": "DepMap_26Q1",
        "alt_status": ALT_POSITIVE_LABEL,
        "sample_type": "cell_line",
        "manifest_has_rna": True,
        "file_backed_rna": True,
        "expression_source": "synthetic",
        "data_type": "log2(TPM+1)",
    }
    return pd.DataFrame([{**base, **r} for r in rows])


def test_verdict_within_cohort_viable_single_cohort():
    rows = (
        [{"sample_id": f"p{i}", "alt_status": ALT_POSITIVE_LABEL}
         for i in range(WITHIN_COHORT_PER_ARM_MIN)]
        + [{"sample_id": f"n{i}", "alt_status": ALT_NEGATIVE_LABEL}
           for i in range(WITHIN_COHORT_PER_ARM_MIN)]
    )
    inv = _inv(rows)
    verdict = decide_aim2_design(inv, assess_confounding_aim2(inv))
    assert verdict["category"] == "within_cohort_viable"
    assert verdict["primary_cohort"] == "DepMap_26Q1"
    assert verdict["primary_n_positive"] >= WITHIN_COHORT_PER_ARM_MIN
    assert verdict["primary_n_negative"] >= WITHIN_COHORT_PER_ARM_MIN


def test_verdict_within_cohort_requires_same_data_type():
    rows = (
        [{"sample_id": f"p{i}", "alt_status": ALT_POSITIVE_LABEL, "data_type": "TPM"}
         for i in range(WITHIN_COHORT_PER_ARM_MIN)]
        + [{"sample_id": f"n{i}", "alt_status": ALT_NEGATIVE_LABEL, "data_type": "microarray"}
           for i in range(WITHIN_COHORT_PER_ARM_MIN)]
    )
    inv = _inv(rows)
    verdict = decide_aim2_design(inv, assess_confounding_aim2(inv))
    assert verdict["category"] != "within_cohort_viable"


def test_verdict_pooled_exploratory():
    rows = [
        {"sample_id": "p1", "cohort": "A", "alt_status": ALT_POSITIVE_LABEL},
        {"sample_id": "p2", "cohort": "B", "alt_status": ALT_POSITIVE_LABEL},
        {"sample_id": "n1", "cohort": "C", "alt_status": ALT_NEGATIVE_LABEL},
        {"sample_id": "n2", "cohort": "D", "alt_status": ALT_NEGATIVE_LABEL},
    ]
    inv = _inv(rows)
    verdict = decide_aim2_design(inv, assess_confounding_aim2(inv))
    assert verdict["category"] == "pooled_exploratory_only"


def test_verdict_not_estimable_zero_negative():
    # Mimics the live data: lots of positives, no file-backed negatives.
    rows = [
        {"sample_id": "p1", "alt_status": ALT_POSITIVE_LABEL},
        {"sample_id": "p2", "alt_status": ALT_POSITIVE_LABEL},
    ]
    inv = _inv(rows)
    verdict = decide_aim2_design(inv, assess_confounding_aim2(inv))
    assert verdict["category"] == "not_estimable"
    assert verdict["fb_n_negative"] == 0


def test_verdict_not_estimable_when_no_file_backed_rna_at_all():
    # Live Aim 2 state in this project: 59 samples in contrast, 0 with RNA.
    rows = [
        {"sample_id": f"x{i}", "cohort": "PCGP_2014",
         "alt_status": ALT_POSITIVE_LABEL if i % 2 else ALT_NEGATIVE_LABEL,
         "file_backed_rna": False}
        for i in range(6)
    ]
    inv = _inv(rows)
    verdict = decide_aim2_design(inv, assess_confounding_aim2(inv))
    assert verdict["category"] == "not_estimable"
    assert verdict["fb_n_positive"] == 0
    assert verdict["fb_n_negative"] == 0


def test_verdict_not_estimable_disjoint_single_cohorts():
    rows = [
        {"sample_id": "p1", "cohort": "A", "alt_status": ALT_POSITIVE_LABEL},
        {"sample_id": "n1", "cohort": "B", "alt_status": ALT_NEGATIVE_LABEL},
    ]
    inv = _inv(rows)
    verdict = decide_aim2_design(inv, assess_confounding_aim2(inv))
    assert verdict["category"] == "not_estimable"


# ---------------------------------------------------------------------------
# assess_confounding_aim2 internals
# ---------------------------------------------------------------------------

def test_assess_confounding_ignores_non_file_backed():
    inv = _inv([
        {"sample_id": "p1", "alt_status": ALT_POSITIVE_LABEL},
        {"sample_id": "n1", "alt_status": ALT_NEGATIVE_LABEL},
        {"sample_id": "ghost", "alt_status": ALT_NEGATIVE_LABEL, "file_backed_rna": False},
    ])
    a = assess_confounding_aim2(inv)
    assert a["n_total"] == 2
    assert a["n_positive"] == 1
    assert a["n_negative"] == 1


def test_assess_confounding_detects_within_cohort_same_dtype():
    inv = _inv([
        {"sample_id": "p1", "cohort": "X", "alt_status": ALT_POSITIVE_LABEL, "data_type": "TPM"},
        {"sample_id": "n1", "cohort": "X", "alt_status": ALT_NEGATIVE_LABEL, "data_type": "TPM"},
        {"sample_id": "p2", "cohort": "Y", "alt_status": ALT_POSITIVE_LABEL, "data_type": "counts"},
        {"sample_id": "n2", "cohort": "Y", "alt_status": ALT_NEGATIVE_LABEL, "data_type": "TPM"},
    ])
    a = assess_confounding_aim2(inv)
    assert a["has_within_cohort_co_occurrence"] is True
    assert ("X", "TPM") in a["cohorts_with_both_arms_same_dtype"]
    assert ("Y", "TPM") not in a["cohorts_with_both_arms_same_dtype"]


def test_assess_confounding_empty():
    inv = pd.DataFrame(columns=[
        "sample_id", "cohort", "alt_status", "sample_type",
        "manifest_has_rna", "file_backed_rna", "expression_source", "data_type",
    ])
    a = assess_confounding_aim2(inv)
    assert a["n_total"] == 0
    assert a["has_within_cohort_co_occurrence"] is False
    assert a["cohort_x_arm"].empty
