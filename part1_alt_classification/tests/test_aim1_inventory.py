"""
Tests for ``src/aim1_inventory.py``.

The file-reading helpers (``depmap_rna_model_ids``, ``target_os_rna_sample_ids``)
are tested with synthetic CSV/TSV-gz fixtures so the suite stays independent
of the real raw data on disk. The decision rule (``decide_aim1_design``) is
the load-bearing piece — it gets the most coverage.
"""
from __future__ import annotations

import gzip
import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from src.aim1_inventory import (  # noqa: E402
    EXPRESSION_SOURCES,
    WITHIN_COHORT_PER_ARM_MIN,
    assess_confounding,
    build_inventory,
    decide_aim1_design,
    depmap_rna_model_ids,
    lookup_file_backed_rna,
    select_pragmatic_contrast_rows,
    target_os_rna_sample_ids,
)


# ---------------------------------------------------------------------------
# Cheap I/O probes
# ---------------------------------------------------------------------------

def test_depmap_rna_model_ids_reads_only_modelid_column(tmp_path):
    p = tmp_path / "rna.csv"
    pd.DataFrame({
        "Unnamed: 0":        [0, 1, 2, 3],
        "SequencingID":      ["s1", "s2", "s3", "s4"],
        "ModelConditionID":  ["mc1", "mc2", "mc3", "mc1"],
        "ModelID":           ["ACH-1", "ACH-2", "ACH-3", "ACH-1"],
        "IsDefaultEntryForMC": ["Yes"] * 4,
        "TP53": [1.0, 2.0, 3.0, 4.0],
    }).to_csv(p, index=False)
    out = depmap_rna_model_ids(p)
    assert out == {"ACH-1", "ACH-2", "ACH-3"}


def test_depmap_rna_model_ids_missing_file_returns_empty(tmp_path):
    assert depmap_rna_model_ids(tmp_path / "nope.csv") == set()


def test_target_os_rna_sample_ids_reads_only_header(tmp_path):
    p = tmp_path / "target.tsv.gz"
    with gzip.open(p, "wt") as fh:
        fh.write("Ensembl_ID\tTARGET-40-A\tTARGET-40-B\tTARGET-40-C\n")
        fh.write("ENSG1\t1\t2\t3\n")
    out = target_os_rna_sample_ids(p)
    assert out == {"TARGET-40-A", "TARGET-40-B", "TARGET-40-C"}


def test_target_os_rna_sample_ids_missing_file_returns_empty(tmp_path):
    assert target_os_rna_sample_ids(tmp_path / "nope.tsv.gz") == set()


# ---------------------------------------------------------------------------
# lookup_file_backed_rna routing
# ---------------------------------------------------------------------------

def test_lookup_depmap_present():
    ok, src, dt = lookup_file_backed_rna(
        "ACH-1", "DepMap_26Q1",
        depmap_models={"ACH-1", "ACH-2"}, target_os_samples=set(),
    )
    assert ok is True
    assert "Expression" in src
    assert "log2" in dt


def test_lookup_depmap_absent():
    ok, _, _ = lookup_file_backed_rna(
        "ACH-99", "DepMap_26Q1",
        depmap_models={"ACH-1"}, target_os_samples=set(),
    )
    assert ok is False


def test_lookup_literature_cohort_always_false():
    for cohort in ("Lawlor_2026_GCCRI", "PCGP_2014",
                   "Taglialatela_2026", "de_Nonneville_2022"):
        ok, src, dt = lookup_file_backed_rna(
            "any-id", cohort, depmap_models=set(), target_os_samples=set(),
        )
        assert ok is False, f"{cohort} should have no file-backed RNA"
        assert "no expression" in src.lower() or "labels only" in src.lower()


def test_lookup_unknown_cohort_returns_false():
    ok, src, dt = lookup_file_backed_rna(
        "id", "TotallyMadeUpCohort",
        depmap_models=set(), target_os_samples=set(),
    )
    assert ok is False
    assert dt == "unknown"


# ---------------------------------------------------------------------------
# Decision rule
# ---------------------------------------------------------------------------

def _inv(rows: list[dict]) -> pd.DataFrame:
    """Build a synthetic inventory DataFrame."""
    base = {
        "sample_id": "",
        "cohort": "DepMap_26Q1",
        "tcab1_overall_status": "disrupted",
        "sample_type": "cell_line",
        "manifest_has_rna": True,
        "file_backed_rna": True,
        "expression_source": "synthetic",
        "data_type": "log2(TPM+1)",
    }
    return pd.DataFrame([{**base, **r} for r in rows])


def test_verdict_within_cohort_viable_single_cohort():
    rows = (
        [{"sample_id": f"d{i}", "tcab1_overall_status": "disrupted"} for i in range(WITHIN_COHORT_PER_ARM_MIN)]
        + [{"sample_id": f"i{i}", "tcab1_overall_status": "intact"}    for i in range(WITHIN_COHORT_PER_ARM_MIN)]
    )
    inv = _inv(rows)
    verdict = decide_aim1_design(inv, assess_confounding(inv))
    assert verdict["category"] == "within_cohort_viable"
    assert verdict["primary_cohort"] == "DepMap_26Q1"
    assert verdict["primary_n_disrupted"] >= WITHIN_COHORT_PER_ARM_MIN
    assert verdict["primary_n_intact"]    >= WITHIN_COHORT_PER_ARM_MIN


def test_verdict_within_cohort_viable_requires_same_data_type():
    rows = (
        [{"sample_id": f"d{i}", "tcab1_overall_status": "disrupted", "data_type": "log2TPM"}
         for i in range(WITHIN_COHORT_PER_ARM_MIN)]
        + [{"sample_id": f"i{i}", "tcab1_overall_status": "intact", "data_type": "microarray"}
           for i in range(WITHIN_COHORT_PER_ARM_MIN)]
    )
    inv = _inv(rows)
    verdict = decide_aim1_design(inv, assess_confounding(inv))
    assert verdict["category"] != "within_cohort_viable"


def test_verdict_pooled_exploratory():
    # Two cohorts on each arm, but no single cohort has both arms.
    rows = [
        {"sample_id": "d1", "cohort": "A", "tcab1_overall_status": "disrupted"},
        {"sample_id": "d2", "cohort": "B", "tcab1_overall_status": "disrupted"},
        {"sample_id": "i1", "cohort": "C", "tcab1_overall_status": "intact"},
        {"sample_id": "i2", "cohort": "D", "tcab1_overall_status": "intact"},
    ]
    inv = _inv(rows)
    verdict = decide_aim1_design(inv, assess_confounding(inv))
    assert verdict["category"] == "pooled_exploratory_only"


def test_verdict_not_estimable_zero_intact():
    # Mimics the live state of this project: all file-backed RNA is one arm.
    rows = [
        {"sample_id": "d1", "tcab1_overall_status": "disrupted"},
        {"sample_id": "d2", "tcab1_overall_status": "disrupted"},
    ]
    inv = _inv(rows)
    verdict = decide_aim1_design(inv, assess_confounding(inv))
    assert verdict["category"] == "not_estimable"
    assert verdict["fb_n_intact"] == 0


def test_verdict_not_estimable_disjoint_single_cohorts():
    # One disrupted from cohort A, one intact from cohort B — disjoint and each
    # arm spans only one cohort, so neither within-cohort nor pooled fits.
    rows = [
        {"sample_id": "d1", "cohort": "A", "tcab1_overall_status": "disrupted"},
        {"sample_id": "i1", "cohort": "B", "tcab1_overall_status": "intact"},
    ]
    inv = _inv(rows)
    verdict = decide_aim1_design(inv, assess_confounding(inv))
    assert verdict["category"] == "not_estimable"


# ---------------------------------------------------------------------------
# Confounding assessment internals
# ---------------------------------------------------------------------------

def test_assess_confounding_ignores_non_file_backed_rows():
    inv = _inv([
        {"sample_id": "d1", "tcab1_overall_status": "disrupted"},
        {"sample_id": "i1", "tcab1_overall_status": "intact"},
        {"sample_id": "ghost", "tcab1_overall_status": "intact",
         "file_backed_rna": False},
    ])
    a = assess_confounding(inv)
    assert a["n_total"] == 2
    assert a["n_disrupted"] == 1
    assert a["n_intact"] == 1


def test_assess_confounding_detects_within_cohort_same_dtype():
    inv = _inv([
        {"sample_id": "d1", "cohort": "X", "tcab1_overall_status": "disrupted", "data_type": "TPM"},
        {"sample_id": "i1", "cohort": "X", "tcab1_overall_status": "intact",    "data_type": "TPM"},
        {"sample_id": "d2", "cohort": "Y", "tcab1_overall_status": "disrupted", "data_type": "counts"},
        {"sample_id": "i2", "cohort": "Y", "tcab1_overall_status": "intact",    "data_type": "TPM"},
    ])
    a = assess_confounding(inv)
    assert a["has_within_cohort_co_occurrence"] is True
    assert ("X", "TPM") in a["cohorts_with_both_arms_same_dtype"]
    # Cohort Y co-occurs but data types differ, so it should NOT be in the
    # same-dtype list.
    assert ("Y", "TPM") not in a["cohorts_with_both_arms_same_dtype"]
    assert ("Y", "counts") not in a["cohorts_with_both_arms_same_dtype"]


def test_assess_confounding_empty():
    inv = pd.DataFrame(columns=[
        "sample_id", "cohort", "tcab1_overall_status", "sample_type",
        "manifest_has_rna", "file_backed_rna", "expression_source", "data_type",
    ])
    a = assess_confounding(inv)
    assert a["n_total"] == 0
    assert a["has_within_cohort_co_occurrence"] is False
    assert a["cohort_x_arm"].empty


# ---------------------------------------------------------------------------
# select_pragmatic_contrast_rows + build_inventory glue
# ---------------------------------------------------------------------------

def test_select_pragmatic_contrast_rows_uses_canonical_classifier():
    # Build a tiny manifest with one disrupted (genomic disrupted_high_conf),
    # one intact (genomic intact), one uncertain (only overall=candidate_disrupted).
    rows = [
        {
            "sample_id": "d1", "cohort": "Lawlor_2026_GCCRI",
            "tcab1_genomic_status": "disrupted_high_conf",
            "tcab1_overall_status": "disrupted_high_conf",
            "tcab1_functional_status": "unknown",
            "tcab1_rna_quantile_within_cohort": "unknown",
            "sv_architecture_17p13_1": "TCAB1-only",
            "include_aim1_primary_strict": True,
            "include_aim1_primary_pragmatic": True,
            "include_aim1_primary_broad": True,
            "include_tp53_sensitivity": False,
            "include_telomerase_negative_secondary": False,
            "include_aim2_strict": False, "include_aim2_pragmatic": False,
            "include_aim2_broad": False, "alt_status": "alt_unknown",
            "tp53_status": "altered", "tert_expression": "unknown",
            "telomerase_status": "unknown", "atrx_status": "unknown",
            "has_rna": False, "sample_type": "pdx",
        },
        {
            "sample_id": "i1", "cohort": "PCGP_2014",
            "tcab1_genomic_status": "intact",
            "tcab1_overall_status": "intact",
            "tcab1_functional_status": "unknown",
            "tcab1_rna_quantile_within_cohort": "unknown",
            "sv_architecture_17p13_1": "Neither",
            "include_aim1_primary_strict": True,
            "include_aim1_primary_pragmatic": True,
            "include_aim1_primary_broad": True,
            "include_tp53_sensitivity": False,
            "include_telomerase_negative_secondary": False,
            "include_aim2_strict": False, "include_aim2_pragmatic": False,
            "include_aim2_broad": False, "alt_status": "alt_unknown",
            "tp53_status": "altered", "tert_expression": "unknown",
            "telomerase_status": "unknown", "atrx_status": "unknown",
            "has_rna": False, "sample_type": "tumor",
        },
        {
            "sample_id": "u1", "cohort": "DepMap_26Q1",
            "tcab1_genomic_status": "unknown",
            "tcab1_overall_status": "candidate_disrupted",
            "tcab1_functional_status": "unknown",
            "tcab1_rna_quantile_within_cohort": "0.10",
            "sv_architecture_17p13_1": "unknown",
            "include_aim1_primary_strict": False,
            "include_aim1_primary_pragmatic": False,
            "include_aim1_primary_broad": True,
            "include_tp53_sensitivity": False,
            "include_telomerase_negative_secondary": False,
            "include_aim2_strict": False, "include_aim2_pragmatic": False,
            "include_aim2_broad": False, "alt_status": "alt_unknown",
            "tp53_status": "altered", "tert_expression": "unknown",
            "telomerase_status": "unknown", "atrx_status": "unknown",
            "has_rna": True, "sample_type": "cell_line",
        },
    ]
    manifest = pd.DataFrame(rows)
    out = select_pragmatic_contrast_rows(manifest)
    assert set(out["arm"]) == {"disrupted", "intact"}
    assert len(out) == 2
    assert "u1" not in set(out["sample_id"])


def test_build_inventory_columns_present(monkeypatch):
    contrast = pd.DataFrame([
        {"sample_id": "ACH-1", "cohort": "DepMap_26Q1", "sample_type": "cell_line",
         "has_rna": True, "arm": "disrupted"},
        {"sample_id": "ACH-99", "cohort": "DepMap_26Q1", "sample_type": "cell_line",
         "has_rna": True, "arm": "intact"},
        {"sample_id": "P-1", "cohort": "PCGP_2014", "sample_type": "tumor",
         "has_rna": False, "arm": "intact"},
    ])
    inv = build_inventory(contrast, depmap_models={"ACH-1"}, target_os_samples=set())
    assert list(inv.columns) == [
        "sample_id", "cohort", "tcab1_overall_status", "sample_type",
        "manifest_has_rna", "file_backed_rna", "expression_source", "data_type",
    ]
    idx = inv.set_index("sample_id")
    assert bool(idx.loc["ACH-1", "file_backed_rna"])     is True
    assert bool(idx.loc["ACH-99", "file_backed_rna"])    is False
    assert bool(idx.loc["P-1",   "file_backed_rna"])     is False


def test_expression_sources_registry_has_all_live_cohorts():
    live_cohorts = {"DepMap_26Q1", "Lawlor_2026_GCCRI", "PCGP_2014",
                    "Taglialatela_2026", "de_Nonneville_2022"}
    assert live_cohorts.issubset(set(EXPRESSION_SOURCES.keys()))
