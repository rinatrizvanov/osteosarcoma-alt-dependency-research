"""Tests for ``src.recovery_inventory.py``."""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent))

from src.recovery_inventory import (  # noqa: E402
    ALT_NEGATIVE_LABEL,
    ALT_POSITIVE_LABEL,
    build_depmap_name_lookup,
    build_route_b_table,
    map_cell_line_to_model_id,
    normalize_cell_line_name,
    route_b_alt_rna_counts,
    select_route_b_candidates,
    summarize_de_nonneville,
    summarize_lawlor,
)


def test_normalize_cell_line_name():
    assert normalize_cell_line_name("Saos-2") == "SAOS2"
    assert normalize_cell_line_name("U-2 OS") == "U2OS"
    assert normalize_cell_line_name("") == ""


def test_map_cell_line_exact_and_prefix(tmp_path):
    model = pd.DataFrame([
        {
            "ModelID": "ACH-000410",
            "CellLineName": "Saos-2",
            "StrippedCellLineName": "SAOS2",
            "CCLEName": "SAOS2_BONE",
            "ModelIDAlias": "",
        },
        {
            "ModelID": "ACH-000082",
            "CellLineName": "G-292, clone A141B1",
            "StrippedCellLineName": "G292CLONEA141B1",
            "CCLEName": "G292CLONEA141B1_BONE",
            "ModelIDAlias": "",
        },
    ])
    lookup = build_depmap_name_lookup(model)
    mid, col = map_cell_line_to_model_id("SAOS2", lookup, model)
    assert mid == "ACH-000410"
    mid2, _ = map_cell_line_to_model_id("G292", lookup, model)
    assert mid2 == "ACH-000082"


def test_select_route_b_candidates_tag_and_altitude():
    manifest = pd.DataFrame([
        {"sample_id": "SAOS2", "cohort": "Taglialatela_2026", "alt_status": "x"},
        {"sample_id": "ACH-000410", "cohort": "DepMap_26Q1", "alt_status": "y"},
        {"sample_id": "OS33", "cohort": "Lawlor_2026_GCCRI", "alt_status": "z"},
    ])
    out = select_route_b_candidates(manifest, altitude_sample_ids={"ACH-000410"})
    assert set(out["sample_id"]) == {"SAOS2", "ACH-000410"}


def test_route_b_alt_rna_counts():
    rb = pd.DataFrame([
        {"in_depmap_rna": True, "alt_status": ALT_POSITIVE_LABEL, "mapped_model_id": "A"},
        {"in_depmap_rna": True, "alt_status": ALT_NEGATIVE_LABEL, "mapped_model_id": "B"},
        {"in_depmap_rna": False, "alt_status": ALT_POSITIVE_LABEL, "mapped_model_id": ""},
    ])
    c = route_b_alt_rna_counts(rb)
    assert c["n_in_depmap_rna"] == 2
    assert c["n_alt_positive_rna"] == 1
    assert c["n_alt_negative_rna"] == 1


def test_summarize_de_nonneville():
    manifest = pd.DataFrame([
        {"sample_id": "a", "cohort": "de_Nonneville_2022", "alt_status": ALT_POSITIVE_LABEL},
        {"sample_id": "b", "cohort": "de_Nonneville_2022", "alt_status": ALT_NEGATIVE_LABEL},
    ])
    s = summarize_de_nonneville(manifest)
    assert s["n_total"] == 2
    assert s["n_alt_positive"] == 1
    assert s["arrayexpress_accession"] == "E-MTAB-11743"


def test_summarize_lawlor_minimal():
    manifest = pd.DataFrame([
        {
            "sample_id": "OS33",
            "cohort": "Lawlor_2026_GCCRI",
            "alt_status": ALT_POSITIVE_LABEL,
            "tcab1_overall_status": "disrupted_high_conf",
            "sv_architecture_17p13_1": "Both-single-SV",
            "tcab1_genomic_status": "unknown",
            "tcab1_functional_status": "unknown",
            "tcab1_rna_quantile_within_cohort": "unknown",
            "include_aim1_primary_pragmatic": False,
            "include_aim1_primary_broad": False,
            "include_aim1_primary_strict": False,
            "include_tp53_sensitivity": False,
            "include_telomerase_negative_secondary": False,
        },
    ])
    s = summarize_lawlor(manifest)
    assert s["n_total_lawlor"] == 1
    assert "OS33" in s["sample_ids_all"]
