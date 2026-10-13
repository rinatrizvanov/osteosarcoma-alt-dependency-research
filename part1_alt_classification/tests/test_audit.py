"""
Tests for ``src/audit.py::classify_aim1_arms`` and ``count_aim1``.

These tests pin down the single source of truth for Aim 1 arm assignment so
the notebook 05 accounting diagnostic cannot drift away from the audit.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from src.audit import (  # noqa: E402
    classify_aim1_arms,
    count_aim1,
    AIM1_DISRUPTED_GENOMIC,
    AIM1_INTACT_GENOMIC,
)


def _toy_manifest(rows: list[dict]) -> pd.DataFrame:
    """Build a minimal manifest with just the columns audit needs."""
    base = {
        "sample_id": "",
        "cohort": "DepMap_26Q1",
        "tcab1_genomic_status": "unknown",
        "tcab1_overall_status": "unknown",
        "tcab1_functional_status": "unknown",
        "tcab1_rna_quantile_within_cohort": "unknown",
        "sv_architecture_17p13_1": "unknown",
        "include_aim1_primary_strict": False,
        "include_aim1_primary_pragmatic": False,
        "include_aim1_primary_broad": False,
        "include_tp53_sensitivity": False,
        "include_telomerase_negative_secondary": False,
    }
    return pd.DataFrame([{**base, **r} for r in rows])


class TestClassifyAim1Arms:
    def test_excluded_when_not_in_universe(self):
        m = _toy_manifest([{
            "sample_id": "A",
            "tcab1_genomic_status": "disrupted_high_conf",
            "include_aim1_primary_broad": False,
        }])
        arms = classify_aim1_arms(m, "broad")
        assert arms.iloc[0] == "excluded"

    def test_disrupted_via_genomic_status(self):
        m = _toy_manifest([{
            "sample_id": "A",
            "tcab1_genomic_status": "disrupted_high_conf",
            "include_aim1_primary_strict": True,
        }])
        assert classify_aim1_arms(m, "strict").iloc[0] == "disrupted"

    def test_disrupted_via_candidate_regulatory(self):
        m = _toy_manifest([{
            "sample_id": "A",
            "tcab1_genomic_status": "candidate_regulatory_disruption",
            "include_aim1_primary_pragmatic": True,
        }])
        assert classify_aim1_arms(m, "pragmatic").iloc[0] == "disrupted"

    def test_disrupted_via_overall_status_only(self):
        """DepMap CN-derived 'disrupted_high_conf' overall without a genomic
        rewrite still counts as disrupted."""
        m = _toy_manifest([{
            "sample_id": "A",
            "tcab1_genomic_status": "unknown",
            "tcab1_overall_status": "disrupted_high_conf",
            "include_aim1_primary_broad": True,
        }])
        assert classify_aim1_arms(m, "broad").iloc[0] == "disrupted"

    def test_intact_via_genomic_status(self):
        m = _toy_manifest([{
            "sample_id": "A",
            "tcab1_genomic_status": "intact",
            "include_aim1_primary_strict": True,
        }])
        assert classify_aim1_arms(m, "strict").iloc[0] == "intact"

    def test_uncertain_when_in_universe_but_no_arm(self):
        """Broad-included DepMap row with RNA-quantile inclusion but
        ``tcab1_genomic_status == 'unknown'`` → uncertain."""
        m = _toy_manifest([{
            "sample_id": "A",
            "tcab1_genomic_status": "unknown",
            "tcab1_overall_status": "unknown",
            "include_aim1_primary_broad": True,
        }])
        assert classify_aim1_arms(m, "broad").iloc[0] == "uncertain"

    def test_disrupted_takes_precedence_over_intact(self):
        """A row simultaneously claiming disruption and intactness must land
        in disrupted (never intact). Should never happen, but pin it."""
        m = _toy_manifest([{
            "sample_id": "A",
            "tcab1_genomic_status": "intact",
            "tcab1_overall_status": "disrupted_high_conf",
            "include_aim1_primary_broad": True,
        }])
        assert classify_aim1_arms(m, "broad").iloc[0] == "disrupted"

    def test_unknown_def_set_returns_all_excluded(self):
        m = _toy_manifest([{
            "sample_id": "A",
            "tcab1_genomic_status": "disrupted_high_conf",
            "include_aim1_primary_strict": True,
        }])
        # missing include_aim1_primary_nonsense column → all excluded
        arms = classify_aim1_arms(m, "nonsense")
        assert (arms == "excluded").all()

    @pytest.mark.parametrize("g", AIM1_DISRUPTED_GENOMIC)
    def test_all_disrupted_states_route_to_disrupted(self, g):
        m = _toy_manifest([{
            "sample_id": "A",
            "tcab1_genomic_status": g,
            "include_aim1_primary_broad": True,
        }])
        assert classify_aim1_arms(m, "broad").iloc[0] == "disrupted"

    @pytest.mark.parametrize("g", AIM1_INTACT_GENOMIC)
    def test_all_intact_states_route_to_intact(self, g):
        m = _toy_manifest([{
            "sample_id": "A",
            "tcab1_genomic_status": g,
            "include_aim1_primary_broad": True,
        }])
        assert classify_aim1_arms(m, "broad").iloc[0] == "intact"


class TestCountAim1Identities:
    """``count_aim1`` is derived from ``classify_aim1_arms``; the relationships
    between the returned counts must be exact."""

    @pytest.fixture
    def manifest(self):
        return _toy_manifest([
            # Disrupted, broad universe
            {"sample_id": "D1", "tcab1_genomic_status": "disrupted_high_conf",
             "include_aim1_primary_broad": True, "include_aim1_primary_strict": True},
            # Intact, all definition sets
            {"sample_id": "I1", "tcab1_genomic_status": "intact",
             "include_aim1_primary_strict": True,
             "include_aim1_primary_pragmatic": True,
             "include_aim1_primary_broad": True},
            # Uncertain, broad only (RNA-proxy with unknown genomic)
            {"sample_id": "U1", "tcab1_genomic_status": "unknown",
             "include_aim1_primary_broad": True},
            # Excluded entirely
            {"sample_id": "X1", "tcab1_genomic_status": "disrupted_high_conf"},
        ])

    def test_universe_equals_arm_sum(self, manifest):
        for def_set in ("strict", "pragmatic", "broad"):
            c = count_aim1(manifest, def_set)
            assert c["aim1_total"] == (
                c["aim1_disrupted"] + c["aim1_intact"] + c["aim1_uncertain"]
            ), f"universe ≠ disrupted+intact+uncertain for {def_set}"

    def test_primary_contrast_equals_disrupted_plus_intact(self, manifest):
        for def_set in ("strict", "pragmatic", "broad"):
            c = count_aim1(manifest, def_set)
            assert c["aim1_primary_contrast"] == c["aim1_disrupted"] + c["aim1_intact"]

    def test_strict_has_no_uncertain_for_this_fixture(self, manifest):
        """No RNA-proxy rows reach strict, so uncertain must be 0."""
        c = count_aim1(manifest, "strict")
        assert c["aim1_uncertain"] == 0
        assert c["aim1_total"] == c["aim1_primary_contrast"] == 2  # D1 + I1

    def test_broad_picks_up_uncertain(self, manifest):
        c = count_aim1(manifest, "broad")
        assert c["aim1_uncertain"] == 1  # U1
        assert c["aim1_total"] == 3       # D1 + I1 + U1
        assert c["aim1_primary_contrast"] == 2

    def test_classify_and_count_agree_per_row(self, manifest):
        """Iterating over the rows classified by classify_aim1_arms must
        recover exactly the same counts that count_aim1 reports."""
        for def_set in ("strict", "pragmatic", "broad"):
            arms = classify_aim1_arms(manifest, def_set)
            c = count_aim1(manifest, def_set)
            assert int((arms == "disrupted").sum()) == c["aim1_disrupted"]
            assert int((arms == "intact").sum()) == c["aim1_intact"]
            assert int((arms == "uncertain").sum()) == c["aim1_uncertain"]


class TestCountAim1OnLiveManifest:
    """Smoke check that the live manifest on disk obeys the universe identity."""

    @pytest.fixture
    def live_manifest(self):
        path = Path(__file__).parent.parent / "data" / "manifests" / "master_manifest.csv"
        if not path.exists():
            pytest.skip("Live manifest not present — run notebooks 01→03 first.")
        return pd.read_csv(path, dtype=str)

    @pytest.mark.parametrize("def_set", ["strict", "pragmatic", "broad"])
    def test_universe_identity_on_live_data(self, live_manifest, def_set):
        """Caller does pre-coercion (defensive style)."""
        flag = f"include_aim1_primary_{def_set}"
        m = live_manifest.copy()
        m[flag] = m[flag].astype(str).str.strip().str.lower().isin(("true", "1"))
        c = count_aim1(m, def_set)
        assert c["aim1_total"] == c["aim1_disrupted"] + c["aim1_intact"] + c["aim1_uncertain"]
        assert c["aim1_primary_contrast"] == c["aim1_disrupted"] + c["aim1_intact"]

    @pytest.mark.parametrize("def_set", ["strict", "pragmatic", "broad"])
    def test_works_on_string_typed_csv_roundtrip(self, live_manifest, def_set):
        """count_aim1 must work when the manifest is loaded with ``dtype=str``
        (the default in notebook 05 / 04 / 03). The string-vs-bool trap:
        ``"True" == True`` is False, so a naïve implementation would silently
        return zeros for every count. This test catches that regression.
        """
        c = count_aim1(live_manifest, def_set)  # no pre-coercion!
        # Universe identity must still hold.
        assert c["aim1_total"] == c["aim1_disrupted"] + c["aim1_intact"] + c["aim1_uncertain"]
        # And at least one row must qualify; if everything is zero we've
        # silently lost the manifest.
        assert c["aim1_total"] > 0, (
            f"count_aim1('{def_set}') returned 0 on string-typed manifest — "
            "the CSV-roundtrip bool coercion has regressed."
        )

    @pytest.mark.parametrize("def_set", ["strict", "pragmatic", "broad"])
    def test_string_and_bool_manifests_give_identical_counts(self, live_manifest, def_set):
        """Coerced manifest and raw string manifest must produce the same counts."""
        m_str = live_manifest.copy()  # strings
        m_bool = live_manifest.copy()
        for col in [c for c in m_bool.columns if c.startswith("include_")]:
            m_bool[col] = m_bool[col].astype(str).str.strip().str.lower().isin(("true", "1"))

        c_str = count_aim1(m_str, def_set)
        c_bool = count_aim1(m_bool, def_set)
        assert c_str == c_bool

    @pytest.mark.parametrize("def_set", ["strict", "pragmatic", "broad"])
    def test_classify_aim1_arms_string_input(self, live_manifest, def_set):
        """classify_aim1_arms must not collapse every row to 'excluded'
        when the flag column is a string."""
        arms = classify_aim1_arms(live_manifest, def_set)
        assert (arms == "excluded").sum() < len(live_manifest), (
            "classify_aim1_arms returned 'excluded' for every row — the "
            "stringified inclusion-flag bug has regressed."
        )
