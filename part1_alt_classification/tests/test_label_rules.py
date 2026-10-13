"""
Tests for ``src/label_rules.py::compute_inclusion_flags`` and the
architecture → genomic / overall status propagation in
``src/manifest_builder.py::refresh_tcab1_status_from_architecture``.

These tests cover the three bugs found in the audit run of 2026-05-14:

1. ``has_rna = bool(row.get('has_rna'))`` silently turned the string
   ``"False"`` into ``True`` (every non-empty string is truthy in Python),
   so the ``has_rna`` gate had been effectively disabled across all flags.

2. ``compute_inclusion_flags`` did not include ``alt_negative_high_conf``
   in any Aim 2 universe, so the 11 high-confidence ALT-negative samples
   were dropped from Pass A's counts (aim2_alt_negative = 0).

3. ``prefill_sv_architecture_from_literature`` / ``prefill_pcgp_sv_architecture``
   wrote ``sv_architecture_17p13_1`` but never updated
   ``tcab1_genomic_status`` or ``tcab1_overall_status``. Literature samples
   carrying architecture calls (PCGP "Neither", "TP53-only", "Both-*";
   Lawlor "Both-*") were invisible to ``count_aim1``.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from src.label_rules import (  # noqa: E402
    compute_inclusion_flags,
    _coerce_bool,
    derive_alt_call_direction,
)
from src.manifest_builder import (  # noqa: E402
    refresh_tcab1_status_from_architecture,
    _genomic_status_from_arch,
    _DISRUPTING_ARCHITECTURES,
    _INTACT_ARCHITECTURES,
)


# ---------------------------------------------------------------------------
# _coerce_bool — the fix for Bug 3 (bool("False") == True)
# ---------------------------------------------------------------------------

class TestCoerceBool:
    @pytest.mark.parametrize(
        "value,expected",
        [
            (True, True),
            (False, False),
            ("True", True),
            ("False", False),  # the critical case
            ("true", True),
            ("false", False),
            ("TRUE", True),
            ("FALSE", False),
            ("1", True),
            ("0", False),
            ("yes", True),
            ("no", False),
            ("", False),
            ("   ", False),
            (None, False),
            (float("nan"), False),
        ],
    )
    def test_roundtrip(self, value, expected):
        assert _coerce_bool(value) is expected

    def test_bare_bool_trap_is_broken(self):
        """Confirms the bug we're protecting against: bool('False') == True."""
        assert bool("False") is True  # the broken behaviour we replaced
        assert _coerce_bool("False") is False  # our fix


# ---------------------------------------------------------------------------
# compute_inclusion_flags — Aim 2 must include negatives
# ---------------------------------------------------------------------------

def _row(**overrides) -> pd.Series:
    """Build a minimal manifest row with sane defaults."""
    base = {
        "has_rna": False,
        "tcab1_genomic_status": "unknown",
        "tcab1_functional_status": "unknown",
        "sv_architecture_17p13_1": "unknown",
        "alt_status": "alt_unknown",
        "tp53_status": "unknown",
        "telomerase_status": "unknown",
        "tcab1_rna_quantile_within_cohort": "unknown",
    }
    base.update(overrides)
    return pd.Series(base)


class TestAim2InclusionIncludesNegatives:
    """Bug 1: alt_negative_high_conf was excluded from every Aim 2 tier."""

    def test_negative_high_conf_in_strict(self):
        flags = compute_inclusion_flags(_row(alt_status="alt_negative_high_conf"))
        assert flags["include_aim2_strict"] is True
        assert flags["include_aim2_pragmatic"] is True
        assert flags["include_aim2_broad"] is True

    def test_gold_positive_still_in_strict(self):
        flags = compute_inclusion_flags(_row(alt_status="alt_gold_positive"))
        assert flags["include_aim2_strict"] is True
        assert flags["include_aim2_pragmatic"] is True
        assert flags["include_aim2_broad"] is True

    def test_silver_only_in_pragmatic_and_broad(self):
        flags = compute_inclusion_flags(_row(alt_status="alt_silver"))
        assert flags["include_aim2_strict"] is False
        assert flags["include_aim2_pragmatic"] is True
        assert flags["include_aim2_broad"] is True

    def test_unknown_excluded_from_all(self):
        flags = compute_inclusion_flags(_row(alt_status="alt_unknown"))
        assert flags["include_aim2_strict"] is False
        assert flags["include_aim2_pragmatic"] is False
        assert flags["include_aim2_broad"] is False

    def test_aim2_no_longer_requires_has_rna(self):
        """A peer-reviewed ALT-negative tumor without DepMap RNA-seq is still
        a valid Aim 2 training sample (negatives anchor the classifier)."""
        flags = compute_inclusion_flags(
            _row(alt_status="alt_negative_high_conf", has_rna="False")
        )
        assert flags["include_aim2_strict"] is True


class TestAim1InclusionUsesFunctionalReadout:
    """Bug 2 + relaxed gate: Aim 1 universe now accepts IHC / architecture
    evidence in lieu of DepMap RNA."""

    def test_lawlor_disrupted_via_architecture_qualifies(self):
        """PDX with Both-single-SV and genomic_status='disrupted_high_conf'
        from refresh_tcab1_status_from_architecture, no RNA."""
        flags = compute_inclusion_flags(_row(
            has_rna="False",
            sv_architecture_17p13_1="Both-single-SV",
            tcab1_genomic_status="disrupted_high_conf",
            tcab1_functional_status="functional_loss_possible",  # Lawlor IHC
        ))
        assert flags["include_aim1_primary_strict"] is True
        assert flags["include_aim1_primary_pragmatic"] is True
        assert flags["include_aim1_primary_broad"] is True

    def test_pcgp_neither_architecture_qualifies_as_intact(self):
        """PCGP tumor with Neither architecture → tcab1_genomic_status='intact'
        (set by refresh_tcab1_status_from_architecture). No RNA."""
        flags = compute_inclusion_flags(_row(
            has_rna="False",
            sv_architecture_17p13_1="Neither",
            tcab1_genomic_status="intact",
        ))
        assert flags["include_aim1_primary_strict"] is True
        assert flags["include_aim1_primary_pragmatic"] is True
        assert flags["include_aim1_primary_broad"] is True

    def test_depmap_disrupted_still_qualifies(self):
        """Backward compatibility: existing DepMap pathway still works."""
        flags = compute_inclusion_flags(_row(
            has_rna=True,
            tcab1_genomic_status="disrupted_high_conf",
            tcab1_functional_status="functional_loss_possible",
        ))
        assert flags["include_aim1_primary_strict"] is True

    def test_literature_without_any_readout_excluded(self):
        """A sample with no RNA, no IHC, no architecture call must not
        accidentally drift into the Aim 1 universe."""
        flags = compute_inclusion_flags(_row(
            has_rna="False",
            tcab1_genomic_status="unknown",
            tcab1_functional_status="unknown",
            sv_architecture_17p13_1="unknown",
        ))
        assert flags["include_aim1_primary_strict"] is False
        assert flags["include_aim1_primary_pragmatic"] is False
        assert flags["include_aim1_primary_broad"] is False


# ---------------------------------------------------------------------------
# refresh_tcab1_status_from_architecture
# ---------------------------------------------------------------------------

def _manifest_seed(rows: list[dict]) -> pd.DataFrame:
    """Build a tiny manifest with only the columns this helper touches."""
    base_cols = {
        "sample_id": "",
        "cohort": "PCGP_2014",
        "sv_architecture_17p13_1": "unknown",
        "sv_architecture_evidence": "",
        "tcab1_genomic_status": "unknown",
        "tcab1_genomic_evidence": "",
        "tcab1_functional_status": "unknown",
        "tcab1_overall_status": "unknown",
        "tcab1_confidence": 0.0,
    }
    full_rows = []
    for r in rows:
        merged = dict(base_cols)
        merged.update(r)
        full_rows.append(merged)
    return pd.DataFrame(full_rows)


class TestArchitectureMapping:
    @pytest.mark.parametrize("arch", sorted(_DISRUPTING_ARCHITECTURES))
    def test_disrupting_architectures_map_to_disrupted(self, arch):
        assert _genomic_status_from_arch(arch) == "disrupted_high_conf"

    @pytest.mark.parametrize("arch", sorted(_INTACT_ARCHITECTURES))
    def test_intact_architectures_map_to_intact(self, arch):
        assert _genomic_status_from_arch(arch) == "intact"

    @pytest.mark.parametrize("arch", ["unknown", "", "garbage", None])
    def test_uninformative_architectures_return_none(self, arch):
        assert _genomic_status_from_arch(arch) is None


class TestRefreshTcab1Status:
    def test_pcgp_neither_becomes_intact(self):
        m = _manifest_seed([{
            "sample_id": "SJOS_A",
            "cohort": "PCGP_2014",
            "sv_architecture_17p13_1": "Neither",
        }])
        out = refresh_tcab1_status_from_architecture(m, verbose=False)
        row = out.iloc[0]
        assert row["tcab1_genomic_status"] == "intact"
        assert "Neither" in row["tcab1_genomic_evidence"]

    def test_pcgp_both_single_sv_becomes_disrupted(self):
        m = _manifest_seed([{
            "sample_id": "SJOS_B",
            "cohort": "PCGP_2014",
            "sv_architecture_17p13_1": "Both-single-SV",
        }])
        out = refresh_tcab1_status_from_architecture(m, verbose=False)
        assert out.iloc[0]["tcab1_genomic_status"] == "disrupted_high_conf"

    def test_lawlor_tcab1_only_becomes_disrupted(self):
        m = _manifest_seed([{
            "sample_id": "OS33",
            "cohort": "Lawlor_2026_GCCRI",
            "sv_architecture_17p13_1": "TCAB1-only",
        }])
        out = refresh_tcab1_status_from_architecture(m, verbose=False)
        assert out.iloc[0]["tcab1_genomic_status"] == "disrupted_high_conf"

    def test_unknown_architecture_leaves_status_unchanged(self):
        m = _manifest_seed([{
            "sample_id": "ACH-000001",
            "cohort": "DepMap_26Q1",
            "sv_architecture_17p13_1": "unknown",
            "tcab1_genomic_status": "unknown",
        }])
        out = refresh_tcab1_status_from_architecture(m, verbose=False)
        assert out.iloc[0]["tcab1_genomic_status"] == "unknown"

    def test_does_not_downgrade_existing_disrupted(self):
        """A DepMap row with CN-derived 'disrupted_high_conf' must not be
        downgraded to 'intact' by a stale 'Neither' architecture call."""
        m = _manifest_seed([{
            "sample_id": "ACH-000410",
            "cohort": "DepMap_26Q1",
            "sv_architecture_17p13_1": "Neither",
            "tcab1_genomic_status": "disrupted_high_conf",
            "tcab1_genomic_evidence": "WRAP53 CN=-1.5 deep deletion",
        }])
        out = refresh_tcab1_status_from_architecture(m, verbose=False)
        assert out.iloc[0]["tcab1_genomic_status"] == "disrupted_high_conf"
        assert "WRAP53 CN" in out.iloc[0]["tcab1_genomic_evidence"]

    def test_overall_status_recomputed_when_unknown(self):
        """An overall_status of 'unknown' should be refreshed after the
        genomic update so audit downstream sees the new evidence."""
        m = _manifest_seed([{
            "sample_id": "SJOS_C",
            "cohort": "PCGP_2014",
            "sv_architecture_17p13_1": "Both-single-SV",
            "tcab1_functional_status": "unknown",
            "tcab1_overall_status": "unknown",
        }])
        out = refresh_tcab1_status_from_architecture(m, verbose=False)
        # gw=1.0 (disrupted_high_conf), fw=None (functional unknown) → disrupted_high_conf
        assert out.iloc[0]["tcab1_overall_status"] == "disrupted_high_conf"

    def test_conflicted_surfaces_even_over_existing_call(self):
        """When architecture says disrupted but functional says intact, the
        resulting 'conflicted' label must be visible (don't bury it under a
        stale 'unknown' or 'candidate_disrupted')."""
        m = _manifest_seed([{
            "sample_id": "SJOS_D",
            "cohort": "PCGP_2014",
            "sv_architecture_17p13_1": "Both-single-SV",
            "tcab1_functional_status": "functional_intact_supported",
            "tcab1_overall_status": "candidate_disrupted",
        }])
        out = refresh_tcab1_status_from_architecture(m, verbose=False)
        assert out.iloc[0]["tcab1_overall_status"] == "conflicted"

    def test_empty_manifest(self):
        m = _manifest_seed([])
        out = refresh_tcab1_status_from_architecture(m, verbose=False)
        assert len(out) == 0

    def test_idempotent(self):
        """Running twice should yield the same result."""
        m = _manifest_seed([
            {"sample_id": "A", "sv_architecture_17p13_1": "Neither"},
            {"sample_id": "B", "sv_architecture_17p13_1": "Both-single-SV"},
            {"sample_id": "C", "sv_architecture_17p13_1": "unknown"},
        ])
        once = refresh_tcab1_status_from_architecture(m, verbose=False)
        twice = refresh_tcab1_status_from_architecture(once, verbose=False)
        pd.testing.assert_frame_equal(once, twice)


# ---------------------------------------------------------------------------
# derive_alt_call_direction
# ---------------------------------------------------------------------------

class TestDeriveAltCallDirection:
    # Gold-positive — direction deterministic, pred_class irrelevant.
    def test_gold_positive_no_pred_class(self):
        assert derive_alt_call_direction("alt_gold_positive") == "positive"

    def test_gold_positive_with_alt_plus(self):
        assert derive_alt_call_direction("alt_gold_positive", "ALT+") == "positive"

    def test_gold_positive_with_alt_minus(self):
        # pred_class is ignored for gold; direction comes from alt_status.
        assert derive_alt_call_direction("alt_gold_positive", "ALT-") == "positive"

    # Gold-negative — direction deterministic.
    def test_negative_high_conf_no_pred_class(self):
        assert derive_alt_call_direction("alt_negative_high_conf") == "negative"

    def test_negative_high_conf_with_alt_minus(self):
        assert derive_alt_call_direction("alt_negative_high_conf", "ALT-") == "negative"

    def test_negative_high_conf_with_alt_plus(self):
        # pred_class is ignored for gold; direction comes from alt_status.
        assert derive_alt_call_direction("alt_negative_high_conf", "ALT+") == "negative"

    # Silver — direction resolved via pred_class.
    def test_silver_alt_plus_gives_positive(self):
        assert derive_alt_call_direction("alt_silver", "ALT+") == "positive"

    def test_silver_alt_minus_gives_negative(self):
        assert derive_alt_call_direction("alt_silver", "ALT-") == "negative"

    def test_silver_no_pred_class_gives_unknown(self):
        assert derive_alt_call_direction("alt_silver") == "unknown"

    def test_silver_none_pred_class_gives_unknown(self):
        assert derive_alt_call_direction("alt_silver", None) == "unknown"

    def test_silver_empty_string_pred_class_gives_unknown(self):
        assert derive_alt_call_direction("alt_silver", "") == "unknown"

    def test_silver_garbage_pred_class_gives_unknown(self):
        assert derive_alt_call_direction("alt_silver", "nan") == "unknown"

    # Unknown / proxy / unrecognised statuses.
    def test_alt_unknown_gives_unknown(self):
        assert derive_alt_call_direction("alt_unknown") == "unknown"

    def test_alt_proxy_gives_unknown(self):
        assert derive_alt_call_direction("alt_proxy", "ALT+") == "unknown"

    def test_empty_status_gives_unknown(self):
        assert derive_alt_call_direction("") == "unknown"

    def test_garbage_status_gives_unknown(self):
        assert derive_alt_call_direction("some_unknown_label", "ALT+") == "unknown"
