"""Tests for the Section 2 inclusion filter and Section 3 label join.

Two layers: synthetic inputs exercising every edge case and abort path, then
real-data assertions pinning the cohort to 16 / 9 / 1.
"""
from __future__ import annotations

import pandas as pd
import pytest

from src import cohort as C


# ---------------------------------------------------------------------------
# Synthetic edge cases
# ---------------------------------------------------------------------------


def _mk_cohort(rows: list[dict]) -> C.Stage1Cohort:
    """Build a Stage1Cohort directly, bypassing file IO."""
    frame = pd.DataFrame(rows)
    pos = frame[frame[C.DIRECTION_COLUMN] == "positive"].reset_index(drop=True)
    neg = frame[frame[C.DIRECTION_COLUMN] == "negative"].reset_index(drop=True)
    return C.Stage1Cohort(
        alt_positive=pos, alt_negative=neg,
        excluded_unknown=pd.DataFrame(columns=frame.columns),
        trace=C.FilterTrace("synthetic"),
    )


def _row(model_id, status, direction):
    return {
        "ModelID": model_id, "StrippedCellLineName": model_id,
        "alt_status": status, C.DIRECTION_COLUMN: direction,
    }


class TestAnchors:
    def test_passes_with_correct_anchors(self):
        c = _mk_cohort([
            _row("ACH-000364", "alt_gold_positive", "positive"),
            _row("ACH-000748", "alt_negative_high_conf", "negative"),
        ])
        C.assert_anchors(c)  # must not raise

    def test_aborts_when_anchor_direction_is_wrong(self):
        """A flipped U2OS label means the label join is broken. Must abort."""
        c = _mk_cohort([
            _row("ACH-000364", "alt_gold_positive", "negative"),   # WRONG
            _row("ACH-000748", "alt_negative_high_conf", "negative"),
        ])
        with pytest.raises(ValueError, match="ACH-000364"):
            C.assert_anchors(c)

    def test_aborts_when_anchor_is_missing(self):
        c = _mk_cohort([_row("ACH-000364", "alt_gold_positive", "positive")])
        with pytest.raises(ValueError, match="ACH-000748.*missing|missing.*ACH-000748"):
            C.assert_anchors(c)


class TestSizeAssertions:
    def test_aborts_on_size_mismatch_and_names_the_discrepancy(self):
        c = _mk_cohort([_row("ACH-000364", "alt_gold_positive", "positive")])
        with pytest.raises(ValueError, match="Cohort size mismatch"):
            C.assert_expected_sizes(c)

    def test_message_shows_both_expected_and_actual(self):
        c = _mk_cohort([_row("X", "alt_gold_positive", "positive")])
        with pytest.raises(ValueError) as exc:
            C.assert_expected_sizes(c)
        assert "expected" in str(exc.value) and "actual" in str(exc.value)
        assert "do not adapt the" in str(exc.value).lower().replace("\n", " ")


class TestGoldOnly:
    def test_drops_silver_only(self):
        frame = pd.DataFrame([
            _row("a", "alt_gold_positive", "positive"),
            _row("b", "alt_silver", "positive"),
            _row("c", "alt_gold_positive", "positive"),
        ])
        gold = C.gold_only(frame)
        assert len(gold) == 2
        assert set(gold["ModelID"]) == {"a", "c"}

    def test_keeps_high_conf_negatives(self):
        frame = pd.DataFrame([
            _row("a", "alt_negative_high_conf", "negative"),
            _row("b", "alt_silver", "negative"),
        ])
        assert list(C.gold_only(frame)["ModelID"]) == ["a"]


# ---------------------------------------------------------------------------
# Real data
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def real_cohort() -> C.Stage1Cohort:
    return C.build_cohort()


class TestRealCohort:
    def test_sizes_are_16_9_1(self, real_cohort):
        assert real_cohort.sizes == {
            "alt_positive": 16, "alt_negative": 9, "excluded_unknown": 1,
        }

    def test_all_assertions_pass(self, real_cohort):
        C.assert_all(real_cohort)

    def test_tier_breakdown(self, real_cohort):
        assert C.tier_counts(real_cohort) == C.EXPECTED_TIER_COUNTS

    def test_gold_only_alt_positive_is_10(self, real_cohort):
        assert len(C.gold_only(real_cohort.alt_positive)) == 10

    def test_filter_trace_records_every_step(self, real_cohort):
        steps = [s for s, _ in real_cohort.trace.steps]
        assert len(steps) == 8
        assert real_cohort.trace.steps[0][1] == 35   # Model.csv osteosarcoma rows
        assert real_cohort.trace.steps[1][1] == 26   # screened

    def test_u2os_and_saos2_are_alt_positive(self, real_cohort):
        names = set(real_cohort.alt_positive["StrippedCellLineName"])
        assert "U2OS" in names and "SAOS2" in names

    def test_no_line_appears_in_two_sets(self, real_cohort):
        pos = set(real_cohort.alt_positive["ModelID"])
        neg = set(real_cohort.alt_negative["ModelID"])
        exc = set(real_cohort.excluded_unknown["ModelID"])
        assert not (pos & neg) and not (pos & exc) and not (neg & exc)

    def test_every_screened_line_is_accounted_for(self, real_cohort):
        """25 labelled + 1 unknown == the 26 screened osteosarcoma lines."""
        assert sum(real_cohort.sizes.values()) == 26
