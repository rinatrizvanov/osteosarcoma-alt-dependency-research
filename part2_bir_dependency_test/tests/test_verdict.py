"""The verdict partition must be exhaustive, mutually exclusive and stable.

The central test is :class:`TestPartitionIsExhaustive`, which enumerates the
whole (direction sign x p-bucket x |delta|-bucket x condition-state) grid and
asserts that every cell reaches exactly one named outcome and that no finite
input ever reaches ``indeterminate``. That is the property
``src.criteria.evaluate`` lacked and the reason Stage 1's own result had no
name.
"""
from __future__ import annotations

import math

import numpy as np
import pytest
from scipy import stats

from src import criteria as K
from src import verdict as V


def _cond(passed: bool, key: str = "c1") -> V.Condition:
    return V.Condition(key, f"condition {key}", passed, "synthetic")


# Grid spanning every boundary: below/at/above the null ceiling, below/at/above
# the support delta, both signs, and exactly zero.
DELTAS = [
    -0.95, -0.60, -0.45, -0.34, -0.33, -0.32, -0.25, -0.21, -0.20, -0.19,
    -0.05, 0.0, 0.05, 0.19, 0.20, 0.21, 0.25, 0.32, 0.33, 0.34, 0.45, 0.60, 0.95,
]
P_VALUES = [0.0, 0.0001, 0.01, 0.0499, 0.05, 0.0501, 0.2, 0.382, 0.9, 1.0]
SIGNS = [-1, 1]
CONDITION_STATES = [
    (),
    (_cond(True),),
    (_cond(False),),
    (_cond(True), _cond(False, "c2")),
    (_cond(True), _cond(True, "c2")),
]


class TestPartitionIsExhaustive:
    def test_every_cell_reaches_exactly_one_named_finite_outcome(self):
        unassigned = []
        seen = set()
        for sign in SIGNS:
            for d in DELTAS:
                for p in P_VALUES:
                    for conds in CONDITION_STATES:
                        v = V.evaluate_iteration(
                            delta=d, p_value=p, predicted_sign=sign,
                            conditions=conds,
                        )
                        if v.outcome not in V.FINITE_OUTCOMES:
                            unassigned.append((sign, d, p, len(conds), v.outcome))
                        seen.add(v.outcome)
        assert not unassigned, (
            f"{len(unassigned)} grid cells reached no named outcome; "
            f"first few: {unassigned[:5]}"
        )
        assert "indeterminate" not in seen

    def test_all_six_finite_outcomes_are_reachable(self):
        seen = set()
        for sign in SIGNS:
            for d in DELTAS:
                for p in P_VALUES:
                    for conds in CONDITION_STATES:
                        seen.add(V.evaluate_iteration(
                            delta=d, p_value=p, predicted_sign=sign,
                            conditions=conds,
                        ).outcome)
        assert seen == set(V.FINITE_OUTCOMES), (
            f"unreachable outcomes: {set(V.FINITE_OUTCOMES) - seen}"
        )

    def test_outcome_is_deterministic(self):
        for _ in range(3):
            v = V.evaluate_iteration(
                delta=0.281, p_value=0.382, predicted_sign=-1,
                conditions=(_cond(False),),
            )
            assert v.outcome == "counter_signal"

    @pytest.mark.parametrize("sign", SIGNS)
    def test_sign_flip_is_symmetric(self, sign):
        """Flipping both the prediction and the effect must give the same verdict."""
        for d in DELTAS:
            for p in P_VALUES:
                a = V.evaluate_iteration(
                    delta=d, p_value=p, predicted_sign=sign,
                    conditions=(_cond(True),),
                )
                b = V.evaluate_iteration(
                    delta=-d, p_value=p, predicted_sign=-sign,
                    conditions=(_cond(True),),
                )
                assert a.outcome == b.outcome
                assert a.aligned_delta == pytest.approx(b.aligned_delta)


class TestTruthTable:
    """The specific cell -> outcome mapping, stated explicitly."""

    @pytest.mark.parametrize("delta,p,conds,expected", [
        # clean null, either sign, any p
        (0.10, 0.9, (_cond(True),), "null"),
        (-0.10, 0.9, (_cond(True),), "null"),
        (0.0, 0.5, (_cond(True),), "null"),
        (0.199, 0.001, (_cond(True),), "null"),
        # predicted direction (sign -1 => negative delta confirms)
        (-0.80, 0.001, (_cond(True),), "supported"),
        (-0.80, 0.001, (_cond(False),), "partial_support"),
        (-0.80, 0.001, (), "partial_support"),
        (-0.80, 0.20, (_cond(True),), "inconclusive_underpowered"),
        (-0.25, 0.60, (_cond(True),), "inconclusive_underpowered"),
        (-0.45, 0.60, (_cond(False),), "inconclusive_underpowered"),
        # counter direction
        (0.80, 0.001, (_cond(True),), "refuted"),
        (0.281, 0.382, (_cond(False),), "counter_signal"),
        (0.25, 0.60, (_cond(True),), "counter_signal"),
        (0.95, 0.90, (_cond(True),), "counter_signal"),
    ])
    def test_cells(self, delta, p, conds, expected):
        v = V.evaluate_iteration(
            delta=delta, p_value=p, predicted_sign=-1, conditions=conds,
        )
        assert v.outcome == expected

    def test_alpha_boundary_is_strict(self):
        """p == alpha is NOT significant."""
        at = V.evaluate_iteration(
            delta=-0.80, p_value=0.05, predicted_sign=-1, conditions=(_cond(True),))
        below = V.evaluate_iteration(
            delta=-0.80, p_value=0.0499, predicted_sign=-1, conditions=(_cond(True),))
        assert at.outcome == "inconclusive_underpowered"
        assert below.outcome == "supported"

    def test_null_ceiling_boundary_is_strict(self):
        """|aligned| == ceiling is NOT a null."""
        at = V.evaluate_iteration(delta=-0.200, p_value=0.9, predicted_sign=-1)
        below = V.evaluate_iteration(delta=-0.199, p_value=0.9, predicted_sign=-1)
        assert at.outcome == "inconclusive_underpowered"
        assert below.outcome == "null"

    def test_only_supported_is_positive(self):
        for outcome_delta, outcome_p, conds, positive in [
            (-0.8, 0.001, (_cond(True),), True),
            (-0.8, 0.001, (_cond(False),), False),
            (-0.8, 0.20, (_cond(True),), False),
            (0.8, 0.001, (_cond(True),), False),
        ]:
            v = V.evaluate_iteration(
                delta=outcome_delta, p_value=outcome_p,
                predicted_sign=-1, conditions=conds,
            )
            assert v.is_positive is positive


class TestStage1HolesAreNowNamed:
    """The three unassigned regions of ``criteria.evaluate``, each now named."""

    def test_hole_1_stage1_actual_result(self):
        """delta=+0.281, p=0.382, prediction negative -> was indeterminate."""
        v = V.evaluate_iteration(
            delta=0.281, p_value=0.382, predicted_sign=-1,
            conditions=(_cond(False, "primary"), _cond(False, "r7_loo")),
        )
        assert v.outcome == "counter_signal"
        assert v.aligned_delta == pytest.approx(-0.281)
        assert not v.is_positive

    def test_hole_2_predicted_direction_mid_band(self):
        """Predicted direction, ns, 0.20 <= |delta| < 0.40 -> was indeterminate."""
        for d in (-0.21, -0.30, -0.399):
            v = V.evaluate_iteration(delta=d, p_value=0.5, predicted_sign=-1)
            assert v.outcome == "inconclusive_underpowered"

    def test_hole_3_large_predicted_effect_with_inconsistent_robustness(self):
        """Robustness failure can no longer drop a result out of the partition."""
        v = V.evaluate_iteration(
            delta=-0.70, p_value=0.5, predicted_sign=-1,
            conditions=(_cond(False, "r7_loo"), _cond(False, "r1_gold_only")),
        )
        assert v.outcome == "inconclusive_underpowered"
        assert v.attributes["failed_conditions"] == ["r7_loo", "r1_gold_only"]

    def test_indeterminate_is_reserved_for_missing_inputs(self):
        for d, p in [(float("nan"), 0.5), (0.3, float("nan")),
                     (float("inf"), 0.5), (None, 0.5)]:
            v = V.evaluate_iteration(delta=d, p_value=p, predicted_sign=-1)
            assert v.outcome == "indeterminate"
            assert v.rule == "missing-input"

    def test_bad_predicted_sign_raises(self):
        with pytest.raises(ValueError):
            V.evaluate_iteration(delta=-0.5, p_value=0.01, predicted_sign=0)


class TestExactMannWhitneyNull:
    @pytest.mark.parametrize("n,m", [(3, 4), (4, 5), (5, 5), (8, 8), (5, 11)])
    def test_two_sided_p_matches_scipy_exact(self, n, m):
        """Every attainable U, against scipy's own exact method."""
        rng = np.random.default_rng(0)
        for _ in range(12):
            vals = rng.permutation(n + m).astype(float)
            x, y = vals[:n], vals[n:]
            res = stats.mannwhitneyu(x, y, alternative="two-sided", method="exact")
            mine = V.exact_mwu_two_sided_p(int(res.statistic), n, m)
            assert mine == pytest.approx(float(res.pvalue), abs=1e-12)

    def test_counts_sum_to_binomial(self):
        from math import comb
        for n, m in [(3, 4), (5, 5), (8, 8)]:
            assert sum(V._u_counts(n, m)) == comb(n + m, n)

    def test_delta_from_u_endpoints(self):
        assert V.delta_from_u(0, 8, 8) == pytest.approx(-1.0)
        assert V.delta_from_u(64, 8, 8) == pytest.approx(1.0)
        assert V.delta_from_u(32, 8, 8) == pytest.approx(0.0)

    @pytest.mark.parametrize("key", list(K.MIN_DETECTABLE_DELTA))
    def test_reproduces_frozen_stage1_table(self, key):
        """Independent derivation must match the hand-computed frozen values."""
        n, m = (int(x) for x in key.split("v"))
        assert V.min_detectable_delta(n, m) == pytest.approx(
            K.MIN_DETECTABLE_DELTA[key], abs=5e-4
        )

    def test_tiny_groups_cannot_reach_significance(self):
        """2v2 has only 6 configurations; min two-sided p is 1/3."""
        assert math.isnan(V.min_detectable_delta(2, 2))
        assert "cannot produce a significant result" in V.power_note(2, 2)

    def test_power_note_quotes_the_boundary(self):
        assert "0.594" in V.power_note(8, 8)


class TestCampaignFDR:
    def test_matches_scipy_false_discovery_control(self):
        ps = [0.001, 0.02, 0.04, 0.2, 0.382, 0.9]
        mine = V.benjamini_hochberg(ps)
        theirs = stats.false_discovery_control(ps, method="bh")
        assert mine == pytest.approx(list(theirs), abs=1e-12)

    def test_nan_passes_through(self):
        q = V.benjamini_hochberg([0.01, float("nan"), 0.5])
        assert math.isnan(q[1])
        assert math.isfinite(q[0]) and math.isfinite(q[2])

    def test_single_iteration_q_equals_p(self):
        assert V.benjamini_hochberg([0.03])[0] == pytest.approx(0.03)

    def test_campaign_threshold_is_the_success_criterion(self):
        rows = V.campaign_fdr(
            [{"iteration": 1, "p_value": 0.04}, {"iteration": 2, "p_value": 0.30},
             {"iteration": 3, "p_value": 0.50}],
            q_threshold=0.100,
        )
        # 0.04 * 3/1 = 0.12 > 0.10: a nominally significant iteration does NOT
        # survive the campaign-level criterion. This is the whole point of 7(a).
        assert rows[0]["q_value_campaign"] == pytest.approx(0.12)
        assert rows[0]["survives_campaign_fdr"] is False

    def test_a_strong_result_does_survive(self):
        rows = V.campaign_fdr(
            [{"p_value": 0.001}, {"p_value": 0.4}, {"p_value": 0.9}],
            q_threshold=0.100,
        )
        assert rows[0]["survives_campaign_fdr"] is True

    def test_stage1_primary_would_not_survive(self):
        rows = V.campaign_fdr([{"p_value": 0.382}], q_threshold=0.100)
        assert rows[0]["survives_campaign_fdr"] is False
