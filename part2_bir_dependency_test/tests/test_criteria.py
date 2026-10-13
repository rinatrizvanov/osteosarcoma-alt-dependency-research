"""Tests for the Section 8/9/10 verdict evaluator.

Every outcome bucket is exercised, including ``indeterminate`` — the evaluator
must refuse to round a result toward the nearest verdict.
"""
from __future__ import annotations

import pandas as pd
import pytest

from src import criteria as K
from src.dependency_stats import BIR_GENES, HR_GENES


def _per_gene(bir_q=0.01, bir_delta=-0.8, hr_q=0.9, hr_delta=-0.05):
    rows = []
    for gene in BIR_GENES:
        rows.append({"name": gene, "gene_set": "BIR", "q_value": bir_q,
                     "p_value": bir_q, "cliffs_delta": bir_delta})
    for gene in HR_GENES:
        rows.append({"name": gene, "gene_set": "HR", "q_value": hr_q,
                     "p_value": hr_q, "cliffs_delta": hr_delta})
    return pd.DataFrame(rows)


def _primary(p, delta):
    return {"p_value": p, "cliffs_delta": delta}


def _interaction(delta_pos=-0.8, delta_neg=-0.1, alt_neg_p=0.8):
    return {
        "delta_alt_positive": delta_pos,
        "delta_alt_negative": delta_neg,
        "alt_negative_p_value": alt_neg_p,
        "selectivity_as_predicted": abs(delta_pos) > abs(delta_neg),
    }


class TestSupported:
    def test_all_conditions_met(self):
        v = K.evaluate(
            primary=_primary(0.002, -0.75), per_gene=_per_gene(),
            r1_gold_only=_primary(0.09, -0.60), r2_median=_primary(0.01, -0.70),
            loo_deltas=[-0.70, -0.72, -0.68], interaction=_interaction(),
        )
        assert v.outcome == "supported"
        assert v.section == "Section 8"
        assert all(c.passed for c in v.conditions)

    def test_fails_if_gold_only_sensitivity_flips_sign(self):
        """R1 is your addition: a silver-driven result must not pass."""
        v = K.evaluate(
            primary=_primary(0.002, -0.75), per_gene=_per_gene(),
            r1_gold_only=_primary(0.5, +0.20), r2_median=_primary(0.01, -0.70),
            loo_deltas=[-0.70], interaction=_interaction(),
        )
        assert v.outcome != "supported"
        assert not next(c for c in v.conditions if c.key == "r1_gold_only").passed

    def test_fails_if_only_one_bir_gene_hits(self):
        pg = _per_gene()
        pg.loc[pg["name"] == "RAD52", "q_value"] = 0.5
        pg.loc[pg["name"] == "SMARCAL1", "q_value"] = 0.5
        v = K.evaluate(
            primary=_primary(0.002, -0.75), per_gene=pg,
            r1_gold_only=_primary(0.09, -0.60), r2_median=_primary(0.01, -0.70),
            loo_deltas=[-0.70], interaction=_interaction(),
        )
        assert not next(c for c in v.conditions if c.key == "bir_genes").passed

    def test_fails_if_an_hr_gene_outperforms_bir(self):
        """Specificity guard: the effect must not be a general HR-gene effect."""
        v = K.evaluate(
            primary=_primary(0.002, -0.75),
            per_gene=_per_gene(bir_delta=-0.5, hr_delta=-0.9),
            r1_gold_only=_primary(0.09, -0.60), r2_median=_primary(0.01, -0.70),
            loo_deltas=[-0.70], interaction=_interaction(),
        )
        cond = next(c for c in v.conditions if c.key == "bir_genes")
        assert not cond.passed
        assert "EXCEEDS" in cond.detail

    def test_fails_if_a_single_line_flips_the_sign(self):
        v = K.evaluate(
            primary=_primary(0.002, -0.75), per_gene=_per_gene(),
            r1_gold_only=_primary(0.09, -0.60), r2_median=_primary(0.01, -0.70),
            loo_deltas=[-0.70, +0.05, -0.68], interaction=_interaction(),
        )
        assert not next(c for c in v.conditions if c.key == "r7_loo").passed


class TestSelectivity:
    def test_passes_when_alt_negative_is_not_significant(self):
        """A noisy but non-significant ALT- delta must NOT block support.

        At the real ALT- group sizes a noise-only stratum reaches |delta|>=0.474
        about 26-29% of the time, so gating on the delta itself would false-fail
        roughly a quarter of true effects.
        """
        v = K.evaluate(
            primary=_primary(0.002, -0.75), per_gene=_per_gene(),
            r1_gold_only=_primary(0.09, -0.60), r2_median=_primary(0.01, -0.70),
            loo_deltas=[-0.70],
            interaction=_interaction(delta_pos=-0.75, delta_neg=-0.80, alt_neg_p=0.42),
        )
        assert next(c for c in v.conditions if c.key == "selectivity").passed
        assert v.outcome == "supported"

    def test_fails_when_alt_negative_shows_a_significant_same_direction_shift(self):
        v = K.evaluate(
            primary=_primary(0.002, -0.75), per_gene=_per_gene(),
            r1_gold_only=_primary(0.09, -0.60), r2_median=_primary(0.01, -0.70),
            loo_deltas=[-0.70],
            interaction=_interaction(delta_pos=-0.75, delta_neg=-0.90, alt_neg_p=0.01),
        )
        assert not next(c for c in v.conditions if c.key == "selectivity").passed
        assert v.outcome != "supported"

    def test_opposite_direction_in_alt_negative_does_not_contradict(self):
        v = K.evaluate(
            primary=_primary(0.002, -0.75), per_gene=_per_gene(),
            r1_gold_only=_primary(0.09, -0.60), r2_median=_primary(0.01, -0.70),
            loo_deltas=[-0.70],
            interaction=_interaction(delta_pos=-0.75, delta_neg=+0.95, alt_neg_p=0.01),
        )
        assert next(c for c in v.conditions if c.key == "selectivity").passed


class TestInconclusive:
    def test_large_effect_but_not_significant(self):
        v = K.evaluate(
            primary=_primary(0.12, -0.55), per_gene=_per_gene(bir_q=0.3),
            r1_gold_only=_primary(0.4, -0.45), r2_median=_primary(0.2, -0.50),
            loo_deltas=[-0.50, -0.58, -0.53], interaction=_interaction(),
        )
        assert v.outcome == "inconclusive"
        assert v.section == "Section 9"
        assert "Stage 2" in v.headline

    def test_requires_directional_consistency(self):
        """A large primary delta with a sign-flipping LOO is not 'inconclusive'."""
        v = K.evaluate(
            primary=_primary(0.12, -0.55), per_gene=_per_gene(bir_q=0.3),
            r1_gold_only=_primary(0.4, +0.45), r2_median=_primary(0.2, -0.50),
            loo_deltas=[-0.50, +0.10], interaction=_interaction(),
        )
        assert v.outcome == "indeterminate"


class TestFalsified:
    def test_null_result(self):
        v = K.evaluate(
            primary=_primary(0.74, -0.05),
            per_gene=_per_gene(bir_q=0.8, bir_delta=-0.05),
            r1_gold_only=_primary(0.8, -0.02), r2_median=_primary(0.7, 0.03),
            loo_deltas=[-0.05, -0.03], interaction=_interaction(delta_pos=-0.05),
        )
        assert v.outcome == "falsified"
        assert v.section == "Section 10"
        assert "do not approach stage 2" in v.headline.lower()

    def test_significant_in_the_opposite_direction(self):
        v = K.evaluate(
            primary=_primary(0.01, +0.80),
            per_gene=_per_gene(bir_q=0.9, bir_delta=+0.8),
            r1_gold_only=_primary(0.05, +0.7), r2_median=_primary(0.02, +0.75),
            loo_deltas=[+0.8], interaction=_interaction(delta_pos=+0.8),
        )
        assert v.outcome == "falsified"
        assert "OPPOSITE" in v.headline

    def test_not_falsified_if_a_bir_gene_still_hits(self):
        """A null composite with a surviving BIR gene is not a clean negative."""
        v = K.evaluate(
            primary=_primary(0.74, -0.05), per_gene=_per_gene(bir_q=0.02),
            r1_gold_only=_primary(0.8, -0.02), r2_median=_primary(0.7, 0.03),
            loo_deltas=[-0.05], interaction=_interaction(delta_pos=-0.05),
        )
        assert v.outcome == "indeterminate"


class TestIndeterminate:
    def test_refuses_to_round_toward_a_verdict(self):
        """p<0.05 and a medium delta, but R1 fails: no bucket fits."""
        v = K.evaluate(
            primary=_primary(0.03, -0.60), per_gene=_per_gene(),
            r1_gold_only=_primary(0.9, -0.10), r2_median=_primary(0.04, -0.55),
            loo_deltas=[-0.60], interaction=_interaction(),
        )
        assert v.outcome == "indeterminate"
        assert "do not round it" in v.headline.lower()
        assert "r1_gold_only" in v.headline


class TestReporting:
    def test_markdown_table_lists_every_condition(self):
        v = K.evaluate(
            primary=_primary(0.002, -0.75), per_gene=_per_gene(),
            r1_gold_only=_primary(0.09, -0.60), r2_median=_primary(0.01, -0.70),
            loo_deltas=[-0.70], interaction=_interaction(),
        )
        md = "\n".join(v.as_markdown())
        for cond in v.conditions:
            assert cond.description in md
        assert md.count("|") >= 5 * 4

    def test_as_dict_is_json_serialisable(self):
        import json
        v = K.evaluate(
            primary=_primary(0.002, -0.75), per_gene=_per_gene(),
            r1_gold_only=_primary(0.09, -0.60), r2_median=_primary(0.01, -0.70),
            loo_deltas=[-0.70], interaction=_interaction(),
        )
        assert json.loads(json.dumps(v.as_dict()))["outcome"] == "supported"


def test_min_detectable_delta_documents_the_real_bar():
    """The p<0.05 requirement, not SUPPORT_DELTA, is the binding constraint."""
    assert K.MIN_DETECTABLE_DELTA["5v11"] > abs(K.SUPPORT_DELTA)
    assert K.MIN_DETECTABLE_DELTA["8v8"] > abs(K.SUPPORT_DELTA)


class TestConcentrationComparison:
    """Section 6.4, pre-registered in prespec v1.1. Descriptive, never gating."""

    def test_flags_concentration_when_tertile_much_stronger(self):
        res = K.concentration_comparison(
            {"cliffs_delta": -0.40}, {"cliffs_delta": -0.90}
        )
        assert res["concentrates_in_extreme_loss"] is True
        assert res["abs_excess_tertile_over_median"] == pytest.approx(0.50)
        assert "concentrates" in res["interpretation"]

    def test_does_not_flag_within_threshold(self):
        res = K.concentration_comparison(
            {"cliffs_delta": -0.50}, {"cliffs_delta": -0.70}
        )
        assert res["concentrates_in_extreme_loss"] is False
        assert "No evidence of concentration" in res["interpretation"]

    def test_does_not_flag_a_stronger_effect_in_the_wrong_direction(self):
        """A large positive tertile delta is not concentration of the predicted effect."""
        res = K.concentration_comparison(
            {"cliffs_delta": -0.30}, {"cliffs_delta": +0.95}
        )
        assert res["concentrates_in_extreme_loss"] is False
        assert res["both_same_direction"] is False

    def test_unavailable_when_tertile_not_run(self):
        assert K.concentration_comparison({"cliffs_delta": -0.4}, None) == {"available": False}

    def test_never_changes_the_verdict(self):
        """The same primary must yield the same outcome whatever the tertile shows."""
        args = dict(
            primary=_primary(0.74, -0.05),
            per_gene=_per_gene(bir_q=0.8, bir_delta=-0.05),
            r1_gold_only=_primary(0.8, -0.02), r2_median=_primary(0.7, 0.03),
            loo_deltas=[-0.05, -0.03], interaction=_interaction(delta_pos=-0.05),
        )
        base = K.evaluate(**args, tertile=None)
        strong = K.evaluate(**args, tertile={"cliffs_delta": -0.95})
        assert base.outcome == strong.outcome == "falsified"
        assert strong.concentration["concentrates_in_extreme_loss"] is True

    def test_appears_in_the_markdown_report(self):
        v = K.evaluate(
            primary=_primary(0.12, -0.45), per_gene=_per_gene(bir_q=0.3),
            r1_gold_only=_primary(0.4, -0.40), r2_median=_primary(0.2, -0.45),
            loo_deltas=[-0.45, -0.43], interaction=_interaction(),
            tertile={"cliffs_delta": -0.90},
        )
        md = "\n".join(v.as_markdown())
        assert "Section 6.4 concentration comparison" in md
        assert "non-gating" in md
