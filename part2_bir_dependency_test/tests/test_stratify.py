"""Tests for the Section 4 split rule, especially the minimum-cell guard.

The guard is the one place where the analysis is allowed to change its own
primary rule, so its boundary behaviour is pinned precisely.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from src import stratify as S


def _series(values, prefix="A"):
    return pd.Series(
        np.asarray(values, dtype=float),
        index=[f"{prefix}{i:02d}" for i in range(len(values))],
    )


class TestTertileCut:
    def test_is_the_one_third_quantile(self):
        values = np.arange(30, dtype=float)
        assert S.tertile_cut(values) == pytest.approx(np.quantile(values, 1 / 3))

    def test_cut_for_zero_to_nine(self):
        assert S.tertile_cut(np.arange(10, dtype=float)) == pytest.approx(3.0)


class TestMinimumCellGuard:
    """The cut is estimated on `reference`, then applied to `alt`."""

    def _run(self, ref_values, alt_values):
        ref = _series(ref_values, "R")
        alt = _series(alt_values, "A")
        return S.stratify(
            alt.index.tolist(), reference_ids=ref.index.tolist(),
            expression=pd.concat([ref, alt]),
        )

    def test_does_not_fire_with_four_in_each_arm(self):
        # cut for 0..29 is 9.667; four alt values below it
        res = self._run(np.arange(30), [0, 1, 2, 3, 20, 21, 22, 23, 24, 25])
        assert res.guard_triggered is False
        assert res.rule == "tertile_os_cohort"
        assert res.sizes == {S.LOW_LABEL: 4, S.INTACT_LABEL: 6}

    def test_fires_with_three_in_the_low_arm(self):
        res = self._run(np.arange(30), [0, 1, 2, 20, 21, 22, 23, 24, 25, 26])
        assert res.guard_triggered is True
        assert res.rule == "median_within_alt_pos"
        assert "n=3" in res.guard_reason

    def test_fires_when_the_intact_arm_is_too_small(self):
        res = self._run(np.arange(30), [0, 1, 2, 3, 4, 5, 6, 7, 8, 25])
        assert res.guard_triggered is True
        assert res.rule == "median_within_alt_pos"

    def test_fallback_produces_a_balanced_split(self):
        res = self._run(np.arange(30), [0, 1, 2, 20, 21, 22, 23, 24, 25, 26])
        assert res.sizes == {S.LOW_LABEL: 5, S.INTACT_LABEL: 5}

    def test_guard_reason_is_recorded_for_the_report(self):
        res = self._run(np.arange(30), [0, 1, 2, 20, 21, 22, 23, 24, 25, 26])
        assert "NOT ESTIMABLE" in res.guard_reason
        assert "pre-specified" in res.guard_reason

    def test_realistic_case_uses_the_tertile_rule(self, synthetic_expression):
        alt_ids = [i for i in synthetic_expression.index if i.startswith("A")]
        ref_ids = [i for i in synthetic_expression.index if i.startswith("R")]
        res = S.stratify(alt_ids, reference_ids=ref_ids, expression=synthetic_expression)
        assert res.rule == "tertile_os_cohort"
        assert res.reference_n == 34
        assert sum(res.sizes.values()) == 16


class TestCutpointIsOutcomeIndependent:
    def test_cut_depends_only_on_the_reference_set(self, synthetic_expression):
        """Changing which lines are ALT+ must not move the cutpoint.

        This is the property that makes the split rule un-tunable.
        """
        ref_ids = [i for i in synthetic_expression.index if i.startswith("R")]
        alt_ids = [i for i in synthetic_expression.index if i.startswith("A")]
        full = S.stratify(alt_ids, reference_ids=ref_ids, expression=synthetic_expression)
        subset = S.stratify(
            alt_ids[:10], reference_ids=ref_ids, expression=synthetic_expression
        )
        assert full.cut == subset.cut


class TestSensitivitySplitters:
    def test_median_split_is_balanced(self):
        e = _series(np.arange(16))
        counts = S.stratify_median(e.index.tolist(), e).value_counts().to_dict()
        assert counts == {S.LOW_LABEL: 8, S.INTACT_LABEL: 8}

    def test_tertile_within_gives_five_eleven(self):
        e = _series(np.arange(16))
        counts = S.stratify_tertile_within(e.index.tolist(), e).value_counts().to_dict()
        assert counts == {S.INTACT_LABEL: 11, S.LOW_LABEL: 5}


class TestCopyNumberLoader:
    """Regression tests for the CN loader.

    `OmicsCNGeneWGS.csv` shares the expression matrix's layout (SequencingID
    rows with ModelID and IsDefaultEntryForModel as columns), not the CRISPR
    matrix's ModelID row index. An earlier version assumed the CRISPR layout,
    which would have silently returned all-NaN copy number.
    """

    def test_returns_values_not_all_nan(self):
        ids = S.osteosarcoma_rna_models()
        cn = S.load_wrap53_cn(ids)
        assert len(cn) == len(ids)
        assert cn.notna().sum() >= 30, "CN came back mostly NaN -- loader layout bug"

    def test_index_is_model_ids(self):
        ids = S.osteosarcoma_rna_models()[:5]
        cn = S.load_wrap53_cn(ids)
        assert list(cn.index) == ids
        assert all(i.startswith("ACH-") for i in cn.index)

    def test_one_row_per_model(self):
        ids = S.osteosarcoma_rna_models()
        cn = S.load_wrap53_cn(ids)
        assert not cn.index.duplicated().any()

    def test_copy_number_concordant_with_expression_at_the_extremes(self):
        """Near-zero WRAP53 CN must coincide with near-zero WRAP53 expression.

        This is the independent check that the expression values are not
        antisense-contaminated: contamination would break the concordance.
        """
        ids = S.osteosarcoma_rna_models()
        cn = S.load_wrap53_cn(ids)
        expr = S.load_expression(ids, [S.TCAB1_GENE])[S.TCAB1_GENE]
        near_zero_cn = cn.index[cn < 0.01]
        assert len(near_zero_cn) >= 1
        assert (expr.loc[near_zero_cn] < 0.5).all()
