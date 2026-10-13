"""Tests for the Section 5/6 statistics.

The effect-size functions are validated against closed-form identities rather
than remembered constants: Cliff's delta is pinned to the Mann-Whitney U
identity ``delta = 2U/(nm) - 1``, which holds exactly in the absence of ties.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from scipy import stats

from src.dependency_stats import (
    BIR_GENES,
    HR_GENES,
    SIX_GENES,
    bh_fdr,
    bootstrap_delta_ci,
    cliffs_delta,
    compare_groups,
    composite_score,
    delta_difference_permutation,
    delta_magnitude,
    hodges_lehmann,
    per_gene_comparisons,
    standardize_within_line,
    strip_gene_symbol,
)


class TestCliffsDelta:
    def test_perfect_separation(self):
        assert cliffs_delta([1, 2, 3], [4, 5, 6]) == -1.0
        assert cliffs_delta([4, 5, 6], [1, 2, 3]) == 1.0

    def test_identical(self):
        assert cliffs_delta([1, 2, 3], [1, 2, 3]) == 0.0

    def test_hand_computed_asymmetric(self):
        # x=[1,4] vs y=[3]: one pair below, one above -> 0
        assert cliffs_delta([1, 4], [3]) == 0.0
        assert cliffs_delta([1, 2], [3]) == -1.0
        # x=[1,5] vs y=[2,3,4]: 3 pairs below, 3 above -> 0
        assert cliffs_delta([1, 5], [2, 3, 4]) == 0.0

    def test_matches_mannwhitney_u_identity(self):
        """delta = 2U/(nm) - 1 exactly, for untied samples."""
        rng = np.random.default_rng(0)
        for _ in range(100):
            n, m = int(rng.integers(3, 12)), int(rng.integers(3, 12))
            x, y = rng.normal(size=n), rng.normal(0.5, 1.0, size=m)
            u = stats.mannwhitneyu(x, y, alternative="two-sided", method="exact").statistic
            assert np.isclose(2 * u / (n * m) - 1, cliffs_delta(x, y))

    def test_empty_group_is_nan(self):
        assert np.isnan(cliffs_delta([], [1, 2]))

    def test_magnitude_labels(self):
        assert delta_magnitude(0.05) == "negligible"
        assert delta_magnitude(0.20) == "small"
        assert delta_magnitude(0.40) == "medium"
        assert delta_magnitude(0.60) == "large"
        # Sign must not affect magnitude.
        assert delta_magnitude(-0.60) == "large"


class TestHodgesLehmann:
    def test_pairwise_median(self):
        # x=[10,12] vs y=[1] -> differences [9, 11] -> median 10
        assert hodges_lehmann([10, 12], [1]) == 10.0

    def test_sign_follows_shift(self):
        assert hodges_lehmann([1, 2, 3], [11, 12, 13]) == -10.0

    def test_in_chronos_units(self):
        """A constant shift of c must return exactly -c."""
        rng = np.random.default_rng(3)
        base = rng.normal(size=9)
        assert np.isclose(hodges_lehmann(base - 0.5, base), -0.5)


class TestCompareGroups:
    def test_p_value_matches_scipy_exact(self):
        rng = np.random.default_rng(5)
        a, b = rng.normal(-1, 0.3, 7), rng.normal(0, 0.3, 8)
        res = compare_groups(a, b)
        expected = stats.mannwhitneyu(a, b, alternative="two-sided", method="exact")
        assert res.p_value == pytest.approx(expected.pvalue)
        assert res.mwu_method == "exact"

    def test_falls_back_to_asymptotic_under_ties(self):
        """Exact p-values are invalid with ties; the fallback must be recorded."""
        res = compare_groups([1.0, 2.0, 3.0], [3.0, 4.0, 5.0])
        assert res.mwu_method == "asymptotic"

    def test_direction_flag_tracks_predicted_direction(self):
        lower = compare_groups([-1, -2, -3], [1, 2, 3])
        assert lower.direction_as_predicted is True
        assert lower.cliffs_delta == -1.0
        higher = compare_groups([1, 2, 3], [-1, -2, -3])
        assert higher.direction_as_predicted is False

    def test_bootstrap_ci_is_deterministic_and_brackets_delta(self):
        rng = np.random.default_rng(9)
        a, b = rng.normal(-0.8, 0.3, 8), rng.normal(0, 0.3, 8)
        first = bootstrap_delta_ci(a, b, n_boot=400)
        second = bootstrap_delta_ci(a, b, n_boot=400)
        assert first == second
        assert first[0] <= cliffs_delta(a, b) <= first[1]


class TestComposite:
    def test_is_bir_mean_minus_hr_mean(self, toy_effects):
        score = composite_score(toy_effects)
        manual = (
            toy_effects[list(BIR_GENES)].mean(axis=1)
            - toy_effects[list(HR_GENES)].mean(axis=1)
        )
        assert np.allclose(score, manual)

    def test_planted_bir_shift_is_negative(self, toy_effects):
        """The fixture plants extra BIR essentiality in the first 8 lines."""
        score = composite_score(toy_effects)
        assert score.iloc[:8].median() < score.iloc[8:].median()
        assert cliffs_delta(score.iloc[:8], score.iloc[8:]) < 0

    def test_missing_gene_raises(self, toy_effects):
        with pytest.raises(KeyError, match="POLD3"):
            composite_score(toy_effects.drop(columns=["POLD3"]))

    def test_insensitive_to_line_level_offset(self, toy_effects):
        """A per-line constant shift cancels in the BIR-HR difference.

        This is the stated reason the composite is the primary endpoint.
        """
        offset = pd.Series(
            np.linspace(-0.5, 0.5, len(toy_effects)), index=toy_effects.index
        )
        shifted = toy_effects.add(offset, axis=0)
        assert np.allclose(composite_score(toy_effects), composite_score(shifted))


class TestPerGene:
    def test_bh_correction_across_six_genes(self, toy_effects):
        ids = list(toy_effects.index)
        frame = per_gene_comparisons(toy_effects, ids[:8], ids[8:])
        assert len(frame) == 6
        assert set(frame["name"]) == set(SIX_GENES)
        assert (frame["q_value"] >= frame["p_value"] - 1e-12).all()
        assert np.allclose(frame["q_value"], bh_fdr(frame["p_value"]))

    def test_gene_set_labels(self, toy_effects):
        ids = list(toy_effects.index)
        frame = per_gene_comparisons(toy_effects, ids[:8], ids[8:]).set_index("name")
        assert all(frame.loc[g, "gene_set"] == "BIR" for g in BIR_GENES)
        assert all(frame.loc[g, "gene_set"] == "HR" for g in HR_GENES)

    def test_planted_signal_is_in_bir_genes_only(self):
        """Constructed, not random: BIR genes separate perfectly, HR genes tie by design."""
        ids = [f"L{i:02d}" for i in range(16)]
        frame = pd.DataFrame(index=ids, columns=list(SIX_GENES), dtype=float)
        for gene in BIR_GENES:
            frame[gene] = [-1.0] * 8 + [-0.1] * 8          # perfect separation
        for gene in HR_GENES:
            frame[gene] = list(np.linspace(-0.5, 0.5, 8)) * 2  # identical in both arms
        res = per_gene_comparisons(frame, ids[:8], ids[8:]).set_index("name")
        for gene in BIR_GENES:
            assert res.loc[gene, "cliffs_delta"] == -1.0
        for gene in HR_GENES:
            assert res.loc[gene, "cliffs_delta"] == 0.0


class TestStandardization:
    def test_zscores_per_row(self):
        frame = pd.DataFrame(
            [[1.0, 2.0, 3.0], [10.0, 20.0, 30.0]], index=["a", "b"], columns=list("xyz")
        )
        z = standardize_within_line(frame)
        assert np.allclose(z.mean(axis=1), 0.0)
        # Both rows are the same shape up to scale, so z-scores match.
        assert np.allclose(z.loc["a"], z.loc["b"])


class TestInteraction:
    def test_detects_effect_present_only_in_alt_positive(self):
        """Constructed so the ALT- stratum has delta exactly 0, not merely small.

        Random noise at the real ALT- group sizes (n=4 vs 5) produces |delta|
        values up to ~0.6 by chance, so a random fixture would be testing the
        RNG rather than the statistic. See
        ``test_alt_negative_delta_is_unstable_at_real_n`` for that property,
        measured deliberately.
        """
        res = delta_difference_permutation(
            [-1.0] * 6, [-0.1] * 10,                      # ALT+: perfect separation
            [-0.5, -0.3, -0.1, 0.1], [-0.5, -0.3, -0.1, 0.1, 0.3],  # ALT-: overlapping
            n_perm=400,
        )
        assert res["delta_alt_positive"] == -1.0
        assert abs(res["delta_alt_negative"]) < 0.3
        assert res["selectivity_as_predicted"] is True

    def test_alt_negative_delta_is_unstable_at_real_n(self):
        """Documents a real limitation of the Section 6.2 negative control.

        With the actual ALT- group sizes, Cliff's delta computed on pure noise
        is wide enough that the Section 8.5 selectivity condition can fail by
        chance alone. This test pins that fact so it cannot be forgotten when
        interpreting a failed selectivity check.
        """
        rng = np.random.default_rng(23)
        deltas = [
            cliffs_delta(rng.normal(size=4), rng.normal(size=5)) for _ in range(2000)
        ]
        assert np.quantile(np.abs(deltas), 0.95) > 0.6

    def test_no_selectivity_when_effect_is_in_both(self):
        rng = np.random.default_rng(17)
        res = delta_difference_permutation(
            rng.normal(-1.0, 0.2, 6), rng.normal(0.0, 0.2, 10),
            rng.normal(-1.0, 0.2, 4), rng.normal(0.0, 0.2, 5),
            n_perm=400,
        )
        assert res["selectivity_as_predicted"] is False

    def test_deterministic(self):
        rng = np.random.default_rng(19)
        args = (rng.normal(-1, .2, 6), rng.normal(0, .2, 10),
                rng.normal(0, .2, 4), rng.normal(0, .2, 5))
        assert (delta_difference_permutation(*args, n_perm=200)["p_value"]
                == delta_difference_permutation(*args, n_perm=200)["p_value"])


def test_strip_gene_symbol():
    assert strip_gene_symbol("POLD3 (10714)") == "POLD3"
    assert strip_gene_symbol("SMARCAL1 (50485)") == "SMARCAL1"
