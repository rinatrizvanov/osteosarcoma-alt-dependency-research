"""Tests for the iteration-03 module-coherence statistic.

The statistic is new code on the critical path of a pre-registered result, so it
is validated the same way ``cliffs_delta`` was in Stage 1: against an independent
implementation, plus the structural properties the pre-specification actually
relies on.

Two of those properties are load-bearing for the freeze and are asserted here
rather than assumed:

* the null distribution is centred on zero and **distribution-free given
  (n_lines, n_genes)** -- this is what lets iteration 03 pre-register its
  detection floor before reading any gene-effect value;
* the module imports nothing that can reach the outcome, so writing it could not
  have leaked the endpoint.
"""
from __future__ import annotations

import importlib.util
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pytest
from scipy import stats

ITER_DIR = (
    Path(__file__).resolve().parents[1]
    / "autoresearch"
    / "iteration_03_bir_module_coherence"
)
MODULE_PATH = ITER_DIR / "src" / "coherence.py"


def _load():
    spec = importlib.util.spec_from_file_location("iter03_coherence", MODULE_PATH)
    module = importlib.util.module_from_spec(spec)
    sys.modules["iter03_coherence"] = module
    spec.loader.exec_module(module)
    return module


coherence = _load()


class TestSpearmanAgreesWithScipy:
    """The pairwise engine must match scipy on real-valued, tied and edge inputs."""

    @pytest.mark.parametrize("seed", range(25))
    def test_matches_scipy_continuous(self, seed):
        rng = np.random.default_rng(seed)
        a = rng.normal(size=16)
        b = rng.normal(size=16)
        assert coherence._spearman(a, b) == pytest.approx(
            stats.spearmanr(a, b).statistic, abs=1e-12
        )

    @pytest.mark.parametrize("seed", range(10))
    def test_matches_scipy_with_ties(self, seed):
        rng = np.random.default_rng(seed)
        a = rng.integers(0, 4, size=16).astype(float)
        b = rng.integers(0, 3, size=16).astype(float)
        assert coherence._spearman(a, b) == pytest.approx(
            stats.spearmanr(a, b).statistic, abs=1e-12
        )

    def test_perfect_monotone_is_one(self):
        x = np.arange(16, dtype=float)
        assert coherence._spearman(x, x**3) == pytest.approx(1.0)
        assert coherence._spearman(x, -(x**3)) == pytest.approx(-1.0)


class TestMeanPairwiseRho:
    def test_is_mean_of_all_pairs(self):
        rng = np.random.default_rng(0)
        values = {g: rng.normal(size=16) for g in ("POLD3", "RAD52", "SMARCAL1")}
        result = coherence.mean_pairwise_rho(values)
        expected = np.mean(
            [
                stats.spearmanr(values[a], values[b]).statistic
                for a, b in combinations(values, 2)
            ]
        )
        assert result.statistic == pytest.approx(expected, abs=1e-12)
        assert len(result.pairwise) == 3
        assert result.n_lines == 16
        assert result.n_genes == 3

    def test_identical_profiles_give_one(self):
        x = np.linspace(-2, 0, 16)
        result = coherence.mean_pairwise_rho({"a": x, "b": x, "c": x})
        assert result.statistic == pytest.approx(1.0)
        assert result.all_pairs_positive
        assert result.min_pairwise == pytest.approx(1.0)

    def test_refuses_to_impute_missing_values(self):
        """Dropping lines is the caller's job, and must be recorded, not silent."""
        x = np.linspace(-1, 1, 8)
        y = x.copy()
        y[3] = np.nan
        with pytest.raises(ValueError, match="non-finite"):
            coherence.mean_pairwise_rho({"a": x, "b": y})

    def test_refuses_ragged_input(self):
        with pytest.raises(ValueError, match="ragged"):
            coherence.mean_pairwise_rho(
                {"a": np.zeros(8) + np.arange(8), "b": np.arange(7.0)}
            )

    def test_needs_two_genes(self):
        with pytest.raises(ValueError, match="at least two"):
            coherence.mean_pairwise_rho({"a": np.arange(8.0)})

    def test_reports_ties_so_the_fallback_can_fire_mechanically(self):
        tied = np.array([0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
        clean = np.arange(8.0)
        result = coherence.mean_pairwise_rho({"tied": tied, "clean": clean})
        assert result.n_ties_per_gene == {"tied": 1, "clean": 0}
        assert result.has_ties
        assert not coherence.mean_pairwise_rho(
            {"a": clean, "b": clean + 1}
        ).has_ties


class TestNullIsCentredAndDistributionFree:
    """The property the pre-registered detection floor rests on."""

    def test_null_is_centred_on_zero(self):
        null = coherence.null_distribution(16, 3, n_permutations=4000, seed=1)
        assert abs(float(null.mean())) < 0.01
        assert float(np.median(null)) == pytest.approx(0.0, abs=0.02)

    def test_statistic_depends_on_values_only_through_ranks(self):
        """The exact property the floor rests on.

        T is a function of within-gene ranks, so any strictly monotone
        transform of a gene's values -- and replacing its values by their ranks
        outright -- leaves T unchanged. This is why a null built from uniformly
        random rank vectors is the correct null whatever the marginals are, and
        hence why the floor could be computed before any value was read.
        """
        rng = np.random.default_rng(7)
        values = {g: rng.normal(size=14) for g in "abc"}
        observed = coherence.mean_pairwise_rho(values).statistic

        as_ranks = {g: coherence._ranks(v) for g, v in values.items()}
        assert coherence.mean_pairwise_rho(as_ranks).statistic == pytest.approx(
            observed, abs=1e-12
        )

        # Different strictly-increasing transform per gene, including one that
        # changes the scale by orders of magnitude.
        warped = {
            "a": np.exp(values["a"]),
            "b": values["b"] * 1e6 - 42.0,
            "c": np.arctan(values["c"]),
        }
        assert coherence.mean_pairwise_rho(warped).statistic == pytest.approx(
            observed, abs=1e-12
        )

    def test_null_quantiles_match_across_marginals(self):
        """Distribution-freeness in its true form: equality in distribution.

        Permuting three real, differently-shaped marginals gives the same null
        as permuting the ranks 1..n that :func:`null_distribution` uses. Not
        value-by-value -- the per-seed values differ, because each marginal
        starts from its own rank order -- but quantile-for-quantile, which is
        all the detection floor needs.
        """
        n, n_perm = 14, 6000
        rng = np.random.default_rng(11)
        marginals = {
            "gaussian": {g: rng.normal(size=n) for g in "abc"},
            "skewed": {g: rng.exponential(size=n) ** 3 for g in "abc"},
            "near-degenerate": {
                g: rng.normal(loc=-1.5, scale=0.01, size=n) for g in "abc"
            },
        }
        reference = coherence.null_distribution(n, 3, n_permutations=n_perm, seed=5)
        probs = [0.05, 0.25, 0.5, 0.75, 0.95]
        reference_q = np.quantile(reference, probs)

        for name, values in marginals.items():
            r = np.random.default_rng(5)
            draws = np.array(
                [
                    coherence.mean_pairwise_rho(
                        {g: r.permutation(v) for g, v in values.items()}
                    ).statistic
                    for _ in range(n_perm)
                ]
            )
            assert np.allclose(
                np.quantile(draws, probs), reference_q, atol=0.02
            ), f"{name} null diverges from the rank null"

    def test_null_spread_shrinks_with_more_lines_and_more_genes(self):
        wide = coherence.null_distribution(9, 3, n_permutations=3000, seed=2).std()
        narrow = coherence.null_distribution(25, 3, n_permutations=3000, seed=2).std()
        more_genes = coherence.null_distribution(
            16, 5, n_permutations=3000, seed=2
        ).std()
        fewer_genes = coherence.null_distribution(
            16, 2, n_permutations=3000, seed=2
        ).std()
        assert narrow < wide
        assert more_genes < fewer_genes

    def test_is_reproducible_under_the_fixed_seed(self):
        a = coherence.null_distribution(12, 3, n_permutations=500, seed=99)
        b = coherence.null_distribution(12, 3, n_permutations=500, seed=99)
        assert np.array_equal(a, b)


class TestPermutationP:
    def test_never_reports_zero(self):
        null = coherence.null_distribution(16, 3, n_permutations=1000, seed=3)
        assert coherence.permutation_p_two_sided(1.0, null) == pytest.approx(
            1 / 1001
        )

    def test_is_two_sided(self):
        null = coherence.null_distribution(16, 3, n_permutations=2000, seed=4)
        assert coherence.permutation_p_two_sided(
            0.45, null
        ) == pytest.approx(coherence.permutation_p_two_sided(-0.45, null))

    def test_zero_effect_is_not_significant(self):
        null = coherence.null_distribution(16, 3, n_permutations=2000, seed=5)
        assert coherence.permutation_p_two_sided(0.0, null) == pytest.approx(1.0)


class TestDetectionFloor:
    def test_frozen_floor_matches_the_prespec(self):
        """The number written into the frozen prespec. If this test fails, the
        prespec's stated floor and the code have diverged."""
        assert coherence.detection_floor(16, 3) == pytest.approx(0.2892, abs=5e-4)

    def test_floor_falls_with_more_lines(self):
        assert coherence.detection_floor(25, 3) < coherence.detection_floor(16, 3)
        assert coherence.detection_floor(16, 3) < coherence.detection_floor(9, 3)

    def test_floor_falls_with_more_genes(self):
        assert coherence.detection_floor(16, 5) < coherence.detection_floor(16, 3)
        assert coherence.detection_floor(16, 3) < coherence.detection_floor(16, 2)

    def test_an_effect_at_the_floor_is_significant_and_below_it_is_not(self):
        null = coherence.null_distribution(16, 3)
        floor = coherence.detection_floor(16, 3)
        assert coherence.permutation_p_two_sided(floor + 1e-6, null) < 0.05
        assert coherence.permutation_p_two_sided(floor - 0.05, null) >= 0.05


class TestFeedsTheCampaignEvaluatorUnmodified:
    """T is on Cliff's-delta scale, so it must need no special-casing."""

    def test_incoherence_maps_to_null(self):
        from src.verdict import evaluate_iteration

        verdict = evaluate_iteration(delta=0.05, p_value=0.72, predicted_sign=+1)
        assert verdict.outcome == "null"

    def test_the_underpowered_band_is_named_not_rounded(self):
        from src.verdict import evaluate_iteration

        # 0.200 <= |T| < 0.289: above the clean-null ceiling, below this n's floor.
        verdict = evaluate_iteration(delta=0.25, p_value=0.11, predicted_sign=+1)
        assert verdict.outcome == "inconclusive_underpowered"

    def test_anticorrelation_is_counter_signal(self):
        from src.verdict import evaluate_iteration

        verdict = evaluate_iteration(delta=-0.30, p_value=0.33, predicted_sign=+1)
        assert verdict.outcome == "counter_signal"


def test_module_cannot_reach_the_outcome():
    """coherence.py must not be able to load gene-effect values.

    Writing the statistic and its floor before the freeze is only legitimate if
    doing so could not have read the endpoint. That is a property of the file, so
    it is tested rather than asserted in prose.
    """
    source = MODULE_PATH.read_text()
    for forbidden in (
        "dependency_stats",
        "CRISPRGeneEffect",
        "read_csv",
        "pandas",
        "open(",
    ):
        assert forbidden not in source, f"coherence.py references {forbidden!r}"
