"""Module-coherence statistic and its distribution-free null, for iteration 03.

Why this module contains no data loading at all
-----------------------------------------------
Iteration 03's outcome variable is the **cross-gene correlation structure** of
Chronos gene effect among POLD3, RAD52 and SMARCAL1. Anti-bias rule (b) of
``autoresearch/PROGRAM.md`` requires that anything which could be tuned against
the outcome be fixed before the outcome is read. The detection floor of the
primary test is exactly such a thing, so it is computed here from rank
permutations alone -- it depends only on the number of lines ``n`` and the
number of genes ``k``, never on a single gene-effect value.

That is a real property of the statistic, not a convenience. ``mean_pairwise_rho``
is Spearman's rho averaged over gene pairs, i.e. a function of within-gene
*ranks*. Permuting each gene's ranks independently across lines destroys
cross-gene association while leaving every marginal distribution untouched, and
the resulting null distribution of the statistic is therefore distribution-free
given ``(n, k)`` and the absence of ties. So the critical value can be, and is,
pre-registered.

Ties would break that guarantee (Spearman on averaged ranks has a
tie-pattern-dependent null). :func:`mean_pairwise_rho` reports the tie count it
saw so the prespec's tie fallback can be triggered mechanically rather than by
judgement.
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations

import numpy as np

#: Permutation count for the primary test. Fixed here, not in the notebook, so
#: it cannot be raised until a p-value crosses a threshold.
N_PERMUTATIONS = 10_000

#: Seed for every permutation and bootstrap draw in this iteration. Fixed so
#: the notebook is byte-reproducible under the sha256 idempotency check.
SEED = 20261012


@dataclass(frozen=True)
class Coherence:
    """Observed module coherence and the bookkeeping needed to judge it."""

    statistic: float
    pairwise: dict[tuple[str, str], float]
    n_lines: int
    n_genes: int
    n_ties_per_gene: dict[str, int]

    @property
    def all_pairs_positive(self) -> bool:
        return bool(self.pairwise) and all(r > 0 for r in self.pairwise.values())

    @property
    def min_pairwise(self) -> float:
        return min(self.pairwise.values()) if self.pairwise else float("nan")

    @property
    def has_ties(self) -> bool:
        return any(c > 0 for c in self.n_ties_per_gene.values())


def _ranks(x: np.ndarray) -> np.ndarray:
    """Average ranks, so the tie case degrades rather than silently misbehaving."""
    order = np.argsort(x, kind="stable")
    ranks = np.empty(len(x), dtype=float)
    ranks[order] = np.arange(1, len(x) + 1, dtype=float)
    # average ranks within tied groups
    _, inverse, counts = np.unique(x, return_inverse=True, return_counts=True)
    for idx in np.flatnonzero(counts > 1):
        mask = inverse == idx
        ranks[mask] = ranks[mask].mean()
    return ranks


def _spearman(a: np.ndarray, b: np.ndarray) -> float:
    ra, rb = _ranks(a), _ranks(b)
    ra = ra - ra.mean()
    rb = rb - rb.mean()
    denom = np.sqrt((ra @ ra) * (rb @ rb))
    if denom == 0:
        return float("nan")
    return float((ra @ rb) / denom)


def mean_pairwise_rho(values: dict[str, np.ndarray]) -> Coherence:
    """Mean Spearman rho over all gene pairs -- the module-coherence statistic.

    ``values`` maps gene symbol to that gene's effect across the *same ordered*
    set of cell lines. Lines with any missing value must be dropped by the
    caller before this is called, and the drop recorded; this function does not
    silently impute.
    """
    genes = list(values)
    if len(genes) < 2:
        raise ValueError("coherence needs at least two genes")
    lengths = {len(v) for v in values.values()}
    if len(lengths) != 1:
        raise ValueError(f"ragged inputs: lengths {lengths}")
    arrays = {g: np.asarray(v, dtype=float) for g, v in values.items()}
    for g, arr in arrays.items():
        if not np.isfinite(arr).all():
            raise ValueError(f"{g} carries non-finite values; drop those lines first")

    pairwise = {
        (a, b): _spearman(arrays[a], arrays[b]) for a, b in combinations(genes, 2)
    }
    ties = {
        g: int(len(arr) - len(np.unique(arr))) for g, arr in arrays.items()
    }
    return Coherence(
        statistic=float(np.mean(list(pairwise.values()))),
        pairwise=pairwise,
        n_lines=int(next(iter(lengths))),
        n_genes=len(genes),
        n_ties_per_gene=ties,
    )


def null_distribution(
    n_lines: int, n_genes: int, *, n_permutations: int = N_PERMUTATIONS, seed: int = SEED
) -> np.ndarray:
    """Null distribution of the statistic under independent rank permutation.

    Distribution-free given ``(n_lines, n_genes)`` in the absence of ties, which
    is what lets the detection floor be pre-registered. Uses ranks ``1..n``
    directly: under the null any gene's rank vector is a uniform random
    permutation, independent across genes.
    """
    rng = np.random.default_rng(seed)
    base = np.arange(1, n_lines + 1, dtype=float)
    out = np.empty(n_permutations, dtype=float)
    for i in range(n_permutations):
        cols = {str(g): rng.permutation(base) for g in range(n_genes)}
        out[i] = mean_pairwise_rho(cols).statistic
    return out


def permutation_p_two_sided(observed: float, null: np.ndarray) -> float:
    """Two-sided permutation p with the conventional +1 correction.

    The +1 keeps the p-value from ever being reported as exactly zero, which
    would overstate the evidence a finite permutation set can carry.
    """
    extreme = int(np.sum(np.abs(null) >= abs(observed)))
    return float((extreme + 1) / (len(null) + 1))


def detection_floor(
    n_lines: int,
    n_genes: int,
    *,
    alpha: float = 0.05,
    n_permutations: int = N_PERMUTATIONS,
    seed: int = SEED,
) -> float:
    """Smallest \\|statistic\\| that can reach two-sided p < alpha at these sizes.

    The analogue of Stage 1 prespec Section 5a, and of
    ``src.verdict.min_detectable_delta`` for the two-group case: the honest
    answer to "what could this test ever have detected?", computed before the
    outcome is read.
    """
    null = null_distribution(
        n_lines, n_genes, n_permutations=n_permutations, seed=seed
    )
    return float(np.quantile(np.abs(null), 1.0 - alpha))
