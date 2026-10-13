"""Iteration 03 analysis: module coherence of the Stage 1 BIR endpoint.

Implements `prespec.md`, frozen at `b24a8ed`. Every threshold, seed, gene set and
tolerance here is quoted from that document; none is chosen in this file.

The coherence statistic and its pre-registered null live in `src/coherence.py`,
which was committed with the draft prespec (`cf73508`) and deliberately contains
no path to the outcome. This module is where the outcome is finally read.
"""
from __future__ import annotations

import importlib.util
import sys
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd

ITER_DIR = Path(__file__).resolve().parent


def _load_coherence():
    """Load the iteration's local coherence module by path.

    Loaded by path rather than imported, because the repository root already has
    a `src` package and the iteration has its own; a plain `import src.coherence`
    resolves to the wrong one.
    """
    spec = importlib.util.spec_from_file_location(
        "iter03_coherence", ITER_DIR / "src" / "coherence.py"
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules["iter03_coherence"] = module
    spec.loader.exec_module(module)
    return module


coherence = _load_coherence()

# --- Gene sets, quoted from the frozen prespec -----------------------------

#: Section 4 primary. The Stage 1 BIR module.
BIR_TRIPLET = ("POLD3", "RAD52", "SMARCAL1")

#: Section 6 primary gate. Contains POLD3, which is the point.
POL_DELTA = ("POLD1", "POLD2", "POLD3")

#: Section 6 secondary control, reported without a gate threshold.
BTR = ("BLM", "TOP3A", "RMI1", "RMI2")

#: Section 6 third module. Audited, never used as a neutral comparator.
HR_MODULE = ("BRCA1", "BRCA2", "RAD51")

#: Section 7 S4. Batenburg's module: MUS81, not SMARCAL1.
BATENBURG_BIR = ("MUS81", "RAD52", "POLD3")

#: Section 7 S4b extended fork-restart set, k = 5.
EXTENDED_RESTART = ("MUS81", "RAD52", "POLD3", "SMARCAL1", "PIF1")

#: The six Stage 1 genes, for the Section 7 S2a variance decomposition.
STAGE1_BIR = ("POLD3", "RAD52", "SMARCAL1")
STAGE1_HR = ("BRCA1", "BRCA2", "RAD51")

# --- Thresholds, quoted from the frozen prespec ----------------------------

GATE_MIN_T = 0.40           # Section 6
MIN_PAIR_RHO = 0.20         # Section 7 support conditions
MATCHED_N_TRIPLETS = 2000   # Section 5
MATCHED_MEAN_TOL = 0.15     # Section 5
MATCHED_SD_TOL = 0.05       # Section 5
MATCHED_MIN_CANDIDATES = 50  # Section 5
MATCHED_PERCENTILE = 95     # Section 7 `exceeds_matched_triplets`
MIN_LINES = 10              # Section 3 minimum-cell guard
N_BOOTSTRAP = 2000          # Section 7 S1
SEED = coherence.SEED


@dataclass(frozen=True)
class ModuleResult:
    """One module's coherence, its p, and the bookkeeping behind both."""

    label: str
    genes: tuple[str, ...]
    n_lines: int
    statistic: float
    p_value: float
    pairwise: dict[tuple[str, str], float]
    floor: float
    has_ties: bool

    @property
    def n_genes(self) -> int:
        return len(self.genes)

    @property
    def all_pairs_positive(self) -> bool:
        return all(r > 0 for r in self.pairwise.values())

    @property
    def min_pair_rho(self) -> float:
        return min(self.pairwise.values())

    def row(self) -> dict:
        return {
            "module": self.label,
            "genes": "+".join(self.genes),
            "k": len(self.genes),
            "n_lines": self.n_lines,
            "T": round(self.statistic, 4),
            "p_perm": round(self.p_value, 4),
            "floor_p05": round(self.floor, 4),
            "reaches_floor": abs(self.statistic) >= self.floor,
            "min_pair_rho": round(min(self.pairwise.values()), 4),
            "max_pair_rho": round(max(self.pairwise.values()), 4),
            "has_ties": self.has_ties,
        }


def module_coherence(
    effects: pd.DataFrame, genes: Sequence[str], label: str
) -> ModuleResult:
    """Coherence of one gene set, with its permutation p and matching floor.

    The null is rebuilt at this set's own (n_lines, k) rather than reused, so a
    4-gene control is never judged against a 3-gene null.
    """
    values = {g: effects[g].to_numpy(dtype=float) for g in genes}
    observed = coherence.mean_pairwise_rho(values)

    if observed.has_ties:
        # Prespec Section 4 tie fallback: the distribution-free guarantee lapses,
        # so permute the observed tied rank vectors instead of 1..n.
        null = _tied_null(values, n_permutations=coherence.N_PERMUTATIONS, seed=SEED)
        floor = float(np.quantile(np.abs(null), 0.95))
    else:
        null = coherence.null_distribution(observed.n_lines, observed.n_genes)
        floor = coherence.detection_floor(observed.n_lines, observed.n_genes)

    return ModuleResult(
        label=label,
        genes=tuple(genes),
        n_lines=observed.n_lines,
        statistic=observed.statistic,
        p_value=coherence.permutation_p_two_sided(observed.statistic, null),
        pairwise=observed.pairwise,
        floor=floor,
        has_ties=observed.has_ties,
    )


def _tied_null(
    values: dict[str, np.ndarray], *, n_permutations: int, seed: int
) -> np.ndarray:
    """Null built by permuting the observed (tied) rank vectors.

    Only reached if some gene has tied values among the analysed lines. Chronos
    effects are continuous floats, so this is expected to be dead code; it exists
    because the prespec requires the fallback to be mechanical.
    """
    rng = np.random.default_rng(seed)
    ranks = {g: coherence._ranks(v) for g, v in values.items()}
    out = np.empty(n_permutations, dtype=float)
    for i in range(n_permutations):
        out[i] = coherence.mean_pairwise_rho(
            {g: rng.permutation(r) for g, r in ranks.items()}
        ).statistic
    return out


def pairwise_bootstrap_ci(
    effects: pd.DataFrame,
    genes: Sequence[str],
    *,
    n_bootstrap: int = N_BOOTSTRAP,
    seed: int = SEED,
) -> pd.DataFrame:
    """Per-pair rho with a percentile bootstrap CI over cell lines (S1).

    Resamples *lines*, which is the unit of observation, keeping each line's gene
    values together.
    """
    rng = np.random.default_rng(seed)
    n = len(effects)
    rows = []
    for a, b in combinations(genes, 2):
        x = effects[a].to_numpy(dtype=float)
        y = effects[b].to_numpy(dtype=float)
        point = coherence._spearman(x, y)
        draws = np.empty(n_bootstrap)
        for i in range(n_bootstrap):
            idx = rng.integers(0, n, size=n)
            # A resample can duplicate lines and so create ties; average ranks
            # handle that, and the CI is a descriptive interval either way.
            draws[i] = coherence._spearman(x[idx], y[idx])
        finite = draws[np.isfinite(draws)]
        rows.append(
            {
                "pair": f"{a}–{b}",
                "rho": round(point, 4),
                "ci_lo": round(float(np.quantile(finite, 0.025)), 4),
                "ci_hi": round(float(np.quantile(finite, 0.975)), 4),
                "n_bootstrap_finite": int(len(finite)),
            }
        )
    return pd.DataFrame(rows)


def pairwise_permutation_p(
    effects: pd.DataFrame, genes: Sequence[str], *, seed: int = SEED
) -> pd.DataFrame:
    """Two-sided permutation p for each pair (S1), at this n and k = 2."""
    null = coherence.null_distribution(len(effects), 2, seed=seed)
    floor = float(np.quantile(np.abs(null), 0.95))
    rows = []
    for a, b in combinations(genes, 2):
        rho = coherence._spearman(
            effects[a].to_numpy(dtype=float), effects[b].to_numpy(dtype=float)
        )
        rows.append(
            {
                "pair": f"{a}–{b}",
                "rho": round(rho, 4),
                "p_perm": round(coherence.permutation_p_two_sided(rho, null), 4),
                "floor_p05_k2": round(floor, 4),
            }
        )
    return pd.DataFrame(rows)


def matched_triplet_null(
    pool: pd.DataFrame,
    target_genes: Sequence[str],
    *,
    n_triplets: int = MATCHED_N_TRIPLETS,
    mean_tol: float = MATCHED_MEAN_TOL,
    sd_tol: float = MATCHED_SD_TOL,
    min_candidates: int = MATCHED_MIN_CANDIDATES,
    seed: int = SEED,
) -> tuple[np.ndarray, pd.DataFrame]:
    """Section 5: T over mean- and SD-matched random triplets.

    `pool` is the genome-wide complete-column gene-effect matrix for the same
    lines. For each target gene, candidate replacements are those whose own mean
    and SD across these lines sit within tolerance of the target's. An unmatched
    random null would be an artefact generator here, because the three target
    genes sit in three different dependency regimes (prespec Section 2.1).
    """
    rng = np.random.default_rng(seed)
    stats = pd.DataFrame(
        {"mean": pool.mean(axis=0), "sd": pool.std(axis=0, ddof=1)}
    )

    candidates: dict[str, np.ndarray] = {}
    diagnostics = []
    for gene in target_genes:
        target_mean = float(pool[gene].mean())
        target_sd = float(pool[gene].std(ddof=1))
        m_tol, s_tol, widenings = mean_tol, sd_tol, 0
        while True:
            hit = stats.index[
                (stats["mean"].sub(target_mean).abs() <= m_tol)
                & (stats["sd"].sub(target_sd).abs() <= s_tol)
            ]
            hit = hit.drop(list(target_genes), errors="ignore")
            if len(hit) >= min_candidates:
                break
            # Stated fallback, not a free parameter (prespec Section 5).
            m_tol += 0.05
            s_tol += 0.02
            widenings += 1
            if widenings > 50:
                raise RuntimeError(f"cannot match {gene} even after 50 widenings")
        candidates[gene] = hit.to_numpy()
        diagnostics.append(
            {
                "gene": gene,
                "target_mean": round(target_mean, 4),
                "target_sd": round(target_sd, 4),
                "mean_tol": round(m_tol, 4),
                "sd_tol": round(s_tol, 4),
                "widenings": widenings,
                "n_candidates": int(len(hit)),
            }
        )

    arrays = {g: pool[list(c)].to_numpy(dtype=float) for g, c in candidates.items()}
    draws = np.empty(n_triplets, dtype=float)
    for i in range(n_triplets):
        picked = {
            g: arrays[g][:, rng.integers(0, arrays[g].shape[1])]
            for g in target_genes
        }
        draws[i] = coherence.mean_pairwise_rho(picked).statistic
    return draws, pd.DataFrame(diagnostics)


def variance_decomposition(effects: pd.DataFrame) -> pd.DataFrame:
    """Section 7 S2a: each Stage 1 gene's contribution to composite variance.

    The composite is `mean(BIR) - mean(HR)`, a fixed linear combination with
    weights +1/3 on each BIR gene and -1/3 on each HR gene. So each gene's
    contribution is decomposed as its own weighted variance plus its share of
    the covariance terms -- the standard decomposition
    `Var(sum w_i X_i) = sum_i w_i cov(X_i, sum_j w_j X_j)`, whose terms sum
    exactly to the composite variance.

    No p-value: this is deterministic given the 16 lines.
    """
    weights = {g: 1 / 3 for g in STAGE1_BIR} | {g: -1 / 3 for g in STAGE1_HR}
    genes = list(weights)
    composite = sum(weights[g] * effects[g] for g in genes)
    total = float(composite.var(ddof=1))

    rows = []
    exact_sum = 0.0
    for g in genes:
        contribution = weights[g] * float(
            np.cov(
                effects[g].to_numpy(dtype=float),
                composite.to_numpy(dtype=float),
                ddof=1,
            )[0, 1]
        )
        exact_sum += contribution
        rows.append(
            {
                "gene": g,
                "arm": "BIR" if g in STAGE1_BIR else "HR",
                "weight": round(weights[g], 4),
                "own_sd": round(float(effects[g].std(ddof=1)), 4),
                "variance_contribution": round(contribution, 5),
                "share_of_total": round(contribution / total, 4) if total else float("nan"),
            }
        )
    frame = pd.DataFrame(rows)
    # The decomposition is exact, so this is a real check on the arithmetic --
    # run on the unrounded contributions, since the rounding written into the
    # frame for display would otherwise break the identity it is checking.
    assert abs(exact_sum - total) < 1e-9, f"{exact_sum} != {total}"
    frame.attrs["composite_variance"] = total
    return frame


def leave_one_out(effects: pd.DataFrame, genes: Sequence[str]) -> pd.DataFrame:
    """Section 7 `loo_no_sign_flip`: refit T dropping each line in turn."""
    rows = []
    for model_id in effects.index:
        subset = effects.drop(index=model_id)
        values = {g: subset[g].to_numpy(dtype=float) for g in genes}
        rows.append(
            {
                "dropped": model_id,
                "n_lines": len(subset),
                "T": round(coherence.mean_pairwise_rho(values).statistic, 4),
            }
        )
    return pd.DataFrame(rows)

