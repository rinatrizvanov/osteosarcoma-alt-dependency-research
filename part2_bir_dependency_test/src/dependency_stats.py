"""Outcome loading and the Section 5/6 statistics.

**This is the only module that reads gene-effect values.** Importing it is the
point at which the outcome variable enters the analysis, which is why notebooks
00-04 do not import it and notebook 02 structurally cannot touch the outcome.

Sign convention (prespec, binding everywhere)
---------------------------------------------
Chronos gene effect: more negative = more essential; 0 is non-essential, -1 is
the median common essential.

    composite_i = mean(POLD3, RAD52, SMARCAL1)_i - mean(BRCA1, BRCA2, RAD51)_i

Under the hypothesis, TCAB1-low lines are *more* dependent on BIR relative to
HR, so their composite is **more negative**:

    median(composite | TCAB1-low) < median(composite | TCAB1-intact)
    cliffs_delta(low, intact) < 0

Tests are two-sided; the success criterion additionally requires this direction.

Determinism: every bootstrap and permutation draws from an explicitly seeded
``numpy`` generator (default :data:`src.utils.RANDOM_SEED`), so reruns are
byte-identical.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd
from scipy import stats

from .utils import CRISPR_PATH, RANDOM_SEED, verify_expected

# ---------------------------------------------------------------------------
# Gene sets (prespec Section 5)
# ---------------------------------------------------------------------------

#: Non-classical break-induced-replication module.
BIR_GENES: tuple[str, ...] = ("POLD3", "RAD52", "SMARCAL1")

#: Canonical homologous-recombination comparators.
HR_GENES: tuple[str, ...] = ("BRCA1", "BRCA2", "RAD51")

#: The six primary-comparison genes, BIR first.
SIX_GENES: tuple[str, ...] = BIR_GENES + HR_GENES

#: Exploratory set, declared now, usable only after the Section 10 stop rule
#: triggers (prespec Section 11.1). Declaring it here prevents it being
#: assembled post-hoc.
BROADER_SET: tuple[str, ...] = (
    "POLD3", "POLD1", "POLD2", "RAD52", "SMARCAL1", "SLX4", "BLM", "WRN",
    "FANCM", "BRCA1", "BRCA2", "RAD51", "RAD51AP1", "MUS81", "GEN1", "NBN",
    "MRE11", "RAD50", "ATRX", "DAXX",
)

N_BOOTSTRAP = 2000
N_PERMUTATIONS = 5000

#: Romano et al. Cliff's delta magnitude thresholds.
DELTA_SMALL, DELTA_MEDIUM, DELTA_LARGE = 0.147, 0.330, 0.474


# ---------------------------------------------------------------------------
# Outcome loading
# ---------------------------------------------------------------------------


def strip_gene_symbol(col: str) -> str:
    """``'POLD3 (10714)'`` -> ``'POLD3'`` (DepMap column convention)."""
    return col.split(" (")[0].strip()


def crispr_gene_columns(
    symbols: Sequence[str], path: Path | None = None
) -> dict[str, str]:
    """Map each gene symbol to its exact ``'SYMBOL (entrez)'`` column name.

    Raises ``KeyError`` naming every symbol absent from the header, rather than
    silently analysing a smaller gene set.
    """
    p = Path(path) if path is not None else CRISPR_PATH
    with open(p) as fh:
        header = next(fh).rstrip("\n").split(",")
    lookup: dict[str, str] = {}
    for col in header:
        lookup.setdefault(strip_gene_symbol(col), col)
    missing = [s for s in symbols if s not in lookup]
    if missing:
        raise KeyError(f"Genes absent from {p.name}: {missing}")
    return {s: lookup[s] for s in symbols}


def load_gene_effects(
    model_ids: Sequence[str],
    symbols: Sequence[str] = SIX_GENES,
    path: Path | None = None,
    *,
    verify: bool = True,
) -> pd.DataFrame:
    """Chronos gene effects for the requested models x genes.

    Reads only the needed columns out of the 440 MB matrix. Returns a frame
    indexed by ModelID in the order given, with plain gene symbols as columns.
    Raises if any requested model or value is missing — Stage 1 does not
    impute the outcome.
    """
    p = Path(path) if path is not None else CRISPR_PATH
    if verify:
        verify_expected(p)
    colmap = crispr_gene_columns(symbols, p)
    index_col = pd.read_csv(p, nrows=0).columns[0]
    frame = pd.read_csv(p, usecols=[index_col, *colmap.values()])
    frame = frame.rename(columns={index_col: "ModelID"})
    frame = frame.set_index(frame["ModelID"].astype(str)).drop(columns=["ModelID"])
    frame = frame.rename(columns={v: k for k, v in colmap.items()})

    wanted = [str(m) for m in model_ids]
    missing = sorted(set(wanted) - set(frame.index))
    if missing:
        raise ValueError(f"{p.name} missing requested models: {missing}")
    out = frame.loc[wanted, list(symbols)].astype(float)
    if out.isna().any().any():
        bad = out.isna().any(axis=1)
        raise ValueError(
            "Missing gene-effect values for "
            f"{out.index[bad].tolist()} — Stage 1 does not impute the outcome."
        )
    out.index.name = "ModelID"
    return out


def load_all_gene_effects_for_line(
    model_ids: Sequence[str], path: Path | None = None
) -> pd.DataFrame:
    """Full gene-effect rows for the given models (all ~18.5k genes).

    Used only by the Section 6.3 within-line standardisation sensitivity, which
    needs each line's genome-wide mean and SD. Streams and keeps matching rows.
    """
    p = Path(path) if path is not None else CRISPR_PATH
    wanted = {str(m) for m in model_ids}
    kept: list[pd.DataFrame] = []
    index_col = None
    for chunk in pd.read_csv(p, chunksize=200):
        if index_col is None:
            index_col = chunk.columns[0]
        sub = chunk[chunk[index_col].astype(str).isin(wanted)]
        if len(sub):
            kept.append(sub)
    if not kept:
        raise ValueError(f"No rows in {p.name} for the {len(wanted)} requested models.")
    frame = pd.concat(kept, ignore_index=True)
    frame = frame.set_index(frame[index_col].astype(str)).drop(columns=[index_col])
    frame.index.name = "ModelID"
    return frame.loc[[str(m) for m in model_ids]]


# ---------------------------------------------------------------------------
# Composite endpoint
# ---------------------------------------------------------------------------


def composite_score(
    effects: pd.DataFrame,
    bir: Sequence[str] = BIR_GENES,
    hr: Sequence[str] = HR_GENES,
) -> pd.Series:
    """mean(BIR) - mean(HR) per line. Lower (more negative) = BIR-shifted."""
    missing = [g for g in (*bir, *hr) if g not in effects.columns]
    if missing:
        raise KeyError(f"composite_score missing gene columns: {missing}")
    score = effects[list(bir)].mean(axis=1) - effects[list(hr)].mean(axis=1)
    score.name = "composite"
    return score


def standardize_within_line(all_effects: pd.DataFrame) -> pd.DataFrame:
    """z-score each line's gene effects across all genes (Section 6.3).

    Controls for a line having globally more negative gene effects (screen
    quality / essentiality scale) in a stronger way than the composite
    difference alone.
    """
    values = all_effects.to_numpy(dtype=float)
    mean = np.nanmean(values, axis=1, keepdims=True)
    sd = np.nanstd(values, axis=1, ddof=1, keepdims=True)
    return pd.DataFrame(
        (values - mean) / sd, index=all_effects.index, columns=all_effects.columns
    )


# ---------------------------------------------------------------------------
# Effect sizes
# ---------------------------------------------------------------------------


def cliffs_delta(x: Sequence[float], y: Sequence[float]) -> float:
    """Cliff's delta = P(x>y) - P(x<y), in [-1, 1].

    ``delta < 0`` means ``x`` tends to sit below ``y``. With ``x`` = TCAB1-low
    and ``y`` = TCAB1-intact composites, that is the predicted direction.
    """
    a = np.asarray(x, dtype=float)
    b = np.asarray(y, dtype=float)
    if a.size == 0 or b.size == 0:
        return float("nan")
    diff = np.sign(a[:, None] - b[None, :])
    return float(diff.sum() / (a.size * b.size))


def hodges_lehmann(x: Sequence[float], y: Sequence[float]) -> float:
    """Median of all pairwise ``x_i - y_j``, in Chronos units."""
    a = np.asarray(x, dtype=float)
    b = np.asarray(y, dtype=float)
    if a.size == 0 or b.size == 0:
        return float("nan")
    return float(np.median((a[:, None] - b[None, :]).ravel()))


def bootstrap_delta_ci(
    x: Sequence[float],
    y: Sequence[float],
    *,
    n_boot: int = N_BOOTSTRAP,
    seed: int = RANDOM_SEED,
    alpha: float = 0.05,
) -> tuple[float, float]:
    """Percentile bootstrap CI for Cliff's delta, resampling within each group."""
    a = np.asarray(x, dtype=float)
    b = np.asarray(y, dtype=float)
    rng = np.random.default_rng(seed)
    draws = np.empty(n_boot, dtype=float)
    for i in range(n_boot):
        draws[i] = cliffs_delta(
            a[rng.integers(0, a.size, a.size)], b[rng.integers(0, b.size, b.size)]
        )
    return (
        float(np.quantile(draws, alpha / 2)),
        float(np.quantile(draws, 1 - alpha / 2)),
    )


def delta_magnitude(delta: float) -> str:
    """Romano et al. label for |delta|."""
    d = abs(delta)
    if np.isnan(d):
        return "undefined"
    if d < DELTA_SMALL:
        return "negligible"
    if d < DELTA_MEDIUM:
        return "small"
    if d < DELTA_LARGE:
        return "medium"
    return "large"


# ---------------------------------------------------------------------------
# Two-group comparison
# ---------------------------------------------------------------------------


@dataclass
class GroupComparison:
    """One Mann-Whitney comparison with both effect sizes."""

    name: str
    n_low: int
    n_intact: int
    median_low: float
    median_intact: float
    u_statistic: float
    p_value: float
    mwu_method: str
    cliffs_delta: float
    delta_ci_low: float
    delta_ci_high: float
    delta_magnitude: str
    hodges_lehmann: float
    direction_as_predicted: bool

    def as_dict(self) -> dict:
        return asdict(self)


def compare_groups(
    low: Sequence[float],
    intact: Sequence[float],
    name: str = "composite",
    *,
    n_boot: int = N_BOOTSTRAP,
    seed: int = RANDOM_SEED,
) -> GroupComparison:
    """Two-sided Mann-Whitney plus Cliff's delta and Hodges-Lehmann.

    Uses the exact null distribution. ``scipy``'s exact method is invalid under
    ties, so if ties are present the method falls back to the asymptotic test
    with continuity correction and records which was used in ``mwu_method`` —
    the choice is never silent.
    """
    a = np.asarray(low, dtype=float)
    b = np.asarray(intact, dtype=float)
    has_ties = len(np.unique(np.concatenate([a, b]))) < a.size + b.size
    method = "asymptotic" if has_ties else "exact"
    res = stats.mannwhitneyu(a, b, alternative="two-sided", method=method)
    delta = cliffs_delta(a, b)
    ci_low, ci_high = bootstrap_delta_ci(a, b, n_boot=n_boot, seed=seed)
    return GroupComparison(
        name=name,
        n_low=int(a.size),
        n_intact=int(b.size),
        median_low=float(np.median(a)),
        median_intact=float(np.median(b)),
        u_statistic=float(res.statistic),
        p_value=float(res.pvalue),
        mwu_method=method,
        cliffs_delta=delta,
        delta_ci_low=ci_low,
        delta_ci_high=ci_high,
        delta_magnitude=delta_magnitude(delta),
        hodges_lehmann=hodges_lehmann(a, b),
        direction_as_predicted=bool(delta < 0),
    )


def per_gene_comparisons(
    effects: pd.DataFrame,
    low_ids: Sequence[str],
    intact_ids: Sequence[str],
    symbols: Sequence[str] = SIX_GENES,
    *,
    seed: int = RANDOM_SEED,
) -> pd.DataFrame:
    """Section 6.1: one comparison per gene, with BH-FDR across exactly these genes.

    For a single gene the predicted direction is the same as for the composite
    for BIR genes (low more essential => more negative) and is *not* predicted
    for HR genes; ``gene_set`` is recorded so the Section 8 criterion can apply
    the right expectation per gene.
    """
    rows = []
    for gene in symbols:
        cmp = compare_groups(
            effects.loc[list(low_ids), gene],
            effects.loc[list(intact_ids), gene],
            name=gene,
            seed=seed,
        )
        row = cmp.as_dict()
        row["gene_set"] = "BIR" if gene in BIR_GENES else "HR"
        rows.append(row)
    frame = pd.DataFrame(rows)
    frame["q_value"] = stats.false_discovery_control(frame["p_value"], method="bh")
    return frame


def bh_fdr(p_values: Sequence[float]) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values (scipy >= 1.11)."""
    return stats.false_discovery_control(np.asarray(p_values, dtype=float), method="bh")


# ---------------------------------------------------------------------------
# Interaction (Section 6.2)
# ---------------------------------------------------------------------------


def delta_difference_permutation(
    alt_pos_low: Sequence[float],
    alt_pos_intact: Sequence[float],
    alt_neg_low: Sequence[float],
    alt_neg_intact: Sequence[float],
    *,
    n_perm: int = N_PERMUTATIONS,
    seed: int = RANDOM_SEED,
) -> dict:
    """Selectivity test: is the TCAB1 effect larger in ALT+ than in ALT-?

    Statistic = delta(ALT+) - delta(ALT-). The null shuffles the TCAB1 label
    **within** each ALT stratum, preserving stratum sizes and the ALT main
    effect, so the p-value speaks only to the interaction.

    Underpowered at the ALT- stratum's N; reported as corroborative, never as
    the primary claim.
    """
    rng = np.random.default_rng(seed)
    pos = np.concatenate([np.asarray(alt_pos_low, float), np.asarray(alt_pos_intact, float)])
    neg = np.concatenate([np.asarray(alt_neg_low, float), np.asarray(alt_neg_intact, float)])
    n_pos_low, n_neg_low = len(alt_pos_low), len(alt_neg_low)

    delta_pos = cliffs_delta(alt_pos_low, alt_pos_intact)
    delta_neg = cliffs_delta(alt_neg_low, alt_neg_intact)
    observed = delta_pos - delta_neg

    # The ALT- stratum's own test. Section 8.5 gates on this rather than on a
    # delta point comparison, because at these group sizes delta alone is far
    # too noisy to gate on (see src/criteria.py S8.5 for the numbers).
    neg_cmp = compare_groups(alt_neg_low, alt_neg_intact, name="alt_negative")

    draws = np.empty(n_perm, dtype=float)
    for i in range(n_perm):
        p = rng.permutation(pos)
        n = rng.permutation(neg)
        draws[i] = (
            cliffs_delta(p[:n_pos_low], p[n_pos_low:])
            - cliffs_delta(n[:n_neg_low], n[n_neg_low:])
        )
    return {
        "delta_alt_positive": delta_pos,
        "delta_alt_negative": delta_neg,
        "delta_difference": observed,
        "p_value": float(np.mean(np.abs(draws) >= abs(observed))),
        "n_perm": n_perm,
        "selectivity_as_predicted": bool(abs(delta_pos) > abs(delta_neg)),
        "alt_negative_p_value": neg_cmp.p_value,
        "alt_negative_mwu_method": neg_cmp.mwu_method,
        "alt_negative_contradicts": bool(neg_cmp.p_value < 0.05 and delta_neg < 0),
    }
