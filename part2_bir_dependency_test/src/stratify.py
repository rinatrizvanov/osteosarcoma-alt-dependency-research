"""Exposure side: expression loading and the Section 4 TCAB1 split rule.

TCAB1/WRAP53 expression is the **exposure (stratifier)**, not an outcome, so
this module and :mod:`src.strand_check` run before the pre-specification is
frozen. Neither imports :mod:`src.dependency_stats`; the outcome stays sealed
until notebook 05.

Section 4 split rule, in order:

1. WRAP53 log2(TPM+1) for **all** osteosarcoma lines with RNA (34), not just
   the 16 ALT+ lines entering the contrast.
2. ``cut`` = 33.3rd percentile of that distribution.
3. Within ALT+: ``TCAB1-low`` if WRAP53 < ``cut``, else ``TCAB1-intact``.
4. Group sizes are an empirical result and are reported, not fixed in advance.

Estimating the cutpoint on 34 lines rather than 16 makes it both better
determined and independent of which lines happen to be ALT+, so it cannot be
tuned toward the outcome.

Minimum-cell guard: if the TCAB1-low arm has fewer than
:data:`MIN_CELL_SIZE` lines the tertile primary is declared not estimable and
the analysis falls back to a median split within the ALT+ set. This rule is
written here **before** the WRAP53 distribution is examined, which is what
keeps the fallback pre-specified rather than post-hoc.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd

from .cohort import EXPRESSION_META_COLS, osteosarcoma_models, rna_model_ids
from .utils import CN_GENE_PATH, RNA_PATH, verify_expected

TCAB1_GENE = "WRAP53"
TP53_GENE = "TP53"

#: Section 4 step 2. Bottom tertile.
TERTILE_Q = 1.0 / 3.0

#: Section 4 minimum-cell guard.
MIN_CELL_SIZE = 4

LOW_LABEL = "TCAB1-low"
INTACT_LABEL = "TCAB1-intact"


# ---------------------------------------------------------------------------
# Expression loading
# ---------------------------------------------------------------------------


def _strip(col: str) -> str:
    return col.split(" (")[0].strip()


def expression_gene_columns(
    symbols: Sequence[str], path: Path | None = None
) -> dict[str, str]:
    p = Path(path) if path is not None else RNA_PATH
    with open(p) as fh:
        header = next(fh).rstrip("\n").split(",")
    lookup: dict[str, str] = {}
    for col in header:
        if col in EXPRESSION_META_COLS:
            continue
        lookup.setdefault(_strip(col), col)
    missing = [s for s in symbols if s not in lookup]
    if missing:
        raise KeyError(f"Genes absent from {p.name}: {missing}")
    return {s: lookup[s] for s in symbols}


def load_expression(
    model_ids: Sequence[str] | None,
    symbols: Sequence[str] = (TCAB1_GENE, TP53_GENE),
    path: Path | None = None,
    *,
    verify: bool = True,
) -> pd.DataFrame:
    """log2(TPM+1) for the requested models x genes.

    The matrix has one row per Sequencing x Model x Condition, so rows are
    reduced to the default model entry before indexing by ModelID — the same
    rule Aim 2's ``aim2_signature.load_expression`` applied.
    ``model_ids=None`` returns every model.
    """
    p = Path(path) if path is not None else RNA_PATH
    if verify:
        verify_expected(p)
    colmap = expression_gene_columns(symbols, p)
    meta = [c for c in ("ModelID", "IsDefaultEntryForModel") if c is not None]
    frame = pd.read_csv(p, usecols=[*meta, *colmap.values()])

    if "IsDefaultEntryForModel" in frame.columns:
        default = frame[frame["IsDefaultEntryForModel"].astype(str) == "Yes"]
        missing = set(frame["ModelID"]) - set(default["ModelID"])
        if missing:
            default = pd.concat([default, frame[frame["ModelID"].isin(missing)]])
        frame = default
    frame = frame.drop_duplicates(subset="ModelID", keep="first")
    frame = frame.set_index(frame["ModelID"].astype(str))
    frame = frame[list(colmap.values())].rename(
        columns={v: k for k, v in colmap.items()}
    )
    frame.index.name = "ModelID"

    if model_ids is None:
        return frame.astype(float)
    wanted = [str(m) for m in model_ids]
    absent = sorted(set(wanted) - set(frame.index))
    if absent:
        raise ValueError(f"{p.name} missing requested models: {absent}")
    return frame.loc[wanted].astype(float)


def osteosarcoma_rna_models() -> list[str]:
    """The osteosarcoma lines with RNA — the 34 defining the Section 4 cutpoint.

    Deliberately *not* restricted to CRISPR-screened lines: the cutpoint is
    estimated on the widest available osteosarcoma expression distribution.
    """
    os_ids = set(osteosarcoma_models()["ModelID"])
    return sorted(os_ids & rna_model_ids())


def load_wrap53_cn(
    model_ids: Sequence[str], path: Path | None = None
) -> pd.Series:
    """WRAP53 gene-level copy number — strand-agnostic corroboration.

    ``OmicsCNGeneWGS.csv`` shares the expression matrix's layout: one row per
    Sequencing x Model x Condition, with ``ModelID`` and
    ``IsDefaultEntryForModel`` as columns rather than a ModelID row index. Rows
    are reduced to the default model entry before indexing, the same rule
    :func:`load_expression` applies. Models absent from the file come back as
    ``NaN`` so the caller can report coverage rather than silently losing them.
    """
    p = Path(path) if path is not None else CN_GENE_PATH
    with open(p) as fh:
        header = next(fh).rstrip("\n").split(",")
    col = next((c for c in header if _strip(c) == TCAB1_GENE), None)
    if col is None:
        raise KeyError(f"{TCAB1_GENE} absent from {p.name}")

    usecols = [c for c in ("ModelID", "IsDefaultEntryForModel") if c in header]
    frame = pd.read_csv(p, usecols=[*usecols, col])
    if "IsDefaultEntryForModel" in frame.columns:
        default = frame[frame["IsDefaultEntryForModel"].astype(str) == "Yes"]
        missing = set(frame["ModelID"]) - set(default["ModelID"])
        if missing:
            default = pd.concat([default, frame[frame["ModelID"].isin(missing)]])
        frame = default
    frame = frame.drop_duplicates(subset="ModelID", keep="first")
    series = frame.set_index(frame["ModelID"].astype(str))[col].astype(float)
    series.index.name = "ModelID"
    series.name = "WRAP53_cn"
    return series.reindex([str(m) for m in model_ids])


# ---------------------------------------------------------------------------
# Section 4 split
# ---------------------------------------------------------------------------


@dataclass
class StratificationResult:
    """Outcome of applying the Section 4 rule, including which rule was used."""

    rule: str                       # 'tertile_os_cohort' or 'median_within_alt_pos'
    cut: float
    reference_n: int                # lines the cutpoint was estimated on
    strata: pd.Series               # ModelID -> LOW_LABEL / INTACT_LABEL
    expression: pd.Series           # ModelID -> WRAP53 log2(TPM+1)
    guard_triggered: bool
    guard_reason: str = ""
    notes: list[str] = field(default_factory=list)

    @property
    def low_ids(self) -> list[str]:
        return self.strata.index[self.strata == LOW_LABEL].tolist()

    @property
    def intact_ids(self) -> list[str]:
        return self.strata.index[self.strata == INTACT_LABEL].tolist()

    @property
    def sizes(self) -> dict[str, int]:
        return {LOW_LABEL: len(self.low_ids), INTACT_LABEL: len(self.intact_ids)}


def tertile_cut(values: Sequence[float]) -> float:
    """33.3rd percentile of the reference distribution (Section 4 step 2)."""
    return float(np.quantile(np.asarray(values, dtype=float), TERTILE_Q))


def stratify(
    alt_positive_ids: Sequence[str],
    *,
    reference_ids: Sequence[str] | None = None,
    expression: pd.Series | None = None,
    min_cell_size: int = MIN_CELL_SIZE,
) -> StratificationResult:
    """Apply the Section 4 rule with its pre-specified minimum-cell fallback.

    Parameters
    ----------
    alt_positive_ids
        The 16 ALT+ lines entering the primary contrast.
    reference_ids
        Lines defining the cutpoint. Defaults to :func:`osteosarcoma_rna_models`
        (the 34 osteosarcoma lines with RNA).
    expression
        Pre-loaded WRAP53 series, for tests. Loaded if omitted.
    """
    ref_ids = list(reference_ids) if reference_ids is not None else osteosarcoma_rna_models()
    if expression is None:
        expression = load_expression(
            sorted(set(ref_ids) | set(map(str, alt_positive_ids))), [TCAB1_GENE]
        )[TCAB1_GENE]

    ref_values = expression.loc[[i for i in ref_ids if i in expression.index]]
    cut = tertile_cut(ref_values)

    alt_expr = expression.loc[[str(m) for m in alt_positive_ids]]
    strata = pd.Series(
        np.where(alt_expr < cut, LOW_LABEL, INTACT_LABEL),
        index=alt_expr.index,
        name="tcab1_stratum",
    )
    n_low = int((strata == LOW_LABEL).sum())

    notes = [
        f"Cutpoint = {TERTILE_Q:.4f} quantile of WRAP53 log2(TPM+1) over "
        f"{len(ref_values)} osteosarcoma lines with RNA = {cut:.6f}",
        f"Applied to {len(alt_expr)} ALT+ lines -> {n_low} {LOW_LABEL} / "
        f"{len(alt_expr) - n_low} {INTACT_LABEL}",
    ]

    if n_low >= min_cell_size and (len(alt_expr) - n_low) >= min_cell_size:
        return StratificationResult(
            rule="tertile_os_cohort",
            cut=cut,
            reference_n=len(ref_values),
            strata=strata,
            expression=alt_expr,
            guard_triggered=False,
            notes=notes,
        )

    # Pre-specified fallback (Section 4 minimum-cell guard).
    reason = (
        f"TCAB1-low arm n={n_low} (intact n={len(alt_expr) - n_low}); "
        f"minimum cell size is {min_cell_size}. Tertile primary declared NOT "
        f"ESTIMABLE; falling back to the median split within the ALT+ set, "
        f"per the pre-specified Section 4 guard."
    )
    median_cut = float(np.median(alt_expr.to_numpy(dtype=float)))
    strata = pd.Series(
        np.where(alt_expr < median_cut, LOW_LABEL, INTACT_LABEL),
        index=alt_expr.index,
        name="tcab1_stratum",
    )
    notes.append(reason)
    notes.append(f"Median cutpoint within ALT+ = {median_cut:.6f}")
    return StratificationResult(
        rule="median_within_alt_pos",
        cut=median_cut,
        reference_n=len(alt_expr),
        strata=strata,
        expression=alt_expr,
        guard_triggered=True,
        guard_reason=reason,
        notes=notes,
    )


def stratify_median(alt_positive_ids: Sequence[str], expression: pd.Series) -> pd.Series:
    """R2 sensitivity: median split within the ALT+ set."""
    expr = expression.loc[[str(m) for m in alt_positive_ids]]
    cut = float(np.median(expr.to_numpy(dtype=float)))
    return pd.Series(
        np.where(expr < cut, LOW_LABEL, INTACT_LABEL), index=expr.index,
        name="tcab1_stratum_median",
    )


def stratify_tertile_within(
    alt_positive_ids: Sequence[str], expression: pd.Series
) -> pd.Series:
    """R3 sensitivity: tertile computed within the ALT+ set only."""
    expr = expression.loc[[str(m) for m in alt_positive_ids]]
    cut = tertile_cut(expr)
    return pd.Series(
        np.where(expr < cut, LOW_LABEL, INTACT_LABEL), index=expr.index,
        name="tcab1_stratum_tertile_within",
    )
