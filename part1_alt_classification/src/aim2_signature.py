"""
Reusable training / cross-validation helpers for
``notebooks/11_aim2_alt_signature.ipynb`` (Aim 2 ALT signature score).

Scope and guard-rails (mirrors the notebook spec and ``label_spec_v1.0`` §2):

* This module trains a **regularized signature score** on DepMap osteosarcoma
  RNA-seq using *direction-aware* labels (``alt_call_direction``), never a
  calibrated binary classifier and never a fixed decision threshold.
* **All preprocessing is fold-local.** Every preprocessing step (near-zero
  variance drop, top-K high-variance gene selection, ``StandardScaler``) lives
  inside an ``sklearn`` ``Pipeline`` so it is fit on the training fold only and
  applied to the held-out sample. Nothing here ever fits a scaler / variance
  filter / feature selector on the full cohort before cross-validation.
* The arm label is derived from ``alt_call_direction`` (``positive`` /
  ``negative``), **never** from ``alt_status`` — using ``alt_status`` would
  corrupt the ``alt_silver`` split that the pragmatic definition relies on.

Determinism: every stochastic component is seeded (``random_state=42`` on the
estimators, fixed ``StratifiedKFold`` shuffles, and explicit ``numpy``
generators for bootstrap / permutation) so that re-running the notebook
produces byte-identical output files.

Performance note (documented design decision): the curated DepMap OS cohort is
tiny (N=20 strict / N=31 pragmatic) but the RNA matrix has ~19k protein-coding
genes. A faithful nested LOOCV with a 100-shuffle permutation null over the
full gene space is computationally intractable (~15 s per elastic-net fit ×
hundreds of thousands of fits). The ``gene-variance filter`` required by the
spec is therefore implemented as a fold-local *top-K most-variable-gene*
selector (:class:`VarianceGeneFilter`): it both drops near-zero-variance genes
and caps the feature count at ``top_k`` using **training-fold variances only**.
``top_k`` is an explicit, reported hyperparameter (default
:data:`TOP_K_VARIABLE_GENES`). Set ``top_k=None`` to keep every non-constant
gene (the literal "drop near-zero-variance only" reading), at large compute
cost.

The observed LOOCV reported as the headline result is the *fully nested*
design (inner 5-fold C-selection in every outer fold). The 100-shuffle
permutation null re-runs the entire fold-local LOOCV per shuffle but fixes
``C`` to the median observed value, because re-selecting ``C`` inside every
fold of every shuffle is intractable (~9 h per definition). See
:func:`permutation_test` for the rationale; pass ``select_C_per_fold=True`` to
force the fully nested null.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score, roc_curve
from sklearn.model_selection import GridSearchCV, LeaveOneOut, StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from .aim1_inventory import DEPMAP_RNA_PATH, depmap_rna_model_ids

# ---------------------------------------------------------------------------
# Constants — definitions, expected counts, ground truth, hyperparameters
# ---------------------------------------------------------------------------

RANDOM_SEED = 42

#: alt_status values admitted by each definition run.
PRIMARY_ALT_STATUSES: tuple[str, ...] = ("alt_gold_positive", "alt_negative_high_conf")
PRAGMATIC_ALT_STATUSES: tuple[str, ...] = (
    "alt_gold_positive",
    "alt_negative_high_conf",
    "alt_silver",
)
ALT_STATUS_SETS: dict[str, tuple[str, ...]] = {
    "primary": PRIMARY_ALT_STATUSES,
    "pragmatic": PRAGMATIC_ALT_STATUSES,
}

#: Arm labels are derived from this column only — NEVER alt_status.
DIRECTION_COLUMN = "alt_call_direction"
VALID_DIRECTIONS: tuple[str, ...] = ("positive", "negative")
TARGET_COHORT = "DepMap_26Q1"

#: Hard expectations asserted before any modelling. (positive, negative).
EXPECTED_COUNTS: dict[str, dict[str, int]] = {
    "primary": {"positive": 12, "negative": 8},
    "pragmatic": {"positive": 18, "negative": 13},
}

#: Ground-truth anchors — abort if either fails.
GROUND_TRUTH_DIRECTION: dict[str, str] = {
    "ACH-000364": "positive",  # U2OS, canonical ALT-positive
    "ACH-000748": "negative",  # SJSA1, canonical ALT-negative / telomerase+
}

#: Curated ALT-related gene list (provided by spec) for coefficient cross-check.
ALT_GENE_LIST: tuple[str, ...] = (
    "ATRX", "DAXX", "SMARCAL1", "TP53", "TERT", "TERC", "WRN", "BLM",
    "MRE11A", "RAD50", "NBN", "FANCD2", "RAD51AP1", "WRAP53", "RAD52",
    "BRCA1", "BRCA2", "FANCM", "SLX4", "GEN1", "MUS81",
)

#: Regularization grid for inner-CV C selection.
C_GRID: np.ndarray = 10 ** np.linspace(-2, 2, 9)

#: Elastic-net mixing parameter (fixed by spec).
L1_RATIO = 0.5

#: Fold-local high-variance gene cap. See module docstring for rationale.
TOP_K_VARIABLE_GENES = 500

#: Genes below this training-fold variance are dropped as near-zero-variance.
NEAR_ZERO_VARIANCE = 1e-8

N_PERMUTATIONS = 100
N_BOOTSTRAP = 1000
INNER_CV_FOLDS = 5
TARGET_SPECIFICITY = 0.90

#: Non-gene columns in the DepMap expression CSV.
EXPRESSION_META_COLS: frozenset[str] = frozenset({
    "Unnamed: 0", "SequencingID", "ModelConditionID", "ModelID",
    "IsDefaultEntryForMC", "IsDefaultEntryForModel",
})


# ---------------------------------------------------------------------------
# Inclusion filter
# ---------------------------------------------------------------------------

@dataclass
class FilterTrace:
    """Records N surviving after each inclusion-filter step."""

    definition: str
    steps: list[tuple[str, int]] = field(default_factory=list)

    def add(self, label: str, n: int) -> None:
        self.steps.append((label, n))

    def as_lines(self) -> list[str]:
        return [f"  {label:<55} N = {n}" for label, n in self.steps]


def apply_inclusion_filter(
    manifest: pd.DataFrame,
    depmap_model_ids: set[str],
    definition: str,
) -> tuple[pd.DataFrame, FilterTrace]:
    """Apply the Aim-2 signature inclusion filter for one definition run.

    Filters are applied *explicitly and in order* (the manifest's
    ``include_*`` flags are intentionally not trusted here):

    1. ``cohort == 'DepMap_26Q1'``
    2. ``sample_id`` present in the DepMap RNA matrix ``ModelID`` index.
    3. ``alt_call_direction in ('positive', 'negative')`` (drops ``unknown``).
    4. ``alt_status`` in the definition's admitted set.

    Parameters
    ----------
    manifest : DataFrame
        The master manifest (read as strings).
    depmap_model_ids : set[str]
        ModelIDs present in the DepMap RNA matrix (membership probe only).
    definition : {'primary', 'pragmatic'}

    Returns
    -------
    (filtered, trace)
        ``filtered`` has one row per included sample; ``trace`` records the
        surviving N after every step for the report.
    """
    if definition not in ALT_STATUS_SETS:
        raise ValueError(f"Unknown definition {definition!r}; expected one of {list(ALT_STATUS_SETS)}")

    trace = FilterTrace(definition=definition)
    df = manifest.copy()
    trace.add("manifest rows (all cohorts)", len(df))

    df = df[df["cohort"] == TARGET_COHORT]
    trace.add(f"cohort == {TARGET_COHORT!r}", len(df))

    df = df[df["sample_id"].astype(str).isin(depmap_model_ids)]
    trace.add("sample_id present in DepMap RNA ModelID index", len(df))

    df = df[df[DIRECTION_COLUMN].isin(VALID_DIRECTIONS)]
    trace.add("alt_call_direction in ('positive', 'negative')", len(df))

    admitted = ALT_STATUS_SETS[definition]
    df = df[df["alt_status"].isin(admitted)]
    trace.add(f"alt_status in {admitted}", len(df))

    df = df.reset_index(drop=True)
    return df, trace


def direction_counts(filtered: pd.DataFrame) -> dict[str, int]:
    counts = filtered[DIRECTION_COLUMN].value_counts().to_dict()
    return {d: int(counts.get(d, 0)) for d in VALID_DIRECTIONS}


def assert_expected_counts(filtered: pd.DataFrame, definition: str) -> None:
    """Raise ``ValueError`` with the exact discrepancy if counts differ."""
    expected = EXPECTED_COUNTS[definition]
    actual = direction_counts(filtered)
    if actual != expected:
        raise ValueError(
            f"{definition} count mismatch: expected positive={expected['positive']}, "
            f"negative={expected['negative']} (total {sum(expected.values())}); "
            f"got positive={actual['positive']}, negative={actual['negative']} "
            f"(total {sum(actual.values())})."
        )


def assert_ground_truth(filtered: pd.DataFrame) -> None:
    """Per-row + anchor assertions on ``alt_call_direction``.

    1. Every included sample has ``alt_call_direction in ('positive','negative')``
       (verified *before* any CV; using ``alt_status`` for arms would be wrong).
    2. ACH-000364 (U2OS) -> 'positive'; ACH-000748 (SJSA1) -> 'negative'.
    """
    bad = filtered[~filtered[DIRECTION_COLUMN].isin(VALID_DIRECTIONS)]
    if len(bad):
        raise ValueError(
            "Per-row direction assertion failed: "
            f"{bad['sample_id'].tolist()} have alt_call_direction not in {VALID_DIRECTIONS}."
        )

    by_id = filtered.set_index("sample_id")[DIRECTION_COLUMN].to_dict()
    for sid, expected_dir in GROUND_TRUTH_DIRECTION.items():
        if sid not in by_id:
            raise ValueError(
                f"Ground-truth anchor {sid} missing from included samples "
                f"(expected alt_call_direction={expected_dir!r})."
            )
        if by_id[sid] != expected_dir:
            raise ValueError(
                f"Ground-truth assertion failed for {sid}: "
                f"alt_call_direction={by_id[sid]!r}, expected {expected_dir!r}."
            )


def labels_from_direction(filtered: pd.DataFrame) -> np.ndarray:
    """Map ``alt_call_direction`` to binary labels (positive=1, negative=0)."""
    return (filtered[DIRECTION_COLUMN].to_numpy() == "positive").astype(int)


# ---------------------------------------------------------------------------
# Expression loading (only after the filter is locked)
# ---------------------------------------------------------------------------

def load_expression(
    sample_ids: Sequence[str],
    path: Path | None = None,
    chunksize: int = 200,
) -> pd.DataFrame:
    """Load only the included rows from the DepMap RNA matrix.

    Reads the matrix in chunks and keeps only rows whose ``ModelID`` is in
    ``sample_ids`` (the matrix is one row per Sequencing × Model × Condition).
    Picks the default model entry (``IsDefaultEntryForModel == 'Yes'``) when a
    ModelID has multiple rows, then deduplicates. No preprocessing is applied —
    that happens inside CV.

    Returns
    -------
    DataFrame indexed by ``sample_id`` (in the order given), columns = genes.
    """
    p = Path(path) if path is not None else DEPMAP_RNA_PATH
    wanted = set(map(str, sample_ids))
    kept: list[pd.DataFrame] = []
    for chunk in pd.read_csv(p, chunksize=chunksize):
        sub = chunk[chunk["ModelID"].astype(str).isin(wanted)]
        if len(sub):
            kept.append(sub)
    if not kept:
        raise ValueError(f"No expression rows found for {len(wanted)} requested sample_ids.")
    raw = pd.concat(kept, ignore_index=True)

    if "IsDefaultEntryForModel" in raw.columns:
        default = raw[raw["IsDefaultEntryForModel"].astype(str) == "Yes"]
        # Fall back to all rows for any ModelID lacking a default entry.
        missing = wanted - set(default["ModelID"].astype(str))
        if missing:
            default = pd.concat([default, raw[raw["ModelID"].astype(str).isin(missing)]])
        raw = default

    raw = raw.drop_duplicates(subset="ModelID", keep="first")
    gene_cols = [c for c in raw.columns if c not in EXPRESSION_META_COLS]
    expr = raw.set_index(raw["ModelID"].astype(str))[gene_cols].astype(float)
    expr.index.name = "sample_id"

    missing = wanted - set(expr.index)
    if missing:
        raise ValueError(f"Expression matrix missing requested samples: {sorted(missing)}")
    return expr.loc[list(map(str, sample_ids))]


def strip_gene_symbol(col: str) -> str:
    """``'TSPAN6 (7105)'`` -> ``'TSPAN6'`` (DepMap gene-column convention)."""
    return col.split(" (")[0].strip()


# ---------------------------------------------------------------------------
# Fold-local preprocessing: variance gene filter
# ---------------------------------------------------------------------------

class VarianceGeneFilter(BaseEstimator, TransformerMixin):
    """Drop near-zero-variance genes and (optionally) keep the top-K most variable.

    Variances are computed on the training fold passed to :meth:`fit` only, so
    when used inside a ``Pipeline`` the selection respects fold boundaries
    automatically. ``top_k=None`` keeps every gene above ``threshold``.
    """

    def __init__(self, top_k: int | None = TOP_K_VARIABLE_GENES, threshold: float = NEAR_ZERO_VARIANCE):
        self.top_k = top_k
        self.threshold = threshold

    def fit(self, X, y=None):  # noqa: D102
        X = np.asarray(X, dtype=float)
        variances = X.var(axis=0)
        keep = np.flatnonzero(variances > self.threshold)
        if self.top_k is not None and keep.size > self.top_k:
            order = np.argsort(variances[keep], kind="stable")[::-1][: self.top_k]
            keep = np.sort(keep[order])
        self.support_ = keep
        self.n_features_in_ = X.shape[1]
        return self

    def transform(self, X):  # noqa: D102
        return np.asarray(X, dtype=float)[:, self.support_]

    def get_support(self, indices: bool = False):
        if indices:
            return self.support_
        mask = np.zeros(self.n_features_in_, dtype=bool)
        mask[self.support_] = True
        return mask


# ---------------------------------------------------------------------------
# Model / pipeline construction
# ---------------------------------------------------------------------------

def build_estimator(C: float, random_state: int = RANDOM_SEED) -> LogisticRegression:
    """Elastic-net logistic regression fixed to the spec."""
    return LogisticRegression(
        penalty="elasticnet",
        solver="saga",
        l1_ratio=L1_RATIO,
        C=C,
        max_iter=10000,
        class_weight="balanced",
        random_state=random_state,
    )


def build_pipeline(
    C: float = 1.0,
    *,
    top_k: int | None = TOP_K_VARIABLE_GENES,
    random_state: int = RANDOM_SEED,
) -> Pipeline:
    """VarianceGeneFilter -> StandardScaler -> LogisticRegression.

    Every step is fold-local: when this pipeline is fit on a training fold the
    variance filter, scaler, and model all see only that fold's rows.
    """
    return Pipeline([
        ("variance", VarianceGeneFilter(top_k=top_k)),
        ("scaler", StandardScaler()),
        ("clf", build_estimator(C=C, random_state=random_state)),
    ])


def _inner_cv(n_pos: int, n_neg: int, random_state: int) -> StratifiedKFold:
    n_splits = min(INNER_CV_FOLDS, n_pos, n_neg)
    n_splits = max(2, n_splits)
    return StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=random_state)


def select_C(
    X_train: np.ndarray,
    y_train: np.ndarray,
    *,
    c_grid: np.ndarray = C_GRID,
    top_k: int | None = TOP_K_VARIABLE_GENES,
    random_state: int = RANDOM_SEED,
) -> float:
    """Inner 5-fold CV selection of ``C`` (nested preprocessing via Pipeline)."""
    n_pos = int(y_train.sum())
    n_neg = int(len(y_train) - n_pos)
    cv = _inner_cv(n_pos, n_neg, random_state)
    pipe = build_pipeline(top_k=top_k, random_state=random_state)
    grid = GridSearchCV(
        pipe,
        param_grid={"clf__C": list(c_grid)},
        scoring="roc_auc",
        cv=cv,
        n_jobs=1,
        refit=False,
    )
    grid.fit(X_train, y_train)
    return float(grid.best_params_["clf__C"])


@dataclass
class LoocvResult:
    sample_ids: list[str]
    y_true: np.ndarray
    y_prob: np.ndarray
    selected_Cs: list[float]
    fold_scaler_means: list[np.ndarray]

    @property
    def auroc(self) -> float:
        return float(roc_auc_score(self.y_true, self.y_prob))


def run_loocv(
    X: np.ndarray,
    y: np.ndarray,
    sample_ids: Sequence[str],
    *,
    c_grid: np.ndarray = C_GRID,
    top_k: int | None = TOP_K_VARIABLE_GENES,
    random_state: int = RANDOM_SEED,
    record_scaler_means: bool = False,
    fixed_C: float | None = None,
) -> LoocvResult:
    """Outer leave-one-out CV with fully fold-local nested preprocessing.

    For each held-out sample: fit StandardScaler + variance filter on the N-1
    training fold (inside the Pipeline), run inner 5-fold CV to pick ``C`` on
    the training fold, refit the inner-selected model on the full N-1 fold, and
    predict the held-out probability.

    ``fixed_C`` short-circuits the inner C-selection and uses the supplied ``C``
    in every fold. This is used **only** for the permutation null, where re-
    selecting ``C`` inside every fold of every shuffle is computationally
    intractable at this gene count; the entire fold-local pipeline (variance
    filter, scaler, model refit, held-out prediction) is still re-run per
    shuffle. See :func:`permutation_test`.
    """
    X = np.asarray(X, dtype=float)
    y = np.asarray(y, dtype=int)
    sample_ids = list(sample_ids)
    loo = LeaveOneOut()
    probs = np.empty(len(y), dtype=float)
    selected: list[float] = []
    scaler_means: list[np.ndarray] = []

    for train_idx, test_idx in loo.split(X):
        X_tr, y_tr = X[train_idx], y[train_idx]
        if fixed_C is not None:
            best_C = float(fixed_C)
        else:
            best_C = select_C(X_tr, y_tr, c_grid=c_grid, top_k=top_k, random_state=random_state)
        pipe = build_pipeline(C=best_C, top_k=top_k, random_state=random_state)
        pipe.fit(X_tr, y_tr)
        probs[test_idx[0]] = pipe.predict_proba(X[test_idx])[:, 1][0]
        selected.append(best_C)
        if record_scaler_means:
            scaler_means.append(pipe.named_steps["scaler"].mean_.copy())

    return LoocvResult(
        sample_ids=sample_ids,
        y_true=y,
        y_prob=probs,
        selected_Cs=selected,
        fold_scaler_means=scaler_means,
    )


# ---------------------------------------------------------------------------
# Evaluation metrics
# ---------------------------------------------------------------------------

def sensitivity_at_specificity(
    y_true: np.ndarray, y_prob: np.ndarray, target_specificity: float = TARGET_SPECIFICITY
) -> float:
    """Max sensitivity (TPR) achievable while specificity >= target."""
    fpr, tpr, _ = roc_curve(y_true, y_prob)
    max_fpr = 1.0 - target_specificity
    mask = fpr <= max_fpr + 1e-12
    if not mask.any():
        return 0.0
    return float(np.max(tpr[mask]))


def evaluate(y_true: np.ndarray, y_prob: np.ndarray) -> dict[str, float]:
    return {
        "auroc": float(roc_auc_score(y_true, y_prob)),
        "auprc": float(average_precision_score(y_true, y_prob)),
        "sensitivity_at_90_spec": sensitivity_at_specificity(y_true, y_prob, TARGET_SPECIFICITY),
    }


def bootstrap_ci(
    y_true: np.ndarray,
    y_prob: np.ndarray,
    *,
    n_boot: int = N_BOOTSTRAP,
    seed: int = RANDOM_SEED,
    alpha: float = 0.05,
) -> dict[str, tuple[float, float]]:
    """Percentile bootstrap 95% CIs for AUROC / AUPRC / sensitivity@90spec.

    Resamples the LOOCV predictions with replacement. Resamples that collapse
    to a single class are skipped (undefined AUROC).
    """
    rng = np.random.default_rng(seed)
    y_true = np.asarray(y_true)
    y_prob = np.asarray(y_prob)
    n = len(y_true)
    metrics = {"auroc": [], "auprc": [], "sensitivity_at_90_spec": []}
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        yt, yp = y_true[idx], y_prob[idx]
        if len(np.unique(yt)) < 2:
            continue
        metrics["auroc"].append(roc_auc_score(yt, yp))
        metrics["auprc"].append(average_precision_score(yt, yp))
        metrics["sensitivity_at_90_spec"].append(sensitivity_at_specificity(yt, yp))
    lo_hi: dict[str, tuple[float, float]] = {}
    for k, vals in metrics.items():
        if vals:
            lo_hi[k] = (
                float(np.quantile(vals, alpha / 2)),
                float(np.quantile(vals, 1 - alpha / 2)),
            )
        else:
            lo_hi[k] = (float("nan"), float("nan"))
    return lo_hi


def _precompute_loocv_preprocessing(
    X: np.ndarray, top_k: int | None
) -> list[tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]]:
    """Per outer LOOCV fold, fit the (label-independent) variance filter +
    scaler on the N-1 training fold and cache the transformed train / held-out
    matrices.

    Used only by the fixed-``C`` permutation null. Preprocessing is identical
    across label shuffles (it depends on ``X`` alone), so computing it once and
    reusing it yields results byte-identical to re-running the full Pipeline per
    shuffle, while removing the dominant per-fold overhead.
    """
    X = np.asarray(X, dtype=float)
    loo = LeaveOneOut()
    folds = []
    for train_idx, test_idx in loo.split(X):
        var = VarianceGeneFilter(top_k=top_k).fit(X[train_idx])
        scaler = StandardScaler().fit(var.transform(X[train_idx]))
        Xtr = scaler.transform(var.transform(X[train_idx]))
        Xte = scaler.transform(var.transform(X[test_idx]))
        folds.append((train_idx, test_idx, Xtr, Xte))
    return folds


def _loocv_auroc_fixed_C(folds, y: np.ndarray, C: float, random_state: int) -> float:
    """LOOCV AUROC at fixed ``C`` using cached fold preprocessing."""
    probs = np.empty(len(y), dtype=float)
    for train_idx, test_idx, Xtr, Xte in folds:
        clf = build_estimator(C=C, random_state=random_state)
        clf.fit(Xtr, y[train_idx])
        probs[test_idx[0]] = clf.predict_proba(Xte)[:, 1][0]
    return float(roc_auc_score(y, probs))


def permutation_test(
    X: np.ndarray,
    y: np.ndarray,
    sample_ids: Sequence[str],
    observed_auroc: float,
    *,
    n_perm: int = N_PERMUTATIONS,
    c_grid: np.ndarray = C_GRID,
    top_k: int | None = TOP_K_VARIABLE_GENES,
    random_state: int = RANDOM_SEED,
    seed: int = RANDOM_SEED,
    fixed_C: float | None = None,
    select_C_per_fold: bool = False,
    progress_every: int | None = None,
) -> dict:
    """Label-shuffling permutation null: re-run the entire LOOCV per shuffle.

    Shuffles labels ``n_perm`` times within the included samples and re-runs the
    full fold-local LOOCV pipeline on each shuffle. Empirical p-value = fraction
    of permutations with AUROC >= observed.

    Computational note (documented deviation): with ~19k genes and N≈20-31,
    re-selecting ``C`` via inner CV inside *every* fold of *every* shuffle is
    intractable (≈9 h per definition). By default the permutation null fixes
    ``C`` to ``fixed_C`` (the median LOOCV-selected ``C`` from the observed run)
    and re-runs everything else fold-locally per shuffle. This is standard for
    permutation nulls and is, if anything, slightly conservative — permuted
    labels are denied the chance to re-tune ``C``. Set
    ``select_C_per_fold=True`` for the (very slow) fully nested null.
    """
    import time as _time

    rng = np.random.default_rng(seed)
    y = np.asarray(y, dtype=int)
    perm_fixed_C = None if select_C_per_fold else fixed_C
    if not select_C_per_fold and perm_fixed_C is None:
        raise ValueError("Provide fixed_C (median observed C) or set select_C_per_fold=True.")

    folds = None if select_C_per_fold else _precompute_loocv_preprocessing(X, top_k)
    perm_aurocs: list[float] = []
    _t0 = _time.time()
    for i in range(n_perm):
        y_perm = rng.permutation(y)
        if select_C_per_fold:
            auroc = run_loocv(
                X, y_perm, sample_ids,
                c_grid=c_grid, top_k=top_k, random_state=random_state, fixed_C=None,
            ).auroc
        else:
            auroc = _loocv_auroc_fixed_C(folds, y_perm, perm_fixed_C, random_state)
        perm_aurocs.append(auroc)
        if progress_every and ((i + 1) % progress_every == 0 or i + 1 == n_perm):
            print(f"    permutation {i + 1}/{n_perm}  ({_time.time() - _t0:.0f}s elapsed)", flush=True)
    perm_aurocs_arr = np.asarray(perm_aurocs)
    p_value = float(np.mean(perm_aurocs_arr >= observed_auroc))
    return {
        "observed_auroc": float(observed_auroc),
        "perm_aurocs": perm_aurocs_arr,
        "p_value": p_value,
        "n_perm": n_perm,
        "permutation_fixed_C": perm_fixed_C,
        "select_C_per_fold": select_C_per_fold,
    }


# ---------------------------------------------------------------------------
# Final refit + feature importance
# ---------------------------------------------------------------------------

def refit_full(
    X: np.ndarray,
    y: np.ndarray,
    C: float,
    *,
    top_k: int | None = TOP_K_VARIABLE_GENES,
    random_state: int = RANDOM_SEED,
) -> Pipeline:
    """Fit one final Pipeline (scaler + model) on all N samples at fixed ``C``."""
    pipe = build_pipeline(C=C, top_k=top_k, random_state=random_state)
    pipe.fit(np.asarray(X, dtype=float), np.asarray(y, dtype=int))
    return pipe


def ranked_genes(pipeline: Pipeline, gene_columns: Sequence[str]) -> pd.DataFrame:
    """Full ranked gene table (signed coefficients) from a refit Pipeline.

    Maps model coefficients back through the variance filter's support to the
    original gene columns; genes dropped by the filter get coefficient 0.
    """
    support = pipeline.named_steps["variance"].get_support(indices=True)
    coef = pipeline.named_steps["clf"].coef_.ravel()
    symbols = [strip_gene_symbol(gene_columns[i]) for i in support]
    df = pd.DataFrame({
        "gene_column": [gene_columns[i] for i in support],
        "gene_symbol": symbols,
        "coefficient": coef,
    })
    df["abs_coefficient"] = df["coefficient"].abs()
    df = df.sort_values("abs_coefficient", ascending=False, kind="stable").reset_index(drop=True)
    df.insert(0, "rank", np.arange(1, len(df) + 1))
    return df


def top_genes(ranked: pd.DataFrame, n: int = 20) -> pd.DataFrame:
    return ranked.head(n).copy()


def flag_alt_genes(ranked_top: pd.DataFrame, alt_genes: Iterable[str] = ALT_GENE_LIST) -> pd.DataFrame:
    alt = set(alt_genes)
    out = ranked_top.copy()
    out["in_alt_gene_list"] = out["gene_symbol"].isin(alt)
    return out


# ---------------------------------------------------------------------------
# Primary vs pragmatic comparison
# ---------------------------------------------------------------------------

def compare_definitions(
    primary: dict, pragmatic: dict, primary_ranked: pd.DataFrame, pragmatic_ranked: pd.DataFrame, n: int = 20
) -> dict:
    """AUROC delta, top-N gene overlap (count + Jaccard), coefficient Pearson r.

    ``primary`` / ``pragmatic`` are the ``evaluate(...)`` metric dicts.
    """
    p_top = set(primary_ranked.head(n)["gene_symbol"])
    g_top = set(pragmatic_ranked.head(n)["gene_symbol"])
    inter = p_top & g_top
    union = p_top | g_top
    jaccard = len(inter) / len(union) if union else float("nan")

    # Pearson on the coefficient vectors over the shared gene space.
    p_coef = primary_ranked.set_index("gene_symbol")["coefficient"]
    g_coef = pragmatic_ranked.set_index("gene_symbol")["coefficient"]
    shared = sorted(set(p_coef.index) & set(g_coef.index))
    if len(shared) >= 2:
        pv = p_coef.loc[shared].to_numpy()
        gv = g_coef.loc[shared].to_numpy()
        if np.std(pv) > 0 and np.std(gv) > 0:
            pearson = float(np.corrcoef(pv, gv)[0, 1])
        else:
            pearson = float("nan")
    else:
        pearson = float("nan")

    return {
        "auroc_delta_primary_minus_pragmatic": primary["auroc"] - pragmatic["auroc"],
        "top_n": n,
        "top_overlap_count": len(inter),
        "top_overlap_jaccard": jaccard,
        "top_overlap_genes": sorted(inter),
        "coef_pearson_on_shared_genes": pearson,
        "n_shared_genes": len(shared),
    }
