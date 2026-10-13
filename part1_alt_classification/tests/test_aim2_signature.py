"""Tests for ``src.aim2_signature`` (Aim 2 ALT signature score helpers).

Covers the inclusion filter (primary vs pragmatic), the direction-based arm
labels, the ground-truth assertions (including a corrupted-manifest abort), and
that the modelling Pipeline fits its scaler fold-locally inside LOOCV.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from src import aim2_signature as a2  # noqa: E402


# ---------------------------------------------------------------------------
# Synthetic manifest fixture covering every relevant edge case
# ---------------------------------------------------------------------------

def _row(sample_id, cohort, alt_status, direction):
    return {
        "sample_id": sample_id,
        "cohort": cohort,
        "alt_status": alt_status,
        "alt_call_direction": direction,
    }


@pytest.fixture
def synthetic_manifest() -> pd.DataFrame:
    """One DepMap row per edge case plus distractors.

    Expected (all DepMap + in RNA matrix unless noted):
      * silver+, silver-       -> pragmatic only
      * gold+, gold-           -> primary AND pragmatic
      * alt_unknown            -> excluded (direction 'unknown')
      * non-OS cohort          -> excluded (wrong cohort)
      * missing-RNA DepMap row -> excluded (not in ModelID index)
    """
    rows = [
        # Ground-truth anchors (gold) — primary + pragmatic, positive/negative.
        _row("ACH-000364", "DepMap_26Q1", "alt_gold_positive", "positive"),
        _row("ACH-000748", "DepMap_26Q1", "alt_negative_high_conf", "negative"),
        # Extra gold pair so primary has >1 per arm.
        _row("ACH-GOLD-P", "DepMap_26Q1", "alt_gold_positive", "positive"),
        _row("ACH-GOLD-N", "DepMap_26Q1", "alt_negative_high_conf", "negative"),
        # Silver pair — pragmatic only, direction drives the arm.
        _row("ACH-SILV-P", "DepMap_26Q1", "alt_silver", "positive"),
        _row("ACH-SILV-N", "DepMap_26Q1", "alt_silver", "negative"),
        # alt_unknown — excluded by the direction filter.
        _row("ACH-UNK", "DepMap_26Q1", "alt_unknown", "unknown"),
        # Wrong cohort — excluded.
        _row("OS-NONDEPMAP", "PCGP_2014", "alt_gold_positive", "positive"),
        # DepMap gold but absent from the RNA ModelID index — excluded.
        _row("ACH-NORNA", "DepMap_26Q1", "alt_gold_positive", "positive"),
    ]
    return pd.DataFrame(rows)


@pytest.fixture
def model_ids() -> set[str]:
    """RNA-matrix ModelIDs: everything except the missing-RNA row."""
    return {
        "ACH-000364", "ACH-000748", "ACH-GOLD-P", "ACH-GOLD-N",
        "ACH-SILV-P", "ACH-SILV-N", "ACH-UNK",
        # ACH-NORNA deliberately absent.
    }


# ---------------------------------------------------------------------------
# Inclusion filter
# ---------------------------------------------------------------------------

def test_primary_filter_selects_expected(synthetic_manifest, model_ids):
    df, trace = a2.apply_inclusion_filter(synthetic_manifest, model_ids, "primary")
    assert set(df["sample_id"]) == {"ACH-000364", "ACH-000748", "ACH-GOLD-P", "ACH-GOLD-N"}
    counts = a2.direction_counts(df)
    assert counts == {"positive": 2, "negative": 2}
    # Trace records a value after every step and is monotonically non-increasing.
    ns = [n for _, n in trace.steps]
    assert ns == sorted(ns, reverse=True)


def test_pragmatic_filter_adds_silver(synthetic_manifest, model_ids):
    df, _ = a2.apply_inclusion_filter(synthetic_manifest, model_ids, "pragmatic")
    assert set(df["sample_id"]) == {
        "ACH-000364", "ACH-000748", "ACH-GOLD-P", "ACH-GOLD-N", "ACH-SILV-P", "ACH-SILV-N",
    }
    assert a2.direction_counts(df) == {"positive": 3, "negative": 3}


def test_filter_excludes_unknown_nondepmap_and_missing_rna(synthetic_manifest, model_ids):
    df, _ = a2.apply_inclusion_filter(synthetic_manifest, model_ids, "pragmatic")
    for excluded in ("ACH-UNK", "OS-NONDEPMAP", "ACH-NORNA"):
        assert excluded not in set(df["sample_id"])


def test_unknown_definition_raises(synthetic_manifest, model_ids):
    with pytest.raises(ValueError):
        a2.apply_inclusion_filter(synthetic_manifest, model_ids, "broad")


# ---------------------------------------------------------------------------
# Arm labels come from alt_call_direction, never alt_status
# ---------------------------------------------------------------------------

def test_labels_from_direction_not_status(synthetic_manifest, model_ids):
    df, _ = a2.apply_inclusion_filter(synthetic_manifest, model_ids, "pragmatic")
    y = a2.labels_from_direction(df)
    # Map back: every positive label must be a 'positive' direction row.
    for sid, label in zip(df["sample_id"], y):
        dirn = df.loc[df["sample_id"] == sid, "alt_call_direction"].iloc[0]
        assert label == (1 if dirn == "positive" else 0)
    # A silver row's arm is its direction, independent of the shared alt_status.
    silver = df[df["alt_status"] == "alt_silver"]
    assert set(a2.labels_from_direction(silver)) == {0, 1}


def test_expected_counts_assertion_message(synthetic_manifest, model_ids):
    # The synthetic cohort does not match the real 12/8 expectation -> ValueError.
    df, _ = a2.apply_inclusion_filter(synthetic_manifest, model_ids, "primary")
    with pytest.raises(ValueError) as exc:
        a2.assert_expected_counts(df, "primary")
    assert "expected positive=12" in str(exc.value)


# ---------------------------------------------------------------------------
# Ground-truth assertions
# ---------------------------------------------------------------------------

def test_ground_truth_passes_on_clean(synthetic_manifest, model_ids):
    df, _ = a2.apply_inclusion_filter(synthetic_manifest, model_ids, "primary")
    a2.assert_ground_truth(df)  # should not raise


def test_corrupted_manifest_triggers_ground_truth_abort(synthetic_manifest, model_ids):
    """U2OS mis-labelled as negative must abort."""
    corrupted = synthetic_manifest.copy()
    corrupted.loc[corrupted["sample_id"] == "ACH-000364", "alt_call_direction"] = "negative"
    df, _ = a2.apply_inclusion_filter(corrupted, model_ids, "primary")
    with pytest.raises(ValueError) as exc:
        a2.assert_ground_truth(df)
    assert "ACH-000364" in str(exc.value)


def test_per_row_direction_assertion(synthetic_manifest, model_ids):
    """An included row with a non positive/negative direction aborts."""
    df, _ = a2.apply_inclusion_filter(synthetic_manifest, model_ids, "primary")
    df.loc[0, "alt_call_direction"] = "unknown"
    with pytest.raises(ValueError):
        a2.assert_ground_truth(df)


# ---------------------------------------------------------------------------
# Pipeline fits the scaler fold-locally inside LOOCV
# ---------------------------------------------------------------------------

def _toy_dataset(n=12, p=8, seed=0):
    rng = np.random.default_rng(seed)
    y = np.array([0, 1] * (n // 2))
    X = rng.normal(size=(n, p)) + y[:, None] * 1.5  # signal so folds are non-degenerate
    return X, y


def test_loocv_fits_scaler_within_each_fold():
    X, y = _toy_dataset()
    # top_k=None keeps all genes so scaler.mean_ vectors are directly comparable.
    res = a2.run_loocv(
        X, y, [f"s{i}" for i in range(len(y))],
        top_k=None, record_scaler_means=True,
    )
    means = res.fold_scaler_means
    assert len(means) == len(y)
    # Each fold leaves out a different row, so the training-fold means must
    # differ across folds — proves the scaler is fit per fold, not globally.
    distinct = {tuple(np.round(m, 10)) for m in means}
    assert len(distinct) == len(y)
    # And none equals the global mean (which would indicate a leak).
    global_mean = X.mean(axis=0)
    for m in means:
        assert not np.allclose(m, global_mean)


def test_pipeline_preprocessing_is_fold_local():
    X, y = _toy_dataset()
    pipe = a2.build_pipeline(C=1.0, top_k=None)
    train_idx = np.arange(2, len(y))
    pipe.fit(X[train_idx], y[train_idx])
    # Scaler mean is over the training fold only, not all rows.
    assert np.allclose(pipe.named_steps["scaler"].mean_, X[train_idx].mean(axis=0))
    assert not np.allclose(pipe.named_steps["scaler"].mean_, X.mean(axis=0))


def test_variance_gene_filter_topk_and_threshold():
    X = np.array([
        [0.0, 1.0, 5.0, 0.0],
        [0.0, 2.0, 1.0, 0.0],
        [0.0, 3.0, 9.0, 0.0],
    ])
    # Columns 0 and 3 are constant (zero variance) -> dropped.
    f = a2.VarianceGeneFilter(top_k=None).fit(X)
    assert set(f.get_support(indices=True)) == {1, 2}
    # top_k=1 keeps only the highest-variance column (col 2).
    f1 = a2.VarianceGeneFilter(top_k=1).fit(X)
    assert list(f1.get_support(indices=True)) == [2]


def test_ranked_genes_maps_through_support():
    X, y = _toy_dataset(p=6)
    gene_cols = [f"G{i} ({100 + i})" for i in range(6)]
    pipe = a2.refit_full(X, y, C=1.0, top_k=4)
    ranked = a2.ranked_genes(pipe, gene_cols)
    assert len(ranked) == 4  # only the 4 kept genes
    assert list(ranked["rank"]) == [1, 2, 3, 4]
    assert (ranked["abs_coefficient"].values[:-1] >= ranked["abs_coefficient"].values[1:]).all()
    assert ranked["gene_symbol"].str.startswith("G").all()
