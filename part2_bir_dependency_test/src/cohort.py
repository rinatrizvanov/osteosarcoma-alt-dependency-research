"""Cohort construction: osteosarcoma inclusion filter and frozen ALT label join.

Pre-specification Sections 2 and 3. This module deliberately handles **no outcome
data**. It probes ``CRISPRGeneEffect.csv`` for *row membership only* (which
ModelIDs were screened) and never reads a gene-effect value; the outcome is
first loaded in :mod:`src.dependency_stats`, used from notebook 05.

Conventions carried over from Aim 2 (``references/aim2_signature.py``):

* Filters are applied **explicitly and in order**, with N logged after each step
  (:class:`FilterTrace`). The manifest's ``include_*`` flags are intentionally
  not trusted.
* Arm labels come from ``alt_call_direction``, **never** ``alt_status`` — using
  ``alt_status`` would corrupt the silver split.
* Ground-truth anchors are asserted before any analysis, and a failure aborts.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

from .utils import (
    CRISPR_PATH,
    MANIFEST_PATH,
    MODEL_PATH,
    RNA_PATH,
    verify_expected,
)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

TARGET_COHORT = "DepMap_26Q1"
DIRECTION_COLUMN = "alt_call_direction"
VALID_DIRECTIONS: tuple[str, ...] = ("positive", "negative")

#: Matches both "Osteosarcoma" and any "...osteosarcoma" subtype, case-safe on
#: the leading character (Aim 2 used the same 'steosarcoma' substring idiom).
OS_SUBSTRING = "steosarcoma"

#: alt_status tiers admitted to Stage 1. Primary ALT+ = gold + silver.
ALT_POSITIVE_STATUSES: tuple[str, ...] = ("alt_gold_positive", "alt_silver")
ALT_NEGATIVE_STATUSES: tuple[str, ...] = ("alt_negative_high_conf", "alt_silver")
ADMITTED_STATUSES: tuple[str, ...] = (
    "alt_gold_positive",
    "alt_silver",
    "alt_negative_high_conf",
)

#: Ground-truth anchors — abort if either fails (Aim 2 GROUND_TRUTH_DIRECTION).
GROUND_TRUTH_DIRECTION: dict[str, str] = {
    "ACH-000364": "positive",  # U2OS, canonical ALT-positive
    "ACH-000748": "negative",  # SJSA1, canonical ALT-negative / telomerase+
}

#: Hard expectations established during planning from membership counts alone.
#: A mismatch means the release or the manifest changed — abort, do not adapt.
EXPECTED_SET_SIZES: dict[str, int] = {
    "alt_positive": 16,   # 10 gold+ and 6 silver+
    "alt_negative": 9,    # 7 gold- and 2 silver-
    "excluded_unknown": 1,
}
EXPECTED_TIER_COUNTS: dict[tuple[str, str], int] = {
    ("alt_gold_positive", "positive"): 10,
    ("alt_silver", "positive"): 6,
    ("alt_negative_high_conf", "negative"): 7,
    ("alt_silver", "negative"): 2,
}

EXPRESSION_META_COLS: frozenset[str] = frozenset({
    "Unnamed: 0", "SequencingID", "ModelConditionID", "ModelID",
    "IsDefaultEntryForMC", "IsDefaultEntryForModel",
})


# ---------------------------------------------------------------------------
# Membership probes (no values read)
# ---------------------------------------------------------------------------


def crispr_model_ids(path: Path | None = None) -> set[str]:
    """ModelIDs with a CRISPR screen — **row index only, no gene effects read**.

    ``CRISPRGeneEffect.csv`` has an unnamed first column holding the ModelID.
    Streaming the first field of each line avoids materialising a 440 MB,
    18,531-column matrix just to ask which lines were screened.
    """
    p = Path(path) if path is not None else CRISPR_PATH
    ids: set[str] = set()
    with open(p) as fh:
        next(fh)  # header
        for line in fh:
            ids.add(line.split(",", 1)[0])
    ids.discard("")
    return ids


def rna_model_ids(path: Path | None = None) -> set[str]:
    """ModelIDs present in the expression matrix (membership probe only)."""
    p = Path(path) if path is not None else RNA_PATH
    with open(p) as fh:
        header = next(fh).rstrip("\n").split(",")
        idx = header.index("ModelID")
        ids = {line.split(",")[idx] for line in fh}
    ids.discard("")
    return ids


def osteosarcoma_models(path: Path | None = None) -> pd.DataFrame:
    """Osteosarcoma rows of ``Model.csv`` with the metadata Stage 1 uses."""
    p = Path(path) if path is not None else MODEL_PATH
    model = pd.read_csv(p, dtype=str).fillna("")
    mask = (
        model["OncotreePrimaryDisease"].str.contains(OS_SUBSTRING, case=True, na=False)
        | model["OncotreeSubtype"].str.contains(OS_SUBSTRING, case=True, na=False)
    )
    keep = [
        "ModelID", "CellLineName", "StrippedCellLineName", "OncotreeLineage",
        "OncotreePrimaryDisease", "OncotreeSubtype", "Sex", "Age",
        "PrimaryOrMetastasis", "CCLEName",
    ]
    return model.loc[mask, keep].reset_index(drop=True)


# ---------------------------------------------------------------------------
# Filter trace (Aim 2 pattern)
# ---------------------------------------------------------------------------


@dataclass
class FilterTrace:
    """Records N surviving after each inclusion-filter step."""

    label: str
    steps: list[tuple[str, int]] = field(default_factory=list)

    def add(self, step: str, n: int) -> None:
        self.steps.append((step, n))

    def as_lines(self) -> list[str]:
        return [f"  {step:<62} N = {n}" for step, n in self.steps]

    def as_markdown(self) -> list[str]:
        out = ["| Step | N |", "|---|---|"]
        out += [f"| {step} | {n} |" for step, n in self.steps]
        return out


# ---------------------------------------------------------------------------
# Cohort construction
# ---------------------------------------------------------------------------


@dataclass
class Stage1Cohort:
    """The three analysis sets plus the trace that produced them."""

    alt_positive: pd.DataFrame
    alt_negative: pd.DataFrame
    excluded_unknown: pd.DataFrame
    trace: FilterTrace

    @property
    def sizes(self) -> dict[str, int]:
        return {
            "alt_positive": len(self.alt_positive),
            "alt_negative": len(self.alt_negative),
            "excluded_unknown": len(self.excluded_unknown),
        }

    def all_screened(self) -> pd.DataFrame:
        """ALT+ and ALT- combined — the 25 lines entering any contrast."""
        return pd.concat([self.alt_positive, self.alt_negative], ignore_index=True)


def build_cohort(
    manifest_path: Path | None = None,
    *,
    verify: bool = True,
) -> Stage1Cohort:
    """Apply the Section 2 inclusion filter and join frozen Section 3 labels.

    Steps, in order and each logged:

    1. ``Model.csv`` osteosarcoma rows
    2. present in ``CRISPRGeneEffect.csv`` (screened)
    3. present in the expression matrix
    4. ``cohort == 'DepMap_26Q1'`` in the frozen manifest
    5. ``alt_status`` in the admitted tiers
    6. split on ``alt_call_direction``

    Parameters
    ----------
    verify
        sha256-verify every input against :data:`src.utils.EXPECTED_SHA256`.
        Leave on except in unit tests with synthetic inputs.
    """
    mpath = Path(manifest_path) if manifest_path is not None else MANIFEST_PATH
    if verify:
        for p in (MODEL_PATH, RNA_PATH, mpath):
            verify_expected(p)

    trace = FilterTrace(label="stage1_cohort")

    os_models = osteosarcoma_models()
    trace.add("Model.csv rows matching 'steosarcoma'", len(os_models))

    screened = crispr_model_ids()
    df = os_models[os_models["ModelID"].isin(screened)].copy()
    trace.add("present in CRISPRGeneEffect.csv (screened)", len(df))

    expressed = rna_model_ids()
    df = df[df["ModelID"].isin(expressed)]
    trace.add("present in expression matrix", len(df))

    manifest = pd.read_csv(mpath, dtype=str).fillna("")
    manifest = manifest[manifest["cohort"] == TARGET_COHORT]
    label_cols = [
        "sample_id", "alt_status", DIRECTION_COLUMN, "alt_source",
        "source_confidence", "atrx_status", "daxx_status", "telomerase_status",
        "tp53_status", "tert_expression", "notes",
    ]
    df = df.merge(
        manifest[label_cols], how="inner", left_on="ModelID", right_on="sample_id"
    )
    trace.add(f"joined to frozen manifest (cohort == {TARGET_COHORT!r})", len(df))

    unknown = df[~df["alt_status"].isin(ADMITTED_STATUSES)].reset_index(drop=True)
    df = df[df["alt_status"].isin(ADMITTED_STATUSES)]
    trace.add(f"alt_status in {ADMITTED_STATUSES}", len(df))

    df = df[df[DIRECTION_COLUMN].isin(VALID_DIRECTIONS)].reset_index(drop=True)
    trace.add(f"{DIRECTION_COLUMN} in {VALID_DIRECTIONS}", len(df))

    pos = df[df[DIRECTION_COLUMN] == "positive"].reset_index(drop=True)
    neg = df[df[DIRECTION_COLUMN] == "negative"].reset_index(drop=True)
    trace.add("ALT+ (primary contrast set)", len(pos))
    trace.add("ALT- (negative-control set)", len(neg))

    return Stage1Cohort(
        alt_positive=pos, alt_negative=neg, excluded_unknown=unknown, trace=trace
    )


# ---------------------------------------------------------------------------
# Assertions — run before any analysis
# ---------------------------------------------------------------------------


def assert_anchors(cohort: Stage1Cohort) -> None:
    """ACH-000364 (U2OS) positive, ACH-000748 (SJSA1) negative. Abort on failure."""
    by_id = (
        cohort.all_screened().set_index("ModelID")[DIRECTION_COLUMN].to_dict()
    )
    for model_id, expected in GROUND_TRUTH_DIRECTION.items():
        if model_id not in by_id:
            raise ValueError(
                f"Ground-truth anchor {model_id} missing from the screened cohort "
                f"(expected {DIRECTION_COLUMN}={expected!r})."
            )
        if by_id[model_id] != expected:
            raise ValueError(
                f"Ground-truth assertion failed for {model_id}: "
                f"{DIRECTION_COLUMN}={by_id[model_id]!r}, expected {expected!r}."
            )


def assert_expected_sizes(cohort: Stage1Cohort) -> None:
    """Abort with the exact discrepancy if any set size differs from plan."""
    actual = cohort.sizes
    if actual != EXPECTED_SET_SIZES:
        raise ValueError(
            f"Cohort size mismatch.\n  expected {EXPECTED_SET_SIZES}\n  actual   {actual}\n"
            "The DepMap release or the frozen manifest changed. Do not adapt the "
            "expectation to the data — investigate the cause first."
        )


def tier_counts(cohort: Stage1Cohort) -> dict[tuple[str, str], int]:
    df = cohort.all_screened()
    return {
        (status, direction): int(n)
        for (status, direction), n in
        df.groupby(["alt_status", DIRECTION_COLUMN]).size().items()
    }


def assert_expected_tiers(cohort: Stage1Cohort) -> None:
    actual = tier_counts(cohort)
    if actual != EXPECTED_TIER_COUNTS:
        raise ValueError(
            f"ALT tier x direction mismatch.\n  expected {EXPECTED_TIER_COUNTS}\n"
            f"  actual   {actual}"
        )


def assert_all(cohort: Stage1Cohort) -> None:
    """Every pre-analysis assertion, in the order the prespec states them."""
    assert_anchors(cohort)
    assert_expected_sizes(cohort)
    assert_expected_tiers(cohort)


def gold_only(frame: pd.DataFrame) -> pd.DataFrame:
    """Gold-tier subset — the R1 sensitivity cohort (prespec Section 7)."""
    return frame[frame["alt_status"] != "alt_silver"].reset_index(drop=True)
