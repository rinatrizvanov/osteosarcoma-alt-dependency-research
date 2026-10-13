"""
Read-only helpers for ``notebooks/08_recovery_route_inventory.ipynb``.

Tests three RNA-recovery routes that notebooks 06/07 did not evaluate because
those notebooks only matched exact ``sample_id`` values against on-disk files.

Route B — cell-line-name → DepMap ``ModelID`` join (Taglialatela cohort and
manifest samples whose ``sample_id`` appears in the ALTitude DepMap panel).

Route A — de Nonneville expression availability (documentation only; no download).

Route C — Lawlor internal RNA (documentation only).

The combined "Route B only" scenario reuses the Aim 1 / Aim 2 verdict machinery
from :mod:`src.aim1_inventory` and :mod:`src.aim2_inventory`.
"""
from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

from .aim1_inventory import (
    DEPMAP_RNA_DATA_TYPE,
    DEPMAP_RNA_PATH,
    WITHIN_COHORT_PER_ARM_MIN,
    decide_aim1_design,
    depmap_rna_model_ids,
    select_pragmatic_contrast_rows,
)
from .aim2_inventory import (
    ALT_NEGATIVE_LABEL,
    ALT_POSITIVE_LABEL,
    decide_aim2_design,
    select_aim2_contrast_rows,
)
from .audit import classify_aim1_arms

DATA_ROOT = Path(__file__).parent.parent / "data"
DEPMAP_DIR = DATA_ROOT / "raw" / "depmap"
MODEL_CSV_PATH = DEPMAP_DIR / "Model.csv"

TAGLIALATELA_COHORT = "Taglialatela_2026"
LAWLOR_COHORT = "Lawlor_2026_GCCRI"
DE_NONNEVILLE_COHORT = "de_Nonneville_2022"

ARRAYEXPRESS_DE_NONNEVILLE = "E-MTAB-11743"

# Columns in Model.csv used for name-based joins (in priority order for reporting).
MODEL_NAME_COLUMNS = (
    "CellLineName",
    "StrippedCellLineName",
    "CCLEName",
    "ModelIDAlias",
)


def normalize_cell_line_name(name: str) -> str:
    """Uppercase alphanumerics only — strips spaces, hyphens, punctuation."""
    if name is None or (isinstance(name, float) and pd.isna(name)):
        return ""
    return re.sub(r"[^A-Z0-9]", "", str(name).upper())


def build_depmap_name_lookup(model_df: pd.DataFrame) -> dict[str, tuple[str, str]]:
    """Map normalized cell-line tokens → (ModelID, matched_column).

  First exact hit wins per normalized key. Multiple ModelIDs sharing a key
  are not stored (ambiguous); callers should treat those as unmapped.
    """
    lookup: dict[str, tuple[str, str]] = {}
    ambiguous: set[str] = set()

    for _, row in model_df.iterrows():
        model_id = str(row["ModelID"])
        for col in MODEL_NAME_COLUMNS:
            if col not in row.index:
                continue
            raw = row[col]
            if pd.isna(raw) or not str(raw).strip():
                continue
            key = normalize_cell_line_name(raw)
            if not key:
                continue
            if key in lookup and lookup[key][0] != model_id:
                ambiguous.add(key)
            else:
                lookup[key] = (model_id, col)

    for key in ambiguous:
        lookup.pop(key, None)

    return lookup


def _prefix_match_model_id(
    query_norm: str,
    model_df: pd.DataFrame,
) -> tuple[str | None, str | None]:
    """Fallback: unique StrippedCellLineName whose normalized form starts with query."""
    if len(query_norm) < 3:
        return None, None
    hits: list[tuple[str, str]] = []
    for _, row in model_df.iterrows():
        stripped = normalize_cell_line_name(row.get("StrippedCellLineName", ""))
        if not stripped:
            continue
        if stripped.startswith(query_norm) or query_norm == stripped[: len(query_norm)]:
            hits.append((str(row["ModelID"]), "StrippedCellLineName(prefix)"))
    if len(hits) == 1:
        return hits[0]
    return None, None


def map_cell_line_to_model_id(
    cell_line_name: str,
    lookup: dict[str, tuple[str, str]],
    model_df: pd.DataFrame,
) -> tuple[str | None, str | None]:
    """Resolve a cell-line string to a DepMap ModelID, or (None, None)."""
    key = normalize_cell_line_name(cell_line_name)
    if not key:
        return None, None
    if key in lookup:
        return lookup[key]
    return _prefix_match_model_id(key, model_df)


def select_route_b_candidates(
    manifest: pd.DataFrame,
    altitude_sample_ids: set[str],
) -> pd.DataFrame:
    """Manifest rows eligible for Route B name → ModelID recovery.

    Includes:
      * all ``Taglialatela_2026`` cohort rows;
      * any row whose ``sample_id`` is in the ALTitude DepMap panel (typically
        DepMap osteosarcoma lines, even if manifest ``alt_source`` does not
        mention ALTitude explicitly).
    """
    is_tag = manifest["cohort"] == TAGLIALATELA_COHORT
    is_altitude = manifest["sample_id"].astype(str).isin(altitude_sample_ids)
    sub = manifest[is_tag | is_altitude].copy()
    sub["route_b_reason"] = "Taglialatela_2026"
    sub.loc[is_altitude & ~is_tag, "route_b_reason"] = "ALTitude_DepMap_panel"
    sub.loc[is_tag & is_altitude, "route_b_reason"] = "Taglialatela_2026+ALTitude"
    return sub.reset_index(drop=True)


def build_route_b_table(
    candidates: pd.DataFrame,
    lookup: dict[str, tuple[str, str]],
    model_df: pd.DataFrame,
    depmap_rna_models: set[str],
) -> pd.DataFrame:
    """Per-sample Route B mapping table."""
    rows = []
    for _, r in candidates.iterrows():
        sid = str(r["sample_id"])
        cohort = str(r["cohort"])
        if cohort == TAGLIALATELA_COHORT or not sid.startswith("ACH-"):
            cell_line = sid
        else:
            # DepMap manifest row: prefer Model.csv name for display.
            mrow = model_df[model_df["ModelID"] == sid]
            if len(mrow):
                cell_line = str(mrow.iloc[0].get("StrippedCellLineName") or sid)
            else:
                cell_line = sid

        if sid.startswith("ACH-") and sid in depmap_rna_models:
            mapped_id = sid
            match_col = "sample_id(exact ModelID)"
        else:
            mapped_id, match_col = map_cell_line_to_model_id(cell_line, lookup, model_df)
            if mapped_id is None and sid.startswith("ACH-"):
                mapped_id = sid
                match_col = "sample_id(fallback)"

        in_rna = bool(mapped_id and mapped_id in depmap_rna_models)
        alt_status = str(r.get("alt_status", "unknown"))
        alt_source = str(r.get("alt_source", "unknown"))

        rows.append({
            "cell_line_name": cell_line,
            "manifest_sample_id": sid,
            "cohort": cohort,
            "mapped_model_id": mapped_id or "",
            "match_method": match_col or "",
            "in_depmap_rna": in_rna,
            "alt_status": alt_status,
            "source": alt_source,
            "route_b_reason": str(r.get("route_b_reason", "")),
        })
    return pd.DataFrame(rows)


def route_b_alt_rna_counts(route_b: pd.DataFrame) -> dict[str, int]:
    """Count ALT arms that gain file-backed RNA via Route B."""
    fb = route_b[route_b["in_depmap_rna"]]
    return {
        "n_mapped_total": int(route_b["mapped_model_id"].astype(bool).sum()),
        "n_in_depmap_rna": len(fb),
        "n_alt_positive_rna": int((fb["alt_status"] == ALT_POSITIVE_LABEL).sum()),
        "n_alt_negative_rna": int((fb["alt_status"] == ALT_NEGATIVE_LABEL).sum()),
        "n_alt_other_rna": int(
            (~fb["alt_status"].isin([ALT_POSITIVE_LABEL, ALT_NEGATIVE_LABEL])).sum()
        ),
    }


def summarize_de_nonneville(manifest: pd.DataFrame) -> dict:
    """Route A — label counts only (no download)."""
    sub = manifest[manifest["cohort"] == DE_NONNEVILLE_COHORT]
    n_pos = int((sub["alt_status"] == ALT_POSITIVE_LABEL).sum())
    n_neg = int((sub["alt_status"] == ALT_NEGATIVE_LABEL).sum())
    n_other = len(sub) - n_pos - n_neg
    return {
        "n_total": len(sub),
        "n_alt_positive": n_pos,
        "n_alt_negative": n_neg,
        "n_other_alt_status": n_other,
        "sample_ids_positive": sorted(
            sub.loc[sub["alt_status"] == ALT_POSITIVE_LABEL, "sample_id"].tolist()
        ),
        "sample_ids_negative": sorted(
            sub.loc[sub["alt_status"] == ALT_NEGATIVE_LABEL, "sample_id"].tolist()
        ),
        "arrayexpress_accession": ARRAYEXPRESS_DE_NONNEVILLE,
    }


def summarize_lawlor(manifest: pd.DataFrame) -> dict:
    """Route C — PDX samples with ALT and/or architecture labels."""
    sub = manifest[manifest["cohort"] == LAWLOR_COHORT].copy()
    has_alt = sub["alt_status"].isin([ALT_POSITIVE_LABEL, ALT_NEGATIVE_LABEL])
    has_arch = sub["sv_architecture_17p13_1"].fillna("unknown").str.lower() != "unknown"
    labeled = sub[has_alt | has_arch]

    arms = classify_aim1_arms(manifest, "pragmatic")
    sub["aim1_arm"] = arms.loc[sub.index].values

    return {
        "n_total_lawlor": len(sub),
        "n_with_alt_or_architecture": len(labeled),
        "n_alt_positive": int((sub["alt_status"] == ALT_POSITIVE_LABEL).sum()),
        "n_alt_negative": int((sub["alt_status"] == ALT_NEGATIVE_LABEL).sum()),
        "n_tcab1_disrupted_high_conf": int(
            (sub["tcab1_overall_status"] == "disrupted_high_conf").sum()
        ),
        "n_aim1_pragmatic_disrupted": int((sub["aim1_arm"] == "disrupted").sum()),
        "n_aim1_pragmatic_intact": int((sub["aim1_arm"] == "intact").sum()),
        "sample_ids_all": sorted(sub["sample_id"].tolist()),
        "sample_ids_alt_positive": sorted(
            sub.loc[sub["alt_status"] == ALT_POSITIVE_LABEL, "sample_id"].tolist()
        ),
        "sample_ids_alt_negative": sorted(
            sub.loc[sub["alt_status"] == ALT_NEGATIVE_LABEL, "sample_id"].tolist()
        ),
        "sample_ids_architecture_labeled": sorted(
            sub.loc[has_arch, "sample_id"].tolist()
        ),
    }


def _route_b_model_ids_for_manifest_row(
    sample_id: str,
    route_b: pd.DataFrame,
) -> str | None:
    hit = route_b[
        (route_b["manifest_sample_id"] == sample_id) & route_b["in_depmap_rna"]
    ]
    if len(hit):
        return str(hit.iloc[0]["mapped_model_id"])
    return None


def build_route_b_only_aim1_inventory(
    manifest: pd.DataFrame,
    route_b: pd.DataFrame,
    *,
    direct_file_backed: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Inventory for Aim 1 pragmatic contrast under Route-B-only RNA recovery.

    ``file_backed_rna`` is True when the sample already had direct DepMap RNA
    (``sample_id`` in expression matrix) OR Route B maps it to an in-matrix
    ModelID. All recovered samples share ``data_type`` = DepMap log2TPM.
    """
    contrast = select_pragmatic_contrast_rows(manifest)
    depmap_rna = depmap_rna_model_ids()
    rows = []
    for _, r in contrast.iterrows():
        sid = str(r["sample_id"])
        direct = sid in depmap_rna
        rb_id = _route_b_model_ids_for_manifest_row(sid, route_b)
        route_b_ok = rb_id is not None
        rows.append({
            "sample_id": sid,
            "cohort": str(r["cohort"]),
            "tcab1_overall_status": r["arm"],
            "sample_type": str(r.get("sample_type", "unknown")),
            "file_backed_rna": bool(direct or route_b_ok),
            "recovery_route": (
                "direct+route_b" if direct and route_b_ok
                else "direct" if direct
                else "route_b" if route_b_ok
                else "none"
            ),
            "mapped_model_id": rb_id or (sid if direct else ""),
            "data_type": DEPMAP_RNA_DATA_TYPE if (direct or route_b_ok) else "unknown",
        })
    return pd.DataFrame(rows)


def build_route_b_only_aim2_inventory(
    manifest: pd.DataFrame,
    route_b: pd.DataFrame,
) -> pd.DataFrame:
    """Inventory for Aim 2 contrast under Route-B-only RNA recovery."""
    contrast = select_aim2_contrast_rows(manifest)
    depmap_rna = depmap_rna_model_ids()
    rows = []
    for _, r in contrast.iterrows():
        sid = str(r["sample_id"])
        direct = sid in depmap_rna
        rb_id = _route_b_model_ids_for_manifest_row(sid, route_b)
        route_b_ok = rb_id is not None
        rows.append({
            "sample_id": sid,
            "cohort": str(r["cohort"]),
            "alt_status": r["arm"],
            "sample_type": str(r.get("sample_type", "unknown")),
            "file_backed_rna": bool(direct or route_b_ok),
            "recovery_route": (
                "direct+route_b" if direct and route_b_ok
                else "direct" if direct
                else "route_b" if route_b_ok
                else "none"
            ),
            "mapped_model_id": rb_id or (sid if direct else ""),
            "data_type": DEPMAP_RNA_DATA_TYPE if (direct or route_b_ok) else "unknown",
        })
    return pd.DataFrame(rows)


def assess_confounding_generic(
    inventory: pd.DataFrame,
    arm_col: str,
    arms: tuple[str, str],
) -> dict:
    """Shared confounding assessment for Route-B-only inventories."""
    fb = inventory[inventory["file_backed_rna"].astype(bool)].copy() if len(inventory) else inventory
    cohort_x_arm = (
        pd.crosstab(fb["cohort"], fb[arm_col]) if len(fb) else pd.DataFrame()
    )
    cohorts_per_arm = {a: [] for a in arms}
    if len(fb):
        for arm in arms:
            cohorts_per_arm[arm] = sorted(
                fb.loc[fb[arm_col] == arm, "cohort"].unique().tolist()
            )
    cohorts_with_both: list[str] = []
    if not cohort_x_arm.empty:
        for cohort in cohort_x_arm.index:
            present = {
                a for a in arms
                if a in cohort_x_arm.columns and cohort_x_arm.at[cohort, a] >= 1
            }
            if set(arms).issubset(present):
                cohorts_with_both.append(cohort)
    same_dtype: list[tuple[str, str]] = []
    for cohort in cohorts_with_both:
        sub = fb[fb["cohort"] == cohort]
        per_dtype = pd.crosstab(sub["data_type"], sub[arm_col])
        for dt in per_dtype.index:
            has_a = all(
                per_dtype.at[dt, a] >= 1 if a in per_dtype.columns else False
                for a in arms
            )
            if has_a:
                same_dtype.append((cohort, dt))
    return {
        "n_total": len(fb),
        "n_arm0": int((fb[arm_col] == arms[0]).sum()) if len(fb) else 0,
        "n_arm1": int((fb[arm_col] == arms[1]).sum()) if len(fb) else 0,
        "cohort_x_arm": cohort_x_arm,
        "cohorts_per_arm": cohorts_per_arm,
        "cohorts_with_both_arms": cohorts_with_both,
        "cohorts_with_both_arms_same_dtype": same_dtype,
        "has_within_cohort_co_occurrence": bool(cohorts_with_both),
        "has_within_cohort_co_occurrence_same_dtype": bool(same_dtype),
    }


def decide_route_b_aim1(
    inventory: pd.DataFrame,
    assessment: dict,
    *,
    per_arm_min: int = WITHIN_COHORT_PER_ARM_MIN,
) -> dict:
    """Verdict for Aim 1 using Route-B-only file-backed inventory."""
    # Reuse aim1_inventory column name expected by decide_aim1_design
    inv = inventory.copy()
    # Build assessment dict shape for decide_aim1_design
    full = {
        "n_total": assessment["n_total"],
        "n_disrupted": assessment["n_arm0"],
        "n_intact": assessment["n_arm1"],
        "cohort_x_arm": assessment["cohort_x_arm"],
        "cohorts_per_arm": {
            "disrupted": assessment["cohorts_per_arm"].get("disrupted", []),
            "intact": assessment["cohorts_per_arm"].get("intact", []),
        },
        "sample_types_per_arm": {"disrupted": [], "intact": []},
        "data_types_per_arm": {"disrupted": [], "intact": []},
        "cohorts_with_both_arms": assessment["cohorts_with_both_arms"],
        "cohorts_with_both_arms_same_dtype": assessment["cohorts_with_both_arms_same_dtype"],
        "has_within_cohort_co_occurrence": assessment["has_within_cohort_co_occurrence"],
        "has_within_cohort_co_occurrence_same_dtype": assessment[
            "has_within_cohort_co_occurrence_same_dtype"
        ],
    }
    verdict = decide_aim1_design(inv, full, per_arm_min=per_arm_min)
    verdict["scenario"] = "route_b_only"
    return verdict


def decide_route_b_aim2(
    inventory: pd.DataFrame,
    assessment: dict,
    *,
    per_arm_min: int = WITHIN_COHORT_PER_ARM_MIN,
) -> dict:
    """Verdict for Aim 2 using Route-B-only file-backed inventory."""
    full = {
        "n_total": assessment["n_total"],
        "n_positive": assessment["n_arm0"],
        "n_negative": assessment["n_arm1"],
        "cohort_x_arm": assessment["cohort_x_arm"],
        "cohorts_per_arm": {
            ALT_POSITIVE_LABEL: assessment["cohorts_per_arm"].get(ALT_POSITIVE_LABEL, []),
            ALT_NEGATIVE_LABEL: assessment["cohorts_per_arm"].get(ALT_NEGATIVE_LABEL, []),
        },
        "sample_types_per_arm": {ALT_POSITIVE_LABEL: [], ALT_NEGATIVE_LABEL: []},
        "data_types_per_arm": {ALT_POSITIVE_LABEL: [], ALT_NEGATIVE_LABEL: []},
        "cohorts_with_both_arms": assessment["cohorts_with_both_arms"],
        "cohorts_with_both_arms_same_dtype": assessment["cohorts_with_both_arms_same_dtype"],
        "has_within_cohort_co_occurrence": assessment["has_within_cohort_co_occurrence"],
        "has_within_cohort_co_occurrence_same_dtype": assessment[
            "has_within_cohort_co_occurrence_same_dtype"
        ],
    }
    verdict = decide_aim2_design(inventory, full, per_arm_min=per_arm_min)
    verdict["scenario"] = "route_b_only"
    return verdict
