"""
Read-only helpers for notebooks/09_altitude_depmap_alt_inventory.ipynb.

Joins ALTitude ALT calls (ALT_pred_class column, Sheet1 of ALTitude-table1.xlsx)
to RNA-backed DepMap ModelIDs from the Route B machinery, applies gold-over-silver
priority, and reports class balance for the Aim 2 contrast.

No manifest mutation. No full expression matrix load.

Column binding
--------------
``Curated.Marker`` in ALTitude-table1 is an ATRX-mutation / chrX SV coordinate
string — NOT the ALT classification. The correct field is ``ALT_pred_class``
("ALT+" / "ALT-" / NaN). ``ALT_pred_prob`` carries the model confidence.

Evidence tier
-------------
A line is ``experimental`` when any of ``ALT source``, ``ALT_assay type``, or
``ALT assay strength`` is non-null / non-blank. Otherwise ``prediction_only``.

Priority merge
--------------
1. Manifest gold (Taglialatela2026 / Lawlor)  →  highest
2. ALTitude experimental (non-null pred_class) →  alt_gold_positive / alt_negative_high_conf
3. ALTitude prediction_only, ALT+ only         →  alt_silver
4. ALTitude prediction_only, ALT-              →  ambiguous_or_missing
   (no lower-tier negative label in system; block negative arm requires high-conf)
5. ALT_pred_class == NaN, any tier             →  ambiguous_or_missing

Sanity guard
------------
``build_altitude_depmap_inventory`` raises ``ValueError`` before returning if
any of four known ground-truth lines maps to the wrong arm.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from .aim1_inventory import WITHIN_COHORT_PER_ARM_MIN, depmap_rna_model_ids  # noqa: F401
from .aim2_inventory import ALT_NEGATIVE_LABEL, ALT_POSITIVE_LABEL

_DATA_ROOT = Path(__file__).parent.parent / "data"
_ALTITUDE_PATH = _DATA_ROOT / "raw" / "literature_tables" / "ALTitude-table1.xlsx"

OS_DISEASE_SUBSTR = "Osteosarcoma"

# Lines in the Taglialatela_2026 cohort that are not osteosarcoma.
# Verified against Model.csv OncotreePrimaryDisease at runtime; listed here
# as a human-readable cross-check for the notebook.
NON_OS_KNOWN = {
    "ACH-000309": "lung_carcinoma (SKLU1)",   # not in DepMap RNA matrix
    "ACH-000592": "glioblastoma (TM31)",
}

# Priority scores for manifest ALT status tiers (higher = more authoritative).
_ALT_PRIORITY: dict[str, int] = {
    "alt_gold_positive": 5,
    "alt_negative_high_conf": 4,
    "alt_silver": 3,
    "alt_proxy": 2,
    "alt_unknown": 1,
}

# Known ground-truth lines for post-build sanity assertions.
# Values: "positive" (any of alt_gold_positive / alt_silver) or "negative".
_SANITY_GROUND_TRUTH: dict[str, str] = {
    "ACH-000364": "positive",   # U2OS  — ALT+ experimental
    "ACH-000410": "positive",   # SAOS2 — ALT+ experimental (also manifest gold)
    "ACH-000748": "negative",   # SJSA1 — ALT- experimental (also manifest gold)
    "ACH-000359": "negative",   # MG63  — ALT- experimental
}


# ---------------------------------------------------------------------------
# Internal: load ALTitude ALT_pred_class from the raw Excel
# ---------------------------------------------------------------------------

def _load_altitude_calls(path: Path | None = None) -> pd.DataFrame:
    """Read ALT_pred_class, ALT_pred_prob, and evidence columns from ALTitude-table1.

    Returns a DataFrame with one row per cell line:
      depmap_id             ACH-XXXXXX (from Achilles column)
      cell_line_name        StripCI
      alt_pred_class        "ALT+" | "ALT-" | None
      alt_pred_prob         float | NaN
      altitude_evidence_tier  "experimental" | "prediction_only"
      altitude_alt_status   mapped status for priority merge (see module docstring)
    """
    p = path or _ALTITUDE_PATH
    raw = pd.read_excel(p, sheet_name="Sheet1", header=0)

    # Locate columns defensively — column names may have trailing spaces.
    def _find(raw_cols: list[str], *needles: str) -> str | None:
        for needle in needles:
            for c in raw_cols:
                if needle.lower() in str(c).lower():
                    return c
        return None

    cols = list(raw.columns)
    achilles_col   = cols[0]                             # always first
    strip_col      = _find(cols, "StripCI", "Strip")
    pred_class_col = _find(cols, "ALT_pred_class")
    pred_prob_col  = _find(cols, "ALT_pred_prob")
    src_col        = _find(cols, "ALT source")
    assay_col      = _find(cols, "ALT_assay type", "assay type")
    strength_col   = _find(cols, "ALT assay strength", "assay strength")

    if pred_class_col is None:
        raise ValueError(
            "ALTitude-table1.xlsx Sheet1: cannot locate 'ALT_pred_class' column. "
            f"Available columns: {cols}"
        )

    def _is_nonnull(val: object) -> bool:
        s = str(val).strip()
        return s not in ("", "nan", "None", "NaN")

    rows = []
    for _, r in raw.iterrows():
        depmap_id = str(r[achilles_col]).strip()
        cell_line = str(r[strip_col]).strip() if strip_col else ""
        pred_class = r[pred_class_col] if pred_class_col else None
        pred_prob  = r[pred_prob_col]  if pred_prob_col  else float("nan")

        # Evidence tier
        experimental = any(
            _is_nonnull(r[c]) for c in (src_col, assay_col, strength_col) if c
        )
        tier = "experimental" if experimental else "prediction_only"

        # Map pred_class + tier → altitude_alt_status
        pc = str(pred_class).strip() if pred_class is not None and str(pred_class) not in ("nan", "", "None") else None

        if pc is None:
            altitude_status = "ambiguous_or_missing"
        elif pc == "ALT+" and tier == "experimental":
            altitude_status = ALT_POSITIVE_LABEL          # alt_gold_positive
        elif pc == "ALT-" and tier == "experimental":
            altitude_status = ALT_NEGATIVE_LABEL          # alt_negative_high_conf
        elif pc == "ALT+" and tier == "prediction_only":
            altitude_status = "alt_silver"
        else:
            # prediction_only ALT-: no lower-tier negative label; exclude from blocks
            altitude_status = "ambiguous_or_missing"

        rows.append({
            "depmap_id":             depmap_id,
            "cell_line_name":        cell_line,
            "alt_pred_class":        pc,
            "alt_pred_prob":         float(pred_prob) if str(pred_prob) not in ("nan", "", "None") else float("nan"),
            "altitude_evidence_tier": tier,
            "altitude_alt_status":   altitude_status,
        })

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Public helpers
# ---------------------------------------------------------------------------

def classify_os_from_model_df(model_id: str, model_df: pd.DataFrame) -> str:
    """Return 'osteosarcoma' | 'non_os' | 'unknown'."""
    hit = model_df[model_df["ModelID"] == model_id]
    if hit.empty:
        return "unknown"
    disease = str(hit.iloc[0].get("OncotreePrimaryDisease", ""))
    return "osteosarcoma" if OS_DISEASE_SUBSTR.lower() in disease.lower() else "non_os"


def _best_manifest_call(model_id: str, manifest: pd.DataFrame) -> tuple[str, str]:
    """Highest-priority ALT call for a model_id across all manifest cohorts.

    A model_id may appear in both DepMap_26Q1 (alt_unknown) and
    Taglialatela_2026 (gold call). Returns the best available status and source.
    """
    rows = manifest[manifest["sample_id"].astype(str) == model_id]
    if rows.empty:
        return "alt_unknown", ""
    best_status = "alt_unknown"
    best_source = ""
    best_score = 0
    for _, row in rows.iterrows():
        status = str(row.get("alt_status", "alt_unknown"))
        score = _ALT_PRIORITY.get(status, 0)
        if score > best_score:
            best_score = score
            best_status = status
            best_source = str(row.get("alt_source", ""))
    return best_status, best_source


def _is_gold_source(alt_source: str) -> bool:
    """True when alt_source traces to a peer-reviewed gold call."""
    s = str(alt_source)
    return "Taglialatela2026" in s or "Lawlor" in s


def build_altitude_depmap_inventory(
    manifest: pd.DataFrame,
    altitude_df: pd.DataFrame,
    model_df: pd.DataFrame,
    rna_model_ids: set[str],
    *,
    altitude_path: Path | None = None,
) -> pd.DataFrame:
    """Per-ModelID ALT inventory for RNA-backed DepMap lines.

    Parameters
    ----------
    manifest      Master manifest CSV loaded as DataFrame.
    altitude_df   Output of ``load_altitude_cellines()``; used only for
                  depmap_id membership.  ALT calls are read fresh from the
                  raw Excel via ``_load_altitude_calls()`` so that
                  ``Curated.Marker`` is never mistaken for the ALT call.
    model_df      Model.csv from DepMap (for OS classification and cell-line name).
    rna_model_ids Set of DepMap ModelIDs with RNA backing (from Route B).
    altitude_path Optional path override for ALTitude-table1.xlsx (tests / CI).

    Output columns
    --------------
    model_id, cell_line_name, os_status, non_os_known_flag,
    manifest_alt_status, manifest_alt_source, is_manifest_gold,
    altitude_alt_status, altitude_evidence_tier, altitude_pred_class,
    altitude_alt_prob, final_alt_status, final_provenance, disagreement_note

    Raises
    ------
    ValueError  If any of the four known ground-truth sanity checks fail.
    """
    # Load ALT_pred_class-based calls from raw Excel
    alt_calls = _load_altitude_calls(altitude_path)
    alt_lookup: dict[str, pd.Series] = {
        str(r["depmap_id"]).strip(): r
        for _, r in alt_calls.iterrows()
        if str(r["depmap_id"]).strip().startswith("ACH-")
    }

    rows = []
    for model_id in sorted(rna_model_ids):
        a_row = alt_lookup.get(model_id)

        os_status = classify_os_from_model_df(model_id, model_df)
        non_os_flag = NON_OS_KNOWN.get(model_id, "")

        # Best manifest call (gold from Taglialatela/Lawlor beats alt_unknown)
        manifest_alt, manifest_source = _best_manifest_call(model_id, manifest)
        is_gold = (
            manifest_alt in (ALT_POSITIVE_LABEL, ALT_NEGATIVE_LABEL)
            and _is_gold_source(manifest_source)
        )

        # ALTitude call (from ALT_pred_class, not Curated.Marker)
        altitude_status  = str(a_row["altitude_alt_status"])  if a_row is not None else "ambiguous_or_missing"
        altitude_tier    = str(a_row["altitude_evidence_tier"]) if a_row is not None else "unknown"
        altitude_pc      = str(a_row["alt_pred_class"])        if a_row is not None and a_row["alt_pred_class"] is not None else ""
        altitude_prob    = float(a_row["alt_pred_prob"])       if a_row is not None and str(a_row["alt_pred_prob"]) not in ("nan", "", "None") else float("nan")

        # Cell line display name from Model.csv
        hit = model_df[model_df["ModelID"] == model_id]
        cell_line = ""
        if not hit.empty:
            cell_line = str(
                hit.iloc[0].get("StrippedCellLineName")
                or hit.iloc[0].get("CellLineName")
                or ""
            )

        # Priority merge: manifest gold > ALTitude experimental > ALTitude prediction_only
        if is_gold:
            final_alt = manifest_alt
            final_prov = "gold_manifest"
        elif altitude_status not in ("ambiguous_or_missing",):
            final_alt = altitude_status
            final_prov = f"altitude_{altitude_tier}"
        else:
            final_alt = "ambiguous_or_missing"
            final_prov = "no_call"

        # Disagreement: manifest gold and ALTitude both have a definitive call but disagree
        disagreement = ""
        if (
            is_gold
            and altitude_status not in ("ambiguous_or_missing",)
            and final_alt != altitude_status
        ):
            disagreement = f"gold={final_alt} vs altitude={altitude_status}"

        rows.append({
            "model_id":               model_id,
            "cell_line_name":         cell_line,
            "os_status":              os_status,
            "non_os_known_flag":      non_os_flag,
            "manifest_alt_status":    manifest_alt,
            "manifest_alt_source":    manifest_source,
            "is_manifest_gold":       is_gold,
            "altitude_alt_status":    altitude_status,
            "altitude_evidence_tier": altitude_tier,
            "altitude_pred_class":    altitude_pc,
            "altitude_alt_prob":      altitude_prob,
            "final_alt_status":       final_alt,
            "final_provenance":       final_prov,
            "disagreement_note":      disagreement,
        })

    inv = pd.DataFrame(rows)

    # ------------------------------------------------------------------
    # Sanity assertions — raise before returning so no partial report
    # is written on a column-binding bug.
    # ------------------------------------------------------------------
    for chk_id, expected_arm in _SANITY_GROUND_TRUTH.items():
        if chk_id not in inv["model_id"].values:
            continue  # line not in rna_model_ids for this run — skip
        actual = str(inv.loc[inv["model_id"] == chk_id, "final_alt_status"].iloc[0])
        if expected_arm == "positive":
            ok = actual in (ALT_POSITIVE_LABEL, "alt_silver")
        else:
            ok = actual == ALT_NEGATIVE_LABEL
        if not ok:
            raise ValueError(
                f"Sanity check FAILED for {chk_id}: expected arm='{expected_arm}' "
                f"but final_alt_status='{actual}'. "
                "Check ALTitude column binding — Curated.Marker is NOT the ALT call."
            )

    return inv


def count_aim2_blocks(inv: pd.DataFrame) -> dict[str, dict]:
    """Arm counts for the three reporting blocks.

    Block 1  OS-only, gold labels only (alt_gold_positive vs alt_negative_high_conf)
    Block 2  OS-only sensitivity (add alt_silver to positive arm)
    Block 3  All lines (OS + non-OS), gold labels only
    """
    os_inv = inv[inv["os_status"] == "osteosarcoma"]

    def _counts(sub: pd.DataFrame, positive_statuses: set[str]) -> dict:
        n_pos   = int(sub["final_alt_status"].isin(positive_statuses).sum())
        n_neg   = int((sub["final_alt_status"] == ALT_NEGATIVE_LABEL).sum())
        n_amb   = int((sub["final_alt_status"] == "ambiguous_or_missing").sum())
        n_silver = int((sub["final_alt_status"] == "alt_silver").sum())
        return {
            "n_positive":             n_pos,
            "n_negative":             n_neg,
            "n_ambiguous_or_missing": n_amb,
            "n_silver_only_positive": n_silver,
            "n_total":                len(sub),
            "min_arm":                min(n_pos, n_neg),
        }

    return {
        "block1_os_gold":        _counts(os_inv, {ALT_POSITIVE_LABEL}),
        "block2_os_sensitivity": _counts(os_inv, {ALT_POSITIVE_LABEL, "alt_silver"}),
        "block3_all_gold":       _counts(inv,    {ALT_POSITIVE_LABEL}),
    }


def apply_altitude_depmap_verdict(
    blocks: dict[str, dict],
    *,
    per_arm_min: int = WITHIN_COHORT_PER_ARM_MIN,
) -> dict:
    """Three-way verdict for the ALTitude-DepMap Aim 2 class balance.

    within_cohort_viable    min(pos, neg) >= per_arm_min using OS-only gold labels
    pooled_exploratory_only only met with silver in positive arm
    not_estimable           neither threshold met
    """
    b1 = blocks["block1_os_gold"]
    b2 = blocks["block2_os_sensitivity"]

    if b1["min_arm"] >= per_arm_min:
        return {
            "category": "within_cohort_viable",
            "block": "block1_os_gold",
            "n_positive": b1["n_positive"],
            "n_negative": b1["n_negative"],
            "justification": (
                f"OS-only gold labels: {b1['n_positive']} ALT-positive and "
                f"{b1['n_negative']} ALT-negative RNA-backed DepMap lines. "
                f"Both arms ≥ {per_arm_min}; all from the same cohort/platform "
                "(DepMap 26Q1 log2TPM). Sufficient for within-cohort Aim 2 analysis."
            ),
        }

    if b2["min_arm"] >= per_arm_min:
        n_silver_added = b2["n_positive"] - b1["n_positive"]
        return {
            "category": "pooled_exploratory_only",
            "block": "block2_os_sensitivity",
            "n_positive": b2["n_positive"],
            "n_negative": b2["n_negative"],
            "n_silver_added": n_silver_added,
            "justification": (
                f"Gold-only is below threshold "
                f"(pos={b1['n_positive']}, neg={b1['n_negative']}, "
                f"min={b1['min_arm']} < {per_arm_min}). "
                f"Including ALTitude silver adds {n_silver_added} positive samples; "
                f"sensitivity gives pos={b2['n_positive']}, neg={b2['n_negative']}. "
                "Silver labels carry higher label noise — flag as exploratory."
            ),
        }

    return {
        "category": "not_estimable",
        "n_positive_gold":   b1["n_positive"],
        "n_negative_gold":   b1["n_negative"],
        "n_positive_silver": b2["n_positive"],
        "n_negative_silver": b2["n_negative"],
        "justification": (
            f"Neither gold (pos={b1['n_positive']}, neg={b1['n_negative']}) nor "
            f"gold+silver (pos={b2['n_positive']}, neg={b2['n_negative']}) "
            f"meets per-arm minimum of {per_arm_min}. "
            "ALTitude-DepMap route alone is not viable for Aim 2 without additional "
            "label sources (Route A or C)."
        ),
    }
