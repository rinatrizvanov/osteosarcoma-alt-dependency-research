"""
Merge helpers for notebooks/10_altitude_depmap_alt_merge.ipynb.

Writes ALTitude ALT calls into the master manifest for RNA-backed osteosarcoma
DepMap cell lines.  This is the first mutating module since the audit freeze.

Priority hierarchy (highest → lowest):
  1. Existing manifest gold  — alt_gold_positive or alt_negative_high_conf already
                               present in the row; never downgraded.
  2. ALTitude experimental   — altitude_evidence_tier = 'experimental' AND
                               altitude_pred_class ∈ {ALT+, ALT-}.
                               Writes alt_gold_positive / alt_negative_high_conf,
                               source_confidence = peer_reviewed.
  3. HOS special case        — ACH-000613 has NaN ALT_pred_class but ALT_assay type
                               = C-circle and ALT assay strength = Negative
                               (PMID: 33344901).  Treated as alt_negative_high_conf,
                               source_confidence = peer_reviewed.
  4. ALTitude prediction-only — altitude_evidence_tier = 'prediction_only' AND
                               altitude_pred_class ∈ {ALT+, ALT-}.
                               Writes alt_silver, source_confidence = computational.
                               Per label_spec_v1.0 §2 — applies to both + and −
                               predictions at this lower-confidence tier.
  5. No call                 — ambiguous_or_missing / unknown; leave manifest
                               unchanged.

Hard invariants (asserted before manifest is written):
  I1. Never downgrade: if current_priority > proposed_priority → skip.
  I2. Never stamp peer_reviewed onto a row that ends up alt_unknown.
  I3. Never overwrite a gold call (priority 2) with silver (priority 1).

Ground-truth assertions (four known literature standards):
  ACH-000364  U2OS   → positive arm (alt_gold_positive)
  ACH-000410  SAOS2  → positive arm (alt_gold_positive)
  ACH-000748  SJSA1  → negative arm (alt_negative_high_conf)
  ACH-000359  MG63   → negative arm (alt_negative_high_conf)
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from .aim2_inventory import ALT_NEGATIVE_LABEL, ALT_POSITIVE_LABEL
from .label_rules import compute_inclusion_flags, derive_alt_call_direction

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# Osteosarcoma lines excluded from merge by disease type.
NON_OS_SKIP: set[str] = {"ACH-000592"}   # TM31 (glioblastoma)

# Special case: no ALT_pred_class but C-circle negative (PMID: 33344901).
HOS_ID = "ACH-000613"
HOS_SOURCE = (
    "ALTitude_2024; PMID:33344901; ALT_assay=C-circle; "
    "ALT_assay_strength=Negative (low C-circle signal)"
)

ALTITUDE_SOURCE_EXPERIMENTAL = "ALTitude_2024"
ALTITUDE_SOURCE_PREDICTED = "ALTitude_2024_predicted"

# Identical to notebook 03 priority table — do not diverge.
ALT_PRIORITY: dict[str, int] = {
    "alt_unknown":            0,
    "unknown":                0,
    "alt_silver":             1,
    "alt_proxy":              1,
    "alt_negative_high_conf": 2,
    "alt_gold_positive":      2,
}

# Ground-truth assertions for four known OS lines.
_GROUND_TRUTH: dict[str, str] = {
    "ACH-000364": "positive",   # U2OS  — experimental ALT+
    "ACH-000410": "positive",   # SAOS2 — experimental ALT+
    "ACH-000748": "negative",   # SJSA1 — experimental ALT−
    "ACH-000359": "negative",   # MG63  — experimental ALT−
}


# ---------------------------------------------------------------------------
# Merge plan
# ---------------------------------------------------------------------------

def _proposed_call(
    model_id: str,
    inv_row: pd.Series | None,
) -> tuple[str, str, str, str]:
    """Return (proposed_alt_status, proposed_alt_source, proposed_source_confidence, proposed_direction).

    Returns ('', '', '', '') when there is no call to write (no_call / skip case).
    The caller is responsible for comparing against the current manifest value.
    """
    # HOS special case: experimental evidence, NaN pred_class → hard-code negative.
    if model_id == HOS_ID:
        return ALT_NEGATIVE_LABEL, HOS_SOURCE, "peer_reviewed", "negative"

    if inv_row is None:
        return "", "", "", ""

    tier = str(inv_row.get("altitude_evidence_tier", "unknown"))
    pc   = str(inv_row.get("altitude_pred_class", ""))

    if tier == "experimental" and pc in ("ALT+", "ALT-"):
        status = ALT_POSITIVE_LABEL if pc == "ALT+" else ALT_NEGATIVE_LABEL
        direction = derive_alt_call_direction(status, pc)
        return status, ALTITUDE_SOURCE_EXPERIMENTAL, "peer_reviewed", direction

    if tier == "prediction_only" and pc in ("ALT+", "ALT-"):
        # Both + and − prediction-only calls go to alt_silver per label_spec_v1.0 §2.
        direction = derive_alt_call_direction("alt_silver", pc)
        return "alt_silver", ALTITUDE_SOURCE_PREDICTED, "computational", direction

    return "", "", "", ""


def build_merge_plan(
    manifest: pd.DataFrame,
    inv: pd.DataFrame,
) -> list[dict]:
    """Return a list of per-row merge action dicts.

    Keys: model_id, action, current_alt, current_src, proposed_alt,
          proposed_src, proposed_conf, reason.

    ``action`` values:
      'write'           — apply the proposed call.
      'skip_gold'       — current call is gold from a non-ALTitude source; keep.
      'skip_downgrade'  — proposed priority < current priority; skip.
      'skip_same'       — proposed call is identical to current; idempotent.
      'skip_no_call'    — ALTitude has no usable call; leave unchanged.
      'skip_non_os'     — model is non-OS; merge does not apply.
      'skip_not_in_manifest' — model_id not found in any DepMap_26Q1 manifest row.
    """
    # Build inventory lookup by model_id (OS lines only).
    inv_lookup: dict[str, pd.Series] = {
        str(r["model_id"]): r
        for _, r in inv.iterrows()
        if str(r.get("os_status", "")) == "osteosarcoma"
    }

    # Add HOS to lookup set even if it resolved ambiguous_or_missing in inv.
    if HOS_ID not in inv_lookup:
        hos_rows = inv[inv["model_id"] == HOS_ID]
        if not hos_rows.empty:
            inv_lookup[HOS_ID] = hos_rows.iloc[0]
        else:
            inv_lookup[HOS_ID] = None   # type: ignore[assignment]

    plan: list[dict] = []
    all_model_ids = set(inv_lookup.keys()) | NON_OS_SKIP

    for model_id in sorted(all_model_ids):
        if model_id in NON_OS_SKIP:
            plan.append({
                "model_id": model_id, "action": "skip_non_os",
                "current_alt": "", "current_src": "",
                "proposed_alt": "", "proposed_src": "", "proposed_conf": "",
                "proposed_direction": "unknown",
                "reason": "non-OS line excluded from OS merge",
            })
            continue

        mrows = manifest[
            (manifest["sample_id"].astype(str) == model_id)
            & (manifest["cohort"] == "DepMap_26Q1")
        ]
        if mrows.empty:
            plan.append({
                "model_id": model_id, "action": "skip_not_in_manifest",
                "current_alt": "", "current_src": "",
                "proposed_alt": "", "proposed_src": "", "proposed_conf": "",
                "proposed_direction": "unknown",
                "reason": "no DepMap_26Q1 row in manifest",
            })
            continue

        current_alt  = str(mrows.iloc[0]["alt_status"])
        current_src  = str(mrows.iloc[0].get("alt_source", ""))
        current_pri  = ALT_PRIORITY.get(current_alt, 0)

        inv_row = inv_lookup.get(model_id)
        proposed_alt, proposed_src, proposed_conf, proposed_dir = _proposed_call(model_id, inv_row)

        if not proposed_alt:
            plan.append({
                "model_id": model_id, "action": "skip_no_call",
                "current_alt": current_alt, "current_src": current_src,
                "proposed_alt": "", "proposed_src": "", "proposed_conf": "",
                "proposed_direction": "unknown",
                "reason": "ALTitude has no usable call for this line",
            })
            continue

        proposed_pri = ALT_PRIORITY.get(proposed_alt, 0)

        if proposed_pri < current_pri:
            plan.append({
                "model_id": model_id, "action": "skip_downgrade",
                "current_alt": current_alt, "current_src": current_src,
                "proposed_alt": proposed_alt, "proposed_src": proposed_src,
                "proposed_conf": proposed_conf,
                "proposed_direction": proposed_dir,
                "reason": f"would downgrade {current_alt!r}→{proposed_alt!r}; skipped",
            })
            continue

        if proposed_alt == current_alt:
            # skip_same: direction is already resolved from inv_row's altitude_pred_class.
            same_direction = derive_alt_call_direction(
                current_alt,
                str(inv_row.get("altitude_pred_class", "")) if inv_row is not None else None,
            )
            plan.append({
                "model_id": model_id, "action": "skip_same",
                "current_alt": current_alt, "current_src": current_src,
                "proposed_alt": proposed_alt, "proposed_src": proposed_src,
                "proposed_conf": proposed_conf,
                "proposed_direction": same_direction,
                "reason": "alt_status already matches proposed call (idempotent re-run)",
            })
            continue

        plan.append({
            "model_id": model_id, "action": "write",
            "current_alt": current_alt, "current_src": current_src,
            "proposed_alt": proposed_alt, "proposed_src": proposed_src,
            "proposed_conf": proposed_conf,
            "proposed_direction": proposed_dir,
            "reason": f"{current_alt!r} → {proposed_alt!r}",
        })

    return plan


# ---------------------------------------------------------------------------
# Apply
# ---------------------------------------------------------------------------

def apply_merge_plan(
    manifest: pd.DataFrame,
    plan: list[dict],
) -> pd.DataFrame:
    """Return a modified copy of manifest with 'write' actions applied.

    Also backfills alt_call_direction for:
    - 'write' action rows (new direction from the merge).
    - 'skip_same' rows whose alt_status is already labeled (idempotent direction stamp).
    - All other rows with alt_gold_positive / alt_negative_high_conf (direction is
      deterministic from alt_status alone regardless of cohort).

    Does not mutate the input DataFrame.
    """
    manifest = manifest.copy()

    # Initialise column so every row has a value before we start writing.
    if "alt_call_direction" not in manifest.columns:
        manifest["alt_call_direction"] = "unknown"

    _DIRECTION_ELIGIBLE = {"alt_gold_positive", "alt_negative_high_conf", "alt_silver"}

    for action in plan:
        model_id = action["model_id"]
        mask = (
            (manifest["sample_id"].astype(str) == model_id)
            & (manifest["cohort"] == "DepMap_26Q1")
        )
        if not mask.any():
            continue

        if action["action"] == "write":
            manifest.loc[mask, "alt_status"]          = action["proposed_alt"]
            manifest.loc[mask, "alt_source"]          = action["proposed_src"]
            manifest.loc[mask, "source_confidence"]   = action["proposed_conf"]
            manifest.loc[mask, "alt_call_direction"]  = action["proposed_direction"]

        elif action["action"] == "skip_same":
            # Backfill direction for already-labeled rows.
            if (
                action["current_alt"] in _DIRECTION_ELIGIBLE
                and action["proposed_direction"] != "unknown"
            ):
                manifest.loc[mask, "alt_call_direction"] = action["proposed_direction"]

    # Vectorised backfill for all cohorts: gold status is direction-deterministic.
    no_dir_mask = manifest["alt_call_direction"] == "unknown"
    manifest.loc[no_dir_mask & (manifest["alt_status"] == "alt_gold_positive"),
                 "alt_call_direction"] = "positive"
    manifest.loc[no_dir_mask & (manifest["alt_status"] == "alt_negative_high_conf"),
                 "alt_call_direction"] = "negative"

    return manifest


# ---------------------------------------------------------------------------
# Invariant assertions
# ---------------------------------------------------------------------------

def assert_merge_invariants(
    manifest_pre: pd.DataFrame,
    manifest_post: pd.DataFrame,
    plan: list[dict],
) -> None:
    """Raise ValueError if any of the four hard invariants is violated.

    I1. No downgrade: every row whose alt_status changed has a new priority
        ≥ the old priority.
    I2. No peer_reviewed stamp on alt_unknown rows.
    I3. No gold→silver overwrite (priority 2 → priority 1).
    I4. Every row whose alt_status is labeled (not alt_unknown / alt_proxy)
        must have alt_call_direction != 'unknown'.
    """
    writes = [a for a in plan if a["action"] == "write"]

    for a in writes:
        mid = a["model_id"]
        old_pri = ALT_PRIORITY.get(a["current_alt"], 0)
        new_pri = ALT_PRIORITY.get(a["proposed_alt"], 0)

        # I1
        if new_pri < old_pri:
            raise ValueError(
                f"Invariant I1 violated for {mid}: "
                f"downgrade {a['current_alt']!r} (pri={old_pri}) → "
                f"{a['proposed_alt']!r} (pri={new_pri})."
            )
        # I3
        if old_pri == 2 and new_pri == 1:
            raise ValueError(
                f"Invariant I3 violated for {mid}: "
                f"gold call {a['current_alt']!r} overwritten by silver "
                f"{a['proposed_alt']!r}."
            )

    # I2: scan post-manifest for peer_reviewed + alt_unknown
    depmap_post = manifest_post[manifest_post["cohort"] == "DepMap_26Q1"]
    bad_i2 = depmap_post[
        (depmap_post["alt_status"].isin(["alt_unknown", "unknown"]))
        & (depmap_post["source_confidence"] == "peer_reviewed")
    ]
    if not bad_i2.empty:
        raise ValueError(
            f"Invariant I2 violated: {len(bad_i2)} DepMap row(s) have "
            f"source_confidence='peer_reviewed' but alt_status='alt_unknown': "
            f"{bad_i2['sample_id'].tolist()}"
        )

    # I4: every labeled row must have alt_call_direction != 'unknown'.
    _LABELED = {"alt_gold_positive", "alt_negative_high_conf", "alt_silver"}
    if "alt_call_direction" not in manifest_post.columns:
        raise ValueError(
            "Invariant I4 violated: alt_call_direction column is missing from manifest_post."
        )
    labeled_rows = manifest_post[manifest_post["alt_status"].isin(_LABELED)]
    bad_i4 = labeled_rows[
        labeled_rows["alt_call_direction"].isna()
        | (labeled_rows["alt_call_direction"].astype(str).isin(["unknown", "", "nan"]))
    ]
    if not bad_i4.empty:
        raise ValueError(
            f"Invariant I4 violated: {len(bad_i4)} row(s) have a labeled alt_status but "
            f"alt_call_direction='unknown': "
            f"{bad_i4[['sample_id', 'cohort', 'alt_status']].to_dict('records')}"
        )


def assert_ground_truth(manifest_post: pd.DataFrame) -> None:
    """Raise ValueError if the four known ground-truth lines map to the wrong arm."""
    for mid, expected_arm in _GROUND_TRUTH.items():
        rows = manifest_post[manifest_post["sample_id"].astype(str) == mid]
        if rows.empty:
            continue  # not in this manifest run
        actual = str(rows.iloc[0]["alt_status"])
        if expected_arm == "positive":
            ok = actual in (ALT_POSITIVE_LABEL, "alt_silver")
        else:
            ok = actual == ALT_NEGATIVE_LABEL
        if not ok:
            raise ValueError(
                f"Ground-truth check FAILED for {mid}: "
                f"expected arm='{expected_arm}' but alt_status='{actual}'. "
                "Abort — do not write the manifest."
            )


# ---------------------------------------------------------------------------
# Include-flag refresh
# ---------------------------------------------------------------------------

def recompute_inclusion_flags(manifest: pd.DataFrame) -> pd.DataFrame:
    """Return manifest copy with all include_* columns recomputed from scratch."""
    manifest = manifest.copy()
    flag_cols = [c for c in manifest.columns if c.startswith("include_")]
    flag_rows = manifest.apply(compute_inclusion_flags, axis=1)
    for col in flag_cols:
        manifest[col] = [r[col] for r in flag_rows]
    return manifest


# ---------------------------------------------------------------------------
# Report helpers
# ---------------------------------------------------------------------------

def _alt_counts(df: pd.DataFrame, cohort: str | None = None) -> pd.Series:
    sub = df[df["cohort"] == cohort] if cohort else df
    return sub["alt_status"].value_counts(dropna=False)


def build_change_summary(
    manifest_pre: pd.DataFrame,
    manifest_post: pd.DataFrame,
    plan: list[dict],
) -> str:
    """Return a markdown string summarising pre→post changes."""
    lines: list[str] = []

    writes = [a for a in plan if a["action"] == "write"]
    skips  = [a for a in plan if a["action"] != "write"]

    lines.append(f"## Merge actions\n")
    lines.append(f"- Rows written:  **{len(writes)}**\n")
    action_counts: dict[str, int] = {}
    for a in skips:
        action_counts[a["action"]] = action_counts.get(a["action"], 0) + 1
    for k, v in sorted(action_counts.items()):
        lines.append(f"- {k}: {v}\n")
    lines.append("\n")

    if writes:
        lines.append("### Written rows (model_id  old → new)\n")
        lines.append("| model_id | old_alt_status | new_alt_status | provenance |\n")
        lines.append("|---|---|---|---|\n")
        for a in writes:
            lines.append(
                f"| {a['model_id']} | {a['current_alt']} "
                f"| {a['proposed_alt']} | {a['proposed_src'][:40]} |\n"
            )
        lines.append("\n")

    lines.append("## alt_status distribution — overall\n")
    lines.append("| alt_status | pre | post |\n")
    lines.append("|---|---|---|\n")
    pre_counts  = _alt_counts(manifest_pre)
    post_counts = _alt_counts(manifest_post)
    all_statuses = sorted(set(pre_counts.index) | set(post_counts.index))
    for s in all_statuses:
        lines.append(f"| {s} | {pre_counts.get(s, 0)} | {post_counts.get(s, 0)} |\n")
    lines.append("\n")

    lines.append("## alt_status distribution — DepMap_26Q1 cohort\n")
    lines.append("| alt_status | pre | post |\n")
    lines.append("|---|---|---|\n")
    pre_d  = _alt_counts(manifest_pre,  "DepMap_26Q1")
    post_d = _alt_counts(manifest_post, "DepMap_26Q1")
    for s in sorted(set(pre_d.index) | set(post_d.index)):
        lines.append(f"| {s} | {pre_d.get(s, 0)} | {post_d.get(s, 0)} |\n")
    lines.append("\n")

    lines.append("## include_aim2 flag counts — before and after\n")
    lines.append("| flag | pre | post |\n")
    lines.append("|---|---|---|\n")
    for col in ("include_aim2_strict", "include_aim2_pragmatic", "include_aim2_broad"):
        if col in manifest_pre.columns and col in manifest_post.columns:
            pre_n  = (manifest_pre[col].astype(str).str.lower() == "true").sum()
            post_n = (manifest_post[col].astype(str).str.lower() == "true").sum()
            lines.append(f"| {col} | {pre_n} | {post_n} |\n")
    lines.append("\n")

    # Silver direction split (only meaningful in post-manifest).
    if "alt_call_direction" in manifest_post.columns:
        lines.append("## alt_silver direction split (DepMap_26Q1)\n")
        silver_rows = manifest_post[
            (manifest_post["cohort"] == "DepMap_26Q1")
            & (manifest_post["alt_status"] == "alt_silver")
        ]
        if silver_rows.empty:
            lines.append("_No alt_silver rows in DepMap_26Q1._\n\n")
        else:
            dir_counts = silver_rows["alt_call_direction"].value_counts(dropna=False)
            n_pos = dir_counts.get("positive", 0)
            n_neg = dir_counts.get("negative", 0)
            n_unk = dir_counts.get("unknown", 0)
            lines.append(f"- positive (ALT+): **{n_pos}**\n")
            lines.append(f"- negative (ALT-): **{n_neg}**\n")
            if n_unk:
                lines.append(f"- unknown (unresolved): **{n_unk}** ← INVESTIGATE\n")
            lines.append("\n")
            lines.append("| sample_id | alt_call_direction |\n")
            lines.append("|---|---|\n")
            for _, r in silver_rows.sort_values("sample_id").iterrows():
                lines.append(f"| {r['sample_id']} | {r['alt_call_direction']} |\n")
            lines.append("\n")

        lines.append("## alt_status × alt_call_direction (all cohorts)\n")
        lines.append("| alt_status | positive | negative | unknown |\n")
        lines.append("|---|---|---|---|\n")
        for status in sorted(manifest_post["alt_status"].unique()):
            sub = manifest_post[manifest_post["alt_status"] == status]
            dc = sub["alt_call_direction"].value_counts(dropna=False)
            lines.append(
                f"| {status} | {dc.get('positive', 0)} "
                f"| {dc.get('negative', 0)} "
                f"| {dc.get('unknown', 0) + dc.get(float('nan'), 0)} |\n"
            )
        lines.append("\n")

    return "".join(lines)
