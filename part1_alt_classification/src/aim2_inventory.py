"""
Read-only helpers for ``notebooks/07_aim2_rna_inventory.ipynb``.

This module mirrors :mod:`src.aim1_inventory` for the Aim 2 (ALT classifier)
contrast. Anything that does not depend on the choice of arm labels —
the expression-source registry, the cheap file-backed-membership probes,
and the per-sample :func:`lookup_file_backed_rna` router — is reused
directly. The Aim 2-specific code lives here and does not modify the
Aim 1 module.

Arm definitions for Aim 2
-------------------------
The Aim 2 primary contrast is the ALT classifier task:

* ``alt_gold_positive``    — ALT-positive by gold-standard assay
* ``alt_negative_high_conf`` — ALT-negative with high-confidence evidence

These are the literal ``alt_status`` values produced by the merge in
``notebooks/03_alt_label_merge.ipynb``. All other ``alt_status`` values
(``alt_unknown``, ``alt_silver``, ``alt_proxy``, ...) are excluded from
the primary contrast.
"""
from __future__ import annotations

from typing import Iterable

import pandas as pd

from .aim1_inventory import (  # noqa: F401 — re-export for notebook convenience
    DEPMAP_RNA_DATA_TYPE,
    DEPMAP_RNA_PATH,
    EXPRESSION_SOURCES,
    TARGET_OS_TPM_DATA_TYPE,
    TARGET_OS_TPM_PATH,
    WITHIN_COHORT_PER_ARM_MIN,
    depmap_rna_model_ids,
    lookup_file_backed_rna,
    target_os_rna_sample_ids,
)

ALT_POSITIVE_LABEL = "alt_gold_positive"
ALT_NEGATIVE_LABEL = "alt_negative_high_conf"
AIM2_ARMS: tuple[str, str] = (ALT_POSITIVE_LABEL, ALT_NEGATIVE_LABEL)


# ---------------------------------------------------------------------------
# Contrast selection
# ---------------------------------------------------------------------------

def select_aim2_contrast_rows(manifest: pd.DataFrame) -> pd.DataFrame:
    """Return manifest rows whose ``alt_status`` is one of the Aim 2 arms.

    Adds an ``arm`` column carrying the literal ``alt_status`` value so the
    downstream helpers share the same shape as the Aim 1 versions.
    """
    sub = manifest[manifest["alt_status"].isin(AIM2_ARMS)].copy()
    sub["arm"] = sub["alt_status"].values
    return sub.reset_index(drop=True)


# ---------------------------------------------------------------------------
# Per-sample inventory
# ---------------------------------------------------------------------------

def build_inventory_aim2(
    contrast_manifest: pd.DataFrame,
    *,
    depmap_models: set[str] | None = None,
    target_os_samples: set[str] | None = None,
) -> pd.DataFrame:
    """Tidy per-sample table for the Aim 2 contrast.

    The output column carrying the arm is named ``alt_status`` (not
    ``tcab1_overall_status`` like the Aim 1 version) so the notebook's
    cross-tabulations read naturally.
    """
    if depmap_models is None:
        depmap_models = depmap_rna_model_ids()
    if target_os_samples is None:
        target_os_samples = target_os_rna_sample_ids()

    rows = []
    for _, r in contrast_manifest.iterrows():
        sid = str(r["sample_id"])
        cohort = str(r["cohort"])
        sample_type = str(r.get("sample_type", "unknown"))
        manifest_has_rna = str(r.get("has_rna", "")).strip().lower() in ("true", "1", "yes")
        file_backed, source, data_type = lookup_file_backed_rna(
            sid, cohort,
            depmap_models=depmap_models,
            target_os_samples=target_os_samples,
        )
        rows.append({
            "sample_id":         sid,
            "cohort":            cohort,
            "alt_status":        str(r["arm"]),
            "sample_type":       sample_type,
            "manifest_has_rna":  bool(manifest_has_rna),
            "file_backed_rna":   bool(file_backed),
            "expression_source": source,
            "data_type":         data_type,
        })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Confounding assessment
# ---------------------------------------------------------------------------

def assess_confounding_aim2(inventory: pd.DataFrame) -> dict:
    """Structured confounding assessment for the Aim 2 file-backed subset.

    Identical shape to :func:`src.aim1_inventory.assess_confounding`, but with
    ALT arm labels instead of disrupted/intact. Returns:

    * ``n_total``, ``n_positive``, ``n_negative`` (file-backed)
    * ``cohort_x_arm``: crosstab cohort × alt_status
    * ``cohorts_per_arm``, ``sample_types_per_arm``, ``data_types_per_arm``
    * ``cohorts_with_both_arms``: list of cohort names with both arms
    * ``cohorts_with_both_arms_same_dtype``: list of (cohort, data_type)
      pairs where the cohort carries both arms under the same data type
    * ``has_within_cohort_co_occurrence`` and
      ``has_within_cohort_co_occurrence_same_dtype`` (bools)
    """
    fb = (
        inventory[inventory["file_backed_rna"].astype(bool)].copy()
        if len(inventory) else inventory
    )

    cohort_x_arm = (
        pd.crosstab(fb["cohort"], fb["alt_status"])
        if len(fb) else pd.DataFrame()
    )

    cohorts_per_arm:      dict[str, list[str]] = {a: [] for a in AIM2_ARMS}
    sample_types_per_arm: dict[str, list[str]] = {a: [] for a in AIM2_ARMS}
    data_types_per_arm:   dict[str, list[str]] = {a: [] for a in AIM2_ARMS}
    if len(fb):
        for arm in AIM2_ARMS:
            sub = fb[fb["alt_status"] == arm]
            cohorts_per_arm[arm]      = sorted(sub["cohort"].unique().tolist())
            sample_types_per_arm[arm] = sorted(sub["sample_type"].unique().tolist())
            data_types_per_arm[arm]   = sorted(sub["data_type"].unique().tolist())

    cohorts_with_both_arms: list[str] = []
    if not cohort_x_arm.empty:
        for cohort in cohort_x_arm.index:
            present = {
                a for a in AIM2_ARMS
                if a in cohort_x_arm.columns and cohort_x_arm.at[cohort, a] >= 1
            }
            if set(AIM2_ARMS).issubset(present):
                cohorts_with_both_arms.append(cohort)

    cohorts_with_both_arms_same_dtype: list[tuple[str, str]] = []
    for cohort in cohorts_with_both_arms:
        sub = fb[fb["cohort"] == cohort]
        per_dtype = pd.crosstab(sub["data_type"], sub["alt_status"])
        for dt in per_dtype.index:
            has_pos = (
                per_dtype.at[dt, ALT_POSITIVE_LABEL] >= 1
                if ALT_POSITIVE_LABEL in per_dtype.columns else False
            )
            has_neg = (
                per_dtype.at[dt, ALT_NEGATIVE_LABEL] >= 1
                if ALT_NEGATIVE_LABEL in per_dtype.columns else False
            )
            if has_pos and has_neg:
                cohorts_with_both_arms_same_dtype.append((cohort, dt))

    n_positive = int((fb["alt_status"] == ALT_POSITIVE_LABEL).sum()) if len(fb) else 0
    n_negative = int((fb["alt_status"] == ALT_NEGATIVE_LABEL).sum()) if len(fb) else 0

    return {
        "n_total":     len(fb),
        "n_positive":  n_positive,
        "n_negative":  n_negative,
        "cohort_x_arm":         cohort_x_arm,
        "cohorts_per_arm":      cohorts_per_arm,
        "sample_types_per_arm": sample_types_per_arm,
        "data_types_per_arm":   data_types_per_arm,
        "cohorts_with_both_arms": cohorts_with_both_arms,
        "cohorts_with_both_arms_same_dtype": cohorts_with_both_arms_same_dtype,
        "has_within_cohort_co_occurrence": bool(cohorts_with_both_arms),
        "has_within_cohort_co_occurrence_same_dtype": bool(cohorts_with_both_arms_same_dtype),
    }


# ---------------------------------------------------------------------------
# Decision rule
# ---------------------------------------------------------------------------

def _within_cohort_qualifies_aim2(
    inventory: pd.DataFrame,
    cohort: str,
    data_type: str,
    *,
    per_arm_min: int = WITHIN_COHORT_PER_ARM_MIN,
) -> tuple[bool, int, int]:
    sub = inventory[
        inventory["file_backed_rna"]
        & (inventory["cohort"] == cohort)
        & (inventory["data_type"] == data_type)
    ]
    n_pos = int((sub["alt_status"] == ALT_POSITIVE_LABEL).sum())
    n_neg = int((sub["alt_status"] == ALT_NEGATIVE_LABEL).sum())
    return (n_pos >= per_arm_min and n_neg >= per_arm_min), n_pos, n_neg


def decide_aim2_design(
    inventory: pd.DataFrame,
    assessment: dict,
    *,
    per_arm_min: int = WITHIN_COHORT_PER_ARM_MIN,
) -> dict:
    """Apply the same three-way decision rule as Aim 1, restated for Aim 2.

    * ``within_cohort_viable``    — ≥1 (cohort × data_type) pair has both
      arms with ≥ ``per_arm_min`` samples each.
    * ``pooled_exploratory_only`` — no single cohort has both arms with
      file-backed RNA, but **each** arm has file-backed samples from ≥2
      distinct cohorts.
    * ``not_estimable``           — anything else (an arm has zero
      file-backed RNA, or the arms come from disjoint single cohorts).
    """
    qualifying: list[tuple[str, str, int, int]] = []
    for cohort, data_type in assessment["cohorts_with_both_arms_same_dtype"]:
        ok, n_pos, n_neg = _within_cohort_qualifies_aim2(
            inventory, cohort, data_type, per_arm_min=per_arm_min,
        )
        if ok:
            qualifying.append((cohort, data_type, n_pos, n_neg))

    if qualifying:
        ranked = sorted(qualifying, key=lambda t: -min(t[2], t[3]))
        primary = ranked[0]
        rest = ranked[1:]
        rest_note = (
            f" {len(rest)} additional cohort(s) qualify for meta-analysis: "
            f"{', '.join(c for c, _, _, _ in rest)}."
            if rest else ""
        )
        return {
            "category": "within_cohort_viable",
            "primary_cohort": primary[0],
            "primary_data_type": primary[1],
            "primary_n_positive": primary[2],
            "primary_n_negative": primary[3],
            "qualifying_pairs": qualifying,
            "justification": (
                f"Cohort '{primary[0]}' provides {primary[2]} ALT-positive and "
                f"{primary[3]} ALT-negative file-backed RNA samples under the same "
                f"data type ('{primary[1]}'), meeting the within-cohort minimum of "
                f"{per_arm_min} per arm. Recommend training the Aim 2 ALT classifier "
                f"within this cohort first to avoid platform/cohort confounding."
                f"{rest_note}"
            ),
        }

    n_pos_cohorts = len(assessment["cohorts_per_arm"][ALT_POSITIVE_LABEL])
    n_neg_cohorts = len(assessment["cohorts_per_arm"][ALT_NEGATIVE_LABEL])

    if n_pos_cohorts >= 2 and n_neg_cohorts >= 2:
        return {
            "category": "pooled_exploratory_only",
            "n_positive_cohorts": n_pos_cohorts,
            "n_negative_cohorts": n_neg_cohorts,
            "justification": (
                f"No single cohort carries both arms with file-backed RNA, but "
                f"the ALT-positive arm spans {n_pos_cohorts} cohorts "
                f"({', '.join(assessment['cohorts_per_arm'][ALT_POSITIVE_LABEL])}) and "
                f"the ALT-negative arm spans {n_neg_cohorts} cohorts "
                f"({', '.join(assessment['cohorts_per_arm'][ALT_NEGATIVE_LABEL])}). "
                f"ALT status is therefore only partially confounded with cohort, so "
                f"a pooled classifier with cohort as a fixed effect / batch covariate "
                f"is feasible as an exploratory analysis. Within-cohort training is "
                f"not possible under the current data."
            ),
        }

    fb_pos = assessment["n_positive"]
    fb_neg = assessment["n_negative"]
    cohorts_pos = assessment["cohorts_per_arm"][ALT_POSITIVE_LABEL] or ["(none)"]
    cohorts_neg = assessment["cohorts_per_arm"][ALT_NEGATIVE_LABEL] or ["(none)"]

    return {
        "category": "not_estimable",
        "fb_n_positive": fb_pos,
        "fb_n_negative": fb_neg,
        "fb_cohorts_positive": cohorts_pos,
        "fb_cohorts_negative": cohorts_neg,
        "justification": (
            f"File-backed RNA is available for {fb_pos} ALT-positive sample(s) "
            f"(from cohort(s) {', '.join(cohorts_pos)}) and {fb_neg} ALT-negative "
            f"sample(s) (from cohort(s) {', '.join(cohorts_neg)}). ALT status is "
            f"perfectly or near-perfectly confounded with cohort / platform / "
            f"sample_type, so the Aim 2 contrast as currently defined is not "
            f"statistically estimable. The contrast must be redefined (e.g. require "
            f"RNA-backed cohorts on both sides) or restricted (e.g. drop arms "
            f"without file-backed RNA) before any ALT classifier can be trained."
        ),
    }
