"""
Audit logic: contingency tables, go/no-go thresholds, report generation.

Pass B (confounding census) must run before Pass A (availability census).
"""

from __future__ import annotations  # noqa: E402

import io
import textwrap
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


# ---------------------------------------------------------------------------
# Contingency tables
# ---------------------------------------------------------------------------

def contingency_table(
    df: pd.DataFrame,
    row_col: str,
    col_col: str,
    normalize: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame, float | None, float | None]:
    """
    Build a count table and optional row-% table with chi-squared p-value.

    Returns (counts, row_pct, chi2_stat, p_value).
    p_value is None if expected counts are too low for a valid chi-squared test.
    """
    counts = pd.crosstab(df[row_col], df[col_col], margins=False)
    row_pct = counts.div(counts.sum(axis=1), axis=0) * 100 if normalize else None

    chi2_stat, p_val = None, None
    if counts.shape[0] >= 2 and counts.shape[1] >= 2:
        observed = counts.values
        expected = np.outer(observed.sum(axis=1), observed.sum(axis=0)) / observed.sum()
        if (expected >= 5).all():
            chi2_stat, p_val, _, _ = stats.chi2_contingency(observed)
        else:
            chi2_stat, p_val = None, None  # expected counts too low

    return counts, row_pct, chi2_stat, p_val


def contingency_to_markdown(
    counts: pd.DataFrame,
    row_pct: pd.DataFrame | None,
    chi2: float | None,
    p_val: float | None,
    title: str,
) -> str:
    lines = [f"### {title}\n"]
    lines.append("**Counts:**\n")
    lines.append(counts.to_markdown())
    if row_pct is not None:
        lines.append("\n\n**Row %:**\n")
        lines.append(row_pct.round(1).to_markdown())
    if chi2 is not None:
        lines.append(f"\n\nchi² = {chi2:.3f}, p = {p_val:.4f}")
    else:
        lines.append("\n\n*(chi-squared not computed — expected counts < 5 or insufficient categories)*")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Pass B: Confounding census
# ---------------------------------------------------------------------------

PASSB_CONTRASTS = [
    ("cohort", "tcab1_genomic_status"),
    ("cohort", "alt_status"),
    ("tcab1_genomic_status", "tp53_status"),
    ("tcab1_genomic_status", "alt_status"),
    ("tcab1_genomic_status", "telomerase_status"),
    ("cohort", "sample_type"),
    ("sv_architecture_17p13_1", "alt_status"),
    ("source_confidence", "alt_status"),
]

COHORT_DOMINANCE_THRESHOLD = 0.70


def run_pass_b(manifest: pd.DataFrame) -> tuple[list[dict], list[str]]:
    """Run all Pass B contingency tables.

    Returns:
      tables  : list of dicts with keys domain, definition_set, title, counts,
                row_pct, chi2, p_val.  domain is 'aim1' for the existing Aim 1
                confounding census and 'aim2_direction' for the new direction-aware
                Aim 2 arm-balance tables.
      warnings: list of human-readable strings.  Blocking issues are prefixed
                "COHORT DOMINANCE:"; informational notes are prefixed "INFO:".
    """
    tables: list[dict] = []
    warnings: list[str] = []

    # ------------------------------------------------------------------
    # Aim 1 confounding census (existing contrasts, Aim 1 universe)
    # ------------------------------------------------------------------
    for def_set in ("strict", "pragmatic", "broad"):
        flag_col = f"include_aim1_primary_{def_set}"
        if flag_col not in manifest.columns:
            continue
        sub = manifest[_coerce_flag_to_bool(manifest[flag_col])].copy()
        if sub.empty:
            warnings.append(f"[{def_set}] No samples pass inclusion — skipping contingency tables")
            continue

        for row_col, col_col in PASSB_CONTRASTS:
            if row_col not in sub.columns or col_col not in sub.columns:
                continue
            counts, row_pct, chi2, p_val = contingency_table(sub, row_col, col_col)

            # Rename the cohort × alt_status table to make its provenance-only
            # role explicit now that direction-aware tables exist for Aim 2.
            if row_col == "cohort" and col_col == "alt_status":
                title = f"[{def_set.upper()}] {row_col} × {col_col} (provenance — Aim 1 universe)"
            else:
                title = f"[{def_set.upper()}] {row_col} × {col_col}"

            tables.append(dict(
                domain="aim1",
                definition_set=def_set,
                title=title,
                row_col=row_col,
                col_col=col_col,
                counts=counts,
                row_pct=row_pct,
                chi2=chi2,
                p_val=p_val,
            ))

            if row_col == "cohort" and col_col == "tcab1_genomic_status":
                _check_cohort_dominance(counts, title, warnings)
            elif row_col == "cohort" and col_col == "alt_status":
                # Label-source coverage — dominance hits are informational only.
                _check_cohort_dominance(
                    counts, title, warnings,
                    informational_classes={"alt_silver", "alt_unknown", "unknown"},
                )

    # ------------------------------------------------------------------
    # Aim 2 direction-aware arm-balance tables (new — one per def_set)
    # ------------------------------------------------------------------
    if "alt_call_direction" in manifest.columns:
        for def_set in ("strict", "pragmatic", "broad"):
            flag_col = f"include_aim2_{def_set}"
            if flag_col not in manifest.columns:
                continue
            sub2 = manifest[_coerce_flag_to_bool(manifest[flag_col])].copy()
            if sub2.empty or "cohort" not in sub2.columns:
                continue

            counts, row_pct, chi2, p_val = contingency_table(sub2, "cohort", "alt_call_direction")
            title = (
                f"[{def_set.upper()}] cohort × alt_call_direction "
                f"(Aim 2 arm balance — {flag_col})"
            )
            tables.append(dict(
                domain="aim2_direction",
                definition_set=def_set,
                title=title,
                row_col="cohort",
                col_col="alt_call_direction",
                counts=counts,
                row_pct=row_pct,
                chi2=chi2,
                p_val=p_val,
            ))
            # Arm-dominance warnings fire on the direction column — blocking.
            _check_cohort_dominance(counts, title, warnings)

    return tables, warnings


def _check_cohort_dominance(
    counts: pd.DataFrame,
    title: str,
    warnings: list[str],
    informational_classes: set[str] | None = None,
) -> None:
    """Append dominance warnings (or informational notes) to ``warnings``.

    Parameters
    ----------
    informational_classes
        Column (class) names for which a dominance hit is downgraded to an
        "INFO:" note rather than a "COHORT DOMINANCE:" blocking warning.
        Use this for label-source coverage columns (e.g. alt_silver, alt_unknown)
        that describe provenance rather than arm imbalance.
    """
    for col in counts.columns:
        col_total = counts[col].sum()
        if col_total == 0:
            continue
        for cohort, count in counts[col].items():
            frac = count / col_total
            if frac > COHORT_DOMINANCE_THRESHOLD:
                is_info = informational_classes and str(col) in informational_classes
                if is_info:
                    warnings.append(
                        f"INFO: {title} — cohort '{cohort}' contributes "
                        f"{frac:.0%} of class '{col}' (label-source coverage, not arm imbalance)."
                    )
                else:
                    warnings.append(
                        f"COHORT DOMINANCE: {title} — cohort '{cohort}' contributes "
                        f"{frac:.0%} of class '{col}' (threshold: {COHORT_DOMINANCE_THRESHOLD:.0%}). "
                        "Do not proceed to Pass A until stratification strategy is agreed."
                    )


def write_pass_b_report(
    tables: list[dict],
    warnings: list[str],
    output_path: Path | str,
) -> None:
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    blocking = [w for w in warnings if not w.startswith("INFO:")]
    informational = [w for w in warnings if w.startswith("INFO:")]

    lines = [
        "# Audit Pass B — Confounding Census",
        "",
        "_Generated by notebook 04_audit_passB_confounding.ipynb_",
        "",
    ]

    if blocking:
        lines += ["## ⚠️ Blocking Warnings", ""]
        for w in blocking:
            lines.append(f"- **{w}**")
        lines.append("")
    else:
        lines += ["## ⚠️ Blocking Warnings", "", "None.", ""]

    if informational:
        lines += ["## ℹ️ Informational Notes (label-source coverage)", ""]
        for w in informational:
            lines.append(f"- {w[len('INFO: '):]}")
        lines.append("")

    # ── Aim 2 direction-aware arm-balance tables ──────────────────────────────
    aim2_tables = [t for t in tables if t.get("domain") == "aim2_direction"]
    if aim2_tables:
        lines += [
            "## Aim 2 — Direction-Aware Arm Balance",
            "",
            "Counts and row-% of cohort × alt_call_direction within each Aim 2 "
            "inclusion universe. This is the modeling-relevant check: does any cohort "
            "dominate the positive or negative arm?",
            "",
        ]
        current_set = None
        for t in aim2_tables:
            if t["definition_set"] != current_set:
                current_set = t["definition_set"]
                lines.append(f"### Definition set: {current_set.upper()}")
                lines.append("")
            md = contingency_to_markdown(
                t["counts"], t["row_pct"], t["chi2"], t["p_val"], t["title"]
            )
            lines.append(md)
            lines.append("")

    # ── Aim 1 confounding census tables ───────────────────────────────────────
    aim1_tables = [t for t in tables if t.get("domain", "aim1") == "aim1"]
    lines += [
        "## Aim 1 — Confounding Census",
        "",
        "Standard confounding contrasts computed on the Aim 1 inclusion universe. "
        "The cohort × alt_status table is a provenance check only; arm balance "
        "for Aim 2 is reported in the direction-aware section above.",
        "",
    ]
    current_set = None
    for t in aim1_tables:
        if t["definition_set"] != current_set:
            current_set = t["definition_set"]
            lines.append(f"### Definition set: {current_set.upper()}")
            lines.append("")
        md = contingency_to_markdown(
            t["counts"], t["row_pct"], t["chi2"], t["p_val"], t["title"]
        )
        lines.append(md)
        lines.append("")

    output_path.write_text("\n".join(lines))
    print(f"Pass B report written to {output_path}")


# ---------------------------------------------------------------------------
# Pass A: Availability census
# ---------------------------------------------------------------------------

GO_NOGO = {
    "aim1_primary_min": 10,
    "aim1_primary_hard_stop": 8,
    "aim1_tp53_sens_min": 8,
    "aim2_min_unpenalized": 20,
    "aim2_min_penalized": 15,
}


# ---------------------------------------------------------------------------
# Aim 1 arm classification (single source of truth)
# ---------------------------------------------------------------------------

# Genomic states that map a row to the **disrupted** arm of the primary contrast.
# `tcab1_overall_status == "disrupted_high_conf"` is treated as equivalent so
# DepMap CN-derived calls (which set overall but not genomic) are preserved.
AIM1_DISRUPTED_GENOMIC = ("disrupted_high_conf", "candidate_regulatory_disruption")
AIM1_INTACT_GENOMIC = ("intact",)


def _coerce_flag_to_bool(s: pd.Series) -> pd.Series:
    """Coerce an inclusion-flag column to a real bool Series.

    Robust against the CSV roundtrip: ``pd.read_csv(..., dtype=str)`` returns
    the strings ``"True"`` / ``"False"``, and ``"False" == True`` is ``False``
    (correct) but ``bool("False")`` is ``True`` (because non-empty strings are
    truthy in Python) — both pitfalls have bitten this project before.
    """
    if s.dtype == bool:
        return s
    return s.astype(str).str.strip().str.lower().isin(("true", "1", "yes"))


def classify_aim1_arms(manifest: pd.DataFrame, def_set: str) -> pd.Series:
    """Assign each manifest row to one of {'disrupted', 'intact', 'uncertain', 'excluded'}.

    Single source of truth for the Aim 1 primary-contrast arm assignment per
    label_spec_v0.9 §1c. ``count_aim1()`` and notebook 05's accounting check
    both call this so they cannot drift.

    Robust against the CSV roundtrip: the include flag is coerced to bool
    internally, so the caller may pass a manifest with either real-bool or
    stringified-bool inclusion columns.

    Parameters
    ----------
    manifest : DataFrame
        Master manifest (after ``compute_inclusion_flags`` has been applied).
    def_set : {'strict', 'pragmatic', 'broad'}
        Definition set; selects the ``include_aim1_primary_{def_set}`` flag.

    Returns
    -------
    pd.Series indexed like ``manifest`` with one of:
      * ``'disrupted'`` — in universe AND tcab1_genomic_status in
        AIM1_DISRUPTED_GENOMIC, OR tcab1_overall_status == 'disrupted_high_conf'.
      * ``'intact'``    — in universe AND tcab1_genomic_status == 'intact'.
      * ``'uncertain'`` — in universe BUT neither arm (e.g. Broad's RNA-quantile
        proxy rows with tcab1_genomic_status == 'unknown').
      * ``'excluded'``  — not in the inclusion universe.
    """
    flag = f"include_aim1_primary_{def_set}"
    if flag not in manifest.columns:
        return pd.Series("excluded", index=manifest.index, name=f"aim1_arm_{def_set}")

    in_universe = _coerce_flag_to_bool(manifest[flag])
    genomic = manifest.get("tcab1_genomic_status", pd.Series("unknown", index=manifest.index))
    overall = manifest.get("tcab1_overall_status", pd.Series("unknown", index=manifest.index))

    disrupted_mask = in_universe & (
        genomic.isin(AIM1_DISRUPTED_GENOMIC)
        | (overall == "disrupted_high_conf")
    )
    intact_mask = in_universe & ~disrupted_mask & genomic.isin(AIM1_INTACT_GENOMIC)
    uncertain_mask = in_universe & ~disrupted_mask & ~intact_mask

    arms = pd.Series("excluded", index=manifest.index, name=f"aim1_arm_{def_set}")
    arms[disrupted_mask] = "disrupted"
    arms[intact_mask] = "intact"
    arms[uncertain_mask] = "uncertain"
    return arms


def count_aim1(manifest: pd.DataFrame, def_set: str) -> dict:
    flag = f"include_aim1_primary_{def_set}"
    if flag not in manifest.columns:
        return {}

    arms = classify_aim1_arms(manifest, def_set)
    n_dis = int((arms == "disrupted").sum())
    n_int = int((arms == "intact").sum())
    n_unc = int((arms == "uncertain").sum())
    n_universe = n_dis + n_int + n_unc

    tp53_flag = "include_tp53_sensitivity"
    telo_flag = "include_telomerase_negative_secondary"

    tp53_sub = (
        manifest[_coerce_flag_to_bool(manifest[tp53_flag])]
        if tp53_flag in manifest.columns else pd.DataFrame()
    )
    telo_sub = (
        manifest[_coerce_flag_to_bool(manifest[telo_flag])]
        if telo_flag in manifest.columns else pd.DataFrame()
    )

    return {
        "definition_set": def_set,
        # Inclusion universe: every row passing include_aim1_primary_{set}.
        "aim1_total": n_universe,
        # Primary-contrast arms (used for go/no-go: min(disrupted, intact)).
        "aim1_disrupted": n_dis,
        "aim1_intact": n_int,
        "aim1_uncertain": n_unc,
        "aim1_primary_contrast": n_dis + n_int,
        "aim1_tp53_total": len(tp53_sub),
        "aim1_telo_neg_total": len(telo_sub),
    }


def count_aim2(manifest: pd.DataFrame, def_set: str) -> dict:
    flag = f"include_aim2_{def_set}"
    if flag not in manifest.columns:
        return {}
    sub = manifest[_coerce_flag_to_bool(manifest[flag])]

    # Arms are defined by direction, not label source.
    # alt_silver lines split across arms depending on their ALT_pred_class;
    # counting by alt_status would lump all silver into the positive arm.
    if "alt_call_direction" in sub.columns:
        alt_pos = int((sub["alt_call_direction"] == "positive").sum())
        alt_neg = int((sub["alt_call_direction"] == "negative").sum())
        n_unknown_dir = int(
            (~sub["alt_call_direction"].isin(["positive", "negative"])).sum()
        )
        # I4 invariant: every included row should have a known direction.
        assert n_unknown_dir == 0, (
            f"count_aim2 [{def_set}]: {n_unknown_dir} included row(s) have "
            f"alt_call_direction not in {{positive, negative}}. "
            "Run assert_merge_invariants (I4) to locate offending rows."
        )
    else:
        # Fallback for manifests built before the alt_call_direction column was added.
        alt_pos = int(sub["alt_status"].isin(["alt_gold_positive", "alt_silver"]).sum())
        alt_neg = int((sub["alt_status"] == "alt_negative_high_conf").sum())

    return {
        "definition_set": def_set,
        "aim2_total": len(sub),
        "aim2_alt_positive": alt_pos,
        "aim2_alt_negative": alt_neg,
    }


def apply_go_nogo(counts: dict) -> dict[str, str]:
    """Apply Section 6 thresholds. Returns dict of aim → 'GO' / 'NO-GO' / 'EXPLORATORY'."""
    results = {}

    n_dis = counts.get("aim1_disrupted", 0)
    n_int = counts.get("aim1_intact", 0)
    min_group = min(n_dis, n_int)

    if min_group >= GO_NOGO["aim1_primary_min"]:
        results["aim1_primary"] = "GO"
    elif min_group >= GO_NOGO["aim1_primary_hard_stop"]:
        results["aim1_primary"] = "MARGINAL"
    else:
        results["aim1_primary"] = "NO-GO"

    n_tp53 = counts.get("aim1_tp53_total", 0)
    results["aim1_tp53_sensitivity"] = (
        "GO" if n_tp53 >= GO_NOGO["aim1_tp53_sens_min"] else "EXPLORATORY"
    )

    n_alt_pos = counts.get("aim2_alt_positive", 0)
    n_alt_neg = counts.get("aim2_alt_negative", 0)
    min_aim2 = min(n_alt_pos, n_alt_neg)

    if min_aim2 >= GO_NOGO["aim2_min_unpenalized"]:
        results["aim2"] = "GO (unpenalized)"
    elif min_aim2 >= GO_NOGO["aim2_min_penalized"]:
        results["aim2"] = "GO (penalized only)"
    else:
        results["aim2"] = "NO-GO (signature scoring only)"

    return results


def write_pass_a_report(
    all_counts: list[dict],
    all_go_nogo: list[dict],
    output_path: Path | str,
    recommended_set: str | None = None,
) -> None:
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    lines = [
        "# Audit Pass A — Availability Census",
        "",
        "_Generated by notebook 05_audit_passA_availability.ipynb_",
        "",
        "## Summary Table",
        "",
    ]

    header = (
        "| Set | Aim1 universe | disrupted | intact | uncertain | primary contrast | "
        "Aim1 TP53 | Aim2 total | ALT+ | ALT- | Aim1 verdict | Aim2 verdict |"
    )
    sep = "|---|---|---|---|---|---|---|---|---|---|---|---|"
    lines += [header, sep]

    for c, g in zip(all_counts, all_go_nogo):
        def_set = c.get("definition_set", "?")
        lines.append(
            f"| {def_set} "
            f"| {c.get('aim1_total', 0)} "
            f"| {c.get('aim1_disrupted', 0)} "
            f"| {c.get('aim1_intact', 0)} "
            f"| {c.get('aim1_uncertain', 0)} "
            f"| {c.get('aim1_primary_contrast', 0)} "
            f"| {c.get('aim1_tp53_total', 0)} "
            f"| {c.get('aim2_total', 0)} "
            f"| {c.get('aim2_alt_positive', 0)} "
            f"| {c.get('aim2_alt_negative', 0)} "
            f"| {g.get('aim1_primary', '?')} "
            f"| {g.get('aim2', '?')} |"
        )

    lines += [
        "",
        "_Aim1 universe_ = all rows with `include_aim1_primary_{set}` (may include RNA-proxy "
        "or other broad inclusion without a genomic arm assignment). "
        "_Primary contrast_ = disrupted + intact only; go/no-go uses "
        "`min(disrupted, intact)`, not universe size.",
        "",
    ]

    if recommended_set:
        lines += [
            f"## Recommended primary definition: **{recommended_set}**",
            "",
            "(Strictest definition that clears all go/no-go thresholds per label_spec_v0.9 §5)",
            "",
        ]

    lines += ["## Per-definition detail", ""]
    for c, g in zip(all_counts, all_go_nogo):
        def_set = c.get("definition_set", "?")
        lines.append(f"### {def_set.upper()}")
        for k, v in c.items():
            if k != "definition_set":
                lines.append(f"- {k}: {v}")
        for k, v in g.items():
            lines.append(f"- verdict_{k}: **{v}**")
        lines.append("")

    lines += [
        "## Next steps",
        "",
        "1. Review both audit reports.",
        "2. Confirm or adjust the primary definition set.",
        "3. Authorize freezing label_spec_v1.0.md.",
        "4. Only then begin DE / classifier / single-cell phases.",
        "",
        "_No DE or modeling code has been written._",
    ]

    output_path.write_text("\n".join(lines))
    print(f"Pass A report written to {output_path}")
