"""
Deterministic label-assignment rules as specified in label_spec_v0.9.md.

All functions follow the rule: absence of evidence → 'unknown'.
Never assign 'intact' or 'not_altered' by default.
"""

from __future__ import annotations

import pandas as pd
import numpy as np


# ---------------------------------------------------------------------------
# TP53 / gene alteration status
# ---------------------------------------------------------------------------

def tp53_status_from_damaging_matrix(
    model_id: str,
    damaging_matrix: pd.DataFrame,
) -> tuple[str, str]:
    """
    Returns (status, evidence) for TP53 from OmicsSomaticMutationsMatrixDamaging.

    The damaging matrix has ModelID as index, gene columns like 'TP53 (7157)'.
    A value > 0 means at least one damaging mutation detected.
    """
    tp53_col = _find_gene_col(damaging_matrix, "TP53", 7157)
    cn_col = None  # CN handled separately

    if tp53_col is None:
        return "unknown", "TP53 column not found in damaging matrix"

    if model_id not in damaging_matrix.index:
        return "unknown", f"{model_id} not in damaging matrix"

    val = damaging_matrix.loc[model_id, tp53_col]
    if pd.isna(val):
        return "unknown", "NaN in damaging matrix"
    if float(val) > 0:
        return "altered", f"OmicsSomaticMutationsMatrixDamaging value={val:.3g}"
    return "not_altered", f"OmicsSomaticMutationsMatrixDamaging value={val:.3g}; no CN check yet"


def gene_status_from_damaging_matrix(
    model_id: str,
    gene_symbol: str,
    entrez_id: int,
    damaging_matrix: pd.DataFrame,
) -> tuple[str, str]:
    col = _find_gene_col(damaging_matrix, gene_symbol, entrez_id)
    if col is None:
        return "unknown", f"{gene_symbol} column not found"
    if model_id not in damaging_matrix.index:
        return "unknown", f"{model_id} not in damaging matrix"
    val = damaging_matrix.loc[model_id, col]
    if pd.isna(val):
        return "unknown", "NaN in damaging matrix"
    if float(val) > 0:
        return "altered", f"{col} damaging value={val:.3g}"
    return "not_altered", f"{col} damaging value={val:.3g}"


def tp53_status_include_deep_deletion(
    model_id: str,
    damaging_matrix: pd.DataFrame,
    cn_gene_matrix: pd.DataFrame,
) -> tuple[str, str]:
    """
    Checks damaging mutations AND deep deletion (CN ≤ -1) in TP53.
    Returns (status, evidence).
    """
    status, evidence = tp53_status_from_damaging_matrix(model_id, damaging_matrix)
    if status == "altered":
        return status, evidence

    tp53_cn_col = _find_gene_col(cn_gene_matrix, "TP53", 7157)
    if tp53_cn_col and model_id in cn_gene_matrix.index:
        cn_val = cn_gene_matrix.loc[model_id, tp53_cn_col]
        if not pd.isna(cn_val) and float(cn_val) <= -1.0:
            return "altered", f"TP53 deep deletion CN={cn_val:.3g}"

    return status, evidence


# ---------------------------------------------------------------------------
# TCAB1 / WRAP53 RNA
# ---------------------------------------------------------------------------

def tcab1_rna_quantile(
    model_id: str,
    expr_matrix: pd.DataFrame,
    cohort_model_ids: list[str] | None = None,
) -> tuple[float | None, float | None, float | None]:
    """
    Returns (global_quantile, rna_log2tpm, within_cohort_quantile).
    quantile is 0-1 rank / (n-1), within cohort if cohort_model_ids provided.
    """
    wrap53_col = _find_gene_col(expr_matrix, "WRAP53", 55135)
    if wrap53_col is None:
        return None, None, None
    if model_id not in expr_matrix.index:
        return None, None, None

    val = expr_matrix.loc[model_id, wrap53_col]
    if pd.isna(val):
        return None, None, None

    all_vals = expr_matrix[wrap53_col].dropna()
    global_q = float((all_vals < val).sum()) / max(len(all_vals) - 1, 1)

    within_q = None
    if cohort_model_ids:
        cohort_present = [m for m in cohort_model_ids if m in expr_matrix.index]
        if cohort_present:
            cohort_vals = expr_matrix.loc[cohort_present, wrap53_col].dropna()
            within_q = float((cohort_vals < val).sum()) / max(len(cohort_vals) - 1, 1)

    return global_q, float(val), within_q


# ---------------------------------------------------------------------------
# TCAB1 functional status
# ---------------------------------------------------------------------------

def tcab1_functional_status(
    rna_quantile_within_cohort: float | None,
    has_corroborating: bool = False,
    genomic_disrupted: bool = False,
) -> tuple[str, str]:
    """
    Apply label_spec_v0.9 Table 1b rules.
    rna_quantile_within_cohort: 0–1 (0.25 = bottom quartile boundary)
    has_corroborating: low protein / methylation evidence
    genomic_disrupted: genomic disruption flag
    """
    if rna_quantile_within_cohort is None:
        return "unknown", "RNA not available"

    low_rna = rna_quantile_within_cohort < 0.25
    corroborating = has_corroborating or genomic_disrupted

    if low_rna and corroborating:
        return "functional_loss_supported", (
            f"RNA quantile within cohort={rna_quantile_within_cohort:.3f} (bottom quartile)"
            + (" + genomic disruption" if genomic_disrupted else " + corroborating evidence")
        )
    if low_rna:
        return "functional_loss_possible", (
            f"RNA quantile within cohort={rna_quantile_within_cohort:.3f} (bottom quartile), "
            "no corroborating evidence"
        )
    if rna_quantile_within_cohort >= 0.25:
        return "functional_intact_supported", (
            f"RNA quantile within cohort={rna_quantile_within_cohort:.3f} (upper three quartiles)"
        )
    return "unknown", "Ambiguous RNA evidence"


# ---------------------------------------------------------------------------
# TCAB1 overall status + confidence
# ---------------------------------------------------------------------------

def tcab1_overall_and_confidence(
    genomic_status: str,
    functional_status: str,
) -> tuple[str, float]:
    """
    Derive tcab1_overall_status and confidence score (0–1).
    Confidence reflects strength of combined genomic + functional evidence.
    """
    genomic_weight = {
        "disrupted_high_conf": 1.0,
        "candidate_regulatory_disruption": 0.5,
        "intact": 0.0,
        "unknown": None,
    }
    functional_weight = {
        "functional_loss_supported": 1.0,
        "functional_loss_possible": 0.5,
        "functional_intact_supported": 0.0,
        "unknown": None,
    }

    gw = genomic_weight.get(genomic_status)
    fw = functional_weight.get(functional_status)

    if gw is None and fw is None:
        return "unknown", 0.0

    if gw == 1.0 and fw not in (0.0,):
        status = "disrupted_high_conf"
    elif gw == 0.0 and fw not in (1.0, 0.5):
        status = "intact"
    elif gw == 1.0 and fw == 0.0:
        status = "conflicted"
    elif gw == 0.5 or fw == 0.5:
        status = "candidate_disrupted"
    else:
        status = "unknown"

    weights = [w for w in (gw, fw) if w is not None]
    confidence = float(np.mean(weights)) if weights else 0.0

    return status, confidence


# ---------------------------------------------------------------------------
# Telomerase status
# ---------------------------------------------------------------------------

def telomerase_status_from_tert(
    model_id: str,
    expr_matrix: pd.DataFrame,
    cohort_model_ids: list[str] | None = None,
    pos_threshold_sd: float = 0.5,
    neg_threshold_sd: float = 0.5,
) -> tuple[str, float | None, str]:
    """
    Returns (status, tert_log2tpm, evidence).
    Provisional rule: positive if TERT > cohort_median + 0.5 SD,
    negative if < cohort_median - 0.5 SD, else indeterminate.
    """
    tert_col = _find_gene_col(expr_matrix, "TERT", 7015)
    if tert_col is None:
        return "unknown", None, "TERT column not found"
    if model_id not in expr_matrix.index:
        return "unknown", None, f"{model_id} not in expression matrix"

    val = expr_matrix.loc[model_id, tert_col]
    if pd.isna(val):
        return "unknown", None, "TERT expression NaN"

    ref_ids = cohort_model_ids if cohort_model_ids else list(expr_matrix.index)
    ref_present = [m for m in ref_ids if m in expr_matrix.index]
    ref_vals = expr_matrix.loc[ref_present, tert_col].dropna()

    median = float(ref_vals.median())
    sd = float(ref_vals.std())
    val_f = float(val)

    if val_f > median + pos_threshold_sd * sd:
        return "positive", val_f, f"TERT log2(TPM+1)={val_f:.3f} > median+{pos_threshold_sd}SD"
    if val_f < median - neg_threshold_sd * sd:
        return "negative", val_f, f"TERT log2(TPM+1)={val_f:.3f} < median-{neg_threshold_sd}SD"
    return "indeterminate", val_f, f"TERT log2(TPM+1)={val_f:.3f} within ±{pos_threshold_sd}SD"


# ---------------------------------------------------------------------------
# ALT tier mapping
# ---------------------------------------------------------------------------

ALTITUDE_ALT_MAP = {
    "ALT+": "alt_silver",
    "ALT-": "alt_negative_high_conf",
    "ALT+ (high confidence)": "alt_gold_positive",
}


def map_altitude_call(curated_marker: str | None) -> str:
    if curated_marker is None or (isinstance(curated_marker, float) and np.isnan(curated_marker)):
        return "alt_unknown"
    return ALTITUDE_ALT_MAP.get(str(curated_marker).strip(), "alt_unknown")


def derive_alt_call_direction(
    alt_status: str,
    altitude_pred_class: str | None = None,
) -> str:
    """Return 'positive' | 'negative' | 'unknown' for an alt_status + optional pred_class.

    Gold calls are direction-deterministic from alt_status alone.
    Silver calls require altitude_pred_class ('ALT+' or 'ALT-') to resolve direction.
    """
    if alt_status == "alt_gold_positive":
        return "positive"
    if alt_status == "alt_negative_high_conf":
        return "negative"
    if alt_status == "alt_silver":
        pc = str(altitude_pred_class or "").strip()
        if pc == "ALT+":
            return "positive"
        if pc == "ALT-":
            return "negative"
    return "unknown"


# ---------------------------------------------------------------------------
# Inclusion flags (definition sets)
# ---------------------------------------------------------------------------

def _coerce_bool(v) -> bool:
    """
    Robust string-aware coercion to bool.

    The manifest is serialised to / loaded from CSV, so boolean columns arrive
    as the strings ``"True"`` / ``"False"`` (or ``"true"``/``"false"``/``"1"``/``"0"``).
    ``bool("False")`` returns ``True`` in Python because non-empty strings are
    truthy — using bare ``bool(...)`` on these values silently turns every
    "False" into "True" and effectively disables any flag gated on it. Use
    this helper instead.
    """
    if isinstance(v, bool):
        return v
    if v is None:
        return False
    try:
        if pd.isna(v):
            return False
    except (TypeError, ValueError):
        pass
    s = str(v).strip().lower()
    return s in ("true", "1", "yes", "t")


# Architecture calls that are treated as informative readouts of the locus.
# Mirrors the mapping in manifest_builder.refresh_tcab1_status_from_architecture.
_INFORMATIVE_ARCH = {
    "Both-single-SV", "Both-independent-SVs", "TCAB1-only", "Complex",
    "Neither", "TP53-only",
}


def compute_inclusion_flags(row: pd.Series) -> dict[str, bool]:
    """
    Compute all include_* columns for one manifest row.
    Returns a dict of bool values.

    Aim 1 universe:
      A sample qualifies if we have a downstream readout of TCAB1 function — either
      DepMap RNA-seq (``has_rna=True``) or a literature IHC-derived
      ``tcab1_functional_status`` (Lawlor PDX-style), or an SV-architecture-derived
      informative call (PCGP "Neither" / "TP53-only" → intact;
      "Both-*" / "TCAB1-only" / "Complex" → disrupted). The earlier ``has_rna``-only
      gate excluded every literature sample without DepMap expression data, even
      when they carried peer-reviewed orthogonal evidence.

    Aim 2 universe:
      A sample qualifies if it carries an informative ALT label — *positive
      or negative*. We deliberately drop the ``has_rna`` gate here: the
      classifier may be trained on cross-cohort features (FISH, IHC,
      WGS-derived) and the audit's job is to count what is available, not to
      pre-impose a single training feature set. The previous logic only
      counted positives, which silently dropped the 11 high-confidence
      ALT-negative samples from the entire feasibility analysis — i.e. the
      classifier was being told only one class exists.
    """
    has_rna = _coerce_bool(row.get("has_rna"))
    genomic = str(row.get("tcab1_genomic_status", "unknown"))
    functional = str(row.get("tcab1_functional_status", "unknown"))
    arch = str(row.get("sv_architecture_17p13_1", "unknown"))
    alt = str(row.get("alt_status", "alt_unknown"))
    tp53 = str(row.get("tp53_status", "unknown"))
    telomerase = str(row.get("telomerase_status", "unknown"))
    rna_q_raw = row.get("tcab1_rna_quantile_within_cohort")
    # rna quantile is stored as a float when known and the string "unknown"
    # (or NaN) when missing. Coerce to float-or-None for downstream checks.
    try:
        rna_q = float(rna_q_raw)
        if pd.isna(rna_q):
            rna_q = None
    except (TypeError, ValueError):
        rna_q = None

    # A sample has *some* TCAB1 readout if any of: DepMap RNA, IHC-derived
    # functional status, or SV-architecture-derived genomic status.
    has_functional_readout = (
        has_rna
        or functional in (
            "functional_loss_supported",
            "functional_loss_possible",
            "functional_intact_supported",
        )
        or arch in _INFORMATIVE_ARCH
    )

    # ------- Aim 1 -------
    # strict: disrupted_high_conf OR intact (no unknown in either)
    strict_disrupted = genomic == "disrupted_high_conf" and functional != "functional_intact_supported"
    strict_intact = genomic == "intact" and functional != "functional_loss_supported"
    aim1_strict = has_functional_readout and (strict_disrupted or strict_intact)

    # pragmatic: + candidate_regulatory_disruption with low RNA
    pragmatic_disrupted = strict_disrupted or (
        genomic == "candidate_regulatory_disruption"
        and functional in ("functional_loss_supported", "functional_loss_possible")
    )
    pragmatic_intact = strict_intact
    aim1_pragmatic = has_functional_readout and (pragmatic_disrupted or pragmatic_intact)

    # broad: any SV in locus OR bottom-quartile RNA OR upper-quartile RNA OR
    # architecture-derived disrupted / intact (handled via genomic_status, which
    # refresh_tcab1_status_from_architecture writes from sv_architecture_17p13_1).
    broad_disrupted = (
        genomic in ("disrupted_high_conf", "candidate_regulatory_disruption")
        or (rna_q is not None and rna_q < 0.25)
    )
    broad_intact = (
        genomic == "intact"
        or (rna_q is not None and rna_q >= 0.75)
    )
    aim1_broad = has_functional_readout and (broad_disrupted or broad_intact)

    # ------- Aim 2 -------
    # Universe = any sample with a confident ALT label, positive OR negative.
    # No has_rna gate (see docstring).
    aim2_strict = alt in ("alt_gold_positive", "alt_negative_high_conf")
    aim2_pragmatic = alt in (
        "alt_gold_positive", "alt_silver", "alt_negative_high_conf",
    )
    aim2_broad = alt in (
        "alt_gold_positive", "alt_silver", "alt_proxy", "alt_negative_high_conf",
    )

    # ------- Sensitivities -------
    tp53_sens = aim1_strict and tp53 == "altered"
    telo_neg = aim1_strict and telomerase in ("negative",)

    return {
        "include_aim1_primary_strict": aim1_strict,
        "include_aim1_primary_pragmatic": aim1_pragmatic,
        "include_aim1_primary_broad": aim1_broad,
        "include_aim2_strict": aim2_strict,
        "include_aim2_pragmatic": aim2_pragmatic,
        "include_aim2_broad": aim2_broad,
        "include_tp53_sensitivity": tp53_sens,
        "include_telomerase_negative_secondary": telo_neg,
    }


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _find_gene_col(df: pd.DataFrame, symbol: str, entrez_id: int) -> str | None:
    """Find the column matching 'SYMBOL (ENTREZ)' format in DepMap matrices."""
    target = f"{symbol} ({entrez_id})"
    if target in df.columns:
        return target
    for col in df.columns:
        if col.startswith(f"{symbol} (") or col == symbol:
            return col
    return None
