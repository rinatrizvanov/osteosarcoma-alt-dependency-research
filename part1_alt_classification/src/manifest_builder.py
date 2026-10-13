"""
Build the master manifest from DepMap + literature sources.

One row per sample. All unknown columns set to 'unknown'; never implies
absence-of-evidence as evidence-of-absence.
"""

from __future__ import annotations  # noqa: E402 — must be first

import pandas as pd
import numpy as np
from pathlib import Path

from .utils import MANIFEST_COLUMNS, GENE_ENTREZ
from .label_rules import (
    tp53_status_include_deep_deletion,
    gene_status_from_damaging_matrix,
    tcab1_rna_quantile,
    tcab1_functional_status,
    tcab1_overall_and_confidence,
    telomerase_status_from_tert,
    compute_inclusion_flags,
    _find_gene_col,
)

DATA_ROOT = Path(__file__).parent.parent / "data"
DEPMAP_DIR = DATA_ROOT / "raw" / "depmap"
LIT_DIR = DATA_ROOT / "raw" / "literature_tables"


# ---------------------------------------------------------------------------
# DepMap loaders
# ---------------------------------------------------------------------------

def load_depmap_model(path: Path | None = None) -> pd.DataFrame:
    p = path or DEPMAP_DIR / "Model.csv"
    df = pd.read_csv(p, low_memory=False)
    return df


def _load_depmap_matrix(p: Path) -> pd.DataFrame:
    """
    Load a DepMap omics matrix CSV with mixed metadata + gene columns.

    DepMap matrices contain leading metadata columns (SequencingID,
    ModelConditionID, ModelID, IsDefaultEntryForMC, IsDefaultEntryForModel)
    followed by per-gene values. A given ModelID can appear in multiple rows
    (one per sequencing run / condition); only the row with
    IsDefaultEntryForModel == 'Yes' is the canonical entry. Without filtering,
    df.loc[model_id, gene_col] would return a Series and break downstream
    truth-value checks ("truth value of a Series is ambiguous").

    This loader:
      1. Reads the file fully.
      2. Filters to IsDefaultEntryForModel == 'Yes' if present.
      3. Drops metadata columns and sets ModelID as the unique index.
    """
    df = pd.read_csv(p, low_memory=False)
    if "IsDefaultEntryForModel" in df.columns:
        df = df[df["IsDefaultEntryForModel"].astype(str).str.upper() == "YES"].copy()
    if "ModelID" not in df.columns:
        raise ValueError(f"{p.name}: no 'ModelID' column; got {list(df.columns[:6])}")
    df = df.drop_duplicates(subset="ModelID", keep="first").set_index("ModelID")
    meta_cols = [
        c for c in (
            "Unnamed: 0",
            "SequencingID",
            "ModelConditionID",
            "IsDefaultEntryForMC",
            "IsDefaultEntryForModel",
        )
        if c in df.columns
    ]
    if meta_cols:
        df = df.drop(columns=meta_cols)
    return df


def load_expression_matrix(path: Path | None = None) -> pd.DataFrame:
    p = path or DEPMAP_DIR / "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"
    return _load_depmap_matrix(p)


def load_damaging_matrix(path: Path | None = None) -> pd.DataFrame:
    p = path or DEPMAP_DIR / "OmicsSomaticMutationsMatrixDamaging.csv"
    return _load_depmap_matrix(p)


def load_cn_gene_matrix(path: Path | None = None) -> pd.DataFrame:
    p = path or DEPMAP_DIR / "OmicsCNGeneWGS.csv"
    return _load_depmap_matrix(p)


def load_cn_segments(path: Path | None = None) -> pd.DataFrame:
    p = path or DEPMAP_DIR / "OmicsCNSegmentsWGS.csv"
    df = pd.read_csv(p, low_memory=False)
    # Keep default-entry rows to avoid duplicates
    if "IsDefaultEntryForModel" in df.columns:
        df = df[df["IsDefaultEntryForModel"].astype(str).str.upper() == "YES"].copy()
    return df


def load_fusions(path: Path | None = None) -> pd.DataFrame:
    p = path or DEPMAP_DIR / "OmicsFusionFiltered.csv"
    df = pd.read_csv(p, low_memory=False)
    if "IsDefaultEntryForModel" in df.columns:
        df = df[df["IsDefaultEntryForModel"].astype(str).str.upper() == "YES"].copy()
    return df


# ---------------------------------------------------------------------------
# Literature loaders
# ---------------------------------------------------------------------------

def load_de_nonneville(path: Path | None = None) -> pd.DataFrame:
    p = path or LIT_DIR / "de_Nonneville_2022_table_EV3.csv"
    raw = pd.read_csv(p, skiprows=2, header=0)
    raw.columns = [c.strip() for c in raw.columns]
    raw = raw.dropna(subset=["N°"])

    df = pd.DataFrame()
    df["sample_id"] = "deNonneville_" + raw["N°"].astype(str).str.strip()
    df["cohort"] = "de_Nonneville_2022"
    df["sample_type"] = "tumor"
    df["alt_status"] = raw["ALT"].str.strip().str.lower().map(
        {"yes": "alt_gold_positive", "no": "alt_negative_high_conf"}
    ).fillna("alt_unknown")
    df["alt_source"] = "de Nonneville et al. EMBO Mol Med 2022"
    df["source_confidence"] = "peer_reviewed"
    df["alt_assay_evidence"] = "Literature-reported: ALT call from Table EV3"

    # TP53
    def _tp53(row):
        parts = []
        if str(row.get("P53 mut", "")).strip() not in ("no", "nd", ""):
            parts.append(f"mut:{row['P53 mut']}")
        if str(row.get("P53 break", "")).strip() == "yes":
            parts.append("break")
        if str(row.get("P53 del", "")).strip() == "yes":
            parts.append("del")
        return ("altered", "; ".join(parts)) if parts else ("unknown", "nd")

    tp53_cols = [c for c in raw.columns if "P53" in c or "p53" in c.lower()]
    if len(tp53_cols) >= 3:
        tp53_res = raw.apply(_tp53, axis=1)
        df["tp53_status"] = [r[0] for r in tp53_res]
        df["tp53_evidence"] = [r[1] for r in tp53_res]
    else:
        df["tp53_status"] = "unknown"
        df["tp53_evidence"] = "TP53 columns not parsed"

    # ATRX
    atrx_col = next((c for c in raw.columns if "ATRX" in c.upper()), None)
    if atrx_col:
        df["atrx_status"] = raw[atrx_col].apply(
            lambda x: "altered" if str(x).strip() not in ("no", "nd", "", "nan") else "unknown"
        )
    else:
        df["atrx_status"] = "unknown"

    return _fill_unknown_cols(df)


def load_lawlor_pdx_labels(path: Path | None = None) -> pd.DataFrame:
    p = path or LIT_DIR / "tcab1_pdx_labels.csv"
    raw = pd.read_csv(p)
    raw.columns = [c.strip() for c in raw.columns]

    df = pd.DataFrame()
    df["sample_id"] = raw["case_id"].str.strip()
    df["cohort"] = "Lawlor_2026_GCCRI"
    df["sample_type"] = "pdx"

    df["alt_status"] = raw["alt_binary"].map(
        {1: "alt_gold_positive", 0: "alt_negative_high_conf"}
    ).fillna("alt_unknown")
    df["alt_source"] = "Lawlor et al. bioRxiv 2026.01.02.691864"
    df["source_confidence"] = "preprint"
    df["alt_assay_evidence"] = raw.apply(
        lambda r: f"TRAP={r.get('trap','?')} ALT={r.get('alt','?')} (TEL-FISH + TRAP)", axis=1
    )

    df["tp53_status"] = raw["tp53_binary"].map(
        {1: "altered", 0: "not_altered"}
    ).fillna("unknown")
    df["tp53_evidence"] = "Lawlor 2026 IHC TP53"

    df["tcab1_functional_status"] = raw["tcab1_binary"].map(
        {0: "functional_loss_possible", 1: "functional_intact_supported"}
    ).fillna("unknown")
    df["tcab1_functional_evidence"] = "Lawlor 2026 IHC TCAB1"

    return _fill_unknown_cols(df)


def load_ccle_svs(path: Path | None = None) -> pd.DataFrame:
    """
    Load the CCLE 2019 SvABA structural-variant table (Ghandi et al. 2019).

    Returns the raw breakpoint-pair frame with columns standardised:
    ``CCLE_name, map_id, bp1, bp2, class, gene1, gene2, site1, site2,
    fusion, multi_sv_fusion, cosmic_fus``. Coordinates in ``bp1``/``bp2``
    are GRCh37 / hg19 (per the file's own ``descriptions`` sheet); they are
    lifted to GRCh38 inside ``sv_annotator.annotate_ccle_svs``.

    Substitutes for the missing ``OmicsStructuralVariants.csv`` in DepMap 26Q1.
    """
    p = path or DATA_ROOT / "raw" / "ccle" / "CCLE_translocations_SvABA_20181221.xlsx"
    df = pd.read_excel(p, sheet_name="CCLE_translocations_SvABA_20181")
    return df


def load_chen_2014_pcgp(
    s4_path: Path | None = None,
    s6_path: Path | None = None,
    s3a_path: Path | None = None,
    s3b_path: Path | None = None,
) -> pd.DataFrame:
    """
    Build a one-row-per-sample manifest fragment for the Chen 2014 PCGP
    osteosarcoma cohort (Cell Reports 2014).

    Sources:
      * Table S3a (discovery cohort)        — sample list, n=22
      * Table S3b (validation cohort)       — sample list, n=16
      * Table S4b (TP53 SV breakpoints)     — TP53 status per sample
      * Table S6  (ATRX IHC + telomere FISH)— atrx_status, alt_status (n=19)

    All upstream coordinates are GRCh37; this loader does NOT lift them. The
    SV architecture call is computed downstream by
    ``prefill_pcgp_sv_architecture`` and ``sv_annotator.classify_chen_sv_architecture``,
    in-frame, with only the categorical result stored in the manifest.

    Maps:
      ATRX status (Table S6 col 'ATRX status') → atrx_status
        - 'wt' (case-insensitive) → not_altered
        - missing / 'Indeterminate' / 'nan' / '' → unknown
        - anything else (Y266*, D1383fs, DEL (CNV), SV ...) → altered

      Telomere FISH (Table S6 col 'Telomere FISH') → alt_status
        - 'Positive' (any variant) → alt_gold_positive
        - 'Negative' → alt_negative_high_conf
        - 'Failed' / 'Not performed' / missing → alt_unknown

    TP53 (presence in Table S4b) → tp53_status='altered',
    evidence quotes count and translocation_type list.
    """
    chen_dir = DATA_ROOT / "raw" / "literature_tables"
    s4_path  = s4_path  or chen_dir / "chen_2014_tableS4_TP53_SVs.xlsx"
    s6_path  = s6_path  or chen_dir / "chen_2014_tableS6_ATRX_ALT.xlsx"
    s3a_path = s3a_path or chen_dir / "chen_2014_tableS3_SVs.xlsx"
    s3b_path = s3b_path or chen_dir / "chen_2014_tableS3_SVs.xlsx"

    # ---- Sample roster: union of discovery (S3a) + validation (S3b) -----
    sample_ids: list[str] = []
    for sheet, src in [("Discovery cohort", s3a_path), ("Validation cohort", s3b_path)]:
        try:
            raw = pd.read_excel(src, sheet_name=sheet, header=None)
            samples = raw.iloc[2:, 0].dropna().astype(str).str.strip()
            samples = [s for s in samples if s and s.lower() != "sample"]
            sample_ids.extend(samples)
        except Exception:
            continue
    sample_ids = sorted(set(sample_ids))

    # ---- TP53 from Table S4b (each row = one TP53 SV breakpoint pair) ----
    # NB: in the source xlsx, row 0 is the title, row 1 is blank, row 2 is the
    # column header line, row 3 is also blank, real data starts at row 4. So
    # ``header=3`` is the correct row index (verified 2026-05-01 by direct
    # inspection of the file).
    s4_lookup: dict[str, dict] = {}
    try:
        s4 = pd.read_excel(s4_path, sheet_name="TP53 SV breakpoints", header=3)
        s4 = s4.rename(columns=lambda c: str(c).strip())
        if "sample" in s4.columns:
            for sample, sub in s4.dropna(subset=["sample"]).groupby("sample"):
                tp53_status_vals = sub.get("TP53 status", pd.Series(dtype=str))
                tp53_status_vals = tp53_status_vals.dropna().astype(str).str.strip().unique()
                trans_types = (
                    sub.get("translocation_type", pd.Series(dtype=str))
                    .dropna().astype(str).str.strip().unique()
                )
                s4_lookup[str(sample).strip()] = {
                    "tp53_status": "altered",
                    "tp53_evidence": (
                        f"{len(sub)} SV breakpoints in Chen 2014 Table S4 "
                        f"(types: {','.join(sorted(trans_types)) or 'unknown'}; "
                        f"TP53 status: {','.join(sorted(tp53_status_vals)) or 'unannotated'})"
                    ),
                }
    except Exception as e:  # pragma: no cover — defensive
        print(f"[load_chen_2014_pcgp] warning: could not load Table S4b: {e}")

    # ---- ATRX + ALT FISH from Table S6 ----------------------------------
    atrx_lookup: dict[str, tuple[str, str, str, str]] = {}
    try:
        s6 = pd.read_excel(s6_path, sheet_name="Sheet1", header=None)
        # Hand-shaped layout (verified 2026-05-01):
        #   row 0  = title
        #   row 1  = blank
        #   row 2  = group header row: NaN | "ATRX IHC" | "Telomere FISH" | NaN
        #   row 3  = column-name row:  "Sample" | NaN | NaN | "ATRX status"
        #   row 4  = blank
        #   row 5+ = data
        # Columns are 0=Sample, 1=ATRX IHC, 2=Telomere FISH, 3=ATRX status
        # (mutation detail). NaN rows are skipped by the pd.isna(sample) check.
        for _, r in s6.iloc[5:, :].iterrows():
            sample = r.iloc[0]
            if pd.isna(sample):
                continue
            sample = str(sample).strip()
            atrx_ihc      = str(r.iloc[1]).strip() if not pd.isna(r.iloc[1]) else ""
            telomere_fish = str(r.iloc[2]).strip() if not pd.isna(r.iloc[2]) else ""
            atrx_mut      = str(r.iloc[3]).strip() if not pd.isna(r.iloc[3]) else ""

            # ATRX status: prefer the mutation column when present
            atrx_marker = atrx_mut or atrx_ihc
            atrx_marker_lower = atrx_marker.lower()
            if not atrx_marker or atrx_marker_lower in {"nan", "indeterminate", ""}:
                atrx_status = "unknown"
            elif atrx_marker_lower == "wt":
                atrx_status = "not_altered"
            else:
                atrx_status = "altered"

            telomere_lower = telomere_fish.lower()
            if "positive" in telomere_lower:
                alt_status = "alt_gold_positive"
            elif "negative" in telomere_lower:
                alt_status = "alt_negative_high_conf"
            elif telomere_lower in {"failed", "not performed", "nd", "nan", ""} \
                    or "fail" in telomere_lower or "not perform" in telomere_lower:
                alt_status = "alt_unknown"
            else:
                alt_status = "alt_unknown"

            atrx_lookup[sample] = (
                atrx_status,
                f"Chen 2014 Table S6: ATRX IHC='{atrx_ihc}'; mutation='{atrx_mut}'",
                alt_status,
                f"Chen 2014 Table S6: telomere FISH='{telomere_fish}'",
            )
    except Exception as e:  # pragma: no cover — defensive
        print(f"[load_chen_2014_pcgp] warning: could not load Table S6: {e}")

    if not sample_ids and (s4_lookup or atrx_lookup):
        sample_ids = sorted(set(s4_lookup) | set(atrx_lookup))

    def _patient_stem(sid: str) -> str:
        """
        S3 uses suffixed sample IDs (SJOS001_M, SJOS002_D, SJOS001103_D1)
        while S6 keys patients by the unsuffixed stem (SJOS001, SJOS002).
        Strip the trailing ``_<letters/digits>`` to recover the patient stem.
        """
        parts = sid.split("_", 1)
        return parts[0] if len(parts) >= 2 else sid

    rows = []
    for sid in sample_ids:
        row = {"sample_id": sid, "cohort": "PCGP_2014", "sample_type": "tumor"}

        s4_hit = s4_lookup.get(sid)
        if s4_hit:
            row["tp53_status"]   = s4_hit["tp53_status"]
            row["tp53_evidence"] = s4_hit["tp53_evidence"]
        else:
            row["tp53_status"]   = "unknown"
            row["tp53_evidence"] = "no TP53 SV in Chen 2014 Table S4 for this sample"

        # S6 keys are unsuffixed patient IDs; resolve via stem if direct miss.
        s6_hit = atrx_lookup.get(sid) or atrx_lookup.get(_patient_stem(sid))
        if s6_hit:
            atrx_s, _atrx_ev, alt_s, alt_ev = s6_hit
            row["atrx_status"]         = atrx_s
            row["alt_status"]          = alt_s
            row["alt_source"]          = "Chen 2014 Table S6 (telomere FISH)"
            row["alt_assay_evidence"]  = alt_ev
            row["source_confidence"]   = "peer_reviewed"
        else:
            row["atrx_status"]         = "unknown"
            row["alt_status"]          = "alt_unknown"
            row["alt_source"]          = "Chen 2014 (no Table S6 entry for this sample)"
            row["alt_assay_evidence"]  = "unknown"
            row["source_confidence"]   = "peer_reviewed"

        # Availability flags — Chen samples are tumor with WGS/WES SV calls
        row["has_wgs"]      = sid in s4_lookup
        row["has_clinical"] = (sid in atrx_lookup) or (_patient_stem(sid) in atrx_lookup)
        rows.append(row)

    df = pd.DataFrame(rows)
    return _fill_unknown_cols(df) if rows else _fill_unknown_cols(
        pd.DataFrame(columns=["sample_id", "cohort", "sample_type"])
    )


def load_chen_2014_tp53_svs(path: Path | None = None) -> pd.DataFrame:
    """
    Load Chen 2014 Table S4b (TP53 SV breakpoints) as a tidy DataFrame.

    Returned columns include ``sample, chra, posa, oria, chrb, posb, orib,
    translocation_type``. Coordinates are GRCh37 / hg19. This is the input
    consumed by ``sv_annotator.classify_chen_sv_architecture``.
    """
    p = path or DATA_ROOT / "raw" / "literature_tables" / "chen_2014_tableS4_TP53_SVs.xlsx"
    # See comment in ``load_chen_2014_pcgp`` for header offset rationale.
    df = pd.read_excel(p, sheet_name="TP53 SV breakpoints", header=3)
    df = df.rename(columns=lambda c: str(c).strip())
    if "sample" in df.columns:
        df = df.dropna(subset=["sample"])
        df["sample"] = df["sample"].astype(str).str.strip()
    return df


def load_taglialatela(path: Path | None = None) -> pd.DataFrame:
    p = path or LIT_DIR / "taglialatela_2026_celllines.csv"
    raw = pd.read_csv(p)
    raw.columns = [c.strip() for c in raw.columns]

    df = pd.DataFrame()
    df["sample_id"] = raw["cell_line"].str.strip()
    df["cohort"] = "Taglialatela_2026"
    df["sample_type"] = "cell_line"
    df["alt_status"] = raw["alt_status"].str.strip().map(
        {"positive": "alt_gold_positive", "negative": "alt_negative_high_conf"}
    ).fillna("alt_unknown")
    df["alt_source"] = raw["source"]
    df["source_confidence"] = "preprint"
    df["alt_assay_evidence"] = raw["evidence"]

    return _fill_unknown_cols(df)


def load_altitude_cellines(path: Path | None = None) -> pd.DataFrame:
    """
    Load ALTitude table1 (Sheet1) and extract curated ALT labels for
    osteosarcoma-relevant cell lines.

    Returns a long-form DataFrame with one row per cell line in ALTitude.
    Columns include ``altitude_pred_class`` (raw ALT_pred_class value: 'ALT+',
    'ALT-', or None) so callers can resolve direction for silver lines.
    """
    p = path or LIT_DIR / "ALTitude-table1.xlsx"
    raw = pd.read_excel(p, sheet_name="Sheet1", header=0)

    def _col(cols: list, *needles: str):
        for needle in needles:
            for c in cols:
                if needle.lower() in str(c).lower():
                    return c
        return None

    cols = list(raw.columns)
    achilles_col      = cols[0]
    strip_col         = _col(cols, "StripCI", "Strip") or cols[1]
    curated_col       = _col(cols, "Curated", "Marker")
    alt_source_col    = _col(cols, "ALT source")
    alt_assay_col     = _col(cols, "ALT_assay", "assay type")
    assay_strength_col = _col(cols, "assay strength")
    pred_class_col    = _col(cols, "ALT_pred_class")
    atrx_class_col    = _col(cols, "ATRX_Class")
    daxx_class_col    = _col(cols, "DAXX_Class")
    smarcal1_dep_col  = _col(cols, "SMARCAL1_dep")

    def _is_nonnull(val) -> bool:
        return str(val).strip() not in ("", "nan", "None", "NaN")

    df = pd.DataFrame()
    df["depmap_id"]      = raw[achilles_col].astype(str).str.strip()
    df["cell_line_name"] = raw[strip_col].astype(str).str.strip() if strip_col else "unknown"
    df["sample_id"]      = df["depmap_id"]
    df["cohort"]         = "ALTitude_DepMap"
    df["sample_type"]    = "cell_line"

    if curated_col:
        df["altitude_curated_marker"] = raw[curated_col].astype(str).str.strip()
    else:
        df["altitude_curated_marker"] = "unknown"

    # Read ALT_pred_class; store as altitude_pred_class (None for missing).
    if pred_class_col:
        def _clean_pc(v):
            s = str(v).strip()
            return s if s not in ("", "nan", "None", "NaN") else None
        df["altitude_pred_class"] = raw[pred_class_col].apply(_clean_pc)
    else:
        df["altitude_pred_class"] = None

    # Evidence tier: experimental if any of alt_source / assay_type / assay_strength is non-null.
    evidence_cols = [c for c in (alt_source_col, alt_assay_col, assay_strength_col) if c]
    if evidence_cols:
        tier_series = raw[evidence_cols].apply(
            lambda row: "experimental" if any(_is_nonnull(row[c]) for c in evidence_cols) else "prediction_only",
            axis=1,
        )
    else:
        tier_series = pd.Series(["prediction_only"] * len(raw), index=raw.index)
    df["altitude_evidence_tier"] = tier_series.values

    from .label_rules import map_altitude_call
    alt_status_list = []
    for i, row_idx in enumerate(raw.index):
        curated = df.loc[i, "altitude_curated_marker"] if i < len(df) else "unknown"
        mapped  = map_altitude_call(curated)
        if mapped == "alt_unknown":
            # Fall back to ALT_pred_class + tier for prediction-only lines.
            pc   = df.loc[i, "altitude_pred_class"]
            tier = df.loc[i, "altitude_evidence_tier"]
            if pc in ("ALT+", "ALT-") and tier == "prediction_only":
                mapped = "alt_silver"
        alt_status_list.append(mapped)
    df["alt_status"] = alt_status_list

    df["alt_source"]        = "ALTitude paper (peer-reviewed), Curated.Marker column"
    df["source_confidence"] = "peer_reviewed"
    df["alt_assay_evidence"] = raw[alt_assay_col].astype(str) if alt_assay_col else "unknown"

    if atrx_class_col:
        df["atrx_status"] = raw[atrx_class_col].apply(
            lambda x: "altered" if str(x).strip() not in ("nan", "", "None", "WT") else "not_altered"
        )
    else:
        df["atrx_status"] = "unknown"

    if daxx_class_col:
        df["daxx_status"] = raw[daxx_class_col].apply(
            lambda x: "altered" if str(x).strip() not in ("nan", "", "None", "WT") else "not_altered"
        )
    else:
        df["daxx_status"] = "unknown"

    if smarcal1_dep_col:
        df["smarcal1_status"] = "see_altitude_smarcal1_dep"
    else:
        df["smarcal1_status"] = "unknown"

    return _fill_unknown_cols(df)


# ---------------------------------------------------------------------------
# DepMap osteosarcoma manifest builder
# ---------------------------------------------------------------------------

def build_depmap_manifest(
    model_df: pd.DataFrame,
    expr_matrix: pd.DataFrame,
    damaging_matrix: pd.DataFrame,
    cn_gene_matrix: pd.DataFrame,
    fusions_df: pd.DataFrame,
    cn_segments_df: pd.DataFrame,
) -> pd.DataFrame:
    """
    Build manifest rows for DepMap osteosarcoma cell lines.
    """
    os_models = model_df[
        model_df["OncotreePrimaryDisease"].str.contains("Osteosarcoma", na=False)
    ].copy()

    cohort_ids = list(os_models["ModelID"].unique())
    rows = []

    for _, m in os_models.iterrows():
        model_id = m["ModelID"]
        row = {c: "unknown" for c in MANIFEST_COLUMNS}
        row["sample_id"] = model_id
        row["cohort"] = "DepMap_26Q1"
        row["sample_type"] = "cell_line"

        # Availability flags
        row["has_rna"] = model_id in expr_matrix.index
        row["has_wes"] = model_id in damaging_matrix.index
        row["has_wgs"] = model_id in cn_segments_df.get("ModelID", pd.Series()).values if "ModelID" in cn_segments_df.columns else False
        row["has_sc"] = False
        row["has_clinical"] = True  # DepMap metadata = clinical equivalent

        # TP53
        tp53_s, tp53_e = tp53_status_include_deep_deletion(
            model_id, damaging_matrix, cn_gene_matrix
        )
        row["tp53_status"] = tp53_s
        row["tp53_evidence"] = tp53_e

        # TCAB1 RNA
        g_q, rna_val, w_q = tcab1_rna_quantile(model_id, expr_matrix, cohort_ids)
        row["tcab1_rna_quantile"] = g_q if g_q is not None else "unknown"
        row["tcab1_rna_value"] = rna_val if rna_val is not None else "unknown"
        row["tcab1_rna_quantile_within_cohort"] = w_q if w_q is not None else "unknown"

        # TCAB1 functional status
        func_s, func_e = tcab1_functional_status(w_q)
        row["tcab1_functional_status"] = func_s
        row["tcab1_functional_evidence"] = func_e

        # TCAB1 genomic status — set unknown here; filled in Phase 3 (SV review)
        row["tcab1_genomic_status"] = "unknown"
        row["tcab1_genomic_evidence"] = "pending_sv_review"

        # TCAB1 overall
        overall, conf = tcab1_overall_and_confidence("unknown", func_s)
        row["tcab1_overall_status"] = overall
        row["tcab1_confidence"] = conf

        # SV architecture — unknown until Phase 3
        row["sv_architecture_17p13_1"] = "unknown"
        row["sv_architecture_evidence"] = "pending_sv_review"

        # ALT — filled by Phase 4 (literature merge) or from ALTitude
        row["alt_status"] = "alt_unknown"
        row["alt_source"] = "pending_literature_merge"
        row["source_confidence"] = "unknown"
        row["alt_assay_evidence"] = "unknown"

        # Telomerase / TERT
        telo_s, tert_v, telo_e = telomerase_status_from_tert(
            model_id, expr_matrix, cohort_ids
        )
        row["telomerase_status"] = telo_s
        row["tert_expression"] = tert_v if tert_v is not None else "unknown"
        row["telomerase_evidence"] = telo_e

        # ATRX, DAXX, SMARCAL1, TOP3A
        for gene, entrez in [
            ("ATRX", GENE_ENTREZ["ATRX"]),
            ("DAXX", GENE_ENTREZ["DAXX"]),
            ("SMARCAL1", GENE_ENTREZ["SMARCAL1"]),
            ("TOP3A", GENE_ENTREZ["TOP3A"]),
        ]:
            s, _ = gene_status_from_damaging_matrix(model_id, gene, entrez, damaging_matrix)
            row[f"{gene.lower()}_status"] = s

        # Purity — not available for cell lines
        row["purity_estimate"] = "not_applicable"
        row["purity_quartile"] = "not_applicable"

        row["notes"] = ""
        rows.append(row)

    df = pd.DataFrame(rows, columns=MANIFEST_COLUMNS)
    df = _apply_inclusion_flags(df)
    return df


def merge_literature_rows(
    depmap_manifest: pd.DataFrame,
    de_nonneville_df: pd.DataFrame,
    lawlor_pdx_df: pd.DataFrame,
    taglialatela_df: pd.DataFrame,
    pcgp_chen_df: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """
    Append literature cohort rows, deduplicating cell lines present in both
    DepMap and a literature table by cohort tag.

    Parameters
    ----------
    pcgp_chen_df : optional
        Chen 2014 PCGP osteosarcoma cohort frame from
        ``load_chen_2014_pcgp``. ``None`` (default) preserves backwards
        compatibility with the pre-PCGP pipeline.
    """
    lit_dfs = [de_nonneville_df, lawlor_pdx_df, taglialatela_df]
    if pcgp_chen_df is not None and len(pcgp_chen_df):
        lit_dfs.append(pcgp_chen_df)
    all_lit = pd.concat([_fill_unknown_cols(d) for d in lit_dfs], ignore_index=True)
    all_lit = all_lit[[c for c in MANIFEST_COLUMNS if c in all_lit.columns]]
    all_lit = _reindex_manifest(all_lit)
    all_lit = _apply_inclusion_flags(all_lit)

    combined = pd.concat([depmap_manifest, all_lit], ignore_index=True)
    return combined


def prefill_pcgp_sv_architecture(
    manifest_df: pd.DataFrame,
    chen_s4_df: pd.DataFrame,
) -> pd.DataFrame:
    """
    Set ``sv_architecture_17p13_1`` for PCGP_2014 samples by classifying
    each sample's Chen 2014 Table S4b rows in-frame (GRCh37) via
    :func:`sv_annotator.classify_chen_sv_architecture`.

    Mirrors :func:`sv_annotator.prefill_sv_architecture_from_literature` for
    Lawlor: only writes when the existing value is ``'unknown'`` (no silent
    overwrite of prior decisions). Updates ``sv_architecture_evidence`` with
    the classifier's reasoning.
    """
    from .sv_annotator import classify_chen_sv_architecture

    df = manifest_df.copy()
    if chen_s4_df is None or len(chen_s4_df) == 0:
        return df

    pcgp_mask = df["cohort"] == "PCGP_2014"
    if not pcgp_mask.any():
        return df

    sample_col = "sample" if "sample" in chen_s4_df.columns else None
    if sample_col is None:
        print("[prefill_pcgp_sv_architecture] no 'sample' column in chen_s4_df — skipped")
        return df

    by_sample = {
        str(s).strip(): sub for s, sub in chen_s4_df.groupby(sample_col)
    }

    n_filled = 0
    for sid in df.loc[pcgp_mask, "sample_id"].unique():
        rows = by_sample.get(str(sid).strip())
        if rows is None or len(rows) == 0:
            call, reason = "Neither", (
                "no TP53 SV in Chen 2014 Table S4 for this sample"
            )
        else:
            call, reason = classify_chen_sv_architecture(rows)

        target = (
            (df["sample_id"] == sid)
            & (df["cohort"] == "PCGP_2014")
            & (df["sv_architecture_17p13_1"] == "unknown")
        )
        if target.any():
            df.loc[target, "sv_architecture_17p13_1"] = call
            df.loc[target, "sv_architecture_evidence"] = (
                f"Chen 2014 Table S4 (in-frame GRCh37): {reason}"
            )
            n_filled += 1

    print(f"[prefill_pcgp_sv_architecture] filled {n_filled} PCGP samples")
    return df


# ---------------------------------------------------------------------------
# Architecture → genomic / overall status propagation
# ---------------------------------------------------------------------------

# WRAP53 (TCAB1) itself is hit by the breakpoint pattern.
_DISRUPTING_ARCHITECTURES = frozenset({
    "Both-single-SV",
    "Both-independent-SVs",
    "TCAB1-only",
    "Complex",
})

# Locus surveyed; WRAP53 itself is *not* hit.
_INTACT_ARCHITECTURES = frozenset({
    "Neither",
    "TP53-only",
})


def _genomic_status_from_arch(arch_call: str) -> str | None:
    """Map an SV architecture call to a ``tcab1_genomic_status`` string.

    Returns ``None`` (= "leave existing value alone") for architectures that
    don't carry locus-level information (``"unknown"`` or anything off-script).
    """
    if arch_call in _DISRUPTING_ARCHITECTURES:
        return "disrupted_high_conf"
    if arch_call in _INTACT_ARCHITECTURES:
        return "intact"
    return None


def refresh_tcab1_status_from_architecture(
    manifest_df: pd.DataFrame,
    *,
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Propagate ``sv_architecture_17p13_1`` into ``tcab1_genomic_status`` and
    ``tcab1_overall_status`` for samples whose architecture has been called
    (typically by ``prefill_sv_architecture_from_literature`` for Lawlor and
    ``prefill_pcgp_sv_architecture`` for the PCGP cohort).

    Mapping
    -------
        Both-single-SV / Both-independent-SVs / TCAB1-only / Complex
            → tcab1_genomic_status = 'disrupted_high_conf'
        Neither / TP53-only
            → tcab1_genomic_status = 'intact'
        unknown or anything else
            → no change

    Safety rails
    ------------
    * Never **downgrades** an already-informative ``tcab1_genomic_status``.
      A row that is already ``'disrupted_high_conf'``, ``'intact'``, or
      ``'candidate_regulatory_disruption'`` keeps its value, even if the
      architecture call would suggest otherwise (we treat the existing call
      as the authoritative DepMap/IHC-grounded evidence).
    * Recomputes ``tcab1_overall_status`` and ``tcab1_confidence`` via
      :func:`tcab1_overall_and_confidence` after updating
      ``tcab1_genomic_status`` so the overall call reflects the refreshed
      genomic evidence.

    This function exists because the prefill helpers only wrote
    ``sv_architecture_17p13_1`` / ``sv_architecture_evidence`` and left
    ``tcab1_genomic_status='unknown'`` on every literature cohort sample.
    Pass A's intact/disrupted counters key off ``tcab1_genomic_status``, so
    without this step PCGP "Neither" (n=23), PCGP "TP53-only" (n=9), and the
    Lawlor architecture-disrupted samples (n=4) were invisible to the
    feasibility audit.
    """
    df = manifest_df.copy()
    n_genomic, n_overall = 0, 0
    for idx in df.index:
        arch = str(df.at[idx, "sv_architecture_17p13_1"])
        derived = _genomic_status_from_arch(arch)
        if derived is None:
            continue

        cur_genomic = str(df.at[idx, "tcab1_genomic_status"])
        if cur_genomic == "unknown":
            df.at[idx, "tcab1_genomic_status"] = derived
            df.at[idx, "tcab1_genomic_evidence"] = (
                f"derived from sv_architecture_17p13_1='{arch}'"
            )
            n_genomic += 1
        # else: keep the prior (DepMap CN-derived) call, never downgrade.

        # Always re-derive overall_status from the (possibly updated) genomic
        # status + existing functional status. tcab1_overall_and_confidence
        # already handles the unknown / intact / disrupted / conflicted cases.
        func = str(df.at[idx, "tcab1_functional_status"])
        overall, conf = tcab1_overall_and_confidence(
            df.at[idx, "tcab1_genomic_status"], func
        )
        # Only commit the overall update if it strengthens or maintains the
        # call. "unknown" is treated as the weakest, "conflicted" as the
        # strongest signal we shouldn't silently overwrite.
        cur_overall = str(df.at[idx, "tcab1_overall_status"])
        if cur_overall == "unknown" and overall != "unknown":
            df.at[idx, "tcab1_overall_status"] = overall
            df.at[idx, "tcab1_confidence"] = conf
            n_overall += 1
        elif overall == "conflicted":
            # Surface genuine conflicts even over an existing call.
            df.at[idx, "tcab1_overall_status"] = overall
            df.at[idx, "tcab1_confidence"] = conf
            n_overall += 1

    if verbose:
        print(
            f"[refresh_tcab1_status_from_architecture] "
            f"updated tcab1_genomic_status for {n_genomic} samples, "
            f"tcab1_overall_status for {n_overall} samples."
        )
    return df


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _fill_unknown_cols(df: pd.DataFrame) -> pd.DataFrame:
    for col in MANIFEST_COLUMNS:
        if col not in df.columns:
            df[col] = "unknown"
    return df


def _reindex_manifest(df: pd.DataFrame) -> pd.DataFrame:
    for col in MANIFEST_COLUMNS:
        if col not in df.columns:
            df[col] = "unknown"
    return df[MANIFEST_COLUMNS].copy()


def _apply_inclusion_flags(df: pd.DataFrame) -> pd.DataFrame:
    flag_cols = [
        "include_aim1_primary_strict",
        "include_aim1_primary_pragmatic",
        "include_aim1_primary_broad",
        "include_aim2_strict",
        "include_aim2_pragmatic",
        "include_aim2_broad",
        "include_tp53_sensitivity",
        "include_telomerase_negative_secondary",
    ]
    flag_rows = df.apply(compute_inclusion_flags, axis=1)
    for col in flag_cols:
        df[col] = [r[col] for r in flag_rows]
    return df
