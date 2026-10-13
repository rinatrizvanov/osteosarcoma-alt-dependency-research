"""
SV/fusion/CN-segment annotator for the chr17:7,560,000-7,710,000 locus (GRCh38).

Annotates whether each breakpoint or segment boundary falls within defined
feature regions around TP53 and WRAP53 (TCAB1).

Feature coordinates are approximate GRCh38 values derived from Ensembl/GENCODE v45.
They will be refined to exact values when the reference GTF is downloaded in notebook 00.

Architecture call rules (deterministic, per label_spec_v0.9):
  WRAP53 coding hit + TP53 coding hit, same record  → Both-single-SV
  Both hits, different records                       → Both-independent-SVs
  Only WRAP53 hits                                   → TCAB1-only
  Only TP53 hits                                     → TP53-only
  ≥5 SVs in locus window                             → Complex (manual review)
  Otherwise                                          → Neither

Reference-genome handling
-------------------------
Two coordinate spaces are tracked:

* **GRCh38** — the project's working frame. DepMap 26Q1 + ALL annotated output
  use these coordinates. Anything ingested in another assembly is lifted before
  use OR classified in its native frame and only the categorical result is kept.

* **GRCh37 / hg19** — needed for two upstream sources:
    - CCLE 2019 SvABA SV calls → lifted to GRCh38 in ``annotate_ccle_svs``.
    - Chen 2014 PCGP osteosarcoma SVs → classified IN-FRAME via
      ``classify_chen_sv_architecture`` and only the architecture call is
      stored. Coordinates never leave GRCh37 in our pipeline.

GRCh37 coordinates were verified 2026-05-01 against:
  - Ensembl GRCh37 ENST00000269305 (TP53-001):   chr17:7,571,720-7,590,868 (-)
  - Ensembl GRCh37 ENSG00000141499 (WRAP53):     chr17:7,589,389-7,606,820 (+)
  - NCBI Gene 7157 / 55135 (GRCh37.p13)
TP53 and WRAP53 share a head-to-head bidirectional promoter; their gene bodies
overlap in approximately chr17:7,589,389-7,590,868 (~1.5 kb) where a single
breakpoint will legitimately flag BOTH genes.
"""

from __future__ import annotations

import pandas as pd
import numpy as np
from typing import Iterable, Literal


# --- GRCh38 (working frame) -------------------------------------------------
FEATURES = [
    # (name, chrom, start, end)  — 0-based half-open [start, end)
    ("TP53_coding",             "chr17", 7_668_402, 7_687_550),
    ("TP53_promoter",           "chr17", 7_687_551, 7_692_550),  # 5 kb upstream of TSS (minus strand)
    ("TP53_5kb_flank",          "chr17", 7_662_402, 7_668_401),  # 5 kb downstream of TES
    ("WRAP53_coding",           "chr17", 7_687_551, 7_725_564),
    ("WRAP53_promoter",         "chr17", 7_682_551, 7_687_550),  # 5 kb upstream of TSS (plus strand)
    ("WRAP53_5kb_flank",        "chr17", 7_725_565, 7_730_564),
    ("intergenic_TP53_WRAP53",  "chr17", 7_687_000, 7_688_000),
]

LOCUS_CHROM = "chr17"
LOCUS_START = 7_560_000
LOCUS_END   = 7_710_000


# --- GRCh37 / hg19 (used for CCLE prefilter + Chen 2014 in-frame call) ------
GRCH37_LOCUS_CHROM    = "chr17"
GRCH37_LOCUS_WINDOW   = (7_540_000, 7_690_000)  # generous neighborhood; spec value
GRCH37_TP53_CODING    = (7_571_720, 7_590_868)  # Ensembl ENST00000269305 (TP53-001), exact
GRCH37_TP53_INTRON1   = (7_577_000, 7_590_000)  # Chen 2014 paper highlights this region
GRCH37_WRAP53_CODING  = (7_589_389, 7_606_820)  # Ensembl ENSG00000141499, exact gene body


COMPLEX_THRESHOLD = 5  # ≥ this many distinct SV records in the window → Complex

ArchitectureCall = Literal[
    "Both-single-SV",
    "Both-independent-SVs",
    "TCAB1-only",
    "TP53-only",
    "Complex",
    "Neither",
    "unknown",
]


def _pos_in_locus(chrom: str, pos: int) -> bool:
    return chrom == LOCUS_CHROM and LOCUS_START <= pos < LOCUS_END


def _features_for_pos(chrom: str, pos: int) -> list[str]:
    if chrom != LOCUS_CHROM:
        return []
    return [
        name for name, fc, fs, fe in FEATURES
        if fc == chrom and fs <= pos < fe
    ]


def _segment_overlaps_feature(
    chrom: str, seg_start: int, seg_end: int, feat_start: int, feat_end: int
) -> bool:
    return chrom == LOCUS_CHROM and seg_start < feat_end and seg_end > feat_start


def _features_for_segment(chrom: str, seg_start: int, seg_end: int) -> list[str]:
    return [
        name for name, fc, fs, fe in FEATURES
        if fc == chrom and seg_start < fe and seg_end > fs
    ]


def annotate_breakpoints(
    sv_df: pd.DataFrame,
    sample_id: str,
    chrom_col: str = "Chrom",
    pos_col: str = "Pos",
    record_id_col: str | None = None,
    extra_cols: list[str] | None = None,
) -> pd.DataFrame:
    """
    Annotate point breakpoints (SNV/indel-style coordinates).

    Parameters
    ----------
    sv_df : rows with at least chrom_col and pos_col
    sample_id : for the output records
    chrom_col, pos_col : column names for chromosome and position
    record_id_col : optional column to carry through as sv_record_id
    extra_cols : additional columns to include in output

    Returns
    -------
    DataFrame with one row per breakpoint that falls in the locus window.
    Columns: sample_id, sv_record_id, chrom, pos, features, in_locus
    """
    rows = []
    for idx, r in sv_df.iterrows():
        chrom = str(r.get(chrom_col, "")).strip()
        try:
            pos = int(r[pos_col])
        except (TypeError, ValueError):
            continue

        in_locus = _pos_in_locus(chrom, pos)
        feats = _features_for_pos(chrom, pos) if in_locus else []

        row = {
            "sample_id": sample_id,
            "sv_record_id": r[record_id_col] if record_id_col else idx,
            "chrom": chrom,
            "pos": pos,
            "features": "|".join(feats) if feats else "",
            "in_locus": in_locus,
        }
        if extra_cols:
            for col in extra_cols:
                row[col] = r.get(col)
        rows.append(row)

    return pd.DataFrame(rows)


def annotate_segments(
    seg_df: pd.DataFrame,
    sample_id: str,
    chrom_col: str = "Chrom",
    start_col: str = "Start",
    end_col: str = "End",
    state_col: str | None = "State",
    record_id_col: str | None = None,
    extra_cols: list[str] | None = None,
) -> pd.DataFrame:
    """
    Annotate CN segments overlapping the locus.

    Returns one row per segment × feature overlap.
    """
    rows = []
    for idx, r in seg_df.iterrows():
        chrom = str(r.get(chrom_col, "")).strip()
        try:
            seg_start = int(r[start_col])
            seg_end = int(r[end_col])
        except (TypeError, ValueError):
            continue

        # Keep only segments that overlap the locus window
        if chrom != LOCUS_CHROM:
            continue
        if seg_start >= LOCUS_END or seg_end <= LOCUS_START:
            continue

        feats = _features_for_segment(chrom, seg_start, seg_end)
        row = {
            "sample_id": sample_id,
            "sv_record_id": r[record_id_col] if record_id_col else idx,
            "chrom": chrom,
            "seg_start": seg_start,
            "seg_end": seg_end,
            "features": "|".join(feats) if feats else "",
            "in_locus": True,
        }
        if state_col and state_col in r.index:
            row["state"] = r[state_col]
        if extra_cols:
            for col in extra_cols:
                row[col] = r.get(col)
        rows.append(row)

    return pd.DataFrame(rows)


def annotate_fusions(
    fusion_df: pd.DataFrame,
    sample_id: str,
    gene1_col: str = "Gene1",
    gene2_col: str = "Gene2",
    canonical_col: str = "CanonicalFusionName",
    record_id_col: str | None = None,
    extra_cols: list[str] | None = None,
) -> pd.DataFrame:
    """
    Flag RNA fusions involving WRAP53 or TP53 as SV-equivalent events.

    Returns one row per fusion involving the target genes.
    """
    target_genes = {"WRAP53", "TP53", "TCAB1"}
    rows = []
    for idx, r in fusion_df.iterrows():
        g1 = str(r.get(gene1_col, "")).strip()
        g2 = str(r.get(gene2_col, "")).strip()
        if not (target_genes & {g1, g2}):
            continue

        involved = sorted(target_genes & {g1, g2})
        feats = []
        if "TP53" in involved:
            feats.append("TP53_coding")
        if "WRAP53" in involved or "TCAB1" in involved:
            feats.append("WRAP53_coding")

        row = {
            "sample_id": sample_id,
            "sv_record_id": r[record_id_col] if record_id_col else idx,
            "event_type": "fusion",
            "canonical_name": r.get(canonical_col, ""),
            "gene1": g1,
            "gene2": g2,
            "genes_involved": "|".join(involved),
            "features": "|".join(feats),
            "in_locus": True,
        }
        if extra_cols:
            for col in extra_cols:
                row[col] = r.get(col)
        rows.append(row)

    return pd.DataFrame(rows)


def annotate_locus_svs(
    sv_df: pd.DataFrame,
    sample_id: str,
    sv_type: Literal["breakpoints", "segments", "fusions"] = "breakpoints",
    **kwargs,
) -> pd.DataFrame:
    """
    Unified entry point.

    sv_type : 'breakpoints' | 'segments' | 'fusions'
    kwargs are forwarded to the type-specific function.
    """
    if sv_type == "breakpoints":
        return annotate_breakpoints(sv_df, sample_id, **kwargs)
    if sv_type == "segments":
        return annotate_segments(sv_df, sample_id, **kwargs)
    if sv_type == "fusions":
        return annotate_fusions(sv_df, sample_id, **kwargs)
    raise ValueError(f"Unknown sv_type: {sv_type!r}")


# ---------------------------------------------------------------------------
# Architecture call
# ---------------------------------------------------------------------------

def architecture_call(annotated_df: pd.DataFrame) -> tuple[ArchitectureCall, str]:
    """
    Derive a deterministic architecture call from the annotated SV records
    for a single sample.

    Returns (call, reasoning_string).
    """
    if annotated_df.empty:
        return "Neither", "No SVs in locus window"

    n_records = annotated_df["sv_record_id"].nunique()
    if n_records >= COMPLEX_THRESHOLD:
        return "Complex", f"{n_records} distinct SV records in locus — manual review required"

    def _has_hit(feature_prefix: str) -> set:
        """Return set of sv_record_ids that hit a given feature prefix."""
        mask = annotated_df["features"].str.contains(feature_prefix, na=False)
        return set(annotated_df.loc[mask, "sv_record_id"].unique())

    wrap53_records = _has_hit("WRAP53_coding")
    tp53_records = _has_hit("TP53_coding")

    both_same = wrap53_records & tp53_records
    if both_same:
        return "Both-single-SV", (
            f"Records hitting both WRAP53_coding and TP53_coding: "
            f"{sorted(both_same)}"
        )

    if wrap53_records and tp53_records:
        return "Both-independent-SVs", (
            f"WRAP53_coding records: {sorted(wrap53_records)}; "
            f"TP53_coding records: {sorted(tp53_records)}"
        )

    if wrap53_records:
        return "TCAB1-only", f"WRAP53_coding records: {sorted(wrap53_records)}"

    if tp53_records:
        return "TP53-only", f"TP53_coding records: {sorted(tp53_records)}"

    locus_records = annotated_df["sv_record_id"].nunique()
    return "Neither", (
        f"{locus_records} SV(s) in locus window but none hit TP53_coding or WRAP53_coding — "
        "flag for review"
    )


# ---------------------------------------------------------------------------
# Pre-fill from literature SV calls (Lawlor_2026_GCCRI)
# ---------------------------------------------------------------------------

LAWLOR_SV_ARCHITECTURE = {
    "CPRIT 1981": "Both-single-SV",
    "0560-000":   "Both-single-SV",
    "OS33":       "Both-independent-SVs",
    "OS37":       "Both-independent-SVs",
    "OS39":       "Both-single-SV",
}

def prefill_sv_architecture_from_literature(
    manifest_df: pd.DataFrame,
) -> pd.DataFrame:
    """
    Set sv_architecture_17p13_1 for Lawlor_2026_GCCRI samples using
    the calls in tcab1_pdx_sv_calls.csv (as specified in Phase 3).

    Writes to sv_architecture_17p13_1 only where the column is currently 'unknown'.
    """
    df = manifest_df.copy()
    for sample_id, call in LAWLOR_SV_ARCHITECTURE.items():
        mask = (
            df["sample_id"].str.contains(sample_id, regex=False)
            & (df["sv_architecture_17p13_1"] == "unknown")
        )
        df.loc[mask, "sv_architecture_17p13_1"] = call
        df.loc[mask, "sv_architecture_evidence"] = (
            "Pre-filled from tcab1_pdx_sv_calls.csv (Lawlor 2026, Figure 2)"
        )
    return df


# ---------------------------------------------------------------------------
# CCLE 2019 SvABA breakpoint annotation (GRCh37 → GRCh38 via liftover)
# ---------------------------------------------------------------------------

def _parse_ccle_bp(bp_str: object) -> tuple[str | None, int | None, str | None]:
    """
    Parse a CCLE bp1/bp2 string of the form ``'CHR:START-END[+-]'``.

    The CCLE SvABA file packs each breakpoint endpoint as a single string;
    typically START == END for point breakpoints. The trailing character is
    the strand (+ / -). Returns ``(chrom_with_chr_prefix, pos, strand)`` or
    ``(None, None, None)`` if the field cannot be parsed.
    """
    if bp_str is None:
        return None, None, None
    s = str(bp_str).strip()
    if not s or s.lower() in {"nan", "none", ""}:
        return None, None, None
    try:
        chrom_part, rest = s.split(":", 1)
    except ValueError:
        return None, None, None

    strand = None
    if rest and rest[-1] in {"+", "-"}:
        strand = rest[-1]
        rest = rest[:-1]

    if "-" in rest:
        start_s, _end_s = rest.split("-", 1)
    else:
        start_s = rest

    try:
        pos = int(start_s)
    except ValueError:
        return None, None, None

    chrom = chrom_part.strip()
    if not chrom.startswith("chr"):
        chrom = f"chr{chrom}"
    return chrom, pos, strand


def annotate_ccle_svs(
    ccle_df: pd.DataFrame,
    ccle_to_model: dict[str, str],
    allowed_model_ids: Iterable[str] | None = None,
    *,
    bp1_col: str = "bp1",
    bp2_col: str = "bp2",
    sample_col: str = "CCLE_name",
    record_id_col: str = "map_id",
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Annotate CCLE 2019 SvABA structural-variant breakpoints at the TP53/WRAP53
    locus, lifting GRCh37 coordinates to GRCh38 so the result is directly
    consumable by ``architecture_call``.

    Pipeline (per spec):
      1. Filter rows whose chr17 endpoint(s) fall in ``GRCH37_LOCUS_WINDOW``.
      2. Lift each kept endpoint hg19 → hg38 via ``utils.liftover_breakpoint``.
         Liftover failures are logged and the failing endpoint is dropped (the
         partner endpoint is still emitted if it lifts successfully).
      3. Map ``CCLE_name`` → DepMap ``ModelID`` via ``ccle_to_model``.
         Rows whose CCLE name is not in the map are dropped (and counted).
      4. If ``allowed_model_ids`` is provided, restrict to those models
         (typically the osteosarcoma cohort).
      5. Emit one row per surviving endpoint with the same schema as
         ``annotate_breakpoints``: ``sample_id, sv_record_id, chrom, pos,
         features, in_locus`` plus ``ccle_name``, ``endpoint`` (a/b),
         ``strand``, ``svclass``, ``gene1``, ``gene2`` for traceability.

    Parameters
    ----------
    ccle_df : DataFrame
        Raw CCLE SvABA table (sheet ``CCLE_translocations_SvABA_20181``).
    ccle_to_model : dict[str, str]
        ``{CCLEName -> ModelID}``. Build with
        ``model_df.set_index('CCLEName')['ModelID'].to_dict()``.
    allowed_model_ids : iterable of str, optional
        If given, results are restricted to these ModelIDs (e.g. the
        DepMap osteosarcoma cohort). Pass ``None`` to keep all mapped rows.
    verbose : bool, default True
        Print a summary of liftover and mapping outcomes.

    Returns
    -------
    DataFrame in the schema of ``annotate_breakpoints`` (one row per
    breakpoint endpoint that survived all filters), pre-annotated with
    ``features`` and ``in_locus`` in GRCh38.
    """
    allowed = set(allowed_model_ids) if allowed_model_ids is not None else None
    win_start, win_end = GRCH37_LOCUS_WINDOW

    n_rows_total = len(ccle_df)
    n_rows_in_window = 0
    n_unmapped_ccle = 0
    n_unmapped_filtered = 0
    n_lift_fail = 0
    out_rows: list[dict] = []
    unmapped_names: set[str] = set()

    for idx, r in ccle_df.iterrows():
        bp1_chrom_37, bp1_pos_37, bp1_strand = _parse_ccle_bp(r.get(bp1_col))
        bp2_chrom_37, bp2_pos_37, bp2_strand = _parse_ccle_bp(r.get(bp2_col))

        endpoints_37 = []
        if bp1_chrom_37 == "chr17" and bp1_pos_37 is not None and win_start <= bp1_pos_37 < win_end:
            endpoints_37.append(("a", bp1_chrom_37, bp1_pos_37, bp1_strand))
        if bp2_chrom_37 == "chr17" and bp2_pos_37 is not None and win_start <= bp2_pos_37 < win_end:
            endpoints_37.append(("b", bp2_chrom_37, bp2_pos_37, bp2_strand))

        if not endpoints_37:
            continue
        n_rows_in_window += 1

        ccle_name = str(r.get(sample_col, "")).strip()
        model_id = ccle_to_model.get(ccle_name)
        if not model_id:
            n_unmapped_ccle += 1
            unmapped_names.add(ccle_name)
            continue
        if allowed is not None and model_id not in allowed:
            n_unmapped_filtered += 1
            continue

        sv_record_id = r.get(record_id_col, idx)

        for endpoint, chrom_37, pos_37, strand in endpoints_37:
            from .utils import liftover_breakpoint
            chrom_38, pos_38 = liftover_breakpoint(chrom_37, pos_37, "hg19", "hg38")
            if chrom_38 is None or pos_38 is None:
                n_lift_fail += 1
                if verbose:
                    print(
                        f"  [annotate_ccle_svs] liftover FAILED for "
                        f"{ccle_name} (model {model_id}, sv {sv_record_id}, "
                        f"endpoint {endpoint}): {chrom_37}:{pos_37}"
                    )
                continue

            in_locus = chrom_38 == LOCUS_CHROM and LOCUS_START <= pos_38 < LOCUS_END
            feats = _features_for_pos(chrom_38, pos_38) if in_locus else []
            out_rows.append({
                "sample_id":      model_id,
                "sv_record_id":   sv_record_id,
                "chrom":          chrom_38,
                "pos":            pos_38,
                "features":       "|".join(feats) if feats else "",
                "in_locus":       in_locus,
                "ccle_name":      ccle_name,
                "endpoint":       endpoint,
                "strand":         strand,
                "chrom_grch37":   chrom_37,
                "pos_grch37":     pos_37,
                "svclass":        r.get("class"),
                "gene1":          r.get("gene1"),
                "gene2":          r.get("gene2"),
            })

    if verbose:
        print(
            f"  [annotate_ccle_svs] CCLE rows scanned: {n_rows_total:,} | "
            f"with chr17 endpoint in GRCh37 locus window: {n_rows_in_window} | "
            f"unmapped CCLE_name → ModelID: {n_unmapped_ccle} | "
            f"filtered out by allowed_model_ids: {n_unmapped_filtered} | "
            f"liftover failures: {n_lift_fail} | "
            f"output rows: {len(out_rows)}"
        )
        if unmapped_names and len(unmapped_names) <= 20:
            print(f"  [annotate_ccle_svs] unmapped CCLE names: "
                  f"{sorted(unmapped_names)}")

    cols = [
        "sample_id", "sv_record_id", "chrom", "pos", "features", "in_locus",
        "ccle_name", "endpoint", "strand",
        "chrom_grch37", "pos_grch37",
        "svclass", "gene1", "gene2",
    ]
    return pd.DataFrame(out_rows, columns=cols)


# ---------------------------------------------------------------------------
# Chen 2014 PCGP architecture classifier (in-frame GRCh37; no liftover)
# ---------------------------------------------------------------------------

def _chen_endpoint_hits(chrom: object, pos: object) -> set[str]:
    """Return the set of feature labels hit by one (chrom, pos) in GRCh37."""
    try:
        chrom_str = str(chrom).strip().replace("chr", "")
        if chrom_str != "17":
            return set()
        p = int(pos)
    except (TypeError, ValueError):
        return set()

    hits: set[str] = set()
    if GRCH37_TP53_CODING[0] <= p < GRCH37_TP53_CODING[1]:
        hits.add("TP53_coding")
    if GRCH37_TP53_INTRON1[0] <= p < GRCH37_TP53_INTRON1[1]:
        hits.add("TP53_intron1")
    if GRCH37_WRAP53_CODING[0] <= p < GRCH37_WRAP53_CODING[1]:
        hits.add("WRAP53_coding")
    return hits


def classify_chen_sv_architecture(
    s4_rows_for_sample: pd.DataFrame,
    *,
    chra_col: str = "chra",
    posa_col: str = "posa",
    chrb_col: str = "chrb",
    posb_col: str = "posb",
) -> tuple[ArchitectureCall, str]:
    """
    Classify ``sv_architecture_17p13_1`` for one PCGP sample using
    Chen 2014 Table S4b (TP53 SV breakpoints).

    Each row in S4b is one validated SV with both endpoints listed
    (``chra/posa`` and ``chrb/posb``) in GRCh37 coordinates. We classify
    in-frame: no liftover is performed, and only the categorical result is
    stored in the manifest.

    Endpoint→feature rule (GRCh37, half-open):
      * TP53 coding   = ``GRCH37_TP53_CODING``  (Ensembl ENST00000269305)
      * TP53 intron 1 = ``GRCH37_TP53_INTRON1`` (Chen 2014 paper focus)
      * WRAP53 coding = ``GRCH37_WRAP53_CODING`` (Ensembl ENSG00000141499)

    A break in the head-to-head bidirectional region (chr17:7,589,389-7,590,868
    in GRCh37, ~1.5 kb) legitimately hits BOTH gene bodies; the classifier
    treats this as a single SV that hits both genes simultaneously.

    Architecture call rule (matches ``architecture_call``):
      * ≥ ``COMPLEX_THRESHOLD`` distinct SVs in window → ``Complex``
      * ≥ 1 SV with both ends touching TP53 and WRAP53 → ``Both-single-SV``
      * Multiple distinct SVs collectively touching both → ``Both-independent-SVs``
      * Only TP53 hit (any SV)                          → ``TP53-only``
      * Only WRAP53 hit (any SV)                        → ``TCAB1-only``
      * Otherwise                                       → ``Neither``

    Returns
    -------
    (call, reasoning) : (str, str)
        ``call`` is one of the ``ArchitectureCall`` values.

    Examples
    --------
    >>> # SJOS001 in Chen 2014 has 3 CTXs into TP53 intron 1 (chr17:7,583,xxx)
    >>> # All hit only TP53 → 'TP53-only'.
    >>> # A single SV with chra=17, posa=7,583,500 (TP53 intron 1)
    >>> #               and chrb=17, posb=7,600,000 (WRAP53 coding)
    >>> # → 'Both-single-SV'.
    """
    if s4_rows_for_sample is None or len(s4_rows_for_sample) == 0:
        return "Neither", "no SVs in Chen 2014 Table S4 for this sample"

    distinct_svs: list[set[str]] = []
    for _, r in s4_rows_for_sample.iterrows():
        hits = _chen_endpoint_hits(r.get(chra_col), r.get(posa_col)) | \
               _chen_endpoint_hits(r.get(chrb_col), r.get(posb_col))
        distinct_svs.append(hits)

    n_in_window = sum(1 for h in distinct_svs if h)

    if n_in_window >= COMPLEX_THRESHOLD:
        return "Complex", (
            f"{n_in_window} distinct SVs touching TP53/WRAP53 in GRCh37 — "
            f"manual review required"
        )

    has_both_in_one = any({"TP53_coding", "WRAP53_coding"} <= h for h in distinct_svs)
    any_tp53   = any("TP53_coding" in h or "TP53_intron1" in h for h in distinct_svs)
    any_wrap53 = any("WRAP53_coding" in h for h in distinct_svs)

    if has_both_in_one:
        return "Both-single-SV", (
            f"{len(distinct_svs)} SV(s) in S4; ≥1 has both endpoints hitting "
            f"TP53 and WRAP53 (head-to-head region)"
        )
    if any_tp53 and any_wrap53:
        return "Both-independent-SVs", (
            f"{len(distinct_svs)} SVs total; TP53 and WRAP53 each hit by "
            f"different SV records"
        )
    if any_tp53:
        return "TP53-only", (
            f"{len(distinct_svs)} SV(s); only TP53 (incl. intron 1) hit, "
            f"no WRAP53 breakpoint"
        )
    if any_wrap53:
        return "TCAB1-only", (
            f"{len(distinct_svs)} SV(s); only WRAP53 hit, no TP53 breakpoint"
        )
    return "Neither", (
        f"{len(distinct_svs)} SV(s) in S4 but none hit TP53/WRAP53 "
        f"in GRCh37 windows"
    )
