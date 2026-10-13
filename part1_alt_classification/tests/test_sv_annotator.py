"""
Tests for src/sv_annotator.py and the GRCh37 helpers / liftover utility.

Synthetic SV inputs covering all architecture call branches plus the
new CCLE liftover and Chen 2014 in-frame classifier behaviour.

Run with: pytest tests/test_sv_annotator.py -v
"""

import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from src.sv_annotator import (
    annotate_breakpoints,
    annotate_segments,
    annotate_fusions,
    architecture_call,
    classify_chen_sv_architecture,
    _parse_ccle_bp,
    LOCUS_CHROM,
    LOCUS_START,
    LOCUS_END,
    COMPLEX_THRESHOLD,
    GRCH37_TP53_CODING,
    GRCH37_WRAP53_CODING,
    GRCH37_TP53_INTRON1,
    GRCH37_LOCUS_WINDOW,
)
from src.utils import liftover_breakpoint


# ---------------------------------------------------------------------------
# Helpers to build synthetic inputs
# ---------------------------------------------------------------------------

def _bp_df(positions: list[tuple[str, int]]) -> pd.DataFrame:
    """positions = [(chrom, pos), ...]"""
    return pd.DataFrame(
        [{"Chrom": c, "Pos": p} for c, p in positions]
    )


def _seg_df(segments: list[tuple[str, int, int, str]]) -> pd.DataFrame:
    """segments = [(chrom, start, end, state), ...]"""
    return pd.DataFrame(
        [{"Chrom": c, "Start": s, "End": e, "State": st} for c, s, e, st in segments]
    )


def _fusion_df(pairs: list[tuple[str, str]]) -> pd.DataFrame:
    return pd.DataFrame(
        [{"Gene1": g1, "Gene2": g2, "CanonicalFusionName": f"{g1}--{g2}"} for g1, g2 in pairs]
    )


# ---------------------------------------------------------------------------
# annotate_breakpoints
# ---------------------------------------------------------------------------

class TestAnnotateBreakpoints:
    def test_in_locus_wrap53(self):
        df = _bp_df([("chr17", 7_700_000)])  # inside WRAP53_coding
        out = annotate_breakpoints(df, "S1")
        assert len(out) == 1
        assert out.iloc[0]["in_locus"] == True
        assert "WRAP53_coding" in out.iloc[0]["features"]

    def test_in_locus_tp53(self):
        df = _bp_df([("chr17", 7_675_000)])  # inside TP53_coding
        out = annotate_breakpoints(df, "S1")
        assert len(out) == 1
        assert "TP53_coding" in out.iloc[0]["features"]

    def test_outside_locus(self):
        df = _bp_df([("chr1", 1_000_000)])
        out = annotate_breakpoints(df, "S1")
        assert len(out) == 1
        assert out.iloc[0]["in_locus"] == False
        assert out.iloc[0]["features"] == ""

    def test_empty_input(self):
        df = pd.DataFrame(columns=["Chrom", "Pos"])
        out = annotate_breakpoints(df, "S1")
        assert len(out) == 0

    def test_sample_id_propagated(self):
        df = _bp_df([("chr17", 7_700_000)])
        out = annotate_breakpoints(df, "MY_SAMPLE")
        assert out.iloc[0]["sample_id"] == "MY_SAMPLE"

    def test_promoter_region(self):
        df = _bp_df([("chr17", 7_684_000)])  # WRAP53_promoter range
        out = annotate_breakpoints(df, "S1")
        assert "WRAP53_promoter" in out.iloc[0]["features"]

    def test_flank_region(self):
        df = _bp_df([("chr17", 7_664_000)])  # TP53_5kb_flank range
        out = annotate_breakpoints(df, "S1")
        assert "TP53_5kb_flank" in out.iloc[0]["features"]


# ---------------------------------------------------------------------------
# annotate_segments
# ---------------------------------------------------------------------------

class TestAnnotateSegments:
    def test_segment_overlaps_wrap53(self):
        df = _seg_df([("chr17", 7_690_000, 7_710_000, "DEL")])
        out = annotate_segments(df, "S1")
        assert len(out) >= 1
        features = " ".join(out["features"].tolist())
        assert "WRAP53" in features

    def test_segment_outside_locus(self):
        df = _seg_df([("chr17", 8_000_000, 9_000_000, "DEL")])
        out = annotate_segments(df, "S1")
        assert len(out) == 0

    def test_segment_wrong_chrom(self):
        df = _seg_df([("chr1", 7_680_000, 7_700_000, "DEL")])
        out = annotate_segments(df, "S1")
        assert len(out) == 0

    def test_segment_spanning_both_genes(self):
        df = _seg_df([("chr17", 7_665_000, 7_730_000, "DEL")])
        out = annotate_segments(df, "S1")
        features = " ".join(out["features"].tolist())
        assert "TP53" in features
        assert "WRAP53" in features


# ---------------------------------------------------------------------------
# annotate_fusions
# ---------------------------------------------------------------------------

class TestAnnotateFusions:
    def test_tp53_fusion(self):
        df = _fusion_df([("TP53", "MDM2")])
        out = annotate_fusions(df, "S1")
        assert len(out) == 1
        assert "TP53_coding" in out.iloc[0]["features"]

    def test_wrap53_fusion(self):
        df = _fusion_df([("WRAP53", "BRCA1")])
        out = annotate_fusions(df, "S1")
        assert len(out) == 1
        assert "WRAP53_coding" in out.iloc[0]["features"]

    def test_unrelated_fusion(self):
        df = _fusion_df([("EGFR", "FGFR1")])
        out = annotate_fusions(df, "S1")
        assert len(out) == 0

    def test_both_genes_in_fusion(self):
        df = _fusion_df([("TP53", "WRAP53")])
        out = annotate_fusions(df, "S1")
        assert len(out) == 1
        features = out.iloc[0]["features"]
        assert "TP53_coding" in features
        assert "WRAP53_coding" in features


# ---------------------------------------------------------------------------
# architecture_call
# ---------------------------------------------------------------------------

class TestArchitectureCall:
    def _make_annotated(self, record_features: list[tuple[int, str]]) -> pd.DataFrame:
        """record_features = [(record_id, features_str), ...]"""
        return pd.DataFrame(
            [{"sv_record_id": rid, "features": feats} for rid, feats in record_features]
        )

    def test_both_single_sv(self):
        df = self._make_annotated([(1, "TP53_coding|WRAP53_coding")])
        call, reason = architecture_call(df)
        assert call == "Both-single-SV"

    def test_both_independent(self):
        df = self._make_annotated([(1, "TP53_coding"), (2, "WRAP53_coding")])
        call, reason = architecture_call(df)
        assert call == "Both-independent-SVs"

    def test_tcab1_only(self):
        df = self._make_annotated([(1, "WRAP53_coding|WRAP53_promoter")])
        call, reason = architecture_call(df)
        assert call == "TCAB1-only"

    def test_tp53_only(self):
        df = self._make_annotated([(1, "TP53_coding")])
        call, reason = architecture_call(df)
        assert call == "TP53-only"

    def test_neither(self):
        df = self._make_annotated([(1, "TP53_promoter"), (2, "WRAP53_5kb_flank")])
        call, reason = architecture_call(df)
        assert call == "Neither"

    def test_complex(self):
        records = [(i, "TP53_5kb_flank") for i in range(COMPLEX_THRESHOLD)]
        df = self._make_annotated(records)
        call, reason = architecture_call(df)
        assert call == "Complex"

    def test_empty_df(self):
        df = pd.DataFrame(columns=["sv_record_id", "features"])
        call, reason = architecture_call(df)
        assert call == "Neither"

    def test_reason_string_nonempty(self):
        df = self._make_annotated([(1, "TP53_coding|WRAP53_coding")])
        call, reason = architecture_call(df)
        assert len(reason) > 0


# ---------------------------------------------------------------------------
# Integration: annotate_breakpoints → architecture_call
# ---------------------------------------------------------------------------

class TestIntegration:
    def test_both_single_from_real_coords(self):
        """
        A single SV with one breakpoint in TP53 coding and another in WRAP53 coding
        should call Both-single-SV when given the same record_id.
        """
        df = pd.DataFrame([
            {"Chrom": "chr17", "Pos": 7_675_000, "record_id": "sv1"},
            {"Chrom": "chr17", "Pos": 7_700_000, "record_id": "sv1"},
        ])
        annotated = annotate_breakpoints(df, "TEST_SAMPLE", record_id_col="record_id")
        call, reason = architecture_call(annotated)
        assert call == "Both-single-SV"

    def test_independent_svs(self):
        df = pd.DataFrame([
            {"Chrom": "chr17", "Pos": 7_675_000, "record_id": "sv1"},
            {"Chrom": "chr17", "Pos": 7_700_000, "record_id": "sv2"},
        ])
        annotated = annotate_breakpoints(df, "TEST_SAMPLE", record_id_col="record_id")
        call, _ = architecture_call(annotated)
        assert call == "Both-independent-SVs"

    def test_segment_focal_deletion_both(self):
        """Single focal deletion spanning both TP53 and WRAP53 → Both-single-SV via segments."""
        df = _seg_df([("chr17", 7_665_000, 7_730_000, "DEL")])
        annotated = annotate_segments(df, "TEST_SAMPLE")
        call, _ = architecture_call(annotated)
        assert call == "Both-single-SV"


# ---------------------------------------------------------------------------
# CCLE breakpoint string parser
# ---------------------------------------------------------------------------

class TestParseCcleBp:
    """The CCLE SvABA file packs breakpoints as 'CHR:START-END[+-]'."""

    @pytest.mark.parametrize("s, expected", [
        ("13:33777850-33777850+", ("chr13", 33777850, "+")),
        ("17:7583762-7583762-",   ("chr17", 7583762, "-")),
        ("17:7583762-7583762",    ("chr17", 7583762, None)),
        ("X:155270000-155270000+", ("chrX", 155270000, "+")),
    ])
    def test_valid_inputs(self, s, expected):
        assert _parse_ccle_bp(s) == expected

    @pytest.mark.parametrize("s", ["", "nan", None, "garbage", "17::"])
    def test_invalid_inputs_return_none_triple(self, s):
        assert _parse_ccle_bp(s) == (None, None, None)


# ---------------------------------------------------------------------------
# liftover_breakpoint  (network-dependent: skips gracefully if unavailable)
# ---------------------------------------------------------------------------

# pyliftover downloads the chain file on first use (~5MB, cached under
# ~/.pyliftover/). These tests fall back to skip if either pyliftover is
# unavailable or the chain file cannot be fetched (offline CI).

@pytest.fixture(scope="module")
def lo_available():
    """Return True if a hg19→hg38 LiftOver can be instantiated."""
    try:
        from pyliftover import LiftOver  # noqa: F401
    except ImportError:
        return False
    try:
        chrom, pos = liftover_breakpoint("chr17", 7_571_720, "hg19", "hg38")
    except Exception:
        return False
    return chrom is not None and pos is not None


class TestLiftover:
    def test_tp53_start_grch37_to_grch38(self, lo_available):
        """TP53 coding start: hg19 chr17:7,571,720 → hg38 chr17:7,668,402."""
        if not lo_available:
            pytest.skip("pyliftover or chain file unavailable")
        chrom, pos = liftover_breakpoint("chr17", 7_571_720, "hg19", "hg38")
        assert chrom == "chr17"
        # UCSC liftover for this exact position is 7,668,402 (single-base mapping)
        assert pos == 7_668_402

    def test_tp53_end_grch37_to_grch38(self, lo_available):
        """TP53 coding end: hg19 chr17:7,590,868 → hg38 chr17:7,687,550 (matches existing TP53_coding)."""
        if not lo_available:
            pytest.skip("pyliftover or chain file unavailable")
        chrom, pos = liftover_breakpoint("chr17", 7_590_868, "hg19", "hg38")
        assert chrom == "chr17"
        assert pos == 7_687_550

    def test_numeric_chrom_input(self, lo_available):
        """Both 17 (int) and 'chr17' (str) inputs should normalise the same way."""
        if not lo_available:
            pytest.skip("pyliftover or chain file unavailable")
        c1, p1 = liftover_breakpoint("chr17", 7_571_720, "hg19", "hg38")
        c2, p2 = liftover_breakpoint(17, 7_571_720, "hg19", "hg38")
        assert (c1, p1) == (c2, p2)

    def test_unmappable_returns_none_pair(self, lo_available):
        """A clearly out-of-range coordinate must return (None, None), not raise."""
        if not lo_available:
            pytest.skip("pyliftover or chain file unavailable")
        chrom, pos = liftover_breakpoint("chrX", 999_999_999, "hg19", "hg38")
        assert chrom is None and pos is None

    def test_invalid_pos_returns_none_pair(self, lo_available):
        """Non-numeric pos values are caught, not raised."""
        if not lo_available:
            pytest.skip("pyliftover or chain file unavailable")
        chrom, pos = liftover_breakpoint("chr17", "not-a-number", "hg19", "hg38")
        assert chrom is None and pos is None


# ---------------------------------------------------------------------------
# classify_chen_sv_architecture (Chen 2014 PCGP, in-frame GRCh37)
# ---------------------------------------------------------------------------

def _chen_row(chra, posa, chrb, posb):
    """Mimic one row of Chen 2014 Table S4b ('TP53 SV breakpoints' sheet)."""
    return {"chra": chra, "posa": posa, "chrb": chrb, "posb": posb}


class TestClassifyChenSvArchitecture:
    """Six possible outputs + boundary cases on the head-to-head region."""

    def test_empty_input_is_neither(self):
        call, reason = classify_chen_sv_architecture(pd.DataFrame())
        assert call == "Neither"
        assert "no svs" in reason.lower()

    def test_tp53_only_three_intron1_ctxs(self):
        # Mirrors SJOS001 in the actual data: 3 CTXs into TP53 intron 1 (chr17:7,583,xxx)
        df = pd.DataFrame([
            _chen_row(5,  34_477_424, 17, 7_583_762),
            _chen_row(17, 7_583_302,  5,  34_478_526),
            _chen_row(17, 7_583_751,  17, 7_583_463),
        ])
        call, _ = classify_chen_sv_architecture(df)
        assert call == "TP53-only"

    def test_tcab1_only_single_wrap53_break(self):
        df = pd.DataFrame([_chen_row(17, 7_600_000, 5, 100)])
        call, _ = classify_chen_sv_architecture(df)
        assert call == "TCAB1-only"

    def test_both_single_sv_one_endpoint_each_gene(self):
        # one SV, one end in TP53 coding, the other in WRAP53 coding
        df = pd.DataFrame([_chen_row(17, 7_580_000, 17, 7_600_000)])
        call, _ = classify_chen_sv_architecture(df)
        assert call == "Both-single-SV"

    def test_both_independent_svs(self):
        df = pd.DataFrame([
            _chen_row(17, 7_580_000, 5, 100),
            _chen_row(17, 7_600_000, 9, 200),
        ])
        call, _ = classify_chen_sv_architecture(df)
        assert call == "Both-independent-SVs"

    def test_neither_when_off_locus(self):
        # chr17 but well outside the TP53/WRAP53 windows
        df = pd.DataFrame([_chen_row(17, 1_000_000, 5, 100)])
        call, _ = classify_chen_sv_architecture(df)
        assert call == "Neither"

    def test_complex_at_threshold(self):
        df = pd.DataFrame([
            _chen_row(17, 7_580_000 + i * 100, 5, 100)
            for i in range(COMPLEX_THRESHOLD)
        ])
        call, _ = classify_chen_sv_architecture(df)
        assert call == "Complex"

    def test_head_to_head_overlap_counts_as_both(self):
        """
        TP53 (-) and WRAP53 (+) gene bodies overlap in
        chr17:7,589,389-7,590,868 (~1.5 kb head-to-head).
        A single break here legitimately hits both genes.
        """
        h2h_pos = 7_590_000
        assert GRCH37_WRAP53_CODING[0] <= h2h_pos < GRCH37_WRAP53_CODING[1]
        assert GRCH37_TP53_CODING[0]   <= h2h_pos < GRCH37_TP53_CODING[1]

        df = pd.DataFrame([_chen_row(17, h2h_pos, 5, 100)])
        call, _ = classify_chen_sv_architecture(df)
        assert call == "Both-single-SV"

    def test_intron1_window_is_subset_of_tp53(self):
        """Sanity: the Chen-paper intron-1 window must lie inside TP53 coding."""
        assert GRCH37_TP53_INTRON1[0] >= GRCH37_TP53_CODING[0]
        assert GRCH37_TP53_INTRON1[1] <= GRCH37_TP53_CODING[1]

    def test_locus_window_contains_both_genes(self):
        """Sanity: the locus filter window contains both TP53 and WRAP53 coding regions."""
        win_start, win_end = GRCH37_LOCUS_WINDOW
        assert win_start <= GRCH37_TP53_CODING[0]
        assert win_end   >= GRCH37_WRAP53_CODING[1]
