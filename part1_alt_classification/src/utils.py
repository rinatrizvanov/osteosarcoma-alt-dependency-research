"""Shared utilities: checksums, version pinning, provenance recording, liftover."""

from __future__ import annotations

import hashlib
import importlib
import json
import os
import time
from pathlib import Path
from typing import Tuple


def sha256_file(path: str | Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def file_mtime_iso(path: str | Path) -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(os.path.getmtime(path)))


def library_versions(names: list[str]) -> dict[str, str]:
    versions = {}
    for name in names:
        try:
            mod = importlib.import_module(name)
            versions[name] = getattr(mod, "__version__", "unknown")
        except ImportError:
            versions[name] = "not_installed"
    return versions


def log_versions(names: list[str] | None = None) -> dict[str, str]:
    default = ["pandas", "numpy", "scipy", "pyranges", "openpyxl", "requests", "GEOparse"]
    vers = library_versions(names or default)
    print("Library versions:")
    for k, v in vers.items():
        print(f"  {k}: {v}")
    return vers


def verify_file(path: str | Path, expected_sha256: str | None = None) -> dict:
    path = Path(path)
    if not path.exists():
        return {"exists": False, "path": str(path)}
    actual = sha256_file(path)
    result = {
        "exists": True,
        "path": str(path),
        "sha256": actual,
        "mtime": file_mtime_iso(path),
        "size_bytes": path.stat().st_size,
    }
    if expected_sha256 is not None:
        result["checksum_match"] = actual == expected_sha256
    return result


def record_provenance(records: list[dict], provenance_path: str | Path) -> None:
    """Append provenance records to a JSON-lines sidecar file."""
    provenance_path = Path(provenance_path)
    provenance_path.parent.mkdir(parents=True, exist_ok=True)
    with open(provenance_path, "a") as f:
        for rec in records:
            f.write(json.dumps(rec) + "\n")


MANIFEST_COLUMNS = [
    "sample_id",
    "cohort",
    "sample_type",
    "has_rna",
    "has_wes",
    "has_wgs",
    "has_sc",
    "has_clinical",
    "tp53_status",
    "tp53_evidence",
    "tcab1_rna_quantile",
    "tcab1_rna_value",
    "tcab1_rna_quantile_within_cohort",
    "tcab1_genomic_status",
    "tcab1_genomic_evidence",
    "tcab1_functional_status",
    "tcab1_functional_evidence",
    "tcab1_overall_status",
    "tcab1_confidence",
    "sv_architecture_17p13_1",
    "sv_architecture_evidence",
    "alt_status",
    "alt_source",
    "source_confidence",
    "alt_assay_evidence",
    "telomerase_status",
    "tert_expression",
    "telomerase_evidence",
    "atrx_status",
    "daxx_status",
    "smarcal1_status",
    "top3a_status",
    "purity_estimate",
    "purity_quartile",
    "include_aim1_primary_strict",
    "include_aim1_primary_pragmatic",
    "include_aim1_primary_broad",
    "include_aim2_strict",
    "include_aim2_pragmatic",
    "include_aim2_broad",
    "include_tp53_sensitivity",
    "include_telomerase_negative_secondary",
    "notes",
]

GENE_ENTREZ = {
    "WRAP53": 55135,
    "TCAB1": 55135,
    "TP53": 7157,
    "TERT": 7015,
    "ATRX": 546,
    "DAXX": 1616,
    "SMARCAL1": 50485,
    "TOP3A": 7156,
}

LOCUS_WINDOW = ("chr17", 7_560_000, 7_710_000)


# ---------------------------------------------------------------------------
# Coordinate liftover (GRCh37/hg19 ↔ GRCh38/hg38)
# ---------------------------------------------------------------------------

def liftover_breakpoint(
    chrom: str | int,
    pos: int,
    source: str = "hg19",
    target: str = "hg38",
) -> Tuple[str | None, int | None]:
    """
    Lift a single chromosomal coordinate between assemblies.

    Returns
    -------
    (chrom, pos)
        Lifted ('chr17', 7_668_402) on success.
        (None, None) if the input coordinate cannot be lifted (deleted region,
        ambiguous mapping, or chain file lookup failure).

    Notes
    -----
    * Uses pyliftover. The first call downloads and caches the UCSC chain file
      under ``~/.pyliftover/`` (~5 MB for hg19↔hg38). Subsequent calls reuse
      the same in-memory ``LiftOver`` instance, attached as a function attribute.
    * Input ``chrom`` is normalised to ``"chr<N>"`` form. Both ``17`` and
      ``"chr17"`` are accepted.
    * Liftover is 0-based half-open in pyliftover, identical to BED conventions.
    * Failures are returned, not raised, so callers can decide policy. Callers
      MUST log dropped breakpoints — silent loss is a project anti-pattern.
    """
    from pyliftover import LiftOver

    cache_key = (source, target)
    cache = getattr(liftover_breakpoint, "_cache", None)
    if cache is None:
        cache = {}
        setattr(liftover_breakpoint, "_cache", cache)
    if cache_key not in cache:
        cache[cache_key] = LiftOver(source, target)
    lo = cache[cache_key]

    chrom_str = str(chrom).strip()
    if not chrom_str.startswith("chr"):
        chrom_str = f"chr{chrom_str}"

    try:
        result = lo.convert_coordinate(chrom_str, int(pos))
    except (TypeError, ValueError):
        return None, None

    if not result:
        return None, None
    new_chrom, new_pos = result[0][0], result[0][1]
    return new_chrom, int(new_pos)
