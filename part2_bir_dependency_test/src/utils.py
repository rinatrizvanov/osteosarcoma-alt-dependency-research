"""Shared utilities: paths, checksums, version pinning, provenance, run reports.

Ported from the Aim 2 project's ``src/utils.py`` (frozen copy at
``references/aim2_utils.py``) and extended with:

* the Stage 1 path constants and their expected sha256 digests, so every
  notebook verifies its inputs instead of trusting the filename;
* :func:`write_report` / :func:`write_sha256_manifest`, which standardise the
  ``notebook_NN_report.md`` + ``notebook_NN_sha256.json`` pair that makes each
  notebook idempotent and auditable (Aim 2 convention).

Nothing in this module reads ``CRISPRGeneEffect.csv``. The outcome variable is
first touched in notebook 05, after the pre-specification is frozen.
"""
from __future__ import annotations

import hashlib
import importlib
import json
import os
import time
from pathlib import Path

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------


def find_root(start: Path | None = None) -> Path:
    """Locate the project root by walking up for a directory containing ``src/``.

    Mirrors the Aim 2 notebook header idiom but works from any depth, so the
    same call succeeds whether a notebook's cwd is ``notebooks/`` or the root.
    """
    here = (start or Path.cwd()).resolve()
    for candidate in (here, *here.parents):
        if (candidate / "src").is_dir() and (candidate / "pyproject.toml").is_file():
            return candidate
    raise RuntimeError(f"Cannot locate project root (no src/ + pyproject.toml) from {here}")


ROOT = find_root(Path(__file__).parent)

REFERENCES_DIR = ROOT / "references"
RAW_DEPMAP_DIR = ROOT / "data" / "raw" / "depmap"
PROCESSED_DIR = ROOT / "data" / "processed"
FIGURES_DIR = ROOT / "figures"
DOCS_DIR = ROOT / "docs"

#: Frozen label input. Stage 1 never re-derives ALT labels (prespec Section 3).
MANIFEST_PATH = REFERENCES_DIR / "master_manifest.csv"

CRISPR_PATH = RAW_DEPMAP_DIR / "CRISPRGeneEffect.csv"
RNA_PATH = RAW_DEPMAP_DIR / "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"
MODEL_PATH = RAW_DEPMAP_DIR / "Model.csv"
CN_GENE_PATH = RAW_DEPMAP_DIR / "OmicsCNGeneWGS.csv"
CN_SEGMENTS_PATH = RAW_DEPMAP_DIR / "OmicsCNSegmentsWGS.csv"

DEPMAP_RELEASE = "26Q1"

#: sha256 recorded when the 26Q1 files were copied from the Aim 2 project
#: (2026-10-11) and verified byte-identical to the originals. Notebook 01
#: re-verifies against these; a mismatch means the release changed under us.
EXPECTED_SHA256: dict[str, str] = {
    "CRISPRGeneEffect.csv":
        "e610a4cefb13a82b5b256b47eb08b63ff14843f8dbd0fb164bc0a32688e5b89e",
    "Model.csv":
        "ea4e0b2a3bc806f81df62689a5ae75f1a100135727a3d7b8a4c7ccc8815183f8",
    "OmicsCNGeneWGS.csv":
        "d2fe65b516802fb2b355b953fcc43dcb908840f275e9166c0dbd7defadca8a1d",
    "OmicsCNSegmentsWGS.csv":
        "690e861b3a35db89b6a6b1122edd842c1adb090ec52ec1b13b0ee91d369a142f",
    "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv":
        "0377be80c525fde98cbd2c6e8b06bdf2a4014a9683eb70182c1f8649d711021a",
    "master_manifest.csv":
        "d51c5f180391ffbce9fd4d5d66ed449532a796859e2bd186a73049c8cbc68b31",
}

#: 17p13.1 locus window (GRCh38), reused verbatim from Aim 2 ``utils.LOCUS_WINDOW``.
LOCUS_WINDOW: tuple[str, int, int] = ("chr17", 7_560_000, 7_710_000)

RANDOM_SEED = 42

# ---------------------------------------------------------------------------
# Checksums and versions (ported from Aim 2)
# ---------------------------------------------------------------------------


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
    default = ["pandas", "numpy", "scipy", "sklearn", "joblib", "matplotlib"]
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


def verify_expected(path: str | Path) -> dict:
    """``verify_file`` against :data:`EXPECTED_SHA256`, keyed on basename.

    Raises
    ------
    FileNotFoundError
        If the file is absent.
    ValueError
        If the digest does not match the recorded value — i.e. the input data
        changed. Analysis must not proceed past that silently.
    """
    path = Path(path)
    expected = EXPECTED_SHA256.get(path.name)
    res = verify_file(path, expected)
    if not res["exists"]:
        raise FileNotFoundError(
            f"{path} is missing. See notebook 01 for the DepMap {DEPMAP_RELEASE} "
            f"archive download fallback."
        )
    if expected is not None and not res["checksum_match"]:
        raise ValueError(
            f"sha256 mismatch for {path.name}:\n"
            f"  expected {expected}\n  actual   {res['sha256']}\n"
            f"Input data changed; re-freeze EXPECTED_SHA256 deliberately or restore the file."
        )
    return res


def record_provenance(records: list[dict], provenance_path: str | Path | None = None) -> None:
    """Append provenance records to a JSON-lines sidecar file."""
    provenance_path = Path(provenance_path or (PROCESSED_DIR / "provenance.jsonl"))
    provenance_path.parent.mkdir(parents=True, exist_ok=True)
    with open(provenance_path, "a") as f:
        for rec in records:
            f.write(json.dumps(rec, sort_keys=True) + "\n")


# ---------------------------------------------------------------------------
# Run reports (Aim 2 notebook_NN_report.md + notebook_NN_sha256.json pattern)
# ---------------------------------------------------------------------------


def run_timestamp() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def write_report(notebook: str, lines: list[str], outdir: Path | None = None) -> Path:
    """Write ``notebook_<nb>_report.md``. Overwrites — notebooks are idempotent."""
    outdir = outdir or PROCESSED_DIR
    outdir.mkdir(parents=True, exist_ok=True)
    path = outdir / f"notebook_{notebook}_report.md"
    path.write_text("\n".join(lines).rstrip() + "\n")
    return path


def write_sha256_manifest(
    notebook: str, paths: list[str | Path], outdir: Path | None = None
) -> Path:
    """Write ``notebook_<nb>_sha256.json`` over this notebook's own outputs.

    Re-running an unchanged notebook must reproduce this file byte-for-byte;
    that equality is the idempotency test (see ``tests/`` and the plan's
    verification section).
    """
    outdir = outdir or PROCESSED_DIR
    outdir.mkdir(parents=True, exist_ok=True)
    digests = {Path(p).name: sha256_file(p) for p in sorted(map(str, paths))}
    path = outdir / f"notebook_{notebook}_sha256.json"
    path.write_text(json.dumps(digests, indent=2, sort_keys=True) + "\n")
    return path


def git_head(root: Path | None = None) -> str:
    """Current git HEAD, used to record the prespec freeze commit in reports."""
    import subprocess

    try:
        out = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=str(root or ROOT), capture_output=True, text=True, check=True,
        )
        return out.stdout.strip()
    except Exception:
        return "unknown"
