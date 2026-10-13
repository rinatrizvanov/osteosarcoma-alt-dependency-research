"""
Read-only helpers for ``notebooks/06_aim1_rna_inventory.ipynb``.

This module exists to keep the notebook tidy. It never writes to disk,
never modifies the manifest, and never loads a full expression matrix
into memory. All membership checks read only the minimum needed
(index columns or header rows) so the notebook can run on a laptop
without paging.

Two pieces of logic worth knowing about:

* ``EXPRESSION_SOURCES`` — the registry of which expression file we
  treat as canonical for each cohort. ``None`` is meaningful: it means
  "no expression file is currently downloaded for this cohort". Callers
  should never silently substitute an alternative; absent means absent.
* ``decide_aim1_design`` — the decision rule documented in notebook 06.
  Mirrors the language in the notebook spec so the verdict is a function,
  not a free-text judgement.
"""
from __future__ import annotations

import gzip
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import pandas as pd

from .audit import classify_aim1_arms

DATA_ROOT = Path(__file__).parent.parent / "data"
DEPMAP_DIR = DATA_ROOT / "raw" / "depmap"
TARGET_OS_DIR = DATA_ROOT / "raw" / "target_os"

DEPMAP_RNA_PATH = DEPMAP_DIR / "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"
DEPMAP_RNA_DATA_TYPE = "log2(TPM+1) protein-coding"

TARGET_OS_TPM_PATH = TARGET_OS_DIR / "TARGET-OS.star_tpm.tsv.gz"
TARGET_OS_TPM_DATA_TYPE = "STAR TPM (UCSC Xena)"


@dataclass(frozen=True)
class ExpressionSource:
    """Where to find the canonical expression matrix for a cohort."""
    path: Path | None      # ``None`` => no expression file available locally.
    data_type: str         # Free-text label of the units (TPM, log2TPM, counts, ...).
    description: str       # Short prose for the report.


# One entry per cohort that appears in the manifest's ``cohort`` column.
# ``None`` path means "literature labels only; no RNA matrix shipped". This is
# the load-bearing observation for the Aim 1 design verdict — do not silently
# substitute a different file.
EXPRESSION_SOURCES: dict[str, ExpressionSource] = {
    "DepMap_26Q1": ExpressionSource(
        path=DEPMAP_RNA_PATH,
        data_type=DEPMAP_RNA_DATA_TYPE,
        description="DepMap 26Q1 RNA-seq; one row per (Sequencing × Model × Condition).",
    ),
    "Lawlor_2026_GCCRI": ExpressionSource(
        path=None,
        data_type="unknown",
        description="Lawlor 2026 ships SV/IHC labels only; no expression matrix in data/raw/.",
    ),
    "PCGP_2014": ExpressionSource(
        path=None,
        data_type="unknown",
        description="Chen 2014 PCGP cohort ships TP53 SV / ATRX / telomere-FISH tables only; "
                    "no expression matrix in data/raw/.",
    ),
    "Taglialatela_2026": ExpressionSource(
        path=None,
        data_type="unknown",
        description="Taglialatela 2026 ships a cell-line list with ALT calls only; "
                    "no expression matrix in data/raw/.",
    ),
    "de_Nonneville_2022": ExpressionSource(
        path=None,
        data_type="unknown",
        description="de Nonneville 2022 ships a tumor ALT-call table (EV3) only; "
                    "no expression matrix in data/raw/.",
    ),
    # TARGET-OS is downloaded but not yet joined to the manifest. Kept here
    # so a future cohort named "TARGET_OS" automatically routes correctly.
    "TARGET_OS": ExpressionSource(
        path=TARGET_OS_TPM_PATH,
        data_type=TARGET_OS_TPM_DATA_TYPE,
        description="UCSC Xena TARGET-OS STAR TPM matrix (gene × sample, gzip TSV).",
    ),
}


# ---------------------------------------------------------------------------
# Cheap file-backed membership probes (no full-matrix loads)
# ---------------------------------------------------------------------------

def depmap_rna_model_ids(path: Path | None = None) -> set[str]:
    """Return the set of ``ModelID`` values present in the DepMap RNA file.

    Reads only the ``ModelID`` column, never the full expression matrix.
    The file may have multiple rows per ModelID (different sequencing runs);
    we deduplicate so the caller gets a clean membership set.
    """
    p = path or DEPMAP_RNA_PATH
    if not p.exists():
        return set()
    col = pd.read_csv(p, usecols=["ModelID"])["ModelID"].dropna().astype(str)
    return set(col.unique())


def target_os_rna_sample_ids(path: Path | None = None) -> set[str]:
    """Return the TARGET-OS sample IDs that appear as columns in the TPM TSV.

    Reads only the gzip header line, never the matrix body.
    """
    p = path or TARGET_OS_TPM_PATH
    if not p.exists():
        return set()
    with gzip.open(p, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
    # First column is the gene ID column; the rest are sample IDs.
    return set(header[1:])


# ---------------------------------------------------------------------------
# Per-cohort sample-level inventory
# ---------------------------------------------------------------------------

def lookup_file_backed_rna(
    sample_id: str,
    cohort: str,
    *,
    depmap_models: set[str],
    target_os_samples: set[str],
) -> tuple[bool, str, str]:
    """Decide whether one manifest row is backed by an actual on-disk RNA matrix.

    Returns ``(file_backed_rna, expression_source, data_type)``.

    Membership is decided by:
      * DepMap_26Q1 → ``sample_id`` ∈ DepMap RNA ``ModelID`` column.
      * TARGET_OS   → ``sample_id`` ∈ TARGET-OS TPM header columns.
      * Anything else → ``False`` (we have no expression file for those cohorts).
    """
    src = EXPRESSION_SOURCES.get(cohort)
    if src is None or src.path is None:
        return (
            False,
            "none (literature labels only)" if src is None else src.description,
            src.data_type if src else "unknown",
        )

    if cohort == "DepMap_26Q1":
        present = sample_id in depmap_models
        return present, src.path.name, src.data_type
    if cohort == "TARGET_OS":
        present = sample_id in target_os_samples
        return present, src.path.name, src.data_type

    # Defensive default: a registered cohort with a path but no membership
    # check above. Treat as not-file-backed until a probe is added.
    return False, src.path.name, src.data_type


def build_inventory(
    contrast_manifest: pd.DataFrame,
    *,
    depmap_models: set[str] | None = None,
    target_os_samples: set[str] | None = None,
) -> pd.DataFrame:
    """Produce the per-sample inventory used by notebook 06.

    Parameters
    ----------
    contrast_manifest : DataFrame
        Manifest rows already restricted to the Aim 1 primary contrast,
        with an ``arm`` column (``'disrupted'`` / ``'intact'``).
    depmap_models, target_os_samples : pre-fetched membership sets.
        If ``None``, this function reads them itself.

    Returns
    -------
    Tidy DataFrame with one row per contrast sample and the columns
    documented in the notebook spec.
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
            "sample_id":             sid,
            "cohort":                cohort,
            "tcab1_overall_status":  r["arm"],
            "sample_type":           sample_type,
            "manifest_has_rna":      bool(manifest_has_rna),
            "file_backed_rna":       bool(file_backed),
            "expression_source":     source,
            "data_type":             data_type,
        })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Confounding assessment + design verdict
# ---------------------------------------------------------------------------

def assess_confounding(inventory: pd.DataFrame) -> dict:
    """Compute the structured confounding assessment described in notebook 06.

    Restricted to rows where ``file_backed_rna == True``. Returns a dict with:
      * ``n_total``, ``n_disrupted``, ``n_intact`` (file-backed)
      * ``cohorts_per_arm`` — {arm: list of cohorts contributing >=1 sample}
      * ``sample_types_per_arm``, ``data_types_per_arm`` — same shape
      * ``cohort_x_arm`` — DataFrame of file-backed counts per cohort × arm
      * ``cohorts_with_both_arms`` — list of cohort names with >=1 of each arm
      * ``cohorts_with_both_arms_same_dtype`` — list of (cohort, data_type) tuples
        where the same cohort has at least one disrupted AND one intact sample
        annotated under the same ``data_type``.
      * ``has_within_cohort_co_occurrence`` (bool)
      * ``has_within_cohort_co_occurrence_same_dtype`` (bool)
    """
    fb = inventory[inventory["file_backed_rna"].astype(bool)].copy() if len(inventory) else inventory

    cohort_x_arm = (
        pd.crosstab(fb["cohort"], fb["tcab1_overall_status"])
          if len(fb) else pd.DataFrame()
    )

    cohorts_per_arm: dict[str, list[str]] = {"disrupted": [], "intact": []}
    sample_types_per_arm: dict[str, list[str]] = {"disrupted": [], "intact": []}
    data_types_per_arm: dict[str, list[str]] = {"disrupted": [], "intact": []}
    if len(fb):
        for arm in ("disrupted", "intact"):
            sub = fb[fb["tcab1_overall_status"] == arm]
            cohorts_per_arm[arm]      = sorted(sub["cohort"].unique().tolist())
            sample_types_per_arm[arm] = sorted(sub["sample_type"].unique().tolist())
            data_types_per_arm[arm]   = sorted(sub["data_type"].unique().tolist())

    cohorts_with_both_arms: list[str] = []
    if not cohort_x_arm.empty:
        for cohort in cohort_x_arm.index:
            arms_present = {a for a in ("disrupted", "intact")
                            if a in cohort_x_arm.columns
                            and cohort_x_arm.at[cohort, a] >= 1}
            if {"disrupted", "intact"}.issubset(arms_present):
                cohorts_with_both_arms.append(cohort)

    cohorts_with_both_arms_same_dtype: list[tuple[str, str]] = []
    for cohort in cohorts_with_both_arms:
        sub = fb[fb["cohort"] == cohort]
        per_dtype = pd.crosstab(sub["data_type"], sub["tcab1_overall_status"])
        for dt in per_dtype.index:
            has_dis = per_dtype.at[dt, "disrupted"] >= 1 if "disrupted" in per_dtype.columns else False
            has_int = per_dtype.at[dt, "intact"]    >= 1 if "intact"    in per_dtype.columns else False
            if has_dis and has_int:
                cohorts_with_both_arms_same_dtype.append((cohort, dt))

    return {
        "n_total":     len(fb),
        "n_disrupted": int((fb["tcab1_overall_status"] == "disrupted").sum()),
        "n_intact":    int((fb["tcab1_overall_status"] == "intact").sum()),
        "cohort_x_arm":       cohort_x_arm,
        "cohorts_per_arm":      cohorts_per_arm,
        "sample_types_per_arm": sample_types_per_arm,
        "data_types_per_arm":   data_types_per_arm,
        "cohorts_with_both_arms": cohorts_with_both_arms,
        "cohorts_with_both_arms_same_dtype": cohorts_with_both_arms_same_dtype,
        "has_within_cohort_co_occurrence": bool(cohorts_with_both_arms),
        "has_within_cohort_co_occurrence_same_dtype": bool(cohorts_with_both_arms_same_dtype),
    }


# Threshold for the "within-cohort viable" verdict, mirrored from the
# notebook 06 spec. Pinned here so the test suite can read it.
WITHIN_COHORT_PER_ARM_MIN = 3


def _within_cohort_qualifies(
    inventory: pd.DataFrame,
    cohort: str,
    data_type: str,
    *,
    per_arm_min: int = WITHIN_COHORT_PER_ARM_MIN,
) -> tuple[bool, int, int]:
    """Helper: does (cohort × data_type) carry ≥ per_arm_min samples per arm?"""
    sub = inventory[
        inventory["file_backed_rna"]
        & (inventory["cohort"] == cohort)
        & (inventory["data_type"] == data_type)
    ]
    n_dis = int((sub["tcab1_overall_status"] == "disrupted").sum())
    n_int = int((sub["tcab1_overall_status"] == "intact").sum())
    return (n_dis >= per_arm_min and n_int >= per_arm_min), n_dis, n_int


def decide_aim1_design(
    inventory: pd.DataFrame,
    assessment: dict,
    *,
    per_arm_min: int = WITHIN_COHORT_PER_ARM_MIN,
) -> dict:
    """Apply the three-way decision rule and return a structured verdict.

    Categories
    ----------
    * ``within_cohort_viable`` — ≥1 (cohort × data_type) pair has both arms
      with ≥ ``per_arm_min`` samples each. Within-cohort DE is the recommended
      primary analysis; if multiple pairs qualify, meta-analyze across them.
    * ``pooled_exploratory_only`` — no single cohort has both arms with
      file-backed RNA, but **each** arm has ≥1 file-backed sample from ≥2
      different cohorts (so cohort and TCAB1 status are not perfectly
      confounded). A pooled model with cohort as a covariate is feasible
      as an exploratory analysis, with caveats.
    * ``not_estimable`` — anything else: one arm has zero file-backed RNA,
      or both arms come from disjoint single cohorts. The Aim 1 primary
      contrast as currently defined cannot be estimated from the available
      data and must be redefined or restricted.

    Returns
    -------
    dict with keys ``category``, ``justification``, and the supporting
    counts. The justification is plain prose suitable for the report.
    """
    qualifying: list[tuple[str, str, int, int]] = []
    for cohort, data_type in assessment["cohorts_with_both_arms_same_dtype"]:
        ok, n_dis, n_int = _within_cohort_qualifies(
            inventory, cohort, data_type, per_arm_min=per_arm_min
        )
        if ok:
            qualifying.append((cohort, data_type, n_dis, n_int))

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
            "primary_n_disrupted": primary[2],
            "primary_n_intact": primary[3],
            "qualifying_pairs": qualifying,
            "justification": (
                f"Cohort '{primary[0]}' provides {primary[2]} disrupted and {primary[3]} "
                f"intact file-backed RNA samples under the same data type "
                f"('{primary[1]}'), meeting the within-cohort minimum of "
                f"{per_arm_min} per arm. Recommend running the Aim 1 primary "
                f"contrast within this cohort first to avoid platform/cohort "
                f"confounding.{rest_note}"
            ),
        }

    n_disrupted_cohorts = len(assessment["cohorts_per_arm"]["disrupted"])
    n_intact_cohorts    = len(assessment["cohorts_per_arm"]["intact"])

    if n_disrupted_cohorts >= 2 and n_intact_cohorts >= 2:
        return {
            "category": "pooled_exploratory_only",
            "n_disrupted_cohorts": n_disrupted_cohorts,
            "n_intact_cohorts": n_intact_cohorts,
            "justification": (
                f"No single cohort carries both arms with file-backed RNA, but "
                f"the disrupted arm spans {n_disrupted_cohorts} cohorts "
                f"({', '.join(assessment['cohorts_per_arm']['disrupted'])}) and "
                f"the intact arm spans {n_intact_cohorts} cohorts "
                f"({', '.join(assessment['cohorts_per_arm']['intact'])}). "
                f"TCAB1 status is therefore only partially confounded with cohort, "
                f"so a pooled model with cohort as a fixed effect is feasible as an "
                f"exploratory analysis. Within-cohort meta-analysis is not possible "
                f"under the current data."
            ),
        }

    fb_dis = assessment["n_disrupted"]
    fb_int = assessment["n_intact"]
    cohorts_dis = assessment["cohorts_per_arm"]["disrupted"] or ["(none)"]
    cohorts_int = assessment["cohorts_per_arm"]["intact"]    or ["(none)"]

    return {
        "category": "not_estimable",
        "fb_n_disrupted": fb_dis,
        "fb_n_intact": fb_int,
        "fb_cohorts_disrupted": cohorts_dis,
        "fb_cohorts_intact": cohorts_int,
        "justification": (
            f"File-backed RNA is available for {fb_dis} disrupted sample(s) "
            f"(from cohort(s) {', '.join(cohorts_dis)}) and {fb_int} intact "
            f"sample(s) (from cohort(s) {', '.join(cohorts_int)}). TCAB1 status "
            f"is perfectly or near-perfectly confounded with cohort/platform/"
            f"sample_type, so the Aim 1 primary contrast as currently defined "
            f"is not statistically estimable. The contrast must be redefined "
            f"(e.g. require RNA-backed cohorts on both sides) or restricted "
            f"(e.g. drop arms without file-backed RNA) before any DE analysis "
            f"can proceed."
        ),
    }


# ---------------------------------------------------------------------------
# Convenience: pragmatic primary-contrast subsetter
# ---------------------------------------------------------------------------

def select_pragmatic_contrast_rows(manifest: pd.DataFrame) -> pd.DataFrame:
    """Return manifest rows assigned to the pragmatic Aim 1 primary contrast.

    Adds an ``arm`` column with values ``'disrupted'`` or ``'intact'``.
    Rows assigned to ``'uncertain'`` or ``'excluded'`` by
    :func:`src.audit.classify_aim1_arms` are dropped. The canonical
    arm-classification logic lives there (see ``label_spec_v0.9 §1c``);
    this helper is a thin presentation wrapper so notebook 06 cannot drift
    from notebook 05.
    """
    arms = classify_aim1_arms(manifest, "pragmatic")
    sub = manifest.copy()
    sub["arm"] = arms.values
    return sub[sub["arm"].isin(("disrupted", "intact"))].reset_index(drop=True)
