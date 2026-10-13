"""Notebook 02: is the WRAP53 expression measurement trustworthy?

The risk. WRAP53 exon 1 carries a ~1,419 bp head-to-head antisense overlap with
TP53. On an **unstranded** library, reads falling in that interval cannot be
assigned to the correct gene, so WRAP53 TPM is contaminated by TP53 — which
would invalidate the Stage 1 stratifier outright.

Why this runs before the prespec is frozen. WRAP53 expression is the
**exposure**, not an outcome. Validating how it was measured is a measurement
question, answerable without touching a single gene-effect value. This module
does not import :mod:`src.dependency_stats`, so notebook 02 structurally cannot
read ``CRISPRGeneEffect.csv``.

Three empirical diagnostics, run regardless of what the documentation says:

1. Spearman rho(WRAP53, TP53) across all lines and within the osteosarcoma
   lines. Contamination inflates this correlation.
2. The same rho against a **null distribution of expression-matched genes**.
   WRAP53 and TP53 are both DNA-damage-associated, so some genuine
   co-regulation is expected; the question is whether the correlation exceeds
   what comparable gene pairs show.
3. Does WRAP53 expression track WRAP53 **copy number** better than it tracks
   TP53 expression? Contamination predicts the reverse, and copy number is
   strand-agnostic.

No single diagnostic is decisive, which is why :func:`strand_decision` reports
all three and the branch they jointly imply rather than thresholding one.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from .cohort import EXPRESSION_META_COLS
from .stratify import (
    TCAB1_GENE,
    TP53_GENE,
    load_expression,
    load_wrap53_cn,
    osteosarcoma_rna_models,
)
from .utils import RANDOM_SEED, RNA_PATH

#: Candidate genes sampled before expression-matching down to N_NULL_GENES.
N_CANDIDATE_GENES = 3000
N_NULL_GENES = 1000

#: rho(WRAP53, TP53) above this percentile of the matched null is treated as
#: evidence of contamination rather than co-regulation.
NULL_PERCENTILE_FLAG = 99.0


@dataclass
class StrandDiagnostics:
    """Everything notebook 02 measured, plus the branch it implies."""

    n_all_lines: int
    n_os_lines: int
    rho_wrap53_tp53_all: float
    p_wrap53_tp53_all: float
    rho_wrap53_tp53_os: float
    p_wrap53_tp53_os: float
    null_n_genes: int
    null_rho_median: float
    null_rho_p95: float
    null_rho_p99: float
    wrap53_percentile_in_null: float
    rho_wrap53_expr_vs_cn_os: float
    p_wrap53_expr_vs_cn_os: float
    exceeds_null: bool
    cn_tracks_better_than_tp53: bool
    branch: str
    rationale: list[str] = field(default_factory=list)

    def as_dict(self) -> dict:
        return asdict(self)


def sample_null_gene_columns(
    *, n_candidates: int = N_CANDIDATE_GENES, seed: int = RANDOM_SEED
) -> list[str]:
    """Randomly sample gene columns from the expression header."""
    with open(RNA_PATH) as fh:
        header = next(fh).rstrip("\n").split(",")
    genes = [c for c in header if c not in EXPRESSION_META_COLS]
    rng = np.random.default_rng(seed)
    n = min(n_candidates, len(genes))
    idx = rng.choice(len(genes), size=n, replace=False)
    return [genes[i] for i in sorted(idx)]


def expression_matched_null(
    wrap53: pd.Series,
    candidates: pd.DataFrame,
    *,
    n_null: int = N_NULL_GENES,
) -> pd.Series:
    """Spearman rho(WRAP53, g) for the ``n_null`` genes closest to WRAP53 in mean expression.

    Matching on mean expression matters: rho is attenuated for low-expressed,
    noisy genes, so an unmatched null would be biased downward and make any
    correlation look anomalous.
    """
    target_mean = float(wrap53.mean())
    means = candidates.mean(axis=0)
    chosen = (means - target_mean).abs().sort_values().head(n_null).index
    rhos = {}
    for gene in chosen:
        values = candidates[gene]
        if values.std(ddof=1) == 0:
            continue
        rhos[gene] = float(stats.spearmanr(wrap53.to_numpy(), values.to_numpy()).statistic)
    return pd.Series(rhos, name="null_rho").dropna()


def strand_decision(
    *, seed: int = RANDOM_SEED, n_candidates: int = N_CANDIDATE_GENES
) -> tuple[StrandDiagnostics, pd.Series]:
    """Run all three diagnostics and return them with the implied branch.

    Returns
    -------
    (diagnostics, null_rhos)
        ``null_rhos`` is kept so the notebook can plot the null distribution
        with the observed rho marked on it.
    """
    pair = load_expression(None, [TCAB1_GENE, TP53_GENE])
    os_ids = [m for m in osteosarcoma_rna_models() if m in pair.index]

    rho_all = stats.spearmanr(pair[TCAB1_GENE], pair[TP53_GENE])
    os_pair = pair.loc[os_ids]
    rho_os = stats.spearmanr(os_pair[TCAB1_GENE], os_pair[TP53_GENE])

    candidate_cols = sample_null_gene_columns(n_candidates=n_candidates, seed=seed)
    candidates = load_expression(
        None, [c.split(" (")[0].strip() for c in candidate_cols], verify=False
    )
    candidates = candidates.loc[pair.index]
    null_rhos = expression_matched_null(pair[TCAB1_GENE], candidates)

    observed = float(rho_all.statistic)
    percentile = float((null_rhos < observed).mean() * 100.0)

    cn = load_wrap53_cn(os_ids)
    valid = cn.notna()
    if valid.sum() >= 4:
        cn_rho = stats.spearmanr(
            os_pair.loc[valid.index[valid], TCAB1_GENE], cn[valid]
        )
        rho_cn, p_cn = float(cn_rho.statistic), float(cn_rho.pvalue)
    else:
        rho_cn, p_cn = float("nan"), float("nan")

    exceeds = percentile >= NULL_PERCENTILE_FLAG
    cn_better = bool(not np.isnan(rho_cn) and abs(rho_cn) > abs(float(rho_os.statistic)))

    rationale = [
        f"rho(WRAP53, TP53) across {len(pair)} lines = {observed:+.4f} "
        f"(p = {rho_all.pvalue:.3g})",
        f"rho(WRAP53, TP53) within {len(os_ids)} osteosarcoma lines = "
        f"{rho_os.statistic:+.4f} (p = {rho_os.pvalue:.3g})",
        f"Expression-matched null over {len(null_rhos)} genes: median "
        f"{null_rhos.median():+.4f}, p95 {np.percentile(null_rhos, 95):+.4f}, "
        f"p99 {np.percentile(null_rhos, 99):+.4f}",
        f"Observed rho sits at the {percentile:.2f}th percentile of that null",
        f"rho(WRAP53 expr, WRAP53 CN) in osteosarcoma lines = {rho_cn:+.4f} "
        f"(p = {p_cn:.3g})",
    ]

    if not exceeds:
        branch = "direct"
        rationale.append(
            f"DIRECT: rho(WRAP53, TP53) is within what expression-matched gene "
            f"pairs show (below the {NULL_PERCENTILE_FLAG:.0f}th percentile). No "
            f"evidence of antisense cross-assignment. Use WRAP53 log2(TPM+1) "
            f"directly as the Section 4 stratifier."
        )
    elif cn_better:
        branch = "direct_cn_corroborated"
        rationale.append(
            "DIRECT (CN-corroborated): rho(WRAP53, TP53) exceeds the null, but "
            "WRAP53 expression tracks WRAP53 copy number more strongly than it "
            "tracks TP53 expression, which is the opposite of what read "
            "cross-assignment would produce. Proceed with direct quantification "
            "and report the CN concordance alongside it."
        )
    else:
        branch = "fallback_cn_conjunction"
        rationale.append(
            "FALLBACK: rho(WRAP53, TP53) exceeds the expression-matched null and "
            "WRAP53 expression does not track its own copy number more strongly. "
            "Exon-aware requantification is impossible on public DepMap data (no "
            "public FASTQ/BAM), so switch to the CN-conjunction stratifier and "
            "run the TP53-residualised split as a sensitivity analysis."
        )

    diagnostics = StrandDiagnostics(
        n_all_lines=len(pair),
        n_os_lines=len(os_ids),
        rho_wrap53_tp53_all=observed,
        p_wrap53_tp53_all=float(rho_all.pvalue),
        rho_wrap53_tp53_os=float(rho_os.statistic),
        p_wrap53_tp53_os=float(rho_os.pvalue),
        null_n_genes=len(null_rhos),
        null_rho_median=float(null_rhos.median()),
        null_rho_p95=float(np.percentile(null_rhos, 95)),
        null_rho_p99=float(np.percentile(null_rhos, 99)),
        wrap53_percentile_in_null=percentile,
        rho_wrap53_expr_vs_cn_os=rho_cn,
        p_wrap53_expr_vs_cn_os=p_cn,
        exceeds_null=exceeds,
        cn_tracks_better_than_tp53=cn_better,
        branch=branch,
        rationale=rationale,
    )
    return diagnostics, null_rhos


def residualize_on_tp53(expression: pd.DataFrame) -> pd.Series:
    """WRAP53 residualised on TP53 (the fallback-ladder sensitivity).

    Removes the shared component, and with it any genuine co-regulation — which
    is why this is a sensitivity analysis and never the primary stratifier.
    """
    x = expression[TP53_GENE].to_numpy(dtype=float)
    y = expression[TCAB1_GENE].to_numpy(dtype=float)
    slope, intercept = np.polyfit(x, y, 1)
    resid = y - (slope * x + intercept)
    return pd.Series(resid, index=expression.index, name="WRAP53_resid_on_TP53")
