"""Figures: per-gene strip plots and the effect-size forest plot (Section 5).

Deliberately plain matplotlib — no seaborn — to keep the dependency set close
to Aim 2's and the output deterministic.
"""
from __future__ import annotations

from pathlib import Path
from typing import Sequence

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from .stratify import INTACT_LABEL, LOW_LABEL  # noqa: E402
from .utils import FIGURES_DIR, RANDOM_SEED  # noqa: E402

LOW_COLOR = "#B2182B"
INTACT_COLOR = "#2166AC"


def _jitter(n: int, width: float = 0.09, seed: int = RANDOM_SEED) -> np.ndarray:
    return np.random.default_rng(seed).uniform(-width, width, n)


def strip_plot(
    effects: pd.DataFrame,
    strata: pd.Series,
    genes: Sequence[str],
    *,
    title: str = "CRISPR gene effect by TCAB1 stratum (ALT+ lines)",
    outpath: Path | None = None,
) -> Path:
    """One panel per gene: TCAB1-low vs TCAB1-intact, points with group medians."""
    fig, axes = plt.subplots(1, len(genes), figsize=(2.05 * len(genes), 3.9), sharey=True)
    axes = np.atleast_1d(axes)
    for ax, gene in zip(axes, genes):
        for i, (label, color) in enumerate(
            ((LOW_LABEL, LOW_COLOR), (INTACT_LABEL, INTACT_COLOR))
        ):
            ids = strata.index[strata == label]
            vals = effects.loc[ids, gene].to_numpy(dtype=float)
            ax.scatter(
                i + _jitter(len(vals)), vals, s=26, color=color,
                alpha=0.85, edgecolor="white", linewidth=0.5, zorder=3,
            )
            if len(vals):
                ax.hlines(np.median(vals), i - 0.26, i + 0.26,
                          color=color, linewidth=2.2, zorder=4)
        ax.axhline(0.0, color="#999999", linewidth=0.7, linestyle=":", zorder=1)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["low", "intact"], fontsize=8)
        ax.set_xlim(-0.55, 1.55)
        ax.set_title(gene, fontsize=10, fontweight="bold")
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel("Chronos gene effect\n(more negative = more essential)", fontsize=9)
    fig.suptitle(title, fontsize=11, y=1.0)
    fig.tight_layout()
    out = Path(outpath or FIGURES_DIR / "per_gene_strip.png")
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return out


def forest_plot(
    per_gene: pd.DataFrame,
    *,
    primary: dict | None = None,
    title: str = "Cliff's delta, TCAB1-low vs TCAB1-intact (ALT+ lines)",
    outpath: Path | None = None,
) -> Path:
    """Effect sizes with bootstrap CIs; negative = predicted direction."""
    rows = per_gene.copy()
    order = list(rows["name"])
    if primary is not None:
        rows = pd.concat([pd.DataFrame([{**primary, "gene_set": "COMPOSITE"}]), rows],
                         ignore_index=True)
        order = ["composite"] + order

    fig, ax = plt.subplots(figsize=(6.4, 0.52 * len(rows) + 1.7))
    colors = {"BIR": LOW_COLOR, "HR": INTACT_COLOR, "COMPOSITE": "#333333"}
    for i, (_, r) in enumerate(rows.iloc[::-1].iterrows()):
        color = colors.get(r["gene_set"], "#333333")
        ax.plot([r["delta_ci_low"], r["delta_ci_high"]], [i, i],
                color=color, linewidth=1.6, alpha=0.75, zorder=2)
        ax.scatter(r["cliffs_delta"], i, s=62 if r["gene_set"] == "COMPOSITE" else 42,
                   color=color, zorder=3, marker="D" if r["gene_set"] == "COMPOSITE" else "o")
    ax.axvline(0.0, color="#999999", linewidth=0.9, linestyle="--", zorder=1)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels(list(rows["name"])[::-1], fontsize=9)
    ax.set_xlabel("Cliff's delta (negative = BIR-shifted in TCAB1-low, as predicted)", fontsize=9)
    ax.set_xlim(-1.05, 1.05)
    ax.set_title(title, fontsize=11)
    ax.spines[["top", "right"]].set_visible(False)

    handles = [
        plt.Line2D([], [], color=colors["BIR"], marker="o", linestyle="", label="BIR module"),
        plt.Line2D([], [], color=colors["HR"], marker="o", linestyle="", label="HR comparator"),
    ]
    if primary is not None:
        handles.insert(0, plt.Line2D([], [], color=colors["COMPOSITE"], marker="D",
                                     linestyle="", label="composite (primary)"))
    ax.legend(handles=handles, fontsize=8, frameon=False, loc="lower right")
    fig.tight_layout()
    out = Path(outpath or FIGURES_DIR / "effect_size_forest.png")
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return out


def null_distribution_plot(
    null_rhos: pd.Series, observed: float, *, outpath: Path | None = None
) -> Path:
    """Notebook 02: observed rho(WRAP53, TP53) against the expression-matched null."""
    fig, ax = plt.subplots(figsize=(6.0, 3.4))
    ax.hist(null_rhos.to_numpy(), bins=40, color="#BBBBBB", edgecolor="white")
    ax.axvline(observed, color=LOW_COLOR, linewidth=2.0,
               label=f"observed rho(WRAP53, TP53) = {observed:+.3f}")
    ax.axvline(float(np.percentile(null_rhos, 99)), color="#444444", linestyle="--",
               linewidth=1.2, label="null 99th percentile")
    ax.set_xlabel("Spearman rho with WRAP53, expression-matched genes", fontsize=9)
    ax.set_ylabel("genes", fontsize=9)
    ax.set_title("Is WRAP53/TP53 correlation anomalous?", fontsize=11)
    ax.legend(fontsize=8, frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    out = Path(outpath or FIGURES_DIR / "strand_check_null.png")
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return out
