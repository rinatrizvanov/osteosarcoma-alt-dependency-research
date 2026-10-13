"""Build the two publication figures from saved notebook outputs.

Reads **only** from ``data/processed/``. Performs no statistics: every delta, CI
and p-value plotted here was computed in notebooks 05 and 06 and frozen. Re-running
this script cannot change a result, only redraw it.

    uv run python scripts/make_publication_figures.py

Design notes (print medium, not screen):

* Two-colour categorical pair ``#B2182B`` / ``#2166AC``, validated for CVD
  separation (worst adjacent deltaE 21.1 protan, 29.8 tritan, 28.7 normal vision,
  contrast >= 3:1 against the surface).
* Identity is **never carried by colour alone** - each arm and each gene set also
  gets its own marker shape, so both figures survive greyscale printing and
  photocopying.
* Direct labels are selective: only OS384 and NRHOS4, the two homozygous
  WRAP53-loss lines, are named. Labelling every point would be noise at n=16.
* Axes are recessive; the zero line and the predicted-direction shading carry the
  argument.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from src.utils import FIGURES_DIR, PROCESSED_DIR, RANDOM_SEED  # noqa: E402

LOW, INTACT = "TCAB1-low", "TCAB1-intact"
C_LOW, C_INTACT = "#B2182B", "#2166AC"
C_BIR, C_HR = "#B2182B", "#2166AC"
C_COMPOSITE = "#333333"
INK, MUTED, GRID = "#1a1a1a", "#666666", "#d8d8d8"

#: The two homozygous WRAP53-loss ALT+ lines (prespec section 4d).
HIGHLIGHT = {"ACH-002834": "OS384", "ACH-003185": "NRHOS4"}

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "axes.edgecolor": MUTED,
    "axes.labelcolor": INK,
    "text.color": INK,
    "xtick.color": MUTED,
    "ytick.color": MUTED,
    "axes.linewidth": 0.8,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "savefig.facecolor": "white",
    "figure.facecolor": "white",
})


def _load():
    comp = pd.read_csv(PROCESSED_DIR / "composite_scores.csv", index_col="ModelID")
    per_gene = pd.read_csv(PROCESSED_DIR / "per_gene_results.csv")
    primary = json.loads((PROCESSED_DIR / "primary_result.json").read_text())
    return comp, per_gene, primary


# ---------------------------------------------------------------------------
# Figure 1 - composite by TCAB1 stratum, both contrasts
# ---------------------------------------------------------------------------

def figure1(comp: pd.DataFrame, primary: dict, outpath: Path) -> Path:
    alt = comp[comp["alt"] == "ALT+"]
    panels = [
        ("arm_8v8", "primary_median_8v8", "a", "Primary: median split (8 v 8)"),
        ("arm_3v13", "sensitivity_tertile_3v13", "b", "Sensitivity: tertile (3 v 13)"),
    ]
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 4.4), sharey=True)
    rng = np.random.default_rng(RANDOM_SEED)

    for ax, (arm_col, key, tag, title) in zip(axes, panels):
        res = primary["contrasts"][key]
        for i, (label, colour, marker) in enumerate(
            ((LOW, C_LOW, "o"), (INTACT, C_INTACT, "s"))
        ):
            ids = alt.index[alt[arm_col] == label]
            vals = alt.loc[ids, "composite"].to_numpy(float)
            jitter = rng.uniform(-0.10, 0.10, len(vals))
            # Median bar sits BEHIND the points: at n=3 a thick bar drawn on top
            # hides a data point entirely.
            if len(vals):
                ax.hlines(np.median(vals), i - 0.28, i + 0.28,
                          color=colour, linewidth=2.4, alpha=0.55, zorder=2)
            ax.scatter(i + jitter, vals, s=38, facecolor=colour, marker=marker,
                       alpha=0.95, edgecolor="white", linewidth=0.8, zorder=4)
            # Selective direct labels: only the two WRAP53-null lines.
            for j, model_id in enumerate(ids):
                if model_id in HIGHLIGHT:
                    ax.annotate(
                        HIGHLIGHT[model_id],
                        xy=(i + jitter[j], alt.loc[model_id, "composite"]),
                        xytext=(i + 0.34, alt.loc[model_id, "composite"]),
                        fontsize=7.5, color=INK, va="center", ha="left",
                        arrowprops=dict(arrowstyle="-", color=MUTED,
                                        linewidth=0.7, shrinkA=0, shrinkB=2),
                        zorder=5,
                    )
        ax.axhline(0.0, color=GRID, linewidth=0.9, linestyle="-", zorder=1)
        ax.set_xticks([0, 1])
        ax.set_xticklabels([f"{LOW}\n(n={res['n_low']})", f"{INTACT}\n(n={res['n_intact']})"],
                           fontsize=8.5, color=INK)
        ax.set_xlim(-0.6, 1.75)
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_title(title, fontsize=9, color=INK, pad=16)
        ax.text(-0.02, 1.085, tag, transform=ax.transAxes,
                fontsize=11, fontweight="bold", va="top", ha="right", color=INK)
        ax.text(0.5, 1.012,
                f"$\\delta$ = {res['cliffs_delta']:+.3f},  p = {res['p_value']:.3f}",
                transform=ax.transAxes, fontsize=8.5, ha="center", color=MUTED)

    ymin, ymax = axes[0].get_ylim()
    axes[0].set_ylabel("BIR − HR composite\n(Chronos gene effect)", fontsize=9, color=INK)
    # Predicted direction: the hypothesis says TCAB1-low sits LOWER.
    axes[0].annotate(
        "", xy=(-0.44, ymin + 0.14 * (ymax - ymin)),
        xytext=(-0.44, ymin + 0.42 * (ymax - ymin)),
        arrowprops=dict(arrowstyle="-|>", color=MUTED, linewidth=1.1),
    )
    axes[0].text(-0.50, ymin + 0.28 * (ymax - ymin), "predicted\ndirection",
                 fontsize=7, color=MUTED, rotation=90, va="center", ha="center")

    handles = [
        Line2D([], [], marker="o", color=C_LOW, linestyle="", markersize=6,
               markeredgecolor="white", label=LOW),
        Line2D([], [], marker="s", color=C_INTACT, linestyle="", markersize=6,
               markeredgecolor="white", label=INTACT),
        Line2D([], [], color=MUTED, linewidth=2.4, label="group median"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False,
               fontsize=8.5, bbox_to_anchor=(0.5, -0.035))
    fig.suptitle(
        "BIR−HR composite dependency by TCAB1/WRAP53 stratum, ALT+ osteosarcoma lines",
        fontsize=10.5, y=1.0,
    )
    fig.tight_layout(rect=[0, 0.045, 1, 0.97])
    outpath.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outpath, dpi=300, bbox_inches="tight")
    fig.savefig(outpath.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    return outpath


# ---------------------------------------------------------------------------
# Figure 2 - per-gene effect sizes
# ---------------------------------------------------------------------------

def figure2(per_gene: pd.DataFrame, primary: dict, outpath: Path) -> Path:
    order = ["composite", "POLD3", "RAD52", "SMARCAL1", "BRCA1", "BRCA2", "RAD51"]
    pri = per_gene[per_gene["contrast"] == "primary_median_8v8"].set_index("name")
    ter = per_gene[per_gene["contrast"] == "sensitivity_tertile_3v13"].set_index("name")
    pc = primary["contrasts"]["primary_median_8v8"]
    tc = primary["contrasts"]["sensitivity_tertile_3v13"]

    def row(gene, src, comp_src):
        if gene == "composite":
            return comp_src["cliffs_delta"], comp_src["delta_ci_low"], comp_src["delta_ci_high"], None
        r = src.loc[gene]
        return r["cliffs_delta"], r["delta_ci_low"], r["delta_ci_high"], r["q_value"]

    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    ax.axvspan(-1.08, 0, color="#f2f2f2", zorder=0)
    ax.axvline(0.0, color=MUTED, linewidth=1.0, linestyle="-", zorder=2)

    ypos = list(range(len(order)))[::-1]
    for y, gene in zip(ypos, order):
        if gene == "composite":
            colour, marker = C_COMPOSITE, "D"
        else:
            is_bir = gene in ("POLD3", "RAD52", "SMARCAL1")
            colour, marker = (C_BIR, "o") if is_bir else (C_HR, "s")

        d, lo, hi, q = row(gene, pri, pc)
        ax.plot([lo, hi], [y + 0.13, y + 0.13], color=colour, linewidth=1.5,
                alpha=0.9, zorder=3, solid_capstyle="round")
        ax.scatter(d, y + 0.13, s=58 if gene == "composite" else 42,
                   color=colour, marker=marker, zorder=4,
                   edgecolor="white", linewidth=0.7)

        dt, lot, hit, _ = row(gene, ter, tc)
        ax.plot([lot, hit], [y - 0.17, y - 0.17], color=colour, linewidth=1.1,
                alpha=0.45, zorder=3, solid_capstyle="round")
        ax.scatter(dt, y - 0.17, s=40 if gene == "composite" else 30,
                   facecolor="white", edgecolor=colour, marker=marker,
                   linewidth=1.1, zorder=4)

        if q is not None:
            ax.text(1.12, y, f"{q:.3f}", fontsize=7.5, va="center",
                    ha="right", color=MUTED)
        else:
            # The composite is a single pre-specified test and carries no BH
            # correction; showing its q would imply one was applied.
            ax.text(1.12, y + 0.13, f"p {pc['p_value']:.3f}*", fontsize=7.5,
                    va="center", ha="right", color=MUTED, style="italic")

    ax.text(1.12, len(order) - 0.35, "q (BH)", fontsize=7.5, va="center",
            ha="right", color=INK, fontweight="bold")
    # Footnote lives below the legend, in axes-relative coords: inside the data
    # area it ran straight through the RAD51 interval.
    ax.text(1.0, -0.265, "*composite is a single pre-specified test: p shown, not BH-adjusted",
            transform=ax.transAxes, fontsize=6.8, va="center", ha="right",
            color=MUTED, style="italic")
    ax.axhline(len(order) - 1.5, color=GRID, linewidth=0.8, zorder=1)

    ax.set_yticks(ypos)
    ax.set_yticklabels(
        ["composite\n(primary endpoint)" if g == "composite" else g for g in order],
        fontsize=8.5, color=INK,
    )
    ax.set_xlim(-1.08, 1.16)
    ax.set_ylim(-0.55, len(order) - 0.25)
    ax.set_xticks([-1.0, -0.5, 0.0, 0.5, 1.0])
    ax.set_xlabel("Cliff's $\\delta$  (TCAB1-low vs TCAB1-intact)", fontsize=9, color=INK)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0)

    ax.text(-0.54, len(order) - 0.55, "predicted direction", fontsize=8,
            color=MUTED, ha="center", style="italic")
    ax.text(0.54, len(order) - 0.55, "observed direction", fontsize=8,
            color=MUTED, ha="center", style="italic")

    handles = [
        Line2D([], [], marker="D", color=C_COMPOSITE, linestyle="", markersize=6,
               markeredgecolor="white", label="composite (primary endpoint)"),
        Line2D([], [], marker="o", color=C_BIR, linestyle="", markersize=6,
               markeredgecolor="white", label="BIR module"),
        Line2D([], [], marker="s", color=C_HR, linestyle="", markersize=6,
               markeredgecolor="white", label="HR comparator"),
        Line2D([], [], marker="o", color=MUTED, linestyle="-", markersize=6,
               markerfacecolor=MUTED, markeredgecolor="white", linewidth=1.5,
               label="8 v 8 primary (filled)"),
        Line2D([], [], marker="o", color=MUTED, linestyle="-", markersize=5.5,
               markerfacecolor="white", markeredgecolor=MUTED, linewidth=1.1,
               alpha=0.6, label="3 v 13 tertile (open)"),
    ]
    # Legend goes BELOW the axes: inside the plotting area it collided with the
    # RAD51 row.
    ax.legend(handles=handles, fontsize=7.8, frameon=False, loc="upper center",
              bbox_to_anchor=(0.5, -0.115), ncol=3, columnspacing=1.6)
    # Title states only what is true: some genes DO have negative point estimates
    # (RAD52 -0.344, BRCA2 -0.250), but every interval spans zero and nothing is
    # significant. "No gene favours the hypothesis" would overstate it.
    ax.set_title(
        "Per-gene effect sizes with 95% bootstrap CIs: every interval spans zero",
        fontsize=10.5, pad=12,
    )
    fig.tight_layout()
    outpath.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outpath, dpi=300, bbox_inches="tight")
    fig.savefig(outpath.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    return outpath


def main() -> None:
    comp, per_gene, primary = _load()
    f1 = figure1(comp, primary, FIGURES_DIR / "figure1_composite_by_stratum.png")
    f2 = figure2(per_gene, primary, FIGURES_DIR / "figure2_per_gene_forest.png")
    for p in (f1, f2):
        print(f"wrote {p.name} and {p.with_suffix('.pdf').name}")


if __name__ == "__main__":
    main()
