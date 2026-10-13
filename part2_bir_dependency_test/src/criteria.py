"""Mechanical evaluation of the Section 8 / 9 / 10 decision rules.

Notebook 07 calls :func:`evaluate` and reports whatever comes back. The point
of implementing the criteria as code is that the verdict is computed from the
rules **as written in the frozen prespec**, not from anyone's reading of the
numbers after seeing them.

Three named outcomes plus an explicit fourth:

* ``supported``    — every Section 8 condition met.
* ``inconclusive`` — Section 9: not significant, but a large effect in the
  predicted direction with consistent sensitivity analyses. At N=16 this is a
  likely outcome, and it is the one that justifies Stage 2.
* ``falsified``    — Section 10: a null, or a significant effect in the wrong
  direction. Stop.
* ``indeterminate`` — matches none of the three. Reported as such rather than
  being forced into the nearest bucket; the per-condition table below shows
  exactly which conditions failed.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Sequence

import pandas as pd

# Section 8 thresholds
ALPHA = 0.05

#: Smallest |Cliff's delta| that can reach two-sided exact p<0.05, by split.
#: Computed exhaustively over rank configurations. These make the real
#: detection boundary explicit: at the expected 5v11 tertile split the p<0.05
#: requirement alone forces |delta| >= 0.673, so SUPPORT_DELTA below is not the
#: binding constraint. Recorded here so no one later mistakes 0.33 for the bar.
MIN_DETECTABLE_DELTA: dict[str, float] = {
    "5v11": 0.673, "8v8": 0.594, "3v7": 0.905, "5v5": 0.840,
    "3v6": 0.889, "4v5": 0.900,
}
SUPPORT_DELTA = -0.330          # Cliff's delta, predicted (negative) direction
R1_MIN_ABS_DELTA = 0.300
BIR_Q_THRESHOLD = 0.100
MIN_BIR_GENES = 2

# Section 9 / 10 thresholds
INCONCLUSIVE_MIN_ABS_DELTA = 0.400
FALSIFIED_MAX_ABS_DELTA = 0.200

#: Section 6.4 concentration comparison (prespec v1.1). If the 3v13 tertile
#: contrast out-deltas the 8v8 median primary by more than this, the signal
#: concentrates in the extreme-WRAP53-loss lines. DESCRIPTIVE AND NON-GATING:
#: it never changes the Section 8/9/10 verdict, it bears on Stage 2
#: justification. Pre-registered before any outcome was read.
CONCENTRATION_DELTA_THRESHOLD = 0.300


@dataclass
class Condition:
    key: str
    description: str
    passed: bool
    detail: str

    def as_dict(self) -> dict:
        return asdict(self)


@dataclass
class Verdict:
    outcome: str
    section: str
    headline: str
    conditions: list[Condition] = field(default_factory=list)
    #: Section 6.4 descriptive comparison. Reported, never gating.
    concentration: dict = field(default_factory=dict)

    def as_dict(self) -> dict:
        return {
            "outcome": self.outcome,
            "section": self.section,
            "headline": self.headline,
            "conditions": [c.as_dict() for c in self.conditions],
            "concentration": self.concentration,
        }

    def as_markdown(self) -> list[str]:
        lines = [
            f"**Outcome: `{self.outcome}`** ({self.section})",
            "",
            self.headline,
            "",
            "| Condition | Met | Detail |",
            "|---|---|---|",
        ]
        for c in self.conditions:
            lines.append(
                f"| {c.description} | {'yes' if c.passed else 'NO'} | {c.detail} |"
            )
        if self.concentration.get("available"):
            lines += [
                "",
                "**Section 6.4 concentration comparison** (descriptive, non-gating):",
                "",
                f"- delta, 8v8 median primary: {self.concentration['delta_median_8v8']:+.3f}",
                f"- delta, 3v13 OS-reference tertile: {self.concentration['delta_tertile_3v13']:+.3f}",
                f"- |delta| excess of tertile over median: "
                f"{self.concentration['abs_excess_tertile_over_median']:+.3f} "
                f"(pre-registered threshold {self.concentration['threshold']})",
                "",
                self.concentration["interpretation"],
            ]
        return lines


def _bir_support(per_gene: pd.DataFrame) -> tuple[bool, str]:
    """>=2 of 3 BIR genes at q<0.10 in the predicted direction, and no HR gene stronger."""
    bir = per_gene[per_gene["gene_set"] == "BIR"]
    hr = per_gene[per_gene["gene_set"] == "HR"]
    hits = bir[(bir["q_value"] < BIR_Q_THRESHOLD) & (bir["cliffs_delta"] < 0)]
    n_hits = len(hits)

    strongest_bir = bir["cliffs_delta"].min() if len(bir) else 0.0
    hr_same_dir = hr[hr["cliffs_delta"] < 0]
    hr_stronger = (
        bool(len(hr_same_dir)) and float(hr_same_dir["cliffs_delta"].min()) < float(strongest_bir)
    )
    detail = (
        f"{n_hits}/3 BIR genes at q<{BIR_Q_THRESHOLD} in predicted direction "
        f"({', '.join(hits['name']) if n_hits else 'none'}); "
        f"strongest BIR delta {strongest_bir:+.3f}; "
    )
    if len(hr_same_dir):
        detail += (
            f"strongest same-direction HR delta "
            f"{float(hr_same_dir['cliffs_delta'].min()):+.3f}"
        )
        if hr_stronger:
            detail += " (HR EXCEEDS BIR - fails specificity)"
    else:
        detail += "no HR gene in the predicted direction"
    return (n_hits >= MIN_BIR_GENES and not hr_stronger), detail


def concentration_comparison(
    primary_median: dict, tertile: dict | None
) -> dict:
    """Section 6.4: does the effect concentrate in the extreme-loss lines?

    Compares |Cliff's delta| for the 3v13 OS-reference tertile contrast against
    the 8v8 median primary. The tertile arm is dominated by the two
    WRAP53-null ALT+ lines, so a markedly stronger tertile effect indicates the
    signal lives in genuine biallelic TCAB1 loss rather than in the diluted
    median split.

    **Descriptive and non-gating.** This never enters the Section 8/9/10
    verdict. It was pre-registered (prespec v1.1) before any gene-effect value
    was read, so it is a planned comparison rather than a post-hoc observation,
    and it speaks to whether Stage 2 is justified.
    """
    if tertile is None:
        return {"available": False}
    d_med = float(primary_median["cliffs_delta"])
    d_ter = float(tertile["cliffs_delta"])
    excess = abs(d_ter) - abs(d_med)
    return {
        "available": True,
        "delta_median_8v8": d_med,
        "delta_tertile_3v13": d_ter,
        "abs_excess_tertile_over_median": excess,
        "threshold": CONCENTRATION_DELTA_THRESHOLD,
        "concentrates_in_extreme_loss": bool(
            excess > CONCENTRATION_DELTA_THRESHOLD and d_ter < 0
        ),
        "both_same_direction": bool(d_med < 0 and d_ter < 0),
        "interpretation": (
            "Signal concentrates in the extreme-WRAP53-loss lines: the 3v13 tertile "
            "contrast exceeds the 8v8 median primary by "
            f"{excess:+.3f} in |delta|, above the pre-registered {CONCENTRATION_DELTA_THRESHOLD} "
            "threshold. Supports the dilution account of the median split and strengthens "
            "the Stage 2 case."
            if excess > CONCENTRATION_DELTA_THRESHOLD and d_ter < 0 else
            "No evidence of concentration in the extreme-loss lines: the tertile contrast "
            f"exceeds the median primary by {excess:+.3f} in |delta|, within the "
            f"pre-registered {CONCENTRATION_DELTA_THRESHOLD} threshold."
        ),
    }


def evaluate(
    primary: dict,
    per_gene: pd.DataFrame,
    r1_gold_only: dict,
    r2_median: dict | None,
    loo_deltas: Sequence[float],
    interaction: dict,
    tertile: dict | None = None,
) -> Verdict:
    """Apply Sections 8, 9 and 10 in order.

    Parameters mirror the artefacts notebooks 05 and 06 produce:
    ``primary`` / ``r1_gold_only`` / ``r2_median`` are
    :meth:`src.dependency_stats.GroupComparison.as_dict` outputs, ``per_gene``
    is :func:`src.dependency_stats.per_gene_comparisons`, ``loo_deltas`` is the
    R7 leave-one-line-out delta series, and ``interaction`` is
    :func:`src.dependency_stats.delta_difference_permutation`.
    """
    p = float(primary["p_value"])
    delta = float(primary["cliffs_delta"])
    r1_delta = float(r1_gold_only["cliffs_delta"])
    loo = [float(d) for d in loo_deltas]

    c1 = Condition(
        "primary", "S8.1 primary p<0.05 and delta <= -0.33",
        p < ALPHA and delta <= SUPPORT_DELTA,
        f"p={p:.4g}, delta={delta:+.3f}",
    )
    c2 = Condition(
        "r1_gold_only", "S8.2 gold-only same direction, |delta| >= 0.30",
        r1_delta < 0 and abs(r1_delta) >= R1_MIN_ABS_DELTA,
        f"gold-only delta={r1_delta:+.3f}",
    )
    bir_ok, bir_detail = _bir_support(per_gene)
    c3 = Condition("bir_genes", "S8.3 >=2/3 BIR genes at q<0.10, no HR gene stronger", bir_ok, bir_detail)
    loo_ok = bool(loo) and all(d < 0 for d in loo)
    c4 = Condition(
        "r7_loo", "S8.4 no sign flip on leave-one-line-out",
        loo_ok,
        f"{len(loo)} refits, delta range [{min(loo):+.3f}, {max(loo):+.3f}]" if loo else "not run",
    )
    # S8.5 is a NON-CONTRADICTION test, not a point comparison of deltas.
    #
    # The original formulation ("ALT- delta closer to zero than ALT+ delta")
    # was replaced before the prespec was frozen because it fails far too often
    # under a true effect: at the real ALT- group sizes (3v6 / 4v5) the exact
    # null distribution of |Cliff's delta| puts 26-29% of its mass at or above
    # 0.474, so a noise-only ALT- stratum routinely out-deltas a genuine ALT+
    # effect. Gating on it would have given a ~20-30% false-failure rate.
    #
    # Requiring instead that the ALT- stratum show no *significant*
    # same-direction shift has a 5% false-failure rate by construction, while
    # still catching the case the criterion exists for: a TCAB1 effect that is
    # general rather than ALT-specific.
    alt_neg_p = float(interaction.get("alt_negative_p_value", float("nan")))
    alt_neg_delta = float(interaction.get("delta_alt_negative", float("nan")))
    contradicted = alt_neg_p < ALPHA and alt_neg_delta < 0
    c5 = Condition(
        "selectivity",
        "S8.5 ALT- stratum shows no significant same-direction shift",
        not contradicted,
        f"ALT- delta={alt_neg_delta:+.3f}, p={alt_neg_p:.4g}; "
        f"ALT+ delta={interaction.get('delta_alt_positive', float('nan')):+.3f} "
        f"(delta difference reported as a diagnostic, not a gate)",
    )
    conditions = [c1, c2, c3, c4, c5]
    conc = concentration_comparison(primary, tertile)

    def _finish(verdict: Verdict) -> Verdict:
        verdict.concentration = conc
        return verdict

    if all(c.passed for c in conditions):
        return _finish(Verdict(
            "supported", "Section 8",
            "Every pre-specified support condition met. TCAB1-low ALT+ lines show "
            "selectively elevated BIR-relative-to-HR dependency.",
            conditions,
        ))

    # Section 9 — inconclusive band.
    r2_ok = r2_median is None or float(r2_median["cliffs_delta"]) < 0
    directionally_consistent = r1_delta < 0 and r2_ok and loo_ok
    if p >= ALPHA and delta <= -INCONCLUSIVE_MIN_ABS_DELTA and directionally_consistent:
        return _finish(Verdict(
            "inconclusive", "Section 9",
            f"Directionally consistent but N-limited: delta={delta:+.3f} in the "
            f"predicted direction with p={p:.4g}. Not a positive result and not "
            f"falsified. This is the outcome that motivates Stage 2, where the PDX "
            f"cohort adds N and replaces the expression proxy with real 17p13.1 SV calls.",
            conditions,
        ))

    # Section 10 — falsified.
    no_bir_hit = not ((per_gene["gene_set"] == "BIR") & (per_gene["q_value"] < BIR_Q_THRESHOLD)).any()
    null_result = p >= ALPHA and abs(delta) < FALSIFIED_MAX_ABS_DELTA
    wrong_direction = p < ALPHA and delta > 0
    if (null_result or wrong_direction) and no_bir_hit:
        why = (
            f"null result (p={p:.4g}, |delta|={abs(delta):.3f} < {FALSIFIED_MAX_ABS_DELTA})"
            if null_result
            else f"significant effect in the OPPOSITE direction (p={p:.4g}, delta={delta:+.3f})"
        )
        return _finish(Verdict(
            "falsified", "Section 10",
            f"Hypothesis falsified: {why}, and no BIR gene reaches q<{BIR_Q_THRESHOLD}. "
            f"Run only the Section 11 secondary set, write the negative-result memo, "
            f"and stop. Do not approach Stage 2 on this result.",
            conditions,
        ))

    failed = [c.key for c in conditions if not c.passed]
    return _finish(Verdict(
        "indeterminate", "none of Sections 8/9/10",
        f"The result matches no pre-specified outcome bucket: p={p:.4g}, "
        f"delta={delta:+.3f}. Failed support conditions: {failed}. Report this "
        f"as-is; do not round it toward the nearest verdict.",
        conditions,
    ))
