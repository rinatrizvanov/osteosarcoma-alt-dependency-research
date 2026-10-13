"""Iteration 04 analysis: does PRIMPOL expression predict SMARCAL1 dependency?

Implements `prespec.md`, frozen at `5f2f039`. Every threshold, seed and gene here
is quoted from that document; none is chosen in this file.

The permutation machinery is reused from iteration 03 rather than reimplemented,
so the null, the seed and the validated Spearman engine are literally the same
code that 58 tests already cover.
"""
from __future__ import annotations

import importlib.util
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd
from scipy import stats

ITER_DIR = Path(__file__).resolve().parent
ITER03_COHERENCE = (
    ITER_DIR.parent / "iteration_03_bir_module_coherence" / "src" / "coherence.py"
)


def _load_coherence():
    spec = importlib.util.spec_from_file_location("iter03_coherence", ITER03_COHERENCE)
    module = importlib.util.module_from_spec(spec)
    sys.modules["iter03_coherence"] = module
    spec.loader.exec_module(module)
    return module


coherence = _load_coherence()

# --- Quoted from the frozen prespec -----------------------------------------

PRIMARY_EXPOSURE = "PRIMPOL"       # §1
NEGATIVE_CONTROL = "ZRANB3"        # §5 S1
DESCRIPTIVE_ONLY = "BLM"           # §5 S4 -- fails the §3 gate, uninterpretable
OUTCOME_GENE = "SMARCAL1"          # §4
PREDICTED_SIGN = -1                # §1

#: §3 exposure-variance gate.
GATE_MIN_SD = 0.50
GATE_MAX_GUARDED_GAP = 1.00
GATE_MIN_SHAPIRO_P = 0.05

MIN_LINES = 12                     # §4 minimum-cell guard
SEED = coherence.SEED
N_PERMUTATIONS = coherence.N_PERMUTATIONS

#: Pre-registered floors from §4, recomputed at run time and asserted to match.
FLOOR_N16_TWO_SIDED = 0.500
CAMPAIGN_BAR_P = 0.033


@dataclass(frozen=True)
class GateResult:
    gene: str
    n: int
    sd: float
    n_zero: int
    guarded_gap: float
    shapiro_p: float

    @property
    def passed(self) -> bool:
        return (
            self.sd >= GATE_MIN_SD
            and self.n_zero == 0
            and self.guarded_gap < GATE_MAX_GUARDED_GAP
            and self.shapiro_p >= GATE_MIN_SHAPIRO_P
        )

    def row(self) -> dict:
        return {
            "gene": self.gene,
            "n": self.n,
            "sd": round(self.sd, 4),
            "n_zero": self.n_zero,
            "guarded_gap": round(self.guarded_gap, 4),
            "shapiro_p": round(self.shapiro_p, 4),
            "sd_ok": self.sd >= GATE_MIN_SD,
            "no_zeros_ok": self.n_zero == 0,
            "gap_ok": self.guarded_gap < GATE_MAX_GUARDED_GAP,
            "shapiro_ok": self.shapiro_p >= GATE_MIN_SHAPIRO_P,
            "gate_passed": self.passed,
        }


def guarded_max_gap(values: Sequence[float], min_cell: int = 4) -> float:
    """Largest interior gap with at least `min_cell` points on both sides.

    The same diagnostic iteration 01 used on WRAP53 and the prespec §3 gate uses
    on the exposures. Returns nan when no guard-respecting split exists.
    """
    v = np.sort(np.asarray(values, dtype=float))
    gaps = np.diff(v)
    eligible = [
        gaps[i] for i in range(len(gaps)) if i + 1 >= min_cell and len(v) - (i + 1) >= min_cell
    ]
    return float(max(eligible)) if eligible else float("nan")


def evaluate_gate(gene: str, values: Sequence[float]) -> GateResult:
    """§3 exposure-variance gate, applied identically to every exposure."""
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    return GateResult(
        gene=gene,
        n=int(len(v)),
        sd=float(v.std(ddof=1)),
        n_zero=int((v == 0).sum()),
        guarded_gap=guarded_max_gap(v),
        shapiro_p=float(stats.shapiro(v).pvalue),
    )


@dataclass(frozen=True)
class Correlation:
    label: str
    exposure: str
    outcome: str
    n: int
    rho: float
    p_value: float
    floor: float

    def row(self) -> dict:
        return {
            "analysis": self.label,
            "exposure": self.exposure,
            "outcome": self.outcome,
            "n": self.n,
            "rho": round(self.rho, 4),
            "p_perm": round(self.p_value, 4),
            "floor_p05": round(self.floor, 4),
            "reaches_floor": abs(self.rho) >= self.floor,
            "aligned": round(self.rho * PREDICTED_SIGN, 4),
        }


def correlate(
    exposure: pd.Series,
    outcome: pd.Series,
    *,
    label: str,
    exposure_name: str,
    outcome_name: str,
) -> Correlation:
    """Spearman rho with the iteration 03 permutation null at this n.

    Both series are aligned on ModelID and any line missing either value is
    dropped by the caller, not here.
    """
    joined = pd.concat([exposure, outcome], axis=1).dropna()
    x = joined.iloc[:, 0].to_numpy(dtype=float)
    y = joined.iloc[:, 1].to_numpy(dtype=float)
    rho = coherence._spearman(x, y)
    null = coherence.null_distribution(len(joined), 2, seed=SEED)
    return Correlation(
        label=label,
        exposure=exposure_name,
        outcome=outcome_name,
        n=int(len(joined)),
        rho=float(rho),
        p_value=coherence.permutation_p_two_sided(rho, null),
        floor=float(np.quantile(np.abs(null), 0.95)),
    )


def leave_one_out(exposure: pd.Series, outcome: pd.Series) -> pd.DataFrame:
    """§5 S5: refit the primary dropping each line in turn."""
    joined = pd.concat([exposure, outcome], axis=1).dropna()
    rows = []
    for model_id in joined.index:
        sub = joined.drop(index=model_id)
        rows.append(
            {
                "dropped": model_id,
                "n": len(sub),
                "rho": round(
                    coherence._spearman(
                        sub.iloc[:, 0].to_numpy(dtype=float),
                        sub.iloc[:, 1].to_numpy(dtype=float),
                    ),
                    4,
                ),
            }
        )
    return pd.DataFrame(rows)
