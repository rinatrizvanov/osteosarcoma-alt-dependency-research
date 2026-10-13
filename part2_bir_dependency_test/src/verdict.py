"""Exhaustive verdict partition for autoresearch iterations.

Why this module exists, and why it is **not** an edit to :mod:`src.criteria`.

``src/criteria.py`` implements the Stage 1 prespec v1.2 Sections 8/9/10 exactly
as frozen at tag ``stage1-v1.0-null``. Editing it would change the recorded
Stage 1 verdict and break notebook 07's idempotency, so it is left untouched.
This module is the generalised replacement used by every autoresearch
iteration.

The defect it fixes
-------------------
``criteria.evaluate`` leaves three regions of the (direction, p, |delta|)
space assigned to no named outcome, so they fall through to
``indeterminate``:

1. counter-direction, ``p >= alpha``, ``|delta| >= 0.20`` — because the
   wrong-direction test required significance and the clean-null test required
   ``|delta| < FALSIFIED_MAX_ABS_DELTA``. *This is exactly Stage 1's own
   result*: delta = +0.281 against a predicted-negative direction, p = 0.382.
2. predicted direction, ``p >= alpha``, ``0.20 <= |delta| < 0.40`` — the gap
   between ``FALSIFIED_MAX_ABS_DELTA`` (0.200) and
   ``INCONCLUSIVE_MIN_ABS_DELTA`` (0.400).
3. predicted direction, ``|delta| >= 0.40``, but the robustness checks
   disagree — Section 9 additionally gated on ``directionally_consistent``.

The fix is structural, not another threshold. :func:`evaluate_iteration`
dispatches through an if/elif chain over the real line that is exhaustive and
mutually exclusive **by construction**, so every finite input reaches exactly
one named outcome. Robustness consistency and secondary support are reported
as *conditions* and can only move a result between ``supported`` and
``partial_support``; they can no longer drop it out of the partition.

``indeterminate`` is reserved for genuinely missing or non-finite inputs. It is
never the answer for an in-range result.

Direction handling
------------------
Hypotheses differ in which sign confirms them, so everything is computed on the
*aligned* effect::

    aligned = delta * predicted_sign

``aligned > 0`` always means "in the direction the hypothesis predicted",
whatever that direction is. Stage 1 predicted ``delta < 0``, i.e.
``predicted_sign = -1``.
"""
from __future__ import annotations

import math
from dataclasses import asdict, dataclass, field
from functools import lru_cache
from typing import Iterable, Sequence

ALPHA = 0.05

#: |aligned delta| below this is a clean null in either direction. Inherited
#: from Stage 1 ``criteria.FALSIFIED_MAX_ABS_DELTA`` so the two evaluators
#: agree on what "no effect" means.
NULL_CEILING = 0.200

#: Reported as a support *condition*, never as a branch gate — gating on it
#: would reopen an unassigned cell at ``p < alpha`` with a small effect.
SUPPORT_MIN_ABS_DELTA = 0.330

OUTCOMES: tuple[str, ...] = (
    "supported",
    "partial_support",
    "inconclusive_underpowered",
    "counter_signal",
    "refuted",
    "null",
    "indeterminate",
)

#: The six outcomes a finite result can reach. ``indeterminate`` is excluded:
#: reaching it from finite inputs is a bug, and the test suite asserts so.
FINITE_OUTCOMES: tuple[str, ...] = OUTCOMES[:-1]


@dataclass
class Condition:
    key: str
    description: str
    passed: bool
    detail: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


@dataclass
class IterationVerdict:
    outcome: str
    rule: str
    headline: str
    aligned_delta: float
    p_value: float
    predicted_sign: int
    conditions: list[Condition] = field(default_factory=list)
    attributes: dict = field(default_factory=dict)

    @property
    def is_positive(self) -> bool:
        """``supported`` only. ``partial_support`` is not a positive result."""
        return self.outcome == "supported"

    def as_dict(self) -> dict:
        return {
            "outcome": self.outcome,
            "rule": self.rule,
            "headline": self.headline,
            "aligned_delta": self.aligned_delta,
            "p_value": self.p_value,
            "predicted_sign": self.predicted_sign,
            "conditions": [c.as_dict() for c in self.conditions],
            "attributes": self.attributes,
        }

    def ledger_line(self) -> str:
        return (
            f"outcome=`{self.outcome}` aligned_delta={self.aligned_delta:+.3f} "
            f"p={self.p_value:.4g} rule={self.rule}"
        )


def _finite(*values: float) -> bool:
    try:
        return all(math.isfinite(float(v)) for v in values)
    except (TypeError, ValueError):
        return False


def evaluate_iteration(
    *,
    delta: float,
    p_value: float,
    predicted_sign: int,
    conditions: Sequence[Condition] = (),
    alpha: float = ALPHA,
    null_ceiling: float = NULL_CEILING,
    attributes: dict | None = None,
) -> IterationVerdict:
    """Map one primary result to exactly one named outcome.

    Parameters
    ----------
    delta
        Observed Cliff's delta (or any signed effect size on the same scale).
    p_value
        Two-sided p for the primary endpoint.
    predicted_sign
        ``-1`` if the hypothesis predicts a negative delta, ``+1`` if positive.
    conditions
        Pre-specified support conditions. They distinguish ``supported`` from
        ``partial_support`` and are reported everywhere; they never leave a
        result unassigned.

    The dispatch below is the whole point: the branches tile the real line, so
    there is no cell without a verdict.
    """
    attributes = dict(attributes or {})
    conditions = list(conditions)

    if predicted_sign not in (-1, 1):
        raise ValueError(f"predicted_sign must be -1 or +1, got {predicted_sign!r}")

    if not _finite(delta, p_value):
        return IterationVerdict(
            "indeterminate", "missing-input",
            "Primary delta or p-value is missing or non-finite, so no criterion "
            "can be applied. This is an input defect, not a result.",
            float("nan"), float("nan"), predicted_sign, conditions, attributes,
        )

    delta = float(delta)
    p = float(p_value)
    aligned = delta * predicted_sign
    all_pass = bool(conditions) and all(c.passed for c in conditions)
    failed = [c.key for c in conditions if not c.passed]
    sig = p < alpha

    def _v(outcome: str, rule: str, headline: str) -> IterationVerdict:
        attributes.setdefault("observed_delta", delta)
        attributes.setdefault("all_conditions_passed", all_pass)
        attributes.setdefault("failed_conditions", failed)
        attributes.setdefault("significant", sig)
        return IterationVerdict(
            outcome, rule, headline, aligned, p, predicted_sign,
            conditions, attributes,
        )

    # --- the partition -----------------------------------------------------
    # Branch 1 catches every |aligned| below the ceiling, including exactly
    # zero, in either direction and at any p. Branches 2 and 3 then split the
    # remainder by sign, and each splits by significance. Exhaustive.
    if abs(aligned) < null_ceiling:
        contradiction = (
            " Note: p is significant while the effect is below the null "
            "ceiling, which at these group sizes should be impossible — check "
            "the inputs." if sig else ""
        )
        return _v(
            "null", "clean-null",
            f"Clean null: |aligned delta| = {abs(aligned):.3f} < {null_ceiling} "
            f"(p = {p:.4g}). No effect of interpretable size in either "
            f"direction.{contradiction}",
        )

    if aligned > 0:
        if sig:
            if all_pass:
                return _v(
                    "supported", "support",
                    f"Supported: aligned delta = {aligned:+.3f} in the predicted "
                    f"direction, p = {p:.4g} < {alpha}, and every pre-specified "
                    f"support condition met.",
                )
            return _v(
                "partial_support", "partial-support",
                f"Partial support: aligned delta = {aligned:+.3f} in the predicted "
                f"direction at p = {p:.4g}, but these pre-specified conditions "
                f"failed: {failed or ['no conditions supplied']}. Not a positive "
                f"result.",
            )
        return _v(
            "inconclusive_underpowered", "inconclusive",
            f"Directionally consistent but N-limited: aligned delta = "
            f"{aligned:+.3f} in the predicted direction with p = {p:.4g}. "
            f"Neither a positive result nor a refutation.",
        )

    # aligned <= -null_ceiling: a counter-directional effect of real size.
    if sig:
        return _v(
            "refuted", "refuted",
            f"Refuted: significant effect OPPOSITE to the prediction "
            f"(aligned delta = {aligned:+.3f}, p = {p:.4g} < {alpha}).",
        )
    return _v(
        "counter_signal", "counter-signal",
        f"Counter-signal: non-significant but substantial effect opposite to "
        f"the prediction (aligned delta = {aligned:+.3f}, |delta| >= "
        f"{null_ceiling}, p = {p:.4g}). Not a clean null — the estimate runs "
        f"the wrong way. This is the cell Stage 1 landed in and that "
        f"criteria.evaluate left unnamed.",
    )


# ---------------------------------------------------------------------------
# Exact Mann-Whitney null, for stating minimum detectable effect in advance
# ---------------------------------------------------------------------------

@lru_cache(maxsize=256)
def _u_counts(n: int, m: int) -> tuple[int, ...]:
    """Number of rank assignments giving each U in 0..n*m. Exact integers."""
    if n < 0 or m < 0:
        raise ValueError("group sizes must be non-negative")
    # table[i][j] = counts for groups of size i and j
    prev_row: list[list[int]] = [[1] for _ in range(m + 1)]
    for i in range(1, n + 1):
        row: list[list[int]] = [[1]]
        for j in range(1, m + 1):
            size = i * j + 1
            acc = [0] * size
            left = row[j - 1]            # f(i, j-1)
            up = prev_row[j]             # f(i-1, j)
            for u, c in enumerate(left):
                if c:
                    acc[u] += c
            for u, c in enumerate(up):
                if c and u + j < size:
                    acc[u + j] += c
            row.append(acc)
        prev_row = row
    return tuple(prev_row[m])


def exact_mwu_two_sided_p(u: int, n: int, m: int) -> float:
    """Exact two-sided Mann-Whitney p for statistic ``u``, no ties.

    Matches ``scipy.stats.mannwhitneyu(..., method="exact")``: twice the
    smaller tail, clipped at 1.
    """
    counts = _u_counts(n, m)
    total = sum(counts)
    if total == 0:
        return float("nan")
    p_le = sum(counts[: u + 1]) / total
    p_ge = sum(counts[u:]) / total
    return min(1.0, 2.0 * min(p_le, p_ge))


def delta_from_u(u: int, n: int, m: int) -> float:
    """Cliff's delta from the Mann-Whitney U: ``delta = 2U/(nm) - 1``."""
    return 2.0 * u / (n * m) - 1.0


def min_detectable_delta(n: int, m: int, alpha: float = ALPHA) -> float:
    """Smallest |Cliff's delta| that can reach two-sided exact p < alpha.

    Computed by exhaustive enumeration of the exact null rather than read from
    a lookup table, so any iteration can state its own detection boundary for
    whatever group sizes it realises. Returns NaN when no configuration can
    reach significance at these sizes — which is itself a finding worth
    pre-registering.
    """
    if n <= 0 or m <= 0:
        return float("nan")
    best = float("inf")
    for u in range(n * m + 1):
        if exact_mwu_two_sided_p(u, n, m) < alpha:
            best = min(best, abs(delta_from_u(u, n, m)))
    return best if math.isfinite(best) else float("nan")


def power_note(n: int, m: int, alpha: float = ALPHA) -> str:
    """One sentence for a prespec's minimum-detectable-effect clause."""
    mde = min_detectable_delta(n, m, alpha)
    if not math.isfinite(mde):
        return (
            f"At {n}v{m} no rank configuration reaches two-sided exact "
            f"p < {alpha}: this comparison cannot produce a significant result "
            f"at any effect size, and a non-significant outcome is therefore "
            f"uninformative."
        )
    return (
        f"At {n}v{m} the smallest |Cliff's delta| able to reach two-sided exact "
        f"p < {alpha} is {mde:.3f}. Any effect smaller than that cannot be "
        f"detected here regardless of the true value."
    )


# ---------------------------------------------------------------------------
# Campaign-wide multiplicity (anti-bias rule 7a)
# ---------------------------------------------------------------------------

def benjamini_hochberg(p_values: Sequence[float]) -> list[float]:
    """BH step-up adjusted q-values, order preserved. NaNs pass through."""
    indexed = [(i, float(p)) for i, p in enumerate(p_values)]
    finite = [(i, p) for i, p in indexed if math.isfinite(p)]
    q: list[float] = [float("nan")] * len(indexed)
    if not finite:
        return q
    finite.sort(key=lambda t: t[1])
    n = len(finite)
    running = 1.0
    for rank in range(n, 0, -1):
        i, p = finite[rank - 1]
        running = min(running, p * n / rank)
        q[i] = min(1.0, running)
    return q


def campaign_fdr(
    records: Iterable[dict],
    *,
    q_threshold: float = 0.100,
    p_key: str = "p_value",
) -> list[dict]:
    """Apply BH across every iteration's *primary* endpoint.

    This is the experiment-level success criterion. Each iteration is
    individually pre-registered, but the hypotheses are chosen in sequence on
    one fixed cohort, so per-iteration p-values understate the campaign error
    rate. A result counts as positive only at ``q < q_threshold`` here.
    """
    rows = [dict(r) for r in records]
    qs = benjamini_hochberg([r.get(p_key, float("nan")) for r in rows])
    for r, q in zip(rows, qs):
        r["q_value_campaign"] = q
        r["survives_campaign_fdr"] = bool(math.isfinite(q) and q < q_threshold)
        r["q_threshold"] = q_threshold
    return rows
