"""Shared fixtures. Adds the project root to ``sys.path`` (Aim 2 test idiom)."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))


@pytest.fixture
def synthetic_expression() -> pd.Series:
    """34 reference lines plus 16 ALT+ lines, WRAP53 log2(TPM+1).

    Deterministic by construction rather than sampled: the reference is an even
    spread over [2, 8] and the ALT+ set is an even spread over the same range,
    so the bottom-tertile cut lands where arithmetic says it does and the tests
    exercise the rule instead of a particular RNG draw.
    """
    ref = pd.Series(
        np.linspace(2.0, 8.0, 34), index=[f"R{i:02d}" for i in range(34)]
    )
    alt = pd.Series(
        np.linspace(2.0, 8.0, 16), index=[f"A{i:02d}" for i in range(16)]
    )
    return pd.concat([ref, alt])


@pytest.fixture
def toy_effects() -> pd.DataFrame:
    """Gene effects with a planted BIR-shift in the first 8 lines."""
    from src.dependency_stats import SIX_GENES

    rng = np.random.default_rng(11)
    ids = [f"L{i:02d}" for i in range(16)]
    # Noise small relative to the planted shift, so the tests below assert on
    # the code's behaviour rather than on a particular RNG draw.
    data = rng.normal(-0.2, 0.05, (16, 6))
    data[:8, :3] -= 0.7  # BIR genes much more essential in the "low" group
    return pd.DataFrame(data, index=ids, columns=list(SIX_GENES))
