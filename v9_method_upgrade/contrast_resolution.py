"""Contrast-level numerical-evidence propagation and decisions."""

from __future__ import annotations

from dataclasses import dataclass, asdict
from statistics import NormalDist
from typing import Sequence

import numpy as np


Array = np.ndarray


@dataclass
class ResolutionDecision:
    mode: str
    estimate: float
    evidence_scale: float
    minimum_effect: float
    confidence: float | None
    conservative_lower_bound: float
    decision_margin: float
    resolved: bool

    def to_dict(self) -> dict:
        return asdict(self)


def deterministic_contrast_bound(
    coefficients: Sequence[float], cell_bounds: Sequence[float]
) -> float:
    """Triangle-inequality bound for a scalar or norm-bounded contrast."""
    c = np.asarray(coefficients, dtype=float)
    e = np.asarray(cell_bounds, dtype=float)
    if c.shape != e.shape:
        raise ValueError("coefficients and cell bounds must have the same shape")
    if np.any(e < 0):
        raise ValueError("cell bounds must be nonnegative")
    return float(np.sum(np.abs(c) * e))


def contrast_standard_deviation(
    coefficients: Sequence[float], covariance: Array
) -> float:
    c = np.asarray(coefficients, dtype=float)
    sigma = np.asarray(covariance, dtype=float)
    if sigma.shape != (len(c), len(c)):
        raise ValueError("covariance shape is incompatible with coefficients")
    variance = float(c @ sigma @ c)
    if variance < -1e-12:
        raise ValueError("covariance produces a negative contrast variance")
    return float(np.sqrt(max(variance, 0.0)))


def deterministic_resolution_decision(
    estimate: float,
    contrast_bound: float,
    minimum_effect: float,
) -> ResolutionDecision:
    """Resolve only if a deterministic lower bound exceeds minimum effect."""
    lower = max(abs(float(estimate)) - float(contrast_bound), 0.0)
    margin = lower - float(minimum_effect)
    return ResolutionDecision(
        mode="deterministic_bound",
        estimate=float(estimate),
        evidence_scale=float(contrast_bound),
        minimum_effect=float(minimum_effect),
        confidence=None,
        conservative_lower_bound=lower,
        decision_margin=margin,
        resolved=bool(margin > 0.0),
    )


def probabilistic_resolution_decision(
    estimate: float,
    contrast_sd: float,
    minimum_effect: float,
    *,
    confidence: float = 0.95,
) -> ResolutionDecision:
    """Normal-model lower confidence bound on absolute contrast magnitude.

    This routine is valid only when the supplied contrast standard deviation
    and approximate normal model are justified.  It is not used to reinterpret
    deterministic mapping sensitivities as probabilities.
    """
    if not 0.5 < confidence < 1.0:
        raise ValueError("confidence must lie between 0.5 and 1")
    z = NormalDist().inv_cdf((1.0 + confidence) / 2.0)
    radius = z * float(contrast_sd)
    lower = max(abs(float(estimate)) - radius, 0.0)
    margin = lower - float(minimum_effect)
    return ResolutionDecision(
        mode="normal_covariance",
        estimate=float(estimate),
        evidence_scale=float(contrast_sd),
        minimum_effect=float(minimum_effect),
        confidence=float(confidence),
        conservative_lower_bound=lower,
        decision_margin=margin,
        resolved=bool(margin > 0.0),
    )


def legacy_ratio_decision(
    estimate: float, floor: float, threshold: float
) -> bool:
    if floor <= 0:
        raise ValueError("floor must be positive")
    return bool(abs(float(estimate)) / float(floor) > float(threshold))

