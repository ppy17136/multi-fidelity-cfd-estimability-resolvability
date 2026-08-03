"""Contrast-level numerical-evidence propagation and three-state decisions."""

from __future__ import annotations

from dataclasses import asdict, dataclass
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
    conservative_upper_bound: float
    lower_decision_margin: float
    upper_decision_margin: float
    state: str

    def to_dict(self) -> dict:
        return asdict(self)

    @property
    def resolved(self) -> bool:
        """Whether either side of the minimum-effect threshold is resolved.

        This compatibility flag must not be interpreted as synonymous with
        ``effect_present``.
        """
        return self.state != "indeterminate"

    @property
    def effect_present(self) -> bool:
        return self.state == "effect_present"


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


def _three_state_interval_decision(
    lower: float, upper: float, minimum_effect: float
) -> str:
    if lower > minimum_effect:
        return "effect_present"
    if upper < minimum_effect:
        return "below_minimum_effect"
    return "indeterminate"


def deterministic_resolution_decision(
    estimate: float,
    contrast_bound: float,
    minimum_effect: float,
) -> ResolutionDecision:
    """Three-state decision from a deterministic interval for |contrast|."""
    if contrast_bound < 0:
        raise ValueError("contrast_bound must be nonnegative")
    magnitude = abs(float(estimate))
    lower = max(magnitude - float(contrast_bound), 0.0)
    upper = magnitude + float(contrast_bound)
    minimum = float(minimum_effect)
    return ResolutionDecision(
        mode="deterministic_bound",
        estimate=float(estimate),
        evidence_scale=float(contrast_bound),
        minimum_effect=minimum,
        confidence=None,
        conservative_lower_bound=lower,
        conservative_upper_bound=upper,
        lower_decision_margin=lower - minimum,
        upper_decision_margin=upper - minimum,
        state=_three_state_interval_decision(lower, upper, minimum),
    )


def probabilistic_resolution_decision(
    estimate: float,
    contrast_sd: float,
    minimum_effect: float,
    *,
    confidence: float = 0.95,
) -> ResolutionDecision:
    """Central normal-model interval for absolute contrast magnitude.

    This routine is valid only when the supplied contrast standard deviation
    and approximate normal model are justified. It is not used to reinterpret
    deterministic mapping sensitivities as probabilities.
    """
    if not 0.5 < confidence < 1.0:
        raise ValueError("confidence must lie between 0.5 and 1")
    if contrast_sd < 0:
        raise ValueError("contrast_sd must be nonnegative")
    z = NormalDist().inv_cdf((1.0 + confidence) / 2.0)
    radius = z * float(contrast_sd)
    magnitude = abs(float(estimate))
    lower = max(magnitude - radius, 0.0)
    upper = magnitude + radius
    minimum = float(minimum_effect)
    return ResolutionDecision(
        mode="normal_covariance",
        estimate=float(estimate),
        evidence_scale=float(contrast_sd),
        minimum_effect=minimum,
        confidence=float(confidence),
        conservative_lower_bound=lower,
        conservative_upper_bound=upper,
        lower_decision_margin=lower - minimum,
        upper_decision_margin=upper - minimum,
        state=_three_state_interval_decision(lower, upper, minimum),
    )


def legacy_ratio_decision(
    estimate: float, floor: float, threshold: float
) -> bool:
    if floor <= 0:
        raise ValueError("floor must be positive")
    return bool(abs(float(estimate)) / float(floor) > float(threshold))
