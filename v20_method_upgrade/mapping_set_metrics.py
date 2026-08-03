"""Dependency-free set functionals for mapping-uncertainty audits."""

from __future__ import annotations

import numpy as np


def set_metrics(vectors: np.ndarray) -> dict:
    """Return the minimum norm S, exact diameter D, and ratio R=S/D."""
    x = np.asarray(vectors, dtype=np.float64)
    gram = x @ x.T
    norm2 = np.maximum(np.diag(gram), 0.0)
    distance2 = (
        norm2[:, None]
        + norm2[None, :]
        - 2.0 * gram
    )
    np.maximum(distance2, 0.0, out=distance2)
    min_index = int(np.argmin(norm2))
    maximum_flat = int(np.argmax(distance2))
    first, second = np.unravel_index(maximum_flat, distance2.shape)
    signal = float(np.linalg.norm(x[min_index]))
    diameter = float(np.linalg.norm(x[first] - x[second]))
    return {
        "set_size": int(len(x)),
        "minimum_norm": signal,
        "diameter": diameter,
        "ratio": signal / diameter if diameter > 0 else float("inf"),
        "three_times_diameter_passed": bool(signal >= 3.0 * diameter),
        "minimum_norm_member": min_index,
        "diameter_pair": [int(first), int(second)],
    }
