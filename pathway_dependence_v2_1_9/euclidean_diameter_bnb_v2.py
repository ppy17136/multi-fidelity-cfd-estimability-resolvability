"""Numerically robust intrinsic-coordinate front end for diameter BnB.

This module supersedes the Gram-eigendecomposition rank extraction in the
initial prototype.  A direct thin SVD avoids squaring the condition number.
"""

from __future__ import annotations

import numpy as np

from euclidean_diameter_bnb import (
    DiameterCertificate,
    brute_force_diameter,
    certified_diameter_bnb,
    ordered_difference_options,
)


def intrinsic_pathway_coordinates(
    contributions: np.ndarray,
    relative_tolerance: float = 1.0e-12,
) -> tuple[np.ndarray, np.ndarray]:
    contributions = np.asarray(contributions, dtype=float)
    if contributions.ndim != 3:
        raise ValueError("contributions must have shape (cells, pathways, dimension)")
    q, p_count, _ = contributions.shape
    generators = np.stack([
        contributions[i, p] - contributions[i, 0]
        for i in range(q)
        for p in range(1, p_count)
    ])
    left_vectors, singular_values, _ = np.linalg.svd(generators, full_matrices=False)
    scale = float(singular_values[0]) if singular_values.size else 0.0
    keep = singular_values > relative_tolerance * scale
    generator_coordinates = left_vectors[:, keep] * singular_values[keep]

    coordinates = np.zeros((q, p_count, int(np.count_nonzero(keep))), dtype=float)
    for i in range(q):
        for p in range(1, p_count):
            coordinates[i, p] = generator_coordinates[i * (p_count - 1) + p - 1]
    return coordinates, singular_values


__all__ = [
    "DiameterCertificate",
    "brute_force_diameter",
    "certified_diameter_bnb",
    "intrinsic_pathway_coordinates",
    "ordered_difference_options",
]

