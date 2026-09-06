"""Distance-preservation tests for the direct-SVD intrinsic coordinates."""

from __future__ import annotations

import itertools
import unittest

import numpy as np

from euclidean_diameter_bnb_v2 import intrinsic_pathway_coordinates


class IntrinsicCoordinateTests(unittest.TestCase):
    def test_all_endpoint_distances_are_preserved(self):
        rng = np.random.default_rng(20260904)
        for q, p_count, ambient in [(2, 3, 7), (3, 4, 11), (4, 3, 5)]:
            contributions = rng.normal(size=(q, p_count, ambient))
            coordinates, _ = intrinsic_pathway_coordinates(contributions)
            choices = list(itertools.product(range(p_count), repeat=q))
            original = np.stack([
                sum(contributions[i, choice[i]] for i in range(q))
                for choice in choices
            ])
            reduced = np.stack([
                sum(coordinates[i, choice[i]] for i in range(q))
                for choice in choices
            ])
            original_delta = original[:, None, :] - original[None, :, :]
            reduced_delta = reduced[:, None, :] - reduced[None, :, :]
            original_distance2 = np.sum(original_delta * original_delta, axis=-1)
            reduced_distance2 = np.sum(reduced_delta * reduced_delta, axis=-1)
            np.testing.assert_allclose(reduced_distance2, original_distance2, rtol=1e-11, atol=1e-11)


if __name__ == "__main__":
    unittest.main()

