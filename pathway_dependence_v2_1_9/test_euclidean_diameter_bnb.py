"""Correctness tests for the certified Euclidean diameter branch-and-bound."""

from __future__ import annotations

import unittest

import numpy as np

from euclidean_diameter_bnb import brute_force_diameter, certified_diameter_bnb


class EuclideanDiameterBranchAndBoundTests(unittest.TestCase):
    def test_random_small_instances_match_complete_enumeration(self):
        rng = np.random.default_rng(20260904)
        for q in range(1, 6):
            for dimension in range(1, 5):
                for _ in range(10):
                    options = []
                    for _cell in range(q):
                        raw = rng.normal(size=(4, dimension))
                        raw[0] = 0.0
                        options.append(raw)
                    expected, _ = brute_force_diameter(options)
                    certificate = certified_diameter_bnb(options)
                    self.assertAlmostEqual(certificate.squared_diameter, expected, places=10)

    def test_zero_instance(self):
        options = [np.zeros((4, 3)), np.zeros((4, 3))]
        certificate = certified_diameter_bnb(options)
        self.assertEqual(certificate.squared_diameter, 0.0)


if __name__ == "__main__":
    unittest.main()
