import unittest

import numpy as np

from contrast_resolution import (
    contrast_standard_deviation,
    deterministic_contrast_bound,
    deterministic_resolution_decision,
)


class ContrastResolutionTests(unittest.TestCase):
    def test_four_cell_triangle_bound(self):
        coefficients = [1, -1, -1, 1]
        bounds = [0.25, 0.25, 0.25, 0.25]
        self.assertAlmostEqual(
            deterministic_contrast_bound(coefficients, bounds), 1.0
        )

    def test_correlated_contrast_sd(self):
        coefficients = np.asarray([1.0, -1.0, -1.0, 1.0])
        covariance = np.eye(4) * 0.25**2
        self.assertAlmostEqual(
            contrast_standard_deviation(coefficients, covariance), 0.5
        )

    def test_conservative_decision(self):
        self.assertFalse(
            deterministic_resolution_decision(3.5, 1.0, 3.0).resolved
        )
        self.assertTrue(
            deterministic_resolution_decision(5.0, 1.0, 3.0).resolved
        )


if __name__ == "__main__":
    unittest.main()

