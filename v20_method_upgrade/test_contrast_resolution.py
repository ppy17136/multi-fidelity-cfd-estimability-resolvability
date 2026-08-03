import unittest

import numpy as np

from contrast_resolution import (
    contrast_standard_deviation,
    deterministic_contrast_bound,
    deterministic_resolution_decision,
    probabilistic_resolution_decision,
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

    def test_deterministic_three_state_decision(self):
        self.assertEqual(
            deterministic_resolution_decision(3.5, 1.0, 3.0).state,
            "indeterminate",
        )
        self.assertEqual(
            deterministic_resolution_decision(5.0, 1.0, 3.0).state,
            "effect_present",
        )
        self.assertEqual(
            deterministic_resolution_decision(1.0, 1.0, 3.0).state,
            "below_minimum_effect",
        )

    def test_two_sided_probabilistic_interval(self):
        decision = probabilistic_resolution_decision(
            3.0, 0.5, 3.0, confidence=0.95
        )
        self.assertAlmostEqual(decision.conservative_lower_bound, 2.020018, 5)
        self.assertAlmostEqual(decision.conservative_upper_bound, 3.979982, 5)
        self.assertEqual(decision.state, "indeterminate")


if __name__ == "__main__":
    unittest.main()
