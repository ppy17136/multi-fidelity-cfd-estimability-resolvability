import unittest

import numpy as np

from support_repair import (
    binary_factorial_design,
    contrast_for_terms,
    exact_minimum_cost_repair,
    greedy_support_repair,
    is_estimable,
    rowspace_residual,
)


class SupportRepairTests(unittest.TestCase):
    def test_missing_2x2_interaction_requires_fourth_corner(self):
        design, _, terms = binary_factorial_design(2)
        contrast = contrast_for_terms(terms, [(0, 1)])
        self.assertFalse(is_estimable(design[:3], contrast))
        self.assertTrue(is_estimable(design, contrast))
        exact = exact_minimum_cost_repair(
            design[:3], design[3:], [7.0], contrast
        )
        self.assertTrue(exact.feasible)
        self.assertEqual(exact.selected, [0])
        self.assertAlmostEqual(exact.added_cost, 7.0)

    def test_multiple_contrasts(self):
        design, _, terms = binary_factorial_design(
            3,
            included_terms=[(), (0,), (1,), (2,), (0, 1), (0, 2)],
        )
        contrast = contrast_for_terms(terms, [(0, 1), (0, 2)])
        observed = design[[0, 1, 2]]
        candidates = design[[3, 4, 5, 6, 7]]
        exact = exact_minimum_cost_repair(
            observed, candidates, [1, 2, 3, 4, 5], contrast
        )
        greedy = greedy_support_repair(
            observed, candidates, [1, 2, 3, 4, 5], contrast
        )
        self.assertTrue(exact.feasible)
        self.assertTrue(greedy.feasible)
        self.assertLessEqual(exact.added_cost, greedy.added_cost)

    def test_projection_residual_is_scale_invariant_above_unit_norm(self):
        x = np.asarray([[1.0, 0.0], [1.0, 1.0]])
        c = np.asarray([0.0, 1.0])
        self.assertAlmostEqual(rowspace_residual(x, c), 0.0)
        self.assertAlmostEqual(rowspace_residual(x, 10.0 * c), 0.0)


if __name__ == "__main__":
    unittest.main()
