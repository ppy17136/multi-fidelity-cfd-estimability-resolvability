"""Cross-check MILP coupling certificates against independent exact solvers."""

from __future__ import annotations

import itertools
import unittest

import numpy as np

from graph_constrained_coupling import (
    exact_decision_relaxation_frontier,
    exact_graph_relaxation_bruteforce,
    pathway_pair_compatibility,
)
from graph_coupling_milp import (
    milp_decision_destroying_relaxation,
    milp_directional_full_width,
)


class GraphCouplingMilpTests(unittest.TestCase):
    def test_directional_milp_matches_bruteforce(self):
        rng = np.random.default_rng(20260907)
        for q in range(2, 8):
            possible = list(itertools.combinations(range(q), 2))
            for _ in range(12):
                z = rng.integers(0, 5, size=(q, 4)).astype(float)
                allowed = pathway_pair_compatibility(z)
                edges = [
                    (left, right, float(rng.integers(1, 8)))
                    for left, right in possible
                    if rng.random() < 0.45
                ]
                brute = exact_graph_relaxation_bruteforce(allowed, edges)
                certificate = milp_directional_full_width(allowed, edges)
                self.assertEqual(certificate.status, "optimal")
                self.assertAlmostEqual(certificate.objective, brute.cost, places=8)
                self.assertLessEqual(certificate.relative_gap, 1.0e-12)

    def test_decision_milp_matches_bruteforce(self):
        rng = np.random.default_rng(20260908)
        checked = 0
        for q in range(2, 7):
            possible = list(itertools.combinations(range(q), 2))
            for _ in range(30):
                z = rng.integers(-4, 5, size=(q, 3)).astype(float)
                edges = [
                    (left, right, float(rng.integers(1, 5)))
                    for left, right in possible
                    if right == left + 1 or rng.random() < 0.5
                ]
                threshold = float(rng.integers(-3, 4))
                brute = exact_decision_relaxation_frontier(z, edges, threshold)
                if brute.coherent_state == "indeterminate":
                    continue
                state, certificate = milp_decision_destroying_relaxation(z, edges, threshold)
                self.assertEqual(state, brute.coherent_state)
                if brute.minimum_relaxation_cost is None:
                    self.assertIsNone(certificate.objective)
                else:
                    self.assertEqual(certificate.status, "optimal")
                    self.assertAlmostEqual(
                        certificate.objective,
                        brute.minimum_relaxation_cost,
                        places=8,
                    )
                checked += 1
        self.assertGreater(checked, 10)

    def test_decision_milp_rejects_disconnected_graph(self):
        z = np.array([[2.0, 0.0], [2.0, 0.0], [2.0, 0.0]])
        with self.assertRaisesRegex(ValueError, "must be connected"):
            milp_decision_destroying_relaxation(z, [(0, 1, 1.0)], 1.0)


if __name__ == "__main__":
    unittest.main()
