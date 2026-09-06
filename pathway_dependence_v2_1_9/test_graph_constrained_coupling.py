"""Tests for graph-constrained coupling relaxation and decision frontiers."""

from __future__ import annotations

import unittest

import numpy as np

from graph_constrained_coupling import (
    exact_decision_relaxation_frontier,
    exact_graph_relaxation_bruteforce,
    pathway_pair_compatibility,
    tree_graph_relaxation,
)


class GraphConstrainedCouplingTests(unittest.TestCase):
    def test_tree_dp_matches_bruteforce_with_ties(self):
        rng = np.random.default_rng(20260904)
        for q in range(2, 8):
            for _ in range(30):
                z = rng.integers(0, 4, size=(q, 4)).astype(float)
                allowed = pathway_pair_compatibility(z)
                edges = [(i - 1, i, float(rng.integers(1, 6))) for i in range(1, q)]
                exact = exact_graph_relaxation_bruteforce(allowed, edges)
                tree = tree_graph_relaxation(allowed, edges)
                self.assertAlmostEqual(tree.cost, exact.cost)

    def test_tree_dp_handles_a_forest(self):
        z = np.array([
            [2.0, 0.0, 1.0],
            [3.0, 0.0, 1.0],
            [0.0, 4.0, 1.0],
            [0.0, 5.0, 1.0],
        ])
        allowed = pathway_pair_compatibility(z)
        edges = [(0, 1, 2.0), (2, 3, 3.0)]
        exact = exact_graph_relaxation_bruteforce(allowed, edges)
        tree = tree_graph_relaxation(allowed, edges)
        self.assertAlmostEqual(tree.cost, exact.cost)
        self.assertEqual(tree.cost, 0.0)

    def test_required_split_has_declared_cost(self):
        z = np.array([[2.0, 0.0], [0.0, 2.0]])
        allowed = pathway_pair_compatibility(z)
        certificate = tree_graph_relaxation(allowed, [(0, 1, 7.5)])
        self.assertEqual(certificate.cost, 7.5)

    def test_decision_destroying_cost_from_above(self):
        # Either common pathway gives value 4 (>3), whereas the mixed choice
        # (1, 0) gives zero and requires cutting the edge.
        z = np.array([[4.0, 0.0], [0.0, 4.0]])
        result = exact_decision_relaxation_frontier(z, [(0, 1, 2.5)], 3.0)
        self.assertEqual(result.coherent_state, "resolved_above")
        self.assertEqual(result.minimum_relaxation_cost, 2.5)
        self.assertLessEqual(result.attained_value, 3.0)
        self.assertEqual(result.lower_frontier, ((0.0, 4.0), (2.5, 0.0)))

    def test_decision_destroying_cost_from_below(self):
        # Common labels both give zero (<3); a mixed choice reaches four.
        z = np.array([[4.0, 0.0], [-4.0, 0.0]])
        result = exact_decision_relaxation_frontier(z, [(0, 1, 1.25)], 3.0)
        self.assertEqual(result.coherent_state, "resolved_below")
        self.assertEqual(result.minimum_relaxation_cost, 1.25)
        self.assertGreaterEqual(result.attained_value, 3.0)
        self.assertEqual(result.upper_frontier, ((0.0, 0.0), (1.25, 4.0)))

    def test_initially_indeterminate_has_zero_additional_cost(self):
        z = np.array([[2.0, 0.0], [2.0, 0.0]])
        result = exact_decision_relaxation_frontier(z, [(0, 1, 1.0)], 1.0)
        self.assertEqual(result.coherent_state, "indeterminate")
        self.assertEqual(result.minimum_relaxation_cost, 0.0)

    def test_decision_frontier_rejects_disconnected_graph(self):
        z = np.array([[2.0, 0.0], [2.0, 0.0], [2.0, 0.0]])
        with self.assertRaisesRegex(ValueError, "must be connected"):
            exact_decision_relaxation_frontier(z, [(0, 1, 1.0)], 1.0)


if __name__ == "__main__":
    unittest.main()
