import unittest

import numpy as np

from graph_constrained_coupling import (
    exact_decision_relaxation_frontier,
    exact_graph_relaxation_bruteforce,
    pathway_pair_compatibility,
    tree_graph_relaxation,
    two_label_graph_relaxation,
)
from graph_coupling_milp import (
    milp_decision_destroying_relaxation,
    milp_directional_full_width,
)


class DegenerateBoundaryCaseTests(unittest.TestCase):
    def test_single_vertex_single_pathway_has_zero_relaxation_cost(self):
        allowed = (((0, 0),),)
        self.assertEqual(exact_graph_relaxation_bruteforce(allowed, []).cost, 0.0)
        self.assertEqual(tree_graph_relaxation(allowed, []).cost, 0.0)
        self.assertEqual(two_label_graph_relaxation(allowed, []).cost, 0.0)
        certificate = milp_directional_full_width(allowed, [])
        self.assertEqual(certificate.status, "optimal")
        self.assertEqual(certificate.objective, 0.0)

    def test_zero_direction_all_ties_still_has_zero_cost(self):
        projections = np.zeros((3, 3))
        allowed = pathway_pair_compatibility(projections)
        self.assertTrue(all(len(labels) == 9 for labels in allowed))
        edges = [(0, 1, 2.0), (1, 2, 3.0)]
        self.assertEqual(exact_graph_relaxation_bruteforce(allowed, edges).cost, 0.0)
        self.assertEqual(tree_graph_relaxation(allowed, edges).cost, 0.0)
        self.assertEqual(milp_directional_full_width(allowed, edges).objective, 0.0)

    def test_disconnected_directional_graph_uses_component_baseline(self):
        # With no retained edges, every vertex is already an independent block.
        allowed = (((0, 0),), ((1, 1),), ((2, 2),))
        self.assertEqual(exact_graph_relaxation_bruteforce(allowed, []).cost, 0.0)
        self.assertEqual(tree_graph_relaxation(allowed, []).cost, 0.0)
        self.assertEqual(milp_directional_full_width(allowed, []).objective, 0.0)

    def test_zero_weight_edge_is_supported_as_code_extension(self):
        allowed = (((0, 0),), ((1, 1),))
        edge = [(0, 1, 0.0)]
        self.assertEqual(exact_graph_relaxation_bruteforce(allowed, edge).cost, 0.0)
        self.assertEqual(tree_graph_relaxation(allowed, edge).cost, 0.0)
        self.assertEqual(milp_directional_full_width(allowed, edge).objective, 0.0)

    def test_threshold_equality_is_indeterminate(self):
        contributions = np.array([[1.0, 0.0], [1.0, 0.0]])
        edges = [(0, 1, 1.0)]
        for threshold in (0.0, 2.0):
            with self.subTest(threshold=threshold):
                direct = exact_decision_relaxation_frontier(
                    contributions, edges, threshold
                )
                state, certificate = milp_decision_destroying_relaxation(
                    contributions, edges, threshold
                )
                self.assertEqual(direct.coherent_state, "indeterminate")
                self.assertEqual(state, "indeterminate")
                self.assertEqual(certificate.status, "initially_indeterminate")

    def test_one_pathway_can_make_crossing_infeasible(self):
        contributions = np.array([[2.0], [3.0]])
        edges = [(0, 1, 1.0)]
        direct = exact_decision_relaxation_frontier(contributions, edges, 4.0)
        state, certificate = milp_decision_destroying_relaxation(
            contributions, edges, 4.0
        )
        self.assertEqual(direct.coherent_state, "resolved_above")
        self.assertIsNone(direct.minimum_relaxation_cost)
        self.assertEqual(state, "resolved_above")
        self.assertIsNone(certificate.objective)

    def test_empty_cell_arrays_are_rejected(self):
        empty = np.empty((0, 2))
        with self.assertRaises(ValueError):
            pathway_pair_compatibility(empty)
        with self.assertRaises(ValueError):
            exact_decision_relaxation_frontier(empty, [], 0.0)
        with self.assertRaises(ValueError):
            milp_decision_destroying_relaxation(empty, [], 0.0)


if __name__ == "__main__":
    unittest.main()
