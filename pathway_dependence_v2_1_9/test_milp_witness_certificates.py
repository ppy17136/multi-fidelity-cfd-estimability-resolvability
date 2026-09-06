import itertools
import unittest

import numpy as np

from graph_constrained_coupling import cut_cost, pathway_pair_compatibility
from graph_coupling_milp import (
    milp_decision_destroying_relaxation,
    milp_directional_full_width,
)


def component_partition(q, labels, edges):
    adjacency = [[] for _ in range(q)]
    for left, right, _ in edges:
        if labels[left] == labels[right]:
            adjacency[left].append(right)
            adjacency[right].append(left)
    unseen = set(range(q))
    components = []
    while unseen:
        root = min(unseen)
        unseen.remove(root)
        stack = [root]
        component = []
        while stack:
            node = stack.pop()
            component.append(node)
            for neighbour in adjacency[node]:
                if neighbour in unseen:
                    unseen.remove(neighbour)
                    stack.append(neighbour)
        components.append(tuple(component))
    return tuple(components)


def directional_width(projections, partition):
    return float(
        sum(
            np.ptp(projections[list(component)].sum(axis=0))
            for component in partition
        )
    )


class MilpWitnessCertificateTests(unittest.TestCase):
    def test_directional_witness_recomputes_cost_and_full_width(self):
        rng = np.random.default_rng(20260911)
        for q in range(2, 7):
            possible = list(itertools.combinations(range(q), 2))
            for _ in range(8):
                projections = rng.integers(0, 5, size=(q, 4)).astype(float)
                edges = [
                    (left, right, float(rng.integers(1, 8)))
                    for left, right in possible
                    if right == left + 1 or rng.random() < 0.35
                ]
                allowed = pathway_pair_compatibility(projections)
                certificate = milp_directional_full_width(allowed, edges)
                self.assertEqual(certificate.status, "optimal")
                self.assertTrue(
                    all(label in allowed[i] for i, label in enumerate(certificate.labels))
                )
                self.assertAlmostEqual(
                    cut_cost(certificate.labels, edges), certificate.objective, places=8
                )
                partition = component_partition(q, certificate.labels, edges)
                self.assertEqual(
                    directional_width(projections, partition),
                    float(np.ptp(projections, axis=1).sum()),
                )

    def test_decision_witness_recomputes_crossing_and_cost(self):
        edges = [(0, 1, 1.0)]
        cases = (
            (np.array([[4.0, 0.0], [0.0, 4.0]]), 3.0, "resolved_above"),
            (np.array([[0.0, 4.0], [4.0, 0.0]]), 5.0, "resolved_below"),
        )
        for contributions, threshold, expected_state in cases:
            with self.subTest(state=expected_state):
                state, certificate = milp_decision_destroying_relaxation(
                    contributions, edges, threshold
                )
                self.assertEqual(state, expected_state)
                self.assertEqual(certificate.status, "optimal")
                self.assertAlmostEqual(
                    cut_cost(certificate.labels, edges), certificate.objective, places=8
                )
                recomputed = float(
                    sum(contributions[i, label] for i, label in enumerate(certificate.labels))
                )
                self.assertEqual(recomputed, certificate.attained_value)
                if state == "resolved_above":
                    self.assertLessEqual(recomputed, threshold)
                else:
                    self.assertGreaterEqual(recomputed, threshold)


if __name__ == "__main__":
    unittest.main()
