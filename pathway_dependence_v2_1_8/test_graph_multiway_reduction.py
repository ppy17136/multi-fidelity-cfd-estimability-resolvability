"""Exhaustive bug check for the three-terminal Multiway-Cut reduction."""

from __future__ import annotations

import itertools
import unittest

import numpy as np

from graph_constrained_coupling import (
    exact_graph_relaxation_bruteforce,
    pathway_pair_compatibility,
)


class MultiwayReductionTests(unittest.TestCase):
    def test_all_simple_graphs_on_five_vertices(self):
        # Vertices 0, 1, 2 are fixed terminals; vertices 3 and 4 may take any
        # terminal label. All 2^10 simple graphs on five labelled vertices are
        # checked against direct multiway-cut enumeration.
        q = 5
        z = np.zeros((q, 4), dtype=float)
        for terminal in range(3):
            z[terminal, 1:] = 0.5
            z[terminal, terminal + 1] = 1.0
        z[3:, 1:] = 1.0
        allowed = pathway_pair_compatibility(z)
        terminal_labels = tuple(allowed[i][0] for i in range(3))
        free_labels = tuple((pathway, 0) for pathway in (1, 2, 3))
        self.assertEqual(tuple(allowed[3]), free_labels)
        self.assertEqual(tuple(allowed[4]), free_labels)

        possible_edges = list(itertools.combinations(range(q), 2))
        for graph_mask in range(1 << len(possible_edges)):
            edges = [
                (left, right, 1.0)
                for bit, (left, right) in enumerate(possible_edges)
                if graph_mask & (1 << bit)
            ]
            coupling = exact_graph_relaxation_bruteforce(allowed, edges)
            direct = min(
                sum(
                    weight
                    for left, right, weight in edges
                    if labels[left] != labels[right]
                )
                for free in itertools.product(free_labels, repeat=2)
                for labels in [terminal_labels + free]
            )
            self.assertEqual(coupling.cost, direct)


if __name__ == "__main__":
    unittest.main()
