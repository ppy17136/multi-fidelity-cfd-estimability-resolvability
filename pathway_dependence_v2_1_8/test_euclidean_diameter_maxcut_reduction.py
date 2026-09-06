"""Executable checks for the Max-Cut reduction to Euclidean diameter."""

from __future__ import annotations

import itertools
import unittest

import numpy as np


def oriented_incidence(vertex_count: int, edges: list[tuple[int, int]]) -> np.ndarray:
    matrix = np.zeros((len(edges), vertex_count), dtype=float)
    for row, (left, right) in enumerate(edges):
        matrix[row, left] = 1.0
        matrix[row, right] = -1.0
    return matrix


def max_cut_size(vertex_count: int, edges: list[tuple[int, int]]) -> int:
    return max(
        sum(labels[i] != labels[j] for i, j in edges)
        for labels in itertools.product((-1, 1), repeat=vertex_count)
    )


def two_pathway_diameter(incidence: np.ndarray) -> float:
    # Cell i contributes either +b_i/2 or -b_i/2, where b_i is column i.
    points = np.stack([
        incidence @ np.asarray(labels, dtype=float) / 2.0
        for labels in itertools.product((-1, 1), repeat=incidence.shape[1])
    ])
    delta = points[:, None, :] - points[None, :, :]
    return float(np.sqrt(np.max(np.sum(delta * delta, axis=-1))))


class MaxCutDiameterReductionTests(unittest.TestCase):
    def assert_graph_identity(self, vertex_count, edges):
        incidence = oriented_incidence(vertex_count, edges)
        diameter_squared = two_pathway_diameter(incidence) ** 2
        self.assertAlmostEqual(diameter_squared, 4.0 * max_cut_size(vertex_count, edges))

    def test_named_graphs(self):
        self.assert_graph_identity(3, [(0, 1), (1, 2), (2, 0)])
        self.assert_graph_identity(4, [(0, 1), (1, 2), (2, 3), (3, 0)])
        self.assert_graph_identity(5, [(0, 1), (0, 2), (0, 3), (0, 4)])

    def test_all_graphs_through_four_vertices(self):
        for vertex_count in range(2, 5):
            possible = list(itertools.combinations(range(vertex_count), 2))
            for edge_mask in range(1 << len(possible)):
                edges = [
                    edge for bit, edge in enumerate(possible)
                    if edge_mask & (1 << bit)
                ]
                self.assert_graph_identity(vertex_count, edges)


if __name__ == "__main__":
    unittest.main()
