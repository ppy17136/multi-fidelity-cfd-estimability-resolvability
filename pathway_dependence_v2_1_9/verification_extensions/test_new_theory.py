"""Exact integer enumeration for the diameter bridge and endpoint example.

These tests verify finite instances, not the proof or approximation complexity.
"""
import itertools
from pathlib import Path
import sys
import unittest
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))


def partitions(q, edges):
    for mask in range(1 << len(edges)):
        parent = list(range(q))
        def root(i):
            while parent[i] != i:
                i = parent[i]
            return i
        for k, (i, j, _) in enumerate(edges):
            if not mask & (1 << k):
                parent[root(i)] = root(j)
        groups = {}
        for i in range(q):
            groups.setdefault(root(i), []).append(i)
        cost = sum(w for k, (_, _, w) in enumerate(edges) if mask & (1 << k))
        yield cost, list(groups.values())


def points(Y, blocks):
    values = np.zeros((1, Y.shape[-1]), dtype=np.int64)
    for block in blocks:
        values = (values[:, None, :] + Y[block].sum(axis=0)[None, :, :]).reshape(-1, Y.shape[-1])
    return np.unique(values, axis=0)


def squared_diameter(values):
    differences = values[:, None, :] - values[None, :, :]
    return int(np.sum(differences**2, axis=2).max())


def label_cost(z, edges):
    lists = [list(itertools.product(np.flatnonzero(r == r.max()),
                                    np.flatnonzero(r == r.min()))) for r in z]
    return min(sum(w for i, j, w in edges if assignment[i] != assignment[j])
               for assignment in itertools.product(*lists))


def bridge_costs(Y, edges):
    independent = points(Y, [[i] for i in range(len(Y))])
    d2 = squared_diameter(independent)
    if d2 == 0:
        raise ValueError("The positive independent-diameter premise is required")
    lhs = min(cost for cost, blocks in partitions(len(Y), edges)
              if squared_diameter(points(Y, blocks)) == d2)
    differences = independent[:, None, :] - independent[None, :, :]
    maximizers = differences[np.sum(differences**2, axis=2) == d2]
    # Integer, unnormalized directions preserve extrema without tie tolerances.
    directions = np.unique(maximizers, axis=0)
    rhs = min(label_cost(Y @ u, edges) for u in directions)
    return lhs, rhs


class NewTheoryTests(unittest.TestCase):
    def setUp(self):
        self.A = np.array([[5, 1, 4], [-1, -4, 0], [-1, 2, 5]], dtype=np.int64)
        self.B = np.array([[-3, 4, 3], [1, 5, 1], [1, 0, -2]], dtype=np.int64)
        self.edges = [(0, 1, 1), (1, 2, 1)]

    def test_same_endpoint_hulls_different_cost(self):
        costs = []
        for z in (self.A, self.B):
            rows = [(c, points(z[:, :, None], b).ravel()) for c, b in partitions(3, self.edges)]
            self.assertEqual((rows[0][1].min(), rows[0][1].max()), (-1, 9))
            self.assertEqual((rows[-1][1].min(), rows[-1][1].max()), (-4, 10))
            costs.append(min(c for c, values in rows if np.ptp(values) == 14))
        self.assertEqual(costs, [2, 1])  # Normalize edge weights by 2: 1 and 1/2.

    def test_equal_hulls_do_not_mean_equal_finite_sets(self):
        a = points(self.A[:, :, None], [[0, 1, 2]]).ravel()
        b = points(self.B[:, :, None], [[0, 1, 2]]).ravel()
        self.assertFalse(np.array_equal(a, b))
        self.assertEqual(a.tolist(), [-1, 3, 9])
        self.assertEqual(b.tolist(), [-1, 2, 9])

    def test_decision_loss_precedes_full_width(self):
        rows = [(c, points(self.A[:, :, None], b).ravel()) for c, b in partitions(3, self.edges)]
        crossing = min(c for c, values in rows if values.min() <= -2)
        full = min(c for c, values in rows if np.ptp(values) == 14)
        self.assertEqual((crossing, full), (1, 2))

    def test_integer_diameter_bridge_on_paths_and_cycles(self):
        rng = np.random.default_rng(2026090604)
        for index in range(30):
            Y = rng.integers(-3, 4, size=(3, 3, 2), dtype=np.int64)
            edges = [(0, 1, 1), (1, 2, 2)]
            if index % 2:
                edges.append((0, 2, 4))
            with self.subTest(instance=index):
                left, right = bridge_costs(Y, edges)
                self.assertEqual(left, right)

    def test_multiple_maximizing_directions(self):
        Y = np.array([[[0, 0], [1, 0]], [[0, 0], [0, 1]]], dtype=np.int64)
        self.assertEqual(bridge_costs(Y, [(0, 1, 3)]), (0, 0))

    def test_zero_diameter_excluded(self):
        with self.assertRaises(ValueError):
            bridge_costs(np.zeros((3, 3, 2), dtype=np.int64), self.edges)

    def test_translation_scaling_and_direction_reversal(self):
        Y = np.stack((self.A, self.B), axis=-1)
        reference = bridge_costs(Y, self.edges)
        shift = np.array([[3, 4], [-2, 7], [5, -8]])[:, None, :]
        self.assertEqual(reference, bridge_costs(3*Y + shift, self.edges))
        z = Y @ np.array([2, -1])
        self.assertEqual(label_cost(z, self.edges), label_cost(-z, self.edges))

    def test_existing_forest_dp_against_integer_labels(self):
        from graph_constrained_coupling import pathway_pair_compatibility, tree_graph_relaxation
        rng = np.random.default_rng(2026090605)
        for _ in range(30):
            z = rng.integers(-2, 3, size=(4, 3))
            edges = [(0, 1, 1), (1, 2, 2), (1, 3, 3)]
            cert = tree_graph_relaxation(pathway_pair_compatibility(z), edges)
            self.assertEqual(cert.cost, label_cost(z, edges))


if __name__ == "__main__":
    unittest.main(verbosity=2)
