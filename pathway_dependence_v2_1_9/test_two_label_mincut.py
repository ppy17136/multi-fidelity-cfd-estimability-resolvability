"""Checks the polynomial two-label special case against complete enumeration."""

from __future__ import annotations

import itertools
import unittest

import numpy as np

from graph_constrained_coupling import (
    exact_graph_relaxation_bruteforce,
    two_label_graph_relaxation,
)


class TwoLabelMinCutTests(unittest.TestCase):
    def test_random_instances_match_bruteforce(self):
        rng = np.random.default_rng(20260910)
        labels = ((0, 0), (1, 1))
        choices = ((labels[0],), (labels[1],), labels)
        for q in range(2, 9):
            possible = list(itertools.combinations(range(q), 2))
            for _ in range(30):
                allowed = tuple(choices[int(rng.integers(0, 3))] for _ in range(q))
                edges = [
                    (left, right, float(rng.integers(1, 8)))
                    for left, right in possible
                    if rng.random() < 0.4
                ]
                exact = exact_graph_relaxation_bruteforce(allowed, edges)
                cut = two_label_graph_relaxation(allowed, edges)
                self.assertAlmostEqual(cut.cost, exact.cost)

    def test_one_label_and_unopposed_forcing_have_zero_cost(self):
        first, second = (0, 0), (1, 1)
        one = two_label_graph_relaxation(((first,), (first,)), [(0, 1, 3.0)])
        self.assertEqual(one.cost, 0.0)
        free = two_label_graph_relaxation(((first,), (first, second)), [(0, 1, 3.0)])
        self.assertEqual(free.cost, 0.0)
        self.assertEqual(free.labels, (first, first))

    def test_three_labels_are_rejected(self):
        with self.assertRaisesRegex(ValueError, "at most two labels"):
            two_label_graph_relaxation(
                (((0, 0),), ((1, 1),), ((2, 2),)),
                [(0, 1, 1.0), (1, 2, 1.0)],
            )


if __name__ == "__main__":
    unittest.main()
