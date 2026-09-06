"""Direct set-level checks of the graph-labelling equivalence theorem."""

from __future__ import annotations

import itertools
import unittest

import numpy as np

from graph_constrained_coupling import (
    exact_graph_relaxation_bruteforce,
    pathway_pair_compatibility,
)


def component_partition(q, retained_edges):
    adjacency = [[] for _ in range(q)]
    for left, right in retained_edges:
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
        components.append(tuple(sorted(component)))
    return tuple(components)


def partition_width(projections, partition):
    width = 0.0
    for component in partition:
        pathway_totals = projections[list(component)].sum(axis=0)
        width += float(pathway_totals.max() - pathway_totals.min())
    return width


def direct_minimum_full_width_cost(projections, edges):
    q = projections.shape[0]
    independent_width = float(np.ptp(projections, axis=1).sum())
    best = np.inf
    for mask in range(1 << len(edges)):
        retained = []
        cost = 0.0
        for edge_index, (left, right, weight) in enumerate(edges):
            if mask & (1 << edge_index):
                cost += weight
            else:
                retained.append((left, right))
        partition = component_partition(q, retained)
        if partition_width(projections, partition) == independent_width:
            best = min(best, cost)
    return float(best)


class CommonPartitionWidthEquivalenceTests(unittest.TestCase):
    def test_direct_partition_enumeration_matches_label_optimum(self):
        rng = np.random.default_rng(20260909)
        for q in range(2, 6):
            possible = list(itertools.combinations(range(q), 2))
            for _ in range(12):
                # Integer projections deliberately create tied extrema while
                # retaining exact comparisons in the direct set calculation.
                projections = rng.integers(0, 4, size=(q, 4)).astype(float)
                edges = [
                    (left, right, float(rng.integers(1, 7)))
                    for left, right in possible
                    if right == left + 1 or rng.random() < 0.35
                ]
                direct = direct_minimum_full_width_cost(projections, edges)
                labels = pathway_pair_compatibility(
                    projections,
                    relative_tolerance=0.0,
                    absolute_tolerance=0.0,
                )
                labelled = exact_graph_relaxation_bruteforce(labels, edges)
                self.assertEqual(labelled.cost, direct)


if __name__ == "__main__":
    unittest.main()
