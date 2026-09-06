"""Independent small-instance tests for the directional coupling-cover identity."""

from __future__ import annotations

import unittest

import numpy as np

from run_directional_coupling_cover_certificate import exact_cover, greedy_cover


def partitions(items):
    items = list(items)
    if not items:
        yield []
        return
    first = items[0]
    for rest in partitions(items[1:]):
        yield [[first], *[list(block) for block in rest]]
        for j in range(len(rest)):
            candidate = [list(block) for block in rest]
            candidate[j] = [first, *candidate[j]]
            yield candidate


def directional_width(z, partition):
    total = 0.0
    for block in partition:
        sums = z[block].sum(axis=0)
        total += float(sums.max() - sums.min())
    return total


def independent_width(z):
    return float(np.sum(z.max(axis=1) - z.min(axis=1)))


def brute_kappa(z):
    target = independent_width(z)
    answer = z.shape[0]
    for partition in partitions(range(z.shape[0])):
        if np.isclose(directional_width(z, partition), target, rtol=1e-11, atol=1e-11):
            answer = min(answer, len(partition))
    return answer


def pathway_pair_masks(z):
    q, p_count = z.shape
    maxima = z.max(axis=1)
    minima = z.min(axis=1)
    masks = {}
    for a in range(p_count):
        for b in range(p_count):
            mask = 0
            for i in range(q):
                if np.isclose(z[i, a], maxima[i]) and np.isclose(z[i, b], minima[i]):
                    mask |= 1 << i
            if mask:
                masks[(a, b)] = mask
    return masks


class DirectionalCouplingCoverTests(unittest.TestCase):
    def assert_cover_identity(self, z):
        q = z.shape[0]
        masks = pathway_pair_masks(z)
        exact = exact_cover(masks, (1 << q) - 1)
        self.assertEqual(len(exact), brute_kappa(z))

    def test_random_unique_extrema(self):
        rng = np.random.default_rng(20260904)
        for q in range(1, 7):
            for p_count in range(2, 5):
                for _ in range(10):
                    self.assert_cover_identity(rng.normal(size=(q, p_count)))

    def test_random_tied_extrema(self):
        rng = np.random.default_rng(20260905)
        for q in range(1, 7):
            for p_count in range(2, 5):
                for _ in range(10):
                    z = rng.integers(0, 4, size=(q, p_count)).astype(float)
                    self.assert_cover_identity(z)

    def test_set_cover_reduction_instance(self):
        subsets = [
            {0, 1},
            {1, 2},
            {2, 3},
            {0, 3},
        ]
        q = 4
        z = np.full((q, len(subsets) + 1), 0.5)
        z[:, 0] = 0.0
        for j, subset in enumerate(subsets, start=1):
            for i in subset:
                z[i, j] = 1.0
        masks = pathway_pair_masks(z)
        exact = exact_cover(masks, (1 << q) - 1)
        self.assertEqual(len(exact), 2)
        self.assertEqual(brute_kappa(z), 2)

    def test_unique_extrema_closed_form(self):
        z = np.array([
            [3.0, 0.0, 1.0],
            [4.0, 0.0, 2.0],
            [0.0, 5.0, 1.0],
            [0.0, 6.0, 2.0],
        ])
        pairs = {
            (int(np.argmax(row)), int(np.argmin(row)))
            for row in z
        }
        masks = pathway_pair_masks(z)
        exact = exact_cover(masks, (1 << z.shape[0]) - 1)
        self.assertEqual(len(exact), len(pairs))

    def test_greedy_never_returns_an_invalid_cover(self):
        rng = np.random.default_rng(20260906)
        for _ in range(100):
            z = rng.integers(0, 5, size=(6, 4)).astype(float)
            masks = pathway_pair_masks(z)
            greedy = greedy_cover(masks, (1 << z.shape[0]) - 1)
            covered = 0
            for pair in greedy:
                covered |= masks[pair]
            self.assertEqual(covered, (1 << z.shape[0]) - 1)


if __name__ == "__main__":
    unittest.main()
