from __future__ import annotations

import unittest

import numpy as np

from mapping_set_metrics import set_metrics


class MappingSetMetricTests(unittest.TestCase):
    def test_minimum_norm_and_exact_diameter(self):
        vectors = np.asarray([[3.0, 4.0], [0.0, 5.0], [6.0, 8.0]])
        result = set_metrics(vectors)
        self.assertAlmostEqual(result["minimum_norm"], 5.0)
        self.assertAlmostEqual(
            result["diameter"], float(np.sqrt(45.0)), places=6
        )
        self.assertAlmostEqual(
            result["ratio"], 5.0 / float(np.sqrt(45.0)), places=6
        )

    def test_superset_cannot_raise_minimum_or_reduce_diameter(self):
        coherent = np.asarray([[1.0, 0.0], [0.0, 1.0]])
        recombined = np.asarray(
            [[1.0, 0.0], [0.0, 1.0], [2.0, 2.0], [-1.0, -1.0]]
        )
        a = set_metrics(coherent)
        b = set_metrics(recombined)
        self.assertLessEqual(b["minimum_norm"], a["minimum_norm"])
        self.assertGreaterEqual(b["diameter"], a["diameter"])


if __name__ == "__main__":
    unittest.main()
