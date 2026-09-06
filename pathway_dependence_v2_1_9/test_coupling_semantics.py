import itertools
import unittest

import numpy as np


class CouplingSemanticsTests(unittest.TestCase):
    def test_budget_union_and_single_partition_diameters_can_differ(self):
        # Three cells on a path, three pathways, and budget one. The two
        # one-edge cuts define different two-block sets. Their robust union has
        # scalar diameter 14, whereas no one common partition exceeds 13.
        z = np.asarray(((5, 1, 4), (-1, -4, 0), (-1, 2, 5)), dtype=float)
        partitions = (((0,), (1, 2)), ((0, 1), (2,)))
        partition_values = []
        for partition in partitions:
            values = []
            for component_labels in itertools.product(range(3), repeat=2):
                labels = [None] * 3
                for component, pathway in zip(partition, component_labels):
                    for cell in component:
                        labels[cell] = pathway
                values.append(sum(z[cell, labels[cell]] for cell in range(3)))
            partition_values.append(values)
        union = partition_values[0] + partition_values[1]
        union_diameter = max(union) - min(union)
        largest_single_partition = max(max(values) - min(values)
                                       for values in partition_values)
        self.assertEqual(union_diameter, 14.0)
        self.assertEqual(largest_single_partition, 13.0)
        self.assertGreater(union_diameter, largest_single_partition)


if __name__ == "__main__":
    unittest.main()
