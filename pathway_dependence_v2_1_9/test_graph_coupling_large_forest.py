import unittest

from graph_constrained_coupling import cut_cost, tree_graph_relaxation


class LargeForestImplementationTests(unittest.TestCase):
    def test_long_chain_does_not_depend_on_python_recursion_depth(self):
        q = 4096
        all_labels = ((0, 0), (0, 1), (1, 0), (1, 1))
        allowed = tuple(
            (((cell // 257) % 2, 0),) if cell % 257 == 0 else all_labels
            for cell in range(q)
        )
        edges = tuple((cell - 1, cell, 1.0) for cell in range(1, q))
        certificate = tree_graph_relaxation(allowed, edges)
        self.assertEqual(certificate.cost, cut_cost(certificate.labels, edges))
        self.assertTrue(all(
            label in cell_labels
            for label, cell_labels in zip(certificate.labels, allowed)
        ))

    def test_cycle_is_rejected(self):
        allowed = (((0, 0),),) * 3
        with self.assertRaisesRegex(ValueError, "requires a forest"):
            tree_graph_relaxation(allowed, ((0, 1, 1.0), (1, 2, 1.0), (2, 0, 1.0)))


if __name__ == "__main__":
    unittest.main()
