import unittest

import numpy as np

from graph_constrained_coupling import pathway_pair_compatibility


class ExactExtremizerListTests(unittest.TestCase):
    def test_default_lists_use_exact_extrema(self):
        projections = np.array([[1.0, 1.0 - 5.0e-13, 0.0]])
        self.assertEqual(pathway_pair_compatibility(projections), (((0, 2),),))

    def test_near_tie_expansion_must_be_explicit(self):
        projections = np.array([[1.0, 1.0 - 5.0e-13, 0.0]])
        labels = pathway_pair_compatibility(
            projections,
            relative_tolerance=1.0e-12,
            absolute_tolerance=0.0,
        )
        self.assertEqual(labels, (((0, 2), (1, 2)),))

    def test_negative_or_nonfinite_tolerances_are_rejected(self):
        projections = np.array([[1.0, 0.0]])
        for relative, absolute in ((-1.0, 0.0), (0.0, -1.0), (np.inf, 0.0)):
            with self.subTest(relative=relative, absolute=absolute):
                with self.assertRaises(ValueError):
                    pathway_pair_compatibility(
                        projections,
                        relative_tolerance=relative,
                        absolute_tolerance=absolute,
                    )


if __name__ == "__main__":
    unittest.main()
