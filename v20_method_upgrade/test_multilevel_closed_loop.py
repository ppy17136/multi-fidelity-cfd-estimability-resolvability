from __future__ import annotations

import runpy
import unittest
from pathlib import Path

import numpy as np


M = runpy.run_path(
    str(Path(__file__).with_name("04b_run_multilevel_closed_loop_acquisition.py")),
    run_name="multilevel_closed_loop_test_target",
)


class MultilevelClosedLoopTests(unittest.TestCase):
    def test_irrelevant_missing_cell_does_not_poison_contrast(self) -> None:
        estimates = np.asarray([10.0, 9.0, 8.5, np.nan, 8.0, 12.0])
        mapping = np.asarray([0.1, 0.1, 0.1, np.inf, 0.1, 0.1])
        repeatability = np.asarray([0.05, 0.05, 0.05, np.inf, 0.05, 0.05])
        self.assertAlmostEqual(M["contrast_estimate"](estimates), 5.0)
        self.assertAlmostEqual(M["contrast_bound"](mapping, repeatability), 0.6)

    def test_declared_bound_prevents_false_present(self) -> None:
        estimates = np.asarray([10.2, 8.8, 8.5, np.nan, 7.8, 8.8])
        mapping = np.asarray([0.2, 0.2, 0.1, np.inf, 0.2, 0.2])
        repeatability = np.asarray([0.1, 0.1, 0.1, np.inf, 0.1, 0.1])
        self.assertNotEqual(
            M["decision_state"](estimates, mapping, repeatability),
            "effect_present",
        )

    def test_policy_names_are_unique_and_duplicate_alias_is_absent(self) -> None:
        policies = M["POLICIES"]
        self.assertEqual(len(policies), len(set(policies)))
        self.assertNotIn("c_optimal", policies)

    def test_shared_scenario_reproduces_a_policy_exactly(self) -> None:
        scenario = M["make_scenario"](np.random.default_rng(11))
        first = M["run_one"]("support_then_c_optimal", 2.0, scenario=scenario)
        second = M["run_one"]("support_then_c_optimal", 2.0, scenario=scenario)
        self.assertEqual(first, second)

    def test_fixed_seed_paths_are_valid(self) -> None:
        for delta in (2.0, 5.0):
            result = M["run_one"](
                "proposed_decision_targeted_switch",
                delta,
                np.random.default_rng(1),
            )
            self.assertTrue(result["valid_decision"])
            self.assertFalse(result["false_present"])
            self.assertFalse(result["false_absent"])


if __name__ == "__main__":
    unittest.main()
