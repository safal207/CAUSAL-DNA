import copy
import json
import unittest
from pathlib import Path

from causal_dna.distributional_robust import (
    DistributionalRobustError,
    DistributionallyRobustAnalyzer,
    _probability_scenarios,
    _reliability_adjusted_multipliers,
)

ROOT = Path(__file__).resolve().parents[1]
PLAN = json.loads((ROOT / "cases" / "CDNA-001.experiment-plan.json").read_text(encoding="utf-8"))
CONFIG = json.loads((ROOT / "cases" / "CDNA-001.distributional-robustness.json").read_text(encoding="utf-8"))
SENSITIVITY = json.loads((ROOT / "cases" / "CDNA-001.sensitivity-analysis.json").read_text(encoding="utf-8"))


class DistributionalRobustTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.analyzer = DistributionallyRobustAnalyzer(PLAN, CONFIG, SENSITIVITY)
        cls.result = cls.analyzer.analyze()

    def test_full_prior_grid_is_used(self):
        self.assertEqual(1140, self.result.prior_points)

    def test_all_experiments_and_stop_are_scored(self):
        ids = {action.action_id for action in self.result.actions}
        self.assertEqual(
            {
                "E1_ARID5B_OCCUPANCY",
                "E2_CUX1_OCCUPANCY",
                "E3_ALLELE_SNATAC",
                "E4_CAPTURE_C",
                "E5_MATCHED_MULTIOMIC",
                "STOP",
            },
            ids,
        )

    def test_probability_uncertainty_preserves_simplex(self):
        scenarios = _probability_scenarios([0.5, 0.5], 0.15)
        self.assertGreater(len(scenarios), 1)
        for scenario in scenarios:
            self.assertAlmostEqual(1.0, sum(scenario))
            self.assertTrue(all(value >= 0 for value in scenario))

    def test_lower_reliability_shrinks_toward_neutral(self):
        original = {"A": 5.0, "B": 0.2}
        weak = _reliability_adjusted_multipliers(original, 0.6)
        strong = _reliability_adjusted_multipliers(original, 1.0)
        self.assertLess(abs(weak["A"] - 1.0), abs(strong["A"] - 1.0))
        self.assertLess(abs(weak["B"] - 1.0), abs(strong["B"] - 1.0))

    def test_minimax_action_really_minimizes_interval_regret(self):
        winner = self.result.action(self.result.minimax_interval_regret_action_id)
        self.assertAlmostEqual(
            winner.max_interval_regret,
            min(action.max_interval_regret for action in self.result.actions),
        )

    def test_maximin_action_really_maximizes_worst_utility(self):
        winner = self.result.action(self.result.maximin_utility_action_id)
        self.assertAlmostEqual(
            winner.worst_case_utility,
            max(action.worst_case_utility for action in self.result.actions),
        )

    def test_result_cannot_claim_causal_status(self):
        self.assertFalse(hasattr(self.result, "cause_found"))
        self.assertFalse(hasattr(self.result, "causal_status"))
        self.assertFalse(hasattr(self.result, "materialized"))

    def test_causal_mutation_field_is_rejected(self):
        config = copy.deepcopy(CONFIG)
        config["cause_found"] = True
        with self.assertRaises(DistributionalRobustError):
            DistributionallyRobustAnalyzer(PLAN, config, SENSITIVITY)

    def test_invalid_uncertainty_bounds_are_rejected(self):
        config = copy.deepcopy(CONFIG)
        config["assay_reliability_min"] = 1.1
        with self.assertRaises(DistributionalRobustError):
            DistributionallyRobustAnalyzer(PLAN, config, SENSITIVITY)


if __name__ == "__main__":
    unittest.main()
