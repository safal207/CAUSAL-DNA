import copy
import json
import unittest
from pathlib import Path

from causal_dna.sensitivity import (
    SensitivityError,
    StrategySensitivityAnalyzer,
    grid_size,
)

ROOT = Path(__file__).resolve().parents[1]
PLAN_PATH = ROOT / "cases" / "CDNA-001.experiment-plan.json"
CONFIG_PATH = ROOT / "cases" / "CDNA-001.sensitivity-analysis.json"


class StrategySensitivityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.plan = json.loads(PLAN_PATH.read_text(encoding="utf-8"))
        cls.config = json.loads(CONFIG_PATH.read_text(encoding="utf-8"))

    def test_grid_contains_more_than_one_thousand_distinct_priors(self):
        self.assertEqual(1140, grid_size(21, 4, 1))
        analyzer = StrategySensitivityAnalyzer(self.plan, self.config)
        priors = list(analyzer.prior_grid())
        self.assertEqual(1140, len(priors))
        self.assertEqual(1140, len({tuple(p.values()) for p in priors}))

    def test_each_grid_prior_preserves_all_hypotheses(self):
        analyzer = StrategySensitivityAnalyzer(self.plan, self.config)
        for prior in analyzer.prior_grid():
            self.assertEqual(set(self.plan["hypotheses"]), set(prior))
            self.assertTrue(all(value > 0 for value in prior.values()))

    def test_analysis_counts_all_grid_points(self):
        result = StrategySensitivityAnalyzer(self.plan, self.config).analyze()
        self.assertEqual(1140, result.grid_points)
        self.assertEqual(1140, sum(result.first_step_counts.values()))
        self.assertEqual(1140, sum(result.policy_counts.values()))

    def test_support_values_are_fractions(self):
        result = StrategySensitivityAnalyzer(self.plan, self.config).analyze()
        for value in (
            result.dominant_first_step_support,
            result.baseline_first_step_support,
            result.dominant_policy_support,
            result.baseline_policy_support,
        ):
            self.assertGreaterEqual(value, 0.0)
            self.assertLessEqual(value, 1.0)

    def test_sensitivity_result_cannot_claim_causal_status(self):
        result = StrategySensitivityAnalyzer(self.plan, self.config).analyze()
        self.assertFalse(hasattr(result, "cause_found"))
        self.assertFalse(hasattr(result, "causal_status"))
        self.assertFalse(hasattr(result, "materialized"))

    def test_too_small_grid_is_rejected(self):
        config = copy.deepcopy(self.config)
        config["simplex_total_units"] = 10
        with self.assertRaises(SensitivityError):
            StrategySensitivityAnalyzer(self.plan, config)

    def test_causal_mutation_field_is_rejected(self):
        config = copy.deepcopy(self.config)
        config["cause_found"] = True
        with self.assertRaises(SensitivityError):
            StrategySensitivityAnalyzer(self.plan, config)


if __name__ == "__main__":
    unittest.main()
