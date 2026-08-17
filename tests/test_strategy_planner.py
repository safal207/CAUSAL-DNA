import json
import unittest
from pathlib import Path

from causal_dna.strategy_planner import StrategyPlanner, StrategyPlannerError

ROOT = Path(__file__).resolve().parents[1]
PLAN_PATH = ROOT / "cases" / "CDNA-001.experiment-plan.json"
CONFIG_PATH = ROOT / "cases" / "CDNA-001.strategy-configs.json"


class StrategyPlannerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.plan = json.loads(PLAN_PATH.read_text(encoding="utf-8"))
        cls.configs = json.loads(CONFIG_PATH.read_text(encoding="utf-8"))
        cls.planner = StrategyPlanner(cls.plan)

    def test_one_step_matches_existing_information_gain_winner(self):
        result = self.planner.best(
            budget=3.0, max_depth=1, objective="information_gain"
        )
        self.assertEqual("E5_MATCHED_MULTIOMIC", result.first_experiment_id)

    def test_every_declared_strategy_respects_path_budget(self):
        for config in self.configs:
            result = self.planner.best(
                budget=config["budget"],
                max_depth=config["max_depth"],
                objective=config["objective"],
            )
            self.assertLessEqual(result.max_path_cost, config["budget"] + 1e-9)
            self.assertGreater(result.expected_information_gain_bits, 0.0)

    def test_two_step_lookahead_is_not_worse_than_one_step_under_same_budget(self):
        one = self.planner.best(
            budget=2.0, max_depth=1, objective="information_gain"
        )
        two = self.planner.best(
            budget=2.0, max_depth=2, objective="information_gain"
        )
        self.assertGreaterEqual(
            two.expected_information_gain_bits + 1e-12,
            one.expected_information_gain_bits,
        )

    def test_strategy_exposes_contingent_next_actions(self):
        result = self.planner.best(
            budget=4.5, max_depth=3, objective="information_gain"
        )
        self.assertTrue(result.next_experiments_by_outcome)
        self.assertTrue(
            all(
                next_id != result.first_experiment_id
                for next_id in result.next_experiments_by_outcome.values()
                if next_id is not None
            )
        )

    def test_gain_per_cost_strategy_is_positive(self):
        result = self.planner.best(
            budget=3.0, max_depth=3, objective="gain_per_cost"
        )
        self.assertGreater(result.gain_per_cost, 0.0)
        self.assertLessEqual(result.max_path_cost, 3.0 + 1e-9)

    def test_strategy_result_cannot_claim_causal_status(self):
        result = self.planner.best(
            budget=2.0, max_depth=2, objective="information_gain"
        )
        self.assertFalse(hasattr(result, "cause_found"))
        self.assertFalse(hasattr(result, "causal_status"))
        self.assertFalse(hasattr(result, "materialized"))

    def test_depth_is_intentionally_bounded(self):
        with self.assertRaises(StrategyPlannerError):
            self.planner.best(budget=10.0, max_depth=4)

    def test_too_small_budget_cannot_fund_experiment(self):
        with self.assertRaises(StrategyPlannerError):
            self.planner.best(budget=0.1, max_depth=2)


if __name__ == "__main__":
    unittest.main()
