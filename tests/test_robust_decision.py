import copy
import json
import unittest
from pathlib import Path

from causal_dna.robust_decision import RobustDecisionError, RobustFirstStepAnalyzer

ROOT = Path(__file__).resolve().parents[1]
PLAN_PATH = ROOT / "cases" / "CDNA-001.experiment-plan.json"
ROBUST_PATH = ROOT / "cases" / "CDNA-001.robust-decision.json"
SENS_PATH = ROOT / "cases" / "CDNA-001.sensitivity-analysis.json"


class RobustFirstStepTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.plan = json.loads(PLAN_PATH.read_text(encoding="utf-8"))
        cls.config = json.loads(ROBUST_PATH.read_text(encoding="utf-8"))
        cls.sensitivity = json.loads(SENS_PATH.read_text(encoding="utf-8"))
        cls.result = RobustFirstStepAnalyzer(
            cls.plan, cls.config, cls.sensitivity
        ).analyze()

    def test_uses_same_full_prior_grid(self):
        self.assertEqual(1140, self.result.grid_points)

    def test_every_affordable_experiment_and_stop_are_scored(self):
        action_ids = {item.action_id for item in self.result.actions}
        self.assertEqual(
            {
                "E1_ARID5B_OCCUPANCY",
                "E2_CUX1_OCCUPANCY",
                "E3_ALLELE_SNATAC",
                "E4_CAPTURE_C",
                "E5_MATCHED_MULTIOMIC",
                "STOP",
            },
            action_ids,
        )

    def test_minimax_really_minimizes_max_regret(self):
        chosen = self.result.action(self.result.minimax_regret_action_id)
        for item in self.result.actions:
            self.assertLessEqual(chosen.max_regret_bits, item.max_regret_bits + 1e-12)

    def test_maximin_really_maximizes_worst_case_gain(self):
        chosen = self.result.action(self.result.maximin_gain_action_id)
        for item in self.result.actions:
            self.assertGreaterEqual(
                chosen.worst_case_gain_bits + 1e-12,
                item.worst_case_gain_bits,
            )

    def test_stop_is_explicit_zero_information_action(self):
        stop = self.result.action("STOP")
        self.assertEqual(0.0, stop.worst_case_gain_bits)
        self.assertEqual(0.0, stop.best_case_gain_bits)
        self.assertEqual(0.0, stop.mean_gain_bits)

    def test_regret_values_are_non_negative(self):
        for item in self.result.actions:
            self.assertGreaterEqual(item.max_regret_bits, 0.0)
            self.assertGreaterEqual(item.mean_regret_bits, 0.0)
            self.assertGreaterEqual(item.optimal_count, 0)
            self.assertLessEqual(item.optimal_count, self.result.grid_points)

    def test_result_cannot_claim_causal_status(self):
        self.assertFalse(hasattr(self.result, "cause_found"))
        self.assertFalse(hasattr(self.result, "causal_status"))
        self.assertFalse(hasattr(self.result, "materialized"))

    def test_causal_mutation_field_is_rejected(self):
        config = copy.deepcopy(self.config)
        config["cause_found"] = True
        with self.assertRaises(RobustDecisionError):
            RobustFirstStepAnalyzer(self.plan, config, self.sensitivity)


if __name__ == "__main__":
    unittest.main()
