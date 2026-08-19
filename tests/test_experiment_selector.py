import copy
import json
import unittest
from pathlib import Path

from causal_dna.experiment_selector import ExperimentSelector, ExperimentSelectorError

ROOT = Path(__file__).resolve().parents[1]
PLAN_PATH = ROOT / "cases" / "CDNA-001.experiment-plan.json"


class ExperimentSelectorTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.plan = json.loads(PLAN_PATH.read_text(encoding="utf-8"))

    def selector(self):
        return ExperimentSelector(copy.deepcopy(self.plan))

    def test_plan_is_valid(self):
        self.selector().assert_valid()

    def test_information_gain_is_non_negative(self):
        for score in self.selector().rank():
            self.assertGreaterEqual(score.expected_information_gain_bits, 0.0)
            self.assertGreaterEqual(score.normalized_information_gain, 0.0)

    def test_matched_multiomic_maximizes_information_gain_under_declared_planning_priors(self):
        best = self.selector().best(objective="information_gain")
        self.assertEqual("E5_MATCHED_MULTIOMIC", best.experiment_id)

    def test_gain_per_cost_can_differ_from_raw_information_gain(self):
        raw = self.selector().best(objective="information_gain")
        efficient = self.selector().best(objective="gain_per_cost")
        self.assertGreater(raw.expected_information_gain_bits, 0.0)
        self.assertGreater(efficient.utility_per_cost, 0.0)

    def test_probabilities_must_sum_to_one(self):
        plan = copy.deepcopy(self.plan)
        plan["experiments"][0]["outcomes"][0]["planning_probability"] = 0.8
        with self.assertRaises(ExperimentSelectorError):
            ExperimentSelector(plan)

    def test_posterior_must_cover_all_hypotheses(self):
        plan = copy.deepcopy(self.plan)
        del plan["experiments"][0]["outcomes"][0]["posterior_weights"]["H4_3D_CONTACT"]
        with self.assertRaises(ExperimentSelectorError):
            ExperimentSelector(plan)

    def test_selector_outputs_scores_only_not_causal_status(self):
        score = self.selector().best()
        self.assertFalse(hasattr(score, "cause_found"))
        self.assertFalse(hasattr(score, "causal_status"))


if __name__ == "__main__":
    unittest.main()
