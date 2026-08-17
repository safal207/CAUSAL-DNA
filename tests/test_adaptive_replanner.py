import copy
import json
import unittest
from pathlib import Path

from causal_dna.adaptive_replanner import AdaptiveReplanner, AdaptiveReplannerError

ROOT = Path(__file__).resolve().parents[1]
PLAN_PATH = ROOT / "cases" / "CDNA-001.experiment-plan.json"
SCENARIOS_PATH = ROOT / "cases" / "CDNA-001.replanning-scenarios.json"


class AdaptiveReplannerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.plan = json.loads(PLAN_PATH.read_text(encoding="utf-8"))
        cls.scenarios = json.loads(SCENARIOS_PATH.read_text(encoding="utf-8"))

    def scenario(self, scenario_id):
        return next(s for s in self.scenarios if s["scenario_id"] == scenario_id)

    def test_negative_arid5b_outcome_lowers_only_planning_weight(self):
        scenario = self.scenario("SIM-E1-NEGATIVE")
        snapshot = AdaptiveReplanner(self.plan, scenario["outcomes"]).replay()
        self.assertEqual(2, snapshot.generation)
        self.assertLess(
            snapshot.planning_weights["H1_ARID5B_REUSE"],
            snapshot.planning_weights["H2_CUX1"],
        )
        self.assertNotEqual("E1_ARID5B_OCCUPANCY", snapshot.best_information_gain.experiment_id)

    def test_positive_arid5b_outcome_raises_planning_weight(self):
        scenario = self.scenario("SIM-E1-POSITIVE")
        snapshot = AdaptiveReplanner(self.plan, scenario["outcomes"]).replay()
        self.assertGreater(
            snapshot.planning_weights["H1_ARID5B_REUSE"],
            snapshot.planning_weights["H2_CUX1"],
        )

    def test_completed_experiment_is_not_selected_again(self):
        scenario = self.scenario("SIM-E1-NEGATIVE")
        snapshot = AdaptiveReplanner(self.plan, scenario["outcomes"]).replay()
        self.assertIn("E1_ARID5B_OCCUPANCY", snapshot.completed_experiments)
        self.assertNotEqual("E1_ARID5B_OCCUPANCY", snapshot.best_gain_per_cost.experiment_id)

    def test_opposite_outcomes_create_different_replanned_states(self):
        negative = AdaptiveReplanner(
            self.plan, self.scenario("SIM-E1-NEGATIVE")["outcomes"]
        ).replay()
        positive = AdaptiveReplanner(
            self.plan, self.scenario("SIM-E1-POSITIVE")["outcomes"]
        ).replay()
        self.assertNotEqual(negative.planning_weights, positive.planning_weights)
        self.assertNotEqual(negative.entropy_bits, positive.entropy_bits)

    def test_planning_outcome_cannot_set_causal_status(self):
        record = copy.deepcopy(self.scenario("SIM-E1-NEGATIVE")["outcomes"][0])
        record["cause_found"] = True
        with self.assertRaises(AdaptiveReplannerError):
            AdaptiveReplanner(self.plan, [record])

    def test_generation_must_advance(self):
        record = copy.deepcopy(self.scenario("SIM-E1-NEGATIVE")["outcomes"][0])
        record["generation"] = 1
        with self.assertRaises(AdaptiveReplannerError):
            AdaptiveReplanner(self.plan, [record])

    def test_simulation_flag_is_mandatory(self):
        record = copy.deepcopy(self.scenario("SIM-E1-NEGATIVE")["outcomes"][0])
        del record["simulation_only"]
        with self.assertRaises(AdaptiveReplannerError):
            AdaptiveReplanner(self.plan, [record])


if __name__ == "__main__":
    unittest.main()
