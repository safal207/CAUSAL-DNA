import copy
import json
import unittest
from pathlib import Path

from causal_dna.adaptive_replanner import AdaptiveReplanner, AdaptiveReplannerError
from causal_dna.processor import CausalProcessor

ROOT = Path(__file__).resolve().parents[1]
PLAN_PATH = ROOT / "cases" / "CDNA-001.experiment-plan.json"
SCENARIOS_PATH = ROOT / "cases" / "CDNA-001.replanning-scenarios.json"
EVENTS_PATH = ROOT / "cases" / "CDNA-001.events.json"


class AdaptiveReplannerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.plan = json.loads(PLAN_PATH.read_text(encoding="utf-8"))
        cls.scenarios = json.loads(SCENARIOS_PATH.read_text(encoding="utf-8"))
        cls.events = json.loads(EVENTS_PATH.read_text(encoding="utf-8"))

    def scenario(self, scenario_id):
        return next(s for s in self.scenarios if s["scenario_id"] == scenario_id)

    def real_event(self, outcome_name="no_allele_specific_occupancy", generation=2):
        return {
            "event_id": "EV-TEST-E1-OUTCOME",
            "generation": generation,
            "event_type": "experiment_outcome",
            "subject": "E1_ARID5B_OCCUPANCY",
            "observer": "experiment",
            "evidence_refs": ["test:synthetic-measurement-fixture"],
            "payload": {
                "experiment_id": "E1_ARID5B_OCCUPANCY",
                "outcome_name": outcome_name,
            },
        }

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

    def test_real_event_advances_generation_without_causal_promotion(self):
        processor = CausalProcessor(self.events)
        before = processor.project()
        result = AdaptiveReplanner(self.plan).replan_from_event(
            processor, self.real_event()
        )
        after = processor.project()
        self.assertEqual(2, result.generation)
        self.assertEqual(2, processor.generation)
        self.assertEqual(before.open_hypotheses, after.open_hypotheses)
        self.assertEqual(before.rejected_hypotheses, after.rejected_hypotheses)
        self.assertEqual(before.superseded_hypotheses, after.superseded_hypotheses)
        self.assertEqual(before.verified_subjects, after.verified_subjects)
        self.assertLess(
            result.planning_weights["H1_ARID5B_REUSE"],
            result.planning_weights["H2_CUX1"],
        )
        self.assertNotEqual("E1_ARID5B_OCCUPANCY", result.best_information_gain.experiment_id)

    def test_real_event_requires_evidence_provenance(self):
        processor = CausalProcessor(self.events)
        event = self.real_event()
        event["evidence_refs"] = []
        with self.assertRaises(AdaptiveReplannerError):
            AdaptiveReplanner(self.plan).replan_from_event(processor, event)

    def test_real_event_must_target_next_processor_generation(self):
        processor = CausalProcessor(self.events)
        with self.assertRaises(AdaptiveReplannerError):
            AdaptiveReplanner(self.plan).replan_from_event(
                processor, self.real_event(generation=3)
            )


if __name__ == "__main__":
    unittest.main()
