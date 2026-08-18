import copy
import json
import unittest
from pathlib import Path

from causal_dna.adversarial_radius import (
    AdversarialRadiusAnalyzer,
    AdversarialRadiusError,
)

ROOT = Path(__file__).resolve().parents[1]


def load(path):
    return json.loads((ROOT / path).read_text(encoding="utf-8"))


class AdversarialRadiusTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.plan = load("cases/CDNA-001.experiment-plan.json")
        cls.joint = load("cases/CDNA-001.joint-assay-frontier.json")
        cls.config = load("cases/CDNA-001.adversarial-radius.json")
        cls.analyzer = AdversarialRadiusAnalyzer(cls.plan, cls.joint, cls.config)
        cls.result = cls.analyzer.analyze()

    def test_production_config_is_valid(self):
        AdversarialRadiusAnalyzer(self.plan, self.joint, self.config)

    def test_reference_winner_is_target(self):
        self.assertEqual(
            self.result.baseline_winner_action_id,
            self.config["target_action_id"],
        )

    def test_radius_one_applies_declared_scales(self):
        active = {axis["axis_id"] for axis in self.config["axes"]}
        settings = self.analyzer._settings_at_radius(1.0, active)
        self.assertAlmostEqual(settings["E1_ARID5B_OCCUPANCY"][0], 0.75)
        self.assertAlmostEqual(settings["E3_ALLELE_SNATAC"][0], 0.75)
        self.assertAlmostEqual(settings["E5_MATCHED_MULTIOMIC"][0], 1.25)
        self.assertAlmostEqual(settings["E5_MATCHED_MULTIOMIC"][1], 0.90)

    def test_coordinated_boundary_if_found_changes_winner(self):
        if not self.result.coordinated_switch_found:
            self.assertIsNone(self.result.coordinated_radius)
            return
        self.assertIsNotNone(self.result.coordinated_radius)
        self.assertNotEqual(
            self.result.coordinated_winner_after_boundary,
            self.config["target_action_id"],
        )
        epsilon = 2 * float(self.config["binary_tolerance"])
        before = max(0.0, float(self.result.coordinated_radius) - epsilon)
        active = {axis["axis_id"] for axis in self.config["axes"]}
        self.assertEqual(
            self.analyzer._cell_at_radius(before, active).winner_action_id,
            self.config["target_action_id"],
        )

    def test_single_axis_thresholds_are_consistent(self):
        self.assertEqual(len(self.result.axis_thresholds), len(self.config["axes"]))
        for threshold in self.result.axis_thresholds:
            if threshold.found:
                self.assertIsNotNone(threshold.normalized_distance)
                self.assertIsNotNone(threshold.boundary_value)
                self.assertNotEqual(
                    threshold.winner_after_boundary,
                    self.config["target_action_id"],
                )
            else:
                self.assertIsNone(threshold.normalized_distance)
                self.assertIsNone(threshold.boundary_value)

    def test_wrong_direction_limit_is_rejected_when_analyzed(self):
        config = copy.deepcopy(self.config)
        config["axes"][0]["limit"] = 1.5
        analyzer = AdversarialRadiusAnalyzer(self.plan, self.joint, config)
        with self.assertRaises(AdversarialRadiusError):
            analyzer.analyze()

    def test_causal_mutation_field_is_rejected(self):
        config = copy.deepcopy(self.config)
        config["cause_found"] = True
        with self.assertRaises(AdversarialRadiusError):
            AdversarialRadiusAnalyzer(self.plan, self.joint, config)

    def test_result_cannot_claim_causal_status(self):
        forbidden = {
            "cause_found",
            "causal_status",
            "edge_status",
            "materialized",
            "verification_status",
        }
        self.assertFalse(forbidden.intersection(vars(self.result)))
        for threshold in self.result.axis_thresholds:
            self.assertFalse(forbidden.intersection(vars(threshold)))


if __name__ == "__main__":
    unittest.main()
