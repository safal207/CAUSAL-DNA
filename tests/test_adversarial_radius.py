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
        fractions = {axis["axis_id"]: 1.0 for axis in self.config["axes"]}
        settings = self.analyzer._settings_at_profile(1.0, fractions)
        self.assertAlmostEqual(settings["E1_ARID5B_OCCUPANCY"][0], 0.75)
        self.assertAlmostEqual(settings["E3_ALLELE_SNATAC"][0], 0.75)
        self.assertAlmostEqual(settings["E5_MATCHED_MULTIOMIC"][0], 1.25)
        self.assertAlmostEqual(settings["E5_MATCHED_MULTIOMIC"][1], 0.90)

    def test_fractional_ray_family_is_normalized(self):
        profiles = self.analyzer._fractional_profiles()
        self.assertGreater(len(profiles), 1)
        for profile in profiles:
            self.assertAlmostEqual(max(profile.values()), 1.0)
        self.assertEqual(self.result.fractional_rays_tested, len(profiles))

    def test_minimum_ray_boundary_if_found_changes_winner(self):
        if not self.result.minimum_ray_switch_found:
            self.assertIsNone(self.result.minimum_ray_radius)
            return
        self.assertIsNotNone(self.result.minimum_ray_radius)
        self.assertIsNotNone(self.result.minimum_ray_fractions)
        self.assertNotEqual(
            self.result.minimum_ray_winner_after_boundary,
            self.config["target_action_id"],
        )
        self.assertIsNotNone(
            self.result.minimum_ray_target_regret_gap_after_boundary
        )

    def test_minimum_ray_is_no_farther_than_coordinated_ray(self):
        if (
            self.result.minimum_ray_switch_found
            and self.result.coordinated_switch_found
        ):
            self.assertLessEqual(
                float(self.result.minimum_ray_radius),
                float(self.result.coordinated_radius) + 1e-12,
            )

    def test_coordinated_boundary_if_found_changes_winner(self):
        if not self.result.coordinated_switch_found:
            self.assertIsNone(self.result.coordinated_radius)
            return
        self.assertNotEqual(
            self.result.coordinated_winner_after_boundary,
            self.config["target_action_id"],
        )
        self.assertIsNotNone(
            self.result.coordinated_target_regret_gap_after_boundary
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
                self.assertIsNotNone(
                    threshold.target_regret_gap_after_boundary
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

    def test_fraction_grid_must_span_zero_to_one(self):
        config = copy.deepcopy(self.config)
        config["fraction_grid"] = [0.0, 0.5]
        with self.assertRaises(AdversarialRadiusError):
            AdversarialRadiusAnalyzer(self.plan, self.joint, config)

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
