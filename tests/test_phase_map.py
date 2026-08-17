import copy
import json
import unittest
from pathlib import Path

from causal_dna.phase_map import PhaseBoundaryAnalyzer, PhaseMapError

ROOT = Path(__file__).resolve().parents[1]


def load(name):
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


class PhaseBoundaryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.plan = load("cases/CDNA-001.experiment-plan.json")
        cls.config = load("cases/CDNA-001.phase-map.json")
        cls.sensitivity = load("cases/CDNA-001.sensitivity-analysis.json")

    def analyzer(self, config=None):
        return PhaseBoundaryAnalyzer(
            self.plan,
            config or self.config,
            self.sensitivity,
        )

    def test_causal_mutation_field_is_rejected(self):
        config = copy.deepcopy(self.config)
        config["cause_found"] = True
        with self.assertRaises(PhaseMapError):
            self.analyzer(config)

    def test_low_target_cost_can_make_multiomic_baseline_winner(self):
        analyzer = self.analyzer()
        utilities = analyzer._utilities(
            self.plan["hypotheses"],
            target_cost_multiplier=0.5,
            target_reliability=1.0,
            max_transfer=0.0,
        )
        self.assertEqual(analyzer._winner(utilities), "E5_MATCHED_MULTIOMIC")

    def test_baseline_target_cost_prefers_cheaper_assay_on_efficiency(self):
        analyzer = self.analyzer()
        utilities = analyzer._utilities(
            self.plan["hypotheses"],
            target_cost_multiplier=1.0,
            target_reliability=1.0,
            max_transfer=0.0,
        )
        self.assertEqual(analyzer._winner(utilities), "E1_ARID5B_OCCUPANCY")

    def test_break_even_separates_target_and_competitor(self):
        analyzer = self.analyzer()
        point = analyzer._break_even(
            self.plan["hypotheses"],
            target_reliability=1.0,
            max_transfer=0.0,
        )
        self.assertEqual(point.competitor_id, "E1_ARID5B_OCCUPANCY")
        self.assertIsNotNone(point.break_even_cost_multiplier)
        threshold = point.break_even_cost_multiplier
        below = analyzer._winner(
            analyzer._utilities(
                self.plan["hypotheses"],
                target_cost_multiplier=threshold * 0.99,
                target_reliability=1.0,
                max_transfer=0.0,
            )
        )
        above = analyzer._winner(
            analyzer._utilities(
                self.plan["hypotheses"],
                target_cost_multiplier=threshold * 1.01,
                target_reliability=1.0,
                max_transfer=0.0,
            )
        )
        self.assertEqual(below, "E5_MATCHED_MULTIOMIC")
        self.assertNotEqual(above, "E5_MATCHED_MULTIOMIC")

    def test_small_surface_still_uses_full_prior_grid(self):
        config = copy.deepcopy(self.config)
        config["target_cost_multipliers"] = [0.5, 1.0]
        config["target_reliabilities"] = [0.8, 1.0]
        config["outcome_probability_max_transfers"] = [0.0]
        result = self.analyzer(config).analyze()
        self.assertEqual(result.prior_points, 1140)
        self.assertEqual(len(result.phase_points), 4)
        self.assertEqual(len(result.break_even_points), 2)

    def test_phase_result_cannot_claim_causal_status(self):
        config = copy.deepcopy(self.config)
        config["target_cost_multipliers"] = [0.5, 1.0]
        config["target_reliabilities"] = [0.8, 1.0]
        config["outcome_probability_max_transfers"] = [0.0]
        result = self.analyzer(config).analyze()
        forbidden = {"cause_found", "causal_status", "materialized", "verification_status"}
        self.assertFalse(forbidden.intersection(vars(result)))


if __name__ == "__main__":
    unittest.main()
