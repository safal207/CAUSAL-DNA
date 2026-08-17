#!/usr/bin/env python3
from __future__ import annotations

import copy
import json
import unittest
from dataclasses import asdict
from pathlib import Path

from causal_dna.phase_boundary import (
    PhaseBoundaryError,
    RobustPhaseBoundaryAnalyzer,
)

ROOT = Path(__file__).resolve().parents[1]


def load(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


class PhaseBoundaryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.plan = load(ROOT / "cases" / "CDNA-001.experiment-plan.json")
        cls.production_config = load(ROOT / "cases" / "CDNA-001.phase-boundary.json")

    def tiny_config(self):
        config = copy.deepcopy(self.production_config)
        config["prior_simplex_total_units"] = 6
        config["target_reliability_grid"] = [0.5, 1.0]
        config["target_cost_multiplier_grid"] = [0.5, 1.5]
        config["outcome_probability_max_transfer_grid"] = [0.0, 0.2]
        return config

    def test_production_config_is_valid(self):
        RobustPhaseBoundaryAnalyzer(self.plan, self.production_config)

    def test_causal_mutation_field_is_rejected(self):
        config = self.tiny_config()
        config["cause_found"] = True
        with self.assertRaises(PhaseBoundaryError):
            RobustPhaseBoundaryAnalyzer(self.plan, config)

    def test_grid_values_must_be_sorted(self):
        config = self.tiny_config()
        config["target_reliability_grid"] = [1.0, 0.5]
        with self.assertRaises(PhaseBoundaryError):
            RobustPhaseBoundaryAnalyzer(self.plan, config)

    def test_tiny_grid_cell_count_is_cartesian_product(self):
        analyzer = RobustPhaseBoundaryAnalyzer(self.plan, self.tiny_config())
        result = analyzer.analyze()
        self.assertEqual(len(result.cells), 8)

    def test_result_cannot_claim_causal_status(self):
        result = RobustPhaseBoundaryAnalyzer(self.plan, self.tiny_config()).analyze()
        payload = asdict(result)
        forbidden = {
            "cause_found",
            "causal_status",
            "edge_status",
            "materialized",
            "verification_status",
        }
        self.assertTrue(forbidden.isdisjoint(payload.keys()))

    def test_target_cost_increase_cannot_improve_target_interval_utility(self):
        analyzer = RobustPhaseBoundaryAnalyzer(self.plan, self.tiny_config())
        prior = next(analyzer.prior_grid())
        target = next(
            item for item in self.plan["experiments"]
            if item["experiment_id"] == "E5_MATCHED_MULTIOMIC"
        )
        cheap = analyzer._experiment_interval(
            prior,
            target,
            target_reliability=0.8,
            target_cost_multiplier=0.5,
            probability_transfer=0.1,
        )
        expensive = analyzer._experiment_interval(
            prior,
            target,
            target_reliability=0.8,
            target_cost_multiplier=1.5,
            probability_transfer=0.1,
        )
        self.assertGreaterEqual(cheap[0], expensive[0])
        self.assertGreaterEqual(cheap[1], expensive[1])

    def test_target_reliability_increase_does_not_reduce_upper_utility_in_reference_cell(self):
        analyzer = RobustPhaseBoundaryAnalyzer(self.plan, self.tiny_config())
        prior = {key: 1.0 for key in self.plan["hypotheses"]}
        target = next(
            item for item in self.plan["experiments"]
            if item["experiment_id"] == "E5_MATCHED_MULTIOMIC"
        )
        low = analyzer._experiment_interval(
            prior,
            target,
            target_reliability=0.5,
            target_cost_multiplier=1.0,
            probability_transfer=0.0,
        )
        high = analyzer._experiment_interval(
            prior,
            target,
            target_reliability=1.0,
            target_cost_multiplier=1.0,
            probability_transfer=0.0,
        )
        self.assertGreaterEqual(high[1], low[1])

    def test_boundaries_only_connect_different_winners(self):
        result = RobustPhaseBoundaryAnalyzer(self.plan, self.tiny_config()).analyze()
        for boundary in result.boundaries:
            self.assertNotEqual(boundary.from_action_id, boundary.to_action_id)


if __name__ == "__main__":
    unittest.main()
