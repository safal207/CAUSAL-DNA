import copy
import json
import unittest
from pathlib import Path

from causal_dna.joint_assay_frontier import (
    JointAssayFrontierAnalyzer,
    JointAssayFrontierError,
)

ROOT = Path(__file__).resolve().parents[1]


def load(path):
    return json.loads((ROOT / path).read_text(encoding="utf-8"))


class JointAssayFrontierTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.plan = load("cases/CDNA-001.experiment-plan.json")
        cls.config = load("cases/CDNA-001.joint-assay-frontier.json")
        cls.analyzer = JointAssayFrontierAnalyzer(cls.plan, cls.config)
        cls.result = cls.analyzer.analyze()

    def test_production_config_is_valid(self):
        JointAssayFrontierAnalyzer(self.plan, self.config)

    def test_full_grid_is_cartesian_product(self):
        states_per_assay = len(self.config["cost_multiplier_grid"]) * len(
            self.config["reliability_grid"]
        )
        expected = states_per_assay ** len(self.config["joint_experiment_ids"])
        self.assertEqual(len(self.result.cells), expected)
        self.assertEqual(expected, 729)

    def test_prior_grid_is_deterministic_165_states(self):
        self.assertEqual(self.result.prior_points_per_cell, 165)

    def test_all_joint_assays_vary_independently(self):
        observed = {
            experiment_id: set()
            for experiment_id in self.config["joint_experiment_ids"]
        }
        for cell in self.result.cells:
            for experiment_id, cost, reliability in cell.settings:
                observed[experiment_id].add((cost, reliability))
        expected_states = {
            (float(cost), float(reliability))
            for cost in self.config["cost_multiplier_grid"]
            for reliability in self.config["reliability_grid"]
        }
        for experiment_id in observed:
            self.assertEqual(observed[experiment_id], expected_states)

    def test_causal_mutation_field_is_rejected(self):
        config = copy.deepcopy(self.config)
        config["cause_found"] = True
        with self.assertRaises(JointAssayFrontierError):
            JointAssayFrontierAnalyzer(self.plan, config)

    def test_boundaries_only_record_distinct_action_pairs(self):
        self.assertTrue(self.result.boundary_pair_counts)
        for pair, count in self.result.boundary_pair_counts.items():
            left, right = pair.split("<->")
            self.assertNotEqual(left, right)
            self.assertGreater(count, 0)

    def test_pareto_frontier_is_nonempty_and_cross_horizon(self):
        self.assertTrue(self.result.pareto_points)
        horizons = {point.horizon for point in self.result.pareto_points}
        self.assertTrue(horizons.issubset({1, 2, 3}))
        self.assertGreaterEqual(len(horizons), 2)

    def test_pareto_points_do_not_dominate_each_other(self):
        points = self.result.pareto_points
        for i, candidate in enumerate(points):
            for j, other in enumerate(points):
                if i == j:
                    continue
                no_worse = (
                    other.expected_information_gain_bits
                    >= candidate.expected_information_gain_bits - 1e-12
                    and other.expected_cost <= candidate.expected_cost + 1e-12
                    and other.first_step_max_regret
                    <= candidate.first_step_max_regret + 1e-12
                )
                strict = (
                    other.expected_information_gain_bits
                    > candidate.expected_information_gain_bits + 1e-12
                    or other.expected_cost < candidate.expected_cost - 1e-12
                    or other.first_step_max_regret
                    < candidate.first_step_max_regret - 1e-12
                )
                self.assertFalse(no_worse and strict)

    def test_result_cannot_claim_causal_status(self):
        forbidden = {
            "cause_found",
            "causal_status",
            "edge_status",
            "materialized",
            "verification_status",
        }
        self.assertFalse(forbidden.intersection(vars(self.result)))
        for point in self.result.pareto_points:
            self.assertFalse(forbidden.intersection(vars(point)))


if __name__ == "__main__":
    unittest.main()
