#!/usr/bin/env python3
"""Joint adversarial assay frontier for CAUSAL-DNA planning.

This layer varies cost and reliability for several experiments at the same time,
then evaluates first-step minimax interval regret across a deterministic prior
grid. It also attaches first-step robust regret to short-horizon policy
frontiers so cost, information, regret, and planning horizon can be inspected
together.

Everything here is planning-only. Reliability, costs, outcome probabilities,
and priors are declared decision-model assumptions, not biological evidence or
causal confidence.
"""
from __future__ import annotations

import copy
import math
from collections import Counter
from dataclasses import dataclass
from itertools import product
from statistics import fmean
from typing import Any, Iterator

from .distributional_robust import (
    _apply,
    _entropy,
    _probability_scenarios,
    _reliability_adjusted_multipliers,
)
from .experiment_selector import ExperimentSelector
from .sensitivity import _positive_compositions, grid_size
from .strategy_planner import Policy, StrategyPlanner


class JointAssayFrontierError(ValueError):
    pass


@dataclass(frozen=True)
class JointActionSummary:
    action_id: str
    max_regret: float
    mean_regret: float
    worst_utility: float
    mean_lower_utility: float


@dataclass(frozen=True)
class JointCell:
    settings: tuple[tuple[str, float, float], ...]
    winner_action_id: str
    winner_max_regret: float
    runner_up_action_id: str
    runner_up_max_regret: float
    regret_margin: float

    def setting(self, experiment_id: str) -> tuple[float, float]:
        for current_id, cost_multiplier, reliability in self.settings:
            if current_id == experiment_id:
                return cost_multiplier, reliability
        raise JointAssayFrontierError(f"setting not found: {experiment_id}")


@dataclass(frozen=True)
class PolicyParetoPoint:
    horizon: int
    first_experiment_id: str
    expected_information_gain_bits: float
    expected_cost: float
    max_path_cost: float
    first_step_max_regret: float


@dataclass(frozen=True)
class JointAssayFrontierResult:
    analysis_id: str
    prior_points_per_cell: int
    cells: tuple[JointCell, ...]
    winner_counts: dict[str, int]
    boundary_pair_counts: dict[str, int]
    boundary_axis_counts: dict[str, int]
    pareto_points: tuple[PolicyParetoPoint, ...]


class JointAssayFrontierAnalyzer:
    FORBIDDEN_CONFIG_FIELDS = {
        "cause_found",
        "causal_status",
        "edge_status",
        "materialized",
        "verification_status",
    }

    def __init__(self, plan: dict[str, Any], config: dict[str, Any]):
        self.plan = copy.deepcopy(plan)
        self.config = copy.deepcopy(config)
        ExperimentSelector(self.plan)
        self.hypothesis_ids = tuple(self.plan["hypotheses"].keys())
        self.experiment_by_id = {
            item["experiment_id"]: item for item in self.plan["experiments"]
        }
        self.joint_ids = tuple(self.config.get("joint_experiment_ids", ()))
        self.assert_valid()
        self._priors = tuple(self.prior_grid())
        self._gain_cache: dict[tuple[int, str, float], tuple[float, float]] = {}

    def assert_valid(self) -> None:
        forbidden = self.FORBIDDEN_CONFIG_FIELDS.intersection(self.config)
        if forbidden:
            raise JointAssayFrontierError(
                "joint-frontier config cannot mutate causal state: "
                + ", ".join(sorted(forbidden))
            )
        if self.config.get("case_id") != self.plan.get("case_id"):
            raise JointAssayFrontierError("case_id must match experiment plan")
        if self.config.get("decision_rule") != "minimax_interval_regret":
            raise JointAssayFrontierError("decision_rule must be minimax_interval_regret")
        if self.config.get("objective") != "information_gain_per_cost":
            raise JointAssayFrontierError("objective must be information_gain_per_cost")
        if self.config.get("prior_model") != "deterministic_simplex_grid":
            raise JointAssayFrontierError("prior_model must be deterministic_simplex_grid")
        if self.config.get("include_stop") is not True:
            raise JointAssayFrontierError("include_stop must be true")
        if len(self.joint_ids) < 2 or len(set(self.joint_ids)) != len(self.joint_ids):
            raise JointAssayFrontierError("joint_experiment_ids must contain unique experiments")
        unknown = set(self.joint_ids) - set(self.experiment_by_id)
        if unknown:
            raise JointAssayFrontierError(
                "unknown joint experiments: " + ", ".join(sorted(unknown))
            )

        self._validate_grid("cost_multiplier_grid", lower=0.0, upper=None, strict=True)
        self._validate_grid("reliability_grid", lower=0.0, upper=1.0, strict=False)
        transfer = self.config.get("outcome_probability_max_transfer")
        if not isinstance(transfer, (int, float)) or not (0 <= transfer < 0.5):
            raise JointAssayFrontierError(
                "outcome_probability_max_transfer must be in [0, 0.5)"
            )

        total = self.config.get("prior_simplex_total_units")
        minimum = self.config.get("prior_min_units_per_hypothesis")
        if not isinstance(total, int) or total <= 0:
            raise JointAssayFrontierError("prior_simplex_total_units must be positive")
        if not isinstance(minimum, int) or minimum <= 0:
            raise JointAssayFrontierError(
                "prior_min_units_per_hypothesis must be positive"
            )
        prior_count = grid_size(total, len(self.hypothesis_ids), minimum)
        if not (10 <= prior_count <= 2000):
            raise JointAssayFrontierError("joint prior grid must contain 10..2000 states")

        states_per_assay = len(self.config["cost_multiplier_grid"]) * len(
            self.config["reliability_grid"]
        )
        cell_count = states_per_assay ** len(self.joint_ids)
        if cell_count > 2500:
            raise JointAssayFrontierError("joint grid is capped at 2500 cells")

        horizons = self.config.get("pareto_horizons")
        if (
            not isinstance(horizons, list)
            or not horizons
            or horizons != sorted(set(horizons))
            or any(not isinstance(value, int) or not (1 <= value <= 3) for value in horizons)
        ):
            raise JointAssayFrontierError(
                "pareto_horizons must be unique ascending integers in [1, 3]"
            )
        budget = self.config.get("pareto_budget")
        if not isinstance(budget, (int, float)) or budget <= 0:
            raise JointAssayFrontierError("pareto_budget must be > 0")

        reference = self.config.get("pareto_reference_settings")
        if not isinstance(reference, dict) or set(reference) != set(self.joint_ids):
            raise JointAssayFrontierError(
                "pareto_reference_settings must cover exactly joint_experiment_ids"
            )
        for experiment_id, settings in reference.items():
            if not isinstance(settings, dict):
                raise JointAssayFrontierError(
                    f"pareto reference for {experiment_id} must be an object"
                )
            cost_multiplier = settings.get("cost_multiplier")
            reliability = settings.get("reliability")
            if not isinstance(cost_multiplier, (int, float)) or cost_multiplier <= 0:
                raise JointAssayFrontierError("reference cost_multiplier must be > 0")
            if not isinstance(reliability, (int, float)) or not (0 <= reliability <= 1):
                raise JointAssayFrontierError("reference reliability must be in [0, 1]")

    def _validate_grid(
        self,
        field: str,
        *,
        lower: float,
        upper: float | None,
        strict: bool,
    ) -> None:
        values = self.config.get(field)
        if not isinstance(values, list) or len(values) < 2:
            raise JointAssayFrontierError(f"{field} must contain at least two values")
        if any(not isinstance(value, (int, float)) for value in values):
            raise JointAssayFrontierError(f"{field} must contain numeric values")
        numeric = [float(value) for value in values]
        if numeric != sorted(set(numeric)):
            raise JointAssayFrontierError(f"{field} must be unique and ascending")
        for value in numeric:
            if (strict and value <= lower) or (not strict and value < lower):
                comparator = ">" if strict else ">="
                raise JointAssayFrontierError(f"{field} values must be {comparator} {lower}")
            if upper is not None and value > upper:
                raise JointAssayFrontierError(f"{field} values must be <= {upper}")

    def prior_grid(self) -> Iterator[dict[str, float]]:
        total = int(self.config["prior_simplex_total_units"])
        minimum = int(self.config["prior_min_units_per_hypothesis"])
        for composition in _positive_compositions(
            total, len(self.hypothesis_ids), minimum
        ):
            yield {
                hypothesis_id: float(units)
                for hypothesis_id, units in zip(self.hypothesis_ids, composition)
            }

    def _gain_interval(
        self,
        prior_index: int,
        experiment_id: str,
        reliability: float,
    ) -> tuple[float, float]:
        key = (prior_index, experiment_id, round(float(reliability), 12))
        cached = self._gain_cache.get(key)
        if cached is not None:
            return cached

        prior = self._priors[prior_index]
        experiment = self.experiment_by_id[experiment_id]
        prior_entropy = _entropy(prior)
        baseline_probabilities = [
            float(outcome["planning_probability"])
            for outcome in experiment["outcomes"]
        ]
        scenarios = _probability_scenarios(
            baseline_probabilities,
            float(self.config["outcome_probability_max_transfer"]),
        )
        gains: list[float] = []
        for probabilities in scenarios:
            expected_entropy = 0.0
            for probability, outcome in zip(probabilities, experiment["outcomes"]):
                multipliers = _reliability_adjusted_multipliers(
                    outcome["posterior_weights"], reliability
                )
                posterior = _apply(prior, multipliers)
                expected_entropy += probability * _entropy(posterior)
            gains.append(max(0.0, prior_entropy - expected_entropy))
        interval = (min(gains), max(gains))
        self._gain_cache[key] = interval
        return interval

    def _cell_summaries(
        self,
        settings: dict[str, tuple[float, float]],
    ) -> tuple[JointActionSummary, ...]:
        action_ids = list(self.experiment_by_id) + ["STOP"]
        regrets: dict[str, list[float]] = {action_id: [] for action_id in action_ids}
        lower_utilities: dict[str, list[float]] = {
            action_id: [] for action_id in action_ids
        }

        for prior_index, _prior in enumerate(self._priors):
            point_lower: dict[str, float] = {}
            point_upper: dict[str, float] = {}
            for experiment_id, experiment in self.experiment_by_id.items():
                if experiment_id in settings:
                    cost_multiplier, reliability = settings[experiment_id]
                else:
                    cost_multiplier, reliability = 1.0, 1.0
                gain_lower, gain_upper = self._gain_interval(
                    prior_index, experiment_id, reliability
                )
                cost = float(experiment["cost"]) * float(cost_multiplier)
                point_lower[experiment_id] = gain_lower / cost
                point_upper[experiment_id] = gain_upper / cost
            point_lower["STOP"] = 0.0
            point_upper["STOP"] = 0.0
            best_possible = max(point_upper.values())
            for action_id in action_ids:
                lower = point_lower[action_id]
                regrets[action_id].append(max(0.0, best_possible - lower))
                lower_utilities[action_id].append(lower)

        return tuple(
            JointActionSummary(
                action_id=action_id,
                max_regret=max(regrets[action_id]),
                mean_regret=fmean(regrets[action_id]),
                worst_utility=min(lower_utilities[action_id]),
                mean_lower_utility=fmean(lower_utilities[action_id]),
            )
            for action_id in action_ids
        )

    @staticmethod
    def _rank_summaries(
        summaries: tuple[JointActionSummary, ...]
    ) -> list[JointActionSummary]:
        return sorted(
            summaries,
            key=lambda item: (
                item.max_regret,
                item.mean_regret,
                -item.worst_utility,
                -item.mean_lower_utility,
                item.action_id,
            ),
        )

    def _cell(self, settings: dict[str, tuple[float, float]]) -> JointCell:
        ranked = self._rank_summaries(self._cell_summaries(settings))
        winner, runner_up = ranked[0], ranked[1]
        packed = tuple(
            (experiment_id, float(settings[experiment_id][0]), float(settings[experiment_id][1]))
            for experiment_id in self.joint_ids
        )
        return JointCell(
            settings=packed,
            winner_action_id=winner.action_id,
            winner_max_regret=winner.max_regret,
            runner_up_action_id=runner_up.action_id,
            runner_up_max_regret=runner_up.max_regret,
            regret_margin=max(0.0, runner_up.max_regret - winner.max_regret),
        )

    def _all_setting_states(self) -> tuple[tuple[float, float], ...]:
        return tuple(
            (float(cost), float(reliability))
            for cost in self.config["cost_multiplier_grid"]
            for reliability in self.config["reliability_grid"]
        )

    def _cell_key(self, cell: JointCell) -> tuple[float, ...]:
        flattened: list[float] = []
        for _experiment_id, cost, reliability in cell.settings:
            flattened.extend((cost, reliability))
        return tuple(flattened)

    def _boundary_counts(
        self,
        cells_by_key: dict[tuple[float, ...], JointCell],
    ) -> tuple[dict[str, int], dict[str, int]]:
        costs = [float(value) for value in self.config["cost_multiplier_grid"]]
        reliabilities = [float(value) for value in self.config["reliability_grid"]]
        axis_values: list[tuple[str, list[float]]] = []
        for experiment_id in self.joint_ids:
            axis_values.append((f"{experiment_id}.cost", costs))
            axis_values.append((f"{experiment_id}.reliability", reliabilities))

        pair_counts: Counter[str] = Counter()
        axis_counts: Counter[str] = Counter()
        for key, cell in cells_by_key.items():
            for axis_index, (axis_name, values) in enumerate(axis_values):
                current = key[axis_index]
                position = values.index(current)
                if position + 1 >= len(values):
                    continue
                neighbor_key = list(key)
                neighbor_key[axis_index] = values[position + 1]
                neighbor = cells_by_key.get(tuple(neighbor_key))
                if neighbor is None or neighbor.winner_action_id == cell.winner_action_id:
                    continue
                pair = "<->".join(sorted((cell.winner_action_id, neighbor.winner_action_id)))
                pair_counts[pair] += 1
                axis_counts[axis_name] += 1
        return dict(sorted(pair_counts.items())), dict(sorted(axis_counts.items()))

    def _reference_settings(self) -> dict[str, tuple[float, float]]:
        return {
            experiment_id: (
                float(self.config["pareto_reference_settings"][experiment_id]["cost_multiplier"]),
                float(self.config["pareto_reference_settings"][experiment_id]["reliability"]),
            )
            for experiment_id in self.joint_ids
        }

    def _modified_plan(
        self, settings: dict[str, tuple[float, float]]
    ) -> dict[str, Any]:
        plan = copy.deepcopy(self.plan)
        for experiment in plan["experiments"]:
            experiment_id = experiment["experiment_id"]
            if experiment_id not in settings:
                continue
            cost_multiplier, reliability = settings[experiment_id]
            experiment["cost"] = float(experiment["cost"]) * cost_multiplier
            for outcome in experiment["outcomes"]:
                outcome["posterior_weights"] = _reliability_adjusted_multipliers(
                    outcome["posterior_weights"], reliability
                )
        return plan

    @staticmethod
    def _policy_point(
        policy: Policy,
        *,
        horizon: int,
        prior_entropy: float,
        risk_by_action: dict[str, float],
    ) -> PolicyParetoPoint:
        first = policy.first_experiment_id or ""
        return PolicyParetoPoint(
            horizon=horizon,
            first_experiment_id=first,
            expected_information_gain_bits=policy.expected_information_gain(prior_entropy),
            expected_cost=policy.expected_cost,
            max_path_cost=policy.max_path_cost,
            first_step_max_regret=risk_by_action[first],
        )

    @staticmethod
    def _pareto_prune(points: list[PolicyParetoPoint]) -> tuple[PolicyParetoPoint, ...]:
        unique: dict[tuple[Any, ...], PolicyParetoPoint] = {}
        for point in points:
            key = (
                point.horizon,
                point.first_experiment_id,
                round(point.expected_information_gain_bits, 12),
                round(point.expected_cost, 12),
                round(point.max_path_cost, 12),
                round(point.first_step_max_regret, 12),
            )
            unique[key] = point
        candidates = list(unique.values())
        kept: list[PolicyParetoPoint] = []
        for candidate in candidates:
            dominated = False
            for other in candidates:
                if other is candidate:
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
                if no_worse and strict:
                    dominated = True
                    break
            if not dominated:
                kept.append(candidate)
        return tuple(
            sorted(
                kept,
                key=lambda item: (
                    item.expected_cost,
                    item.first_step_max_regret,
                    -item.expected_information_gain_bits,
                    item.horizon,
                    item.first_experiment_id,
                ),
            )
        )

    def pareto_frontier(self) -> tuple[PolicyParetoPoint, ...]:
        settings = self._reference_settings()
        summaries = self._cell_summaries(settings)
        risk_by_action = {item.action_id: item.max_regret for item in summaries}
        modified_plan = self._modified_plan(settings)
        planner = StrategyPlanner(modified_plan)
        prior_entropy = _entropy(modified_plan["hypotheses"])
        budget = float(self.config["pareto_budget"])
        points: list[PolicyParetoPoint] = []
        for horizon in self.config["pareto_horizons"]:
            for policy in planner.frontier(budget=budget, max_depth=int(horizon)):
                if not policy.first_experiment_id:
                    continue
                points.append(
                    self._policy_point(
                        policy,
                        horizon=int(horizon),
                        prior_entropy=prior_entropy,
                        risk_by_action=risk_by_action,
                    )
                )
        return self._pareto_prune(points)

    def analyze(self) -> JointAssayFrontierResult:
        states = self._all_setting_states()
        cells: list[JointCell] = []
        for combination in product(states, repeat=len(self.joint_ids)):
            settings = {
                experiment_id: state
                for experiment_id, state in zip(self.joint_ids, combination)
            }
            cells.append(self._cell(settings))

        cells_by_key = {self._cell_key(cell): cell for cell in cells}
        pair_counts, axis_counts = self._boundary_counts(cells_by_key)
        winner_counts = dict(
            sorted(Counter(cell.winner_action_id for cell in cells).items())
        )
        return JointAssayFrontierResult(
            analysis_id=self.config["analysis_id"],
            prior_points_per_cell=len(self._priors),
            cells=tuple(cells),
            winner_counts=winner_counts,
            boundary_pair_counts=pair_counts,
            boundary_axis_counts=axis_counts,
            pareto_points=self.pareto_frontier(),
        )
