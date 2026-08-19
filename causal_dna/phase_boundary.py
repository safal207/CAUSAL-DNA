#!/usr/bin/env python3
"""Robust phase-boundary mapping for CAUSAL-DNA research decisions.

The analyzer maps where a *planning* decision changes as declared assumptions
about one target experiment change. It is intentionally epistemic: the output
is a decision boundary over planning assumptions, not a biological phase
transition, causal confidence, or evidence.

For each grid cell we vary target-experiment reliability and cost, allow a
bounded outcome-probability transfer, sweep a deterministic simplex of planning
priors, and choose the first action that minimizes worst-case interval regret.
Competitor experiments remain at their declared baseline reliability and cost.
"""
from __future__ import annotations

import copy
from collections import Counter
from dataclasses import dataclass
from statistics import fmean
from typing import Any, Iterator

from .distributional_robust import (
    DistributionalRobustError,
    _apply,
    _entropy,
    _probability_scenarios,
    _reliability_adjusted_multipliers,
)
from .experiment_selector import ExperimentSelector
from .sensitivity import _positive_compositions, grid_size


class PhaseBoundaryError(ValueError):
    pass


@dataclass(frozen=True)
class CellActionSummary:
    action_id: str
    max_regret: float
    mean_regret: float
    worst_utility: float
    mean_lower_utility: float


@dataclass(frozen=True)
class PhaseCell:
    target_reliability: float
    target_cost_multiplier: float
    outcome_probability_max_transfer: float
    winner_action_id: str
    winner_max_regret: float
    runner_up_action_id: str
    runner_up_max_regret: float
    regret_margin: float
    target_action_max_regret: float


@dataclass(frozen=True)
class PhaseBoundary:
    axis: str
    left_value: float
    right_value: float
    fixed_coordinates: tuple[tuple[str, float], ...]
    from_action_id: str
    to_action_id: str


@dataclass(frozen=True)
class PhaseBoundaryResult:
    analysis_id: str
    target_experiment_id: str
    prior_points_per_cell: int
    cells: tuple[PhaseCell, ...]
    boundaries: tuple[PhaseBoundary, ...]
    winner_counts: dict[str, int]

    def cell(
        self,
        *,
        reliability: float,
        cost_multiplier: float,
        probability_transfer: float,
    ) -> PhaseCell:
        for item in self.cells:
            if (
                abs(item.target_reliability - reliability) <= 1e-12
                and abs(item.target_cost_multiplier - cost_multiplier) <= 1e-12
                and abs(item.outcome_probability_max_transfer - probability_transfer) <= 1e-12
            ):
                return item
        raise PhaseBoundaryError("requested phase cell is not present")

    def target_cost_ceiling(
        self,
        *,
        reliability: float,
        probability_transfer: float,
    ) -> float | None:
        values = [
            item.target_cost_multiplier
            for item in self.cells
            if abs(item.target_reliability - reliability) <= 1e-12
            and abs(item.outcome_probability_max_transfer - probability_transfer) <= 1e-12
            and item.winner_action_id == self.target_experiment_id
        ]
        return max(values) if values else None


class RobustPhaseBoundaryAnalyzer:
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
        self.target_experiment_id = self.config.get("target_experiment_id")
        self.assert_valid()

    def assert_valid(self) -> None:
        forbidden = self.FORBIDDEN_CONFIG_FIELDS.intersection(self.config)
        if forbidden:
            raise PhaseBoundaryError(
                "phase-boundary config cannot mutate causal state: "
                + ", ".join(sorted(forbidden))
            )
        if self.config.get("case_id") != self.plan.get("case_id"):
            raise PhaseBoundaryError("case_id must match experiment plan")
        if self.config.get("decision_rule") != "minimax_interval_regret":
            raise PhaseBoundaryError("decision_rule must be minimax_interval_regret")
        if self.config.get("prior_model") != "deterministic_simplex_grid":
            raise PhaseBoundaryError("prior_model must be deterministic_simplex_grid")
        if self.config.get("include_stop") is not True:
            raise PhaseBoundaryError("include_stop must be true")

        experiment_ids = {
            item["experiment_id"] for item in self.plan.get("experiments", [])
        }
        if self.target_experiment_id not in experiment_ids:
            raise PhaseBoundaryError("target_experiment_id is not in experiment plan")

        total = self.config.get("prior_simplex_total_units")
        minimum = self.config.get("prior_min_units_per_hypothesis")
        if not isinstance(total, int) or total <= 0:
            raise PhaseBoundaryError("prior_simplex_total_units must be positive")
        if not isinstance(minimum, int) or minimum <= 0:
            raise PhaseBoundaryError("prior_min_units_per_hypothesis must be positive")
        size = grid_size(total, len(self.hypothesis_ids), minimum)
        if size < 10:
            raise PhaseBoundaryError("phase prior grid must contain at least 10 states")
        if size > 2000:
            raise PhaseBoundaryError("phase prior grid is capped at 2000 states")

        self._validate_grid("target_reliability_grid", lower=0.0, upper=1.0)
        self._validate_grid("target_cost_multiplier_grid", lower=0.0, strict_lower=True)
        self._validate_grid(
            "outcome_probability_max_transfer_grid", lower=0.0, upper=0.49
        )

    def _validate_grid(
        self,
        field: str,
        *,
        lower: float,
        upper: float | None = None,
        strict_lower: bool = False,
    ) -> None:
        values = self.config.get(field)
        if not isinstance(values, list) or len(values) < 2:
            raise PhaseBoundaryError(f"{field} must contain at least two values")
        if any(not isinstance(value, (int, float)) for value in values):
            raise PhaseBoundaryError(f"{field} values must be numeric")
        numeric = [float(value) for value in values]
        if numeric != sorted(set(numeric)):
            raise PhaseBoundaryError(f"{field} must be unique and ascending")
        for value in numeric:
            if strict_lower:
                if value <= lower:
                    raise PhaseBoundaryError(f"{field} values must be > {lower}")
            elif value < lower:
                raise PhaseBoundaryError(f"{field} values must be >= {lower}")
            if upper is not None and value > upper:
                raise PhaseBoundaryError(f"{field} values must be <= {upper}")

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

    def _experiment_interval(
        self,
        prior: dict[str, float],
        experiment: dict[str, Any],
        *,
        target_reliability: float,
        target_cost_multiplier: float,
        probability_transfer: float,
    ) -> tuple[float, float]:
        is_target = experiment["experiment_id"] == self.target_experiment_id
        reliability = target_reliability if is_target else 1.0
        cost_multiplier = target_cost_multiplier if is_target else 1.0
        cost = float(experiment["cost"]) * cost_multiplier
        prior_entropy = _entropy(prior)
        baseline_probabilities = [
            float(outcome["planning_probability"])
            for outcome in experiment["outcomes"]
        ]
        scenarios = _probability_scenarios(
            baseline_probabilities, probability_transfer
        )
        utilities: list[float] = []
        for probabilities in scenarios:
            expected_entropy = 0.0
            for probability, outcome in zip(probabilities, experiment["outcomes"]):
                multipliers = _reliability_adjusted_multipliers(
                    outcome["posterior_weights"], reliability
                )
                posterior = _apply(prior, multipliers)
                expected_entropy += probability * _entropy(posterior)
            gain = max(0.0, prior_entropy - expected_entropy)
            utilities.append(gain / cost)
        return min(utilities), max(utilities)

    def _cell(
        self,
        *,
        target_reliability: float,
        target_cost_multiplier: float,
        probability_transfer: float,
    ) -> PhaseCell:
        experiments = self.plan["experiments"]
        action_ids = [item["experiment_id"] for item in experiments] + ["STOP"]
        regrets: dict[str, list[float]] = {action_id: [] for action_id in action_ids}
        lower_utilities: dict[str, list[float]] = {
            action_id: [] for action_id in action_ids
        }

        for prior in self.prior_grid():
            intervals = {
                experiment["experiment_id"]: self._experiment_interval(
                    prior,
                    experiment,
                    target_reliability=target_reliability,
                    target_cost_multiplier=target_cost_multiplier,
                    probability_transfer=probability_transfer,
                )
                for experiment in experiments
            }
            point_lower = {
                action_id: interval[0] for action_id, interval in intervals.items()
            }
            point_upper = {
                action_id: interval[1] for action_id, interval in intervals.items()
            }
            point_lower["STOP"] = 0.0
            point_upper["STOP"] = 0.0
            best_possible = max(point_upper.values())
            for action_id in action_ids:
                lower = point_lower[action_id]
                regrets[action_id].append(max(0.0, best_possible - lower))
                lower_utilities[action_id].append(lower)

        summaries = [
            CellActionSummary(
                action_id=action_id,
                max_regret=max(regrets[action_id]),
                mean_regret=fmean(regrets[action_id]),
                worst_utility=min(lower_utilities[action_id]),
                mean_lower_utility=fmean(lower_utilities[action_id]),
            )
            for action_id in action_ids
        ]
        ranked = sorted(
            summaries,
            key=lambda item: (
                item.max_regret,
                item.mean_regret,
                -item.worst_utility,
                -item.mean_lower_utility,
                item.action_id,
            ),
        )
        winner, runner_up = ranked[0], ranked[1]
        target_summary = next(
            item for item in summaries if item.action_id == self.target_experiment_id
        )
        return PhaseCell(
            target_reliability=target_reliability,
            target_cost_multiplier=target_cost_multiplier,
            outcome_probability_max_transfer=probability_transfer,
            winner_action_id=winner.action_id,
            winner_max_regret=winner.max_regret,
            runner_up_action_id=runner_up.action_id,
            runner_up_max_regret=runner_up.max_regret,
            regret_margin=max(0.0, runner_up.max_regret - winner.max_regret),
            target_action_max_regret=target_summary.max_regret,
        )

    @staticmethod
    def _neighbors(
        cells: dict[tuple[float, float, float], PhaseCell],
        reliabilities: list[float],
        costs: list[float],
        transfers: list[float],
    ) -> tuple[PhaseBoundary, ...]:
        boundaries: list[PhaseBoundary] = []
        axes = (
            ("target_reliability", reliabilities, 0),
            ("target_cost_multiplier", costs, 1),
            ("outcome_probability_max_transfer", transfers, 2),
        )
        for axis_name, values, index in axes:
            for left, right in zip(values, values[1:]):
                for key, cell in cells.items():
                    if abs(key[index] - left) > 1e-12:
                        continue
                    other_key = list(key)
                    other_key[index] = right
                    other = cells.get(tuple(other_key))
                    if other is None or other.winner_action_id == cell.winner_action_id:
                        continue
                    fixed = []
                    coordinate_names = (
                        "target_reliability",
                        "target_cost_multiplier",
                        "outcome_probability_max_transfer",
                    )
                    for coordinate_index, coordinate_name in enumerate(coordinate_names):
                        if coordinate_index != index:
                            fixed.append((coordinate_name, key[coordinate_index]))
                    boundaries.append(
                        PhaseBoundary(
                            axis=axis_name,
                            left_value=left,
                            right_value=right,
                            fixed_coordinates=tuple(fixed),
                            from_action_id=cell.winner_action_id,
                            to_action_id=other.winner_action_id,
                        )
                    )
        return tuple(boundaries)

    def analyze(self) -> PhaseBoundaryResult:
        reliabilities = [float(v) for v in self.config["target_reliability_grid"]]
        costs = [float(v) for v in self.config["target_cost_multiplier_grid"]]
        transfers = [
            float(v) for v in self.config["outcome_probability_max_transfer_grid"]
        ]
        cells: list[PhaseCell] = []
        by_key: dict[tuple[float, float, float], PhaseCell] = {}
        counts: Counter[str] = Counter()
        for reliability in reliabilities:
            for cost in costs:
                for transfer in transfers:
                    cell = self._cell(
                        target_reliability=reliability,
                        target_cost_multiplier=cost,
                        probability_transfer=transfer,
                    )
                    cells.append(cell)
                    by_key[(reliability, cost, transfer)] = cell
                    counts[cell.winner_action_id] += 1

        prior_points = grid_size(
            int(self.config["prior_simplex_total_units"]),
            len(self.hypothesis_ids),
            int(self.config["prior_min_units_per_hypothesis"]),
        )
        return PhaseBoundaryResult(
            analysis_id=self.config["analysis_id"],
            target_experiment_id=self.target_experiment_id,
            prior_points_per_cell=prior_points,
            cells=tuple(cells),
            boundaries=self._neighbors(by_key, reliabilities, costs, transfers),
            winner_counts=dict(sorted(counts.items())),
        )
