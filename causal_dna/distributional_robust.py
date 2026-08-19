#!/usr/bin/env python3
"""Distributionally robust first-step planning for CAUSAL-DNA.

This module stress-tests experiment choice under a rectangular uncertainty set:
planning priors, outcome probabilities, assay reliability, and experiment cost
may all vary inside declared bounds. The result is a conservative planning
robustness measure only. It cannot modify causal status or create biological
evidence.
"""
from __future__ import annotations

import copy
import math
from dataclasses import dataclass
from statistics import fmean
from typing import Any, Iterable

from .experiment_selector import ExperimentSelector
from .sensitivity import StrategySensitivityAnalyzer


class DistributionalRobustError(ValueError):
    pass


def _entropy(weights: dict[str, float]) -> float:
    total = sum(weights.values())
    if total <= 0:
        raise DistributionalRobustError("hypothesis weights must sum to > 0")
    value = 0.0
    for weight in weights.values():
        if weight < 0:
            raise DistributionalRobustError("hypothesis weights must be non-negative")
        if weight == 0:
            continue
        p = weight / total
        value -= p * math.log2(p)
    return value


def _normalize(values: Iterable[float]) -> list[float]:
    values = [float(v) for v in values]
    total = sum(values)
    if total <= 0:
        raise DistributionalRobustError("probability vector must sum to > 0")
    return [value / total for value in values]


def _probability_scenarios(
    baseline: list[float], max_transfer: float
) -> tuple[tuple[float, ...], ...]:
    """Return baseline plus bounded pairwise probability-mass transfers."""
    base = tuple(_normalize(baseline))
    scenarios = {base}
    n = len(base)
    for source in range(n):
        for target in range(n):
            if source == target:
                continue
            amount = min(max_transfer, base[source])
            if amount <= 0:
                continue
            shifted = list(base)
            shifted[source] -= amount
            shifted[target] += amount
            scenarios.add(tuple(round(value, 12) for value in shifted))
    return tuple(sorted(scenarios))


def _reliability_adjusted_multipliers(
    multipliers: dict[str, float], reliability: float
) -> dict[str, float]:
    """Shrink likelihood-like multipliers toward neutral 1.0 as reliability falls."""
    if not (0 <= reliability <= 1):
        raise DistributionalRobustError("reliability must be in [0, 1]")
    adjusted: dict[str, float] = {}
    for key, multiplier in multipliers.items():
        value = float(multiplier)
        if value < 0:
            raise DistributionalRobustError("posterior multiplier must be non-negative")
        adjusted[key] = 1.0 + reliability * (value - 1.0)
    return adjusted


def _apply(prior: dict[str, float], multipliers: dict[str, float]) -> dict[str, float]:
    if set(prior) != set(multipliers):
        raise DistributionalRobustError("multipliers must cover all hypotheses")
    updated = {key: prior[key] * multipliers[key] for key in prior}
    if sum(updated.values()) <= 0:
        raise DistributionalRobustError("stressed outcome eliminates all hypotheses")
    return updated


@dataclass(frozen=True)
class UtilityInterval:
    action_id: str
    lower_utility: float
    upper_utility: float
    lower_information_gain_bits: float
    upper_information_gain_bits: float
    min_cost: float
    max_cost: float


@dataclass(frozen=True)
class DistributionalActionSummary:
    action_id: str
    worst_case_utility: float
    best_case_utility: float
    mean_lower_utility: float
    max_interval_regret: float
    mean_interval_regret: float
    robust_optimal_count: int


@dataclass(frozen=True)
class DistributionalRobustResult:
    analysis_id: str
    prior_points: int
    uncertainty_model: str
    minimax_interval_regret_action_id: str
    maximin_utility_action_id: str
    actions: tuple[DistributionalActionSummary, ...]

    def action(self, action_id: str) -> DistributionalActionSummary:
        for action in self.actions:
            if action.action_id == action_id:
                return action
        raise DistributionalRobustError(f"unknown action: {action_id}")


class DistributionallyRobustAnalyzer:
    FORBIDDEN_CONFIG_FIELDS = {
        "cause_found",
        "causal_status",
        "edge_status",
        "materialized",
        "verification_status",
    }

    def __init__(
        self,
        plan: dict[str, Any],
        config: dict[str, Any],
        sensitivity_config: dict[str, Any],
    ):
        self.plan = copy.deepcopy(plan)
        self.config = copy.deepcopy(config)
        self.sensitivity_config = copy.deepcopy(sensitivity_config)
        ExperimentSelector(self.plan)
        self.assert_valid()

    def assert_valid(self) -> None:
        forbidden = self.FORBIDDEN_CONFIG_FIELDS.intersection(self.config)
        if forbidden:
            raise DistributionalRobustError(
                "distributional config cannot mutate causal state: "
                + ", ".join(sorted(forbidden))
            )
        if self.config.get("case_id") != self.plan.get("case_id"):
            raise DistributionalRobustError("case_id must match experiment plan")
        if self.config.get("sensitivity_analysis_ref") != self.sensitivity_config.get("analysis_id"):
            raise DistributionalRobustError("sensitivity_analysis_ref mismatch")
        if self.config.get("objective") != "information_gain_per_cost":
            raise DistributionalRobustError("only information_gain_per_cost is supported")
        if self.config.get("uncertainty_model") != "rectangular_interval":
            raise DistributionalRobustError("uncertainty_model must be rectangular_interval")

        transfer = self.config.get("outcome_probability_max_transfer")
        if not isinstance(transfer, (int, float)) or not (0 <= transfer < 0.5):
            raise DistributionalRobustError("outcome_probability_max_transfer must be in [0, 0.5)")

        rmin = self.config.get("assay_reliability_min")
        rmax = self.config.get("assay_reliability_max")
        if not all(isinstance(v, (int, float)) for v in (rmin, rmax)):
            raise DistributionalRobustError("assay reliability bounds must be numeric")
        if not (0 <= rmin <= rmax <= 1):
            raise DistributionalRobustError("assay reliability bounds must satisfy 0 <= min <= max <= 1")

        cmin = self.config.get("cost_multiplier_min")
        cmax = self.config.get("cost_multiplier_max")
        if not all(isinstance(v, (int, float)) for v in (cmin, cmax)):
            raise DistributionalRobustError("cost multiplier bounds must be numeric")
        if not (0 < cmin <= cmax):
            raise DistributionalRobustError("cost multiplier bounds must satisfy 0 < min <= max")

        max_cost = self.config.get("max_baseline_first_step_cost")
        if not isinstance(max_cost, (int, float)) or max_cost <= 0:
            raise DistributionalRobustError("max_baseline_first_step_cost must be > 0")
        if self.config.get("include_stop") is not True:
            raise DistributionalRobustError("include_stop must be true")

    def _actions(self) -> list[dict[str, Any]]:
        max_cost = float(self.config["max_baseline_first_step_cost"])
        actions = [
            experiment
            for experiment in self.plan["experiments"]
            if float(experiment["cost"]) <= max_cost + 1e-12
        ]
        if not actions:
            raise DistributionalRobustError("no affordable experiments")
        return actions

    def interval_for_experiment(
        self, prior: dict[str, float], experiment: dict[str, Any]
    ) -> UtilityInterval:
        prior_entropy = _entropy(prior)
        baseline_probabilities = [
            float(outcome["planning_probability"]) for outcome in experiment["outcomes"]
        ]
        probability_scenarios = _probability_scenarios(
            baseline_probabilities,
            float(self.config["outcome_probability_max_transfer"]),
        )
        reliabilities = sorted(
            {
                float(self.config["assay_reliability_min"]),
                float(self.config["assay_reliability_max"]),
            }
        )
        cost_multipliers = sorted(
            {
                float(self.config["cost_multiplier_min"]),
                float(self.config["cost_multiplier_max"]),
            }
        )

        gains: list[float] = []
        utilities: list[float] = []
        costs: list[float] = []
        baseline_cost = float(experiment["cost"])

        for probabilities in probability_scenarios:
            for reliability in reliabilities:
                expected_entropy = 0.0
                for probability, outcome in zip(probabilities, experiment["outcomes"]):
                    multipliers = _reliability_adjusted_multipliers(
                        outcome["posterior_weights"], reliability
                    )
                    posterior = _apply(prior, multipliers)
                    expected_entropy += probability * _entropy(posterior)
                gain = max(0.0, prior_entropy - expected_entropy)
                gains.append(gain)
                for cost_multiplier in cost_multipliers:
                    cost = baseline_cost * cost_multiplier
                    costs.append(cost)
                    utilities.append(gain / cost)

        return UtilityInterval(
            action_id=experiment["experiment_id"],
            lower_utility=min(utilities),
            upper_utility=max(utilities),
            lower_information_gain_bits=min(gains),
            upper_information_gain_bits=max(gains),
            min_cost=min(costs),
            max_cost=max(costs),
        )

    def analyze(self) -> DistributionalRobustResult:
        prior_grid = StrategySensitivityAnalyzer(self.plan, self.sensitivity_config)
        experiments = self._actions()
        action_ids = [experiment["experiment_id"] for experiment in experiments] + ["STOP"]

        lower_utilities: dict[str, list[float]] = {action_id: [] for action_id in action_ids}
        upper_utilities: dict[str, list[float]] = {action_id: [] for action_id in action_ids}
        regrets: dict[str, list[float]] = {action_id: [] for action_id in action_ids}
        robust_optimal_counts = {action_id: 0 for action_id in action_ids}
        points = 0

        for prior in prior_grid.prior_grid():
            intervals = {
                experiment["experiment_id"]: self.interval_for_experiment(prior, experiment)
                for experiment in experiments
            }
            point_lower = {action_id: intervals[action_id].lower_utility for action_id in intervals}
            point_upper = {action_id: intervals[action_id].upper_utility for action_id in intervals}
            point_lower["STOP"] = 0.0
            point_upper["STOP"] = 0.0
            best_possible = max(point_upper.values())

            for action_id in action_ids:
                lower = point_lower[action_id]
                upper = point_upper[action_id]
                regret = max(0.0, best_possible - lower)
                lower_utilities[action_id].append(lower)
                upper_utilities[action_id].append(upper)
                regrets[action_id].append(regret)
                competitor_upper = max(
                    value for other, value in point_upper.items() if other != action_id
                )
                if lower >= competitor_upper - 1e-12:
                    robust_optimal_counts[action_id] += 1
            points += 1

        if points == 0:
            raise DistributionalRobustError("prior grid is empty")

        summaries: list[DistributionalActionSummary] = []
        for action_id in action_ids:
            summaries.append(
                DistributionalActionSummary(
                    action_id=action_id,
                    worst_case_utility=min(lower_utilities[action_id]),
                    best_case_utility=max(upper_utilities[action_id]),
                    mean_lower_utility=fmean(lower_utilities[action_id]),
                    max_interval_regret=max(regrets[action_id]),
                    mean_interval_regret=fmean(regrets[action_id]),
                    robust_optimal_count=robust_optimal_counts[action_id],
                )
            )

        minimax = min(
            summaries,
            key=lambda item: (
                item.max_interval_regret,
                item.mean_interval_regret,
                -item.worst_case_utility,
                item.action_id,
            ),
        )
        maximin = max(
            summaries,
            key=lambda item: (
                item.worst_case_utility,
                -item.max_interval_regret,
                item.mean_lower_utility,
                item.action_id,
            ),
        )

        return DistributionalRobustResult(
            analysis_id=self.config["analysis_id"],
            prior_points=points,
            uncertainty_model=self.config["uncertainty_model"],
            minimax_interval_regret_action_id=minimax.action_id,
            maximin_utility_action_id=maximin.action_id,
            actions=tuple(sorted(summaries, key=lambda item: item.action_id)),
        )
