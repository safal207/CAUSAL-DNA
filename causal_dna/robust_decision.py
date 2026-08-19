#!/usr/bin/env python3
"""Robust first-step experiment choice for CAUSAL-DNA.

This layer evaluates each candidate first experiment across the same explicit
planning-prior grid used by sensitivity analysis. It reports minimax regret and
maximin information gain. These are planning robustness measures only; they are
not biological probabilities, evidence, or causal confidence.
"""
from __future__ import annotations

import copy
from dataclasses import dataclass
from statistics import fmean
from typing import Any

from .experiment_selector import ExperimentSelector
from .sensitivity import StrategySensitivityAnalyzer


class RobustDecisionError(ValueError):
    pass


@dataclass(frozen=True)
class ActionRobustness:
    action_id: str
    worst_case_gain_bits: float
    best_case_gain_bits: float
    mean_gain_bits: float
    max_regret_bits: float
    mean_regret_bits: float
    optimal_count: int


@dataclass(frozen=True)
class RobustDecisionResult:
    analysis_id: str
    grid_points: int
    regret_metric: str
    minimax_regret_action_id: str
    maximin_gain_action_id: str
    actions: tuple[ActionRobustness, ...]

    def action(self, action_id: str) -> ActionRobustness:
        for item in self.actions:
            if item.action_id == action_id:
                return item
        raise RobustDecisionError(f"unknown robust action: {action_id}")


class RobustFirstStepAnalyzer:
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
        self.assert_valid()

    def assert_valid(self) -> None:
        forbidden = self.FORBIDDEN_CONFIG_FIELDS.intersection(self.config)
        if forbidden:
            raise RobustDecisionError(
                "robust-decision config cannot mutate causal state: "
                + ", ".join(sorted(forbidden))
            )
        if self.config.get("case_id") != self.plan.get("case_id"):
            raise RobustDecisionError("robust-decision case_id must match experiment plan")
        if self.config.get("case_id") != self.sensitivity_config.get("case_id"):
            raise RobustDecisionError("robust-decision case_id must match sensitivity config")
        if self.config.get("sensitivity_analysis_ref") != self.sensitivity_config.get("analysis_id"):
            raise RobustDecisionError("sensitivity_analysis_ref does not match sensitivity config")
        if self.config.get("regret_metric") != "expected_information_gain_bits":
            raise RobustDecisionError("only expected_information_gain_bits regret is supported")
        max_cost = self.config.get("max_first_step_cost")
        if not isinstance(max_cost, (int, float)) or max_cost <= 0:
            raise RobustDecisionError("max_first_step_cost must be > 0")
        if self.config.get("include_stop") is not True:
            raise RobustDecisionError("include_stop must be true so doing nothing remains an explicit action")

    def _candidate_action_ids(self) -> list[str]:
        max_cost = float(self.config["max_first_step_cost"])
        actions = [
            exp["experiment_id"]
            for exp in self.plan.get("experiments", [])
            if float(exp.get("cost", 1.0)) <= max_cost + 1e-12
        ]
        if not actions:
            raise RobustDecisionError("no experiment is affordable under max_first_step_cost")
        return sorted(actions) + ["STOP"]

    def analyze(self) -> RobustDecisionResult:
        grid = StrategySensitivityAnalyzer(self.plan, self.sensitivity_config)
        action_ids = self._candidate_action_ids()
        gains: dict[str, list[float]] = {action_id: [] for action_id in action_ids}
        regrets: dict[str, list[float]] = {action_id: [] for action_id in action_ids}
        optimal_counts = {action_id: 0 for action_id in action_ids}
        points = 0

        for weights in grid.prior_grid():
            scenario_plan = copy.deepcopy(self.plan)
            scenario_plan["hypotheses"] = copy.deepcopy(weights)
            selector = ExperimentSelector(scenario_plan)
            point_gains = {
                score.experiment_id: score.expected_information_gain_bits
                for score in selector.rank(objective="information_gain")
                if score.experiment_id in action_ids
            }
            point_gains["STOP"] = 0.0
            best_gain = max(point_gains.values())

            for action_id in action_ids:
                gain = point_gains[action_id]
                regret = max(0.0, best_gain - gain)
                gains[action_id].append(gain)
                regrets[action_id].append(regret)
                if abs(regret) <= 1e-12:
                    optimal_counts[action_id] += 1
            points += 1

        if points == 0:
            raise RobustDecisionError("planning-prior grid is empty")

        summaries = []
        for action_id in action_ids:
            action_gains = gains[action_id]
            action_regrets = regrets[action_id]
            summaries.append(
                ActionRobustness(
                    action_id=action_id,
                    worst_case_gain_bits=min(action_gains),
                    best_case_gain_bits=max(action_gains),
                    mean_gain_bits=fmean(action_gains),
                    max_regret_bits=max(action_regrets),
                    mean_regret_bits=fmean(action_regrets),
                    optimal_count=optimal_counts[action_id],
                )
            )

        minimax = min(
            summaries,
            key=lambda item: (
                item.max_regret_bits,
                item.mean_regret_bits,
                -item.worst_case_gain_bits,
                item.action_id,
            ),
        )
        maximin = max(
            summaries,
            key=lambda item: (
                item.worst_case_gain_bits,
                -item.max_regret_bits,
                item.mean_gain_bits,
                item.action_id,
            ),
        )

        return RobustDecisionResult(
            analysis_id=self.config["analysis_id"],
            grid_points=points,
            regret_metric=self.config["regret_metric"],
            minimax_regret_action_id=minimax.action_id,
            maximin_gain_action_id=maximin.action_id,
            actions=tuple(sorted(summaries, key=lambda item: item.action_id)),
        )
