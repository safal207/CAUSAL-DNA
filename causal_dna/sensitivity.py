#!/usr/bin/env python3
"""Planning-prior sensitivity analysis for CAUSAL-DNA research strategies.

This module asks whether a selected research policy survives across a broad,
explicit grid of alternative *planning priors*. These priors remain advisory;
they are not biological probabilities, evidence, or causal confidence.
"""
from __future__ import annotations

import copy
from collections import Counter
from dataclasses import dataclass
from math import comb
from typing import Any, Iterator

from .adaptive_replanner import entropy_bits
from .strategy_planner import Policy, StrategyPlanner, StrategyResult


class SensitivityError(ValueError):
    pass


def _positive_compositions(total: int, parts: int, minimum: int) -> Iterator[tuple[int, ...]]:
    remaining = total - minimum * parts
    if remaining < 0:
        return

    def walk(index: int, left: int, prefix: list[int]) -> Iterator[tuple[int, ...]]:
        if index == parts - 1:
            yield tuple(prefix + [left + minimum])
            return
        for value in range(left + 1):
            yield from walk(index + 1, left - value, prefix + [value + minimum])

    yield from walk(0, remaining, [])


def grid_size(total_units: int, parts: int, min_units: int = 1) -> int:
    remaining = total_units - min_units * parts
    if remaining < 0:
        return 0
    return comb(remaining + parts - 1, parts - 1)


def _signature(first_experiment_id: str, branches: dict[str, str | None]) -> str:
    if first_experiment_id == "STOP":
        return "STOP"
    branch_text = ",".join(
        f"{name}->{next_id or 'STOP'}"
        for name, next_id in sorted(branches.items())
    )
    return f"{first_experiment_id}|{branch_text}"


def policy_signature(result: StrategyResult) -> str:
    return _signature(result.first_experiment_id, result.next_experiments_by_outcome)


@dataclass(frozen=True)
class SensitivityResult:
    analysis_id: str
    grid_points: int
    budget: float
    max_depth: int
    objective: str
    baseline_first_experiment_id: str
    baseline_policy_signature: str
    dominant_first_experiment_id: str
    dominant_first_step_support: float
    baseline_first_step_support: float
    dominant_policy_signature: str
    dominant_policy_support: float
    baseline_policy_support: float
    first_step_counts: dict[str, int]
    policy_counts: dict[str, int]
    first_step_robust: bool
    policy_signature_robust: bool


class StrategySensitivityAnalyzer:
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
        self.hypothesis_ids = tuple(self.plan.get("hypotheses", {}).keys())
        self.assert_valid()

    def assert_valid(self) -> None:
        forbidden = self.FORBIDDEN_CONFIG_FIELDS.intersection(self.config)
        if forbidden:
            raise SensitivityError(
                "sensitivity config cannot mutate causal state: " + ", ".join(sorted(forbidden))
            )
        if len(self.hypothesis_ids) < 2:
            raise SensitivityError("at least two hypotheses are required")
        total = self.config.get("simplex_total_units")
        minimum = self.config.get("min_units_per_hypothesis")
        if not isinstance(total, int) or total <= 0:
            raise SensitivityError("simplex_total_units must be a positive integer")
        if not isinstance(minimum, int) or minimum <= 0:
            raise SensitivityError("min_units_per_hypothesis must be a positive integer")
        size = grid_size(total, len(self.hypothesis_ids), minimum)
        if size < 1000:
            raise SensitivityError("sensitivity grid must contain at least 1000 prior states")
        if size > 10000:
            raise SensitivityError("sensitivity grid is intentionally capped at 10000 states")
        budget = self.config.get("budget")
        depth = self.config.get("max_depth")
        objective = self.config.get("objective")
        if not isinstance(budget, (int, float)) or budget <= 0:
            raise SensitivityError("budget must be > 0")
        if not isinstance(depth, int) or not (1 <= depth <= 3):
            raise SensitivityError("max_depth must be in [1, 3]")
        if objective not in {"information_gain", "gain_per_cost"}:
            raise SensitivityError("unsupported objective")
        for field in ("first_step_robust_threshold", "policy_signature_robust_threshold"):
            value = self.config.get(field)
            if not isinstance(value, (int, float)) or not (0 <= value <= 1):
                raise SensitivityError(f"{field} must be in [0, 1]")

    def prior_grid(self) -> Iterator[dict[str, float]]:
        total = self.config["simplex_total_units"]
        minimum = self.config["min_units_per_hypothesis"]
        for composition in _positive_compositions(total, len(self.hypothesis_ids), minimum):
            yield {
                hypothesis_id: float(units)
                for hypothesis_id, units in zip(self.hypothesis_ids, composition)
            }

    def _choose_policy(self, weights: dict[str, float]) -> tuple[str, str]:
        plan = copy.deepcopy(self.plan)
        plan["hypotheses"] = copy.deepcopy(weights)
        planner = StrategyPlanner(plan)
        prior_entropy = entropy_bits(weights)
        frontier = planner.frontier(
            budget=float(self.config["budget"]),
            max_depth=int(self.config["max_depth"]),
        )

        objective = self.config["objective"]
        if objective == "information_gain":
            best = min(
                frontier,
                key=lambda p: (
                    p.expected_terminal_entropy_bits,
                    p.expected_cost,
                    p.first_experiment_id or "",
                ),
            )
        else:
            best = max(
                frontier,
                key=lambda p: (
                    p.gain_per_cost(prior_entropy),
                    p.expected_information_gain(prior_entropy),
                    -p.expected_cost,
                    p.first_experiment_id or "",
                ),
            )

        if best.first_experiment_id is None:
            return "STOP", "STOP"
        branches = {
            outcome: child.first_experiment_id if child is not None else None
            for outcome, child in best.branches
        }
        return best.first_experiment_id, _signature(best.first_experiment_id, branches)

    def analyze(self) -> SensitivityResult:
        baseline = StrategyPlanner(self.plan).best(
            budget=float(self.config["budget"]),
            max_depth=int(self.config["max_depth"]),
            objective=self.config["objective"],
        )
        baseline_signature = policy_signature(baseline)

        first_steps: Counter[str] = Counter()
        policies: Counter[str] = Counter()
        points = 0
        for weights in self.prior_grid():
            first_step, signature = self._choose_policy(weights)
            first_steps[first_step] += 1
            policies[signature] += 1
            points += 1

        if points == 0:
            raise SensitivityError("sensitivity grid is empty")

        dominant_first, dominant_first_count = first_steps.most_common(1)[0]
        dominant_policy, dominant_policy_count = policies.most_common(1)[0]
        baseline_first_count = first_steps[baseline.first_experiment_id]
        baseline_policy_count = policies[baseline_signature]

        dominant_first_support = dominant_first_count / points
        dominant_policy_support = dominant_policy_count / points
        baseline_first_support = baseline_first_count / points
        baseline_policy_support = baseline_policy_count / points

        return SensitivityResult(
            analysis_id=self.config["analysis_id"],
            grid_points=points,
            budget=float(self.config["budget"]),
            max_depth=int(self.config["max_depth"]),
            objective=self.config["objective"],
            baseline_first_experiment_id=baseline.first_experiment_id,
            baseline_policy_signature=baseline_signature,
            dominant_first_experiment_id=dominant_first,
            dominant_first_step_support=dominant_first_support,
            baseline_first_step_support=baseline_first_support,
            dominant_policy_signature=dominant_policy,
            dominant_policy_support=dominant_policy_support,
            baseline_policy_support=baseline_policy_support,
            first_step_counts=dict(sorted(first_steps.items())),
            policy_counts=dict(sorted(policies.items())),
            first_step_robust=baseline_first_support >= float(self.config["first_step_robust_threshold"]),
            policy_signature_robust=baseline_policy_support >= float(self.config["policy_signature_robust_threshold"]),
        )
