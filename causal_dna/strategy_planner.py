#!/usr/bin/env python3
"""Budgeted multi-step research strategy planner for CAUSAL-DNA.

The planner searches adaptive experiment policies over a short horizon. Each
branch may choose a different next experiment after observing a different
outcome. Planning weights and outcome probabilities remain advisory: this
module cannot promote causal status, materialize graph states, or verify a
claim.
"""
from __future__ import annotations

import copy
from dataclasses import dataclass
from itertools import product
from typing import Any

from .adaptive_replanner import entropy_bits
from .experiment_selector import ExperimentSelector


class StrategyPlannerError(ValueError):
    pass


@dataclass(frozen=True)
class Policy:
    first_experiment_id: str | None
    expected_terminal_entropy_bits: float
    expected_cost: float
    max_path_cost: float
    branches: tuple[tuple[str, "Policy | None"], ...] = ()

    def expected_information_gain(self, prior_entropy_bits: float) -> float:
        return max(0.0, prior_entropy_bits - self.expected_terminal_entropy_bits)

    def gain_per_cost(self, prior_entropy_bits: float) -> float:
        if self.expected_cost <= 0:
            return 0.0
        return self.expected_information_gain(prior_entropy_bits) / self.expected_cost

    def next_experiments(self) -> dict[str, str | None]:
        return {
            outcome: child.first_experiment_id if child is not None else None
            for outcome, child in self.branches
        }


@dataclass(frozen=True)
class StrategyResult:
    objective: str
    budget: float
    max_depth: int
    prior_entropy_bits: float
    first_experiment_id: str
    expected_terminal_entropy_bits: float
    expected_information_gain_bits: float
    expected_cost: float
    max_path_cost: float
    gain_per_cost: float
    next_experiments_by_outcome: dict[str, str | None]


class StrategyPlanner:
    """Enumerate and Pareto-prune adaptive policies up to depth 3."""

    FORBIDDEN_CONFIG_FIELDS = {
        "cause_found",
        "causal_status",
        "edge_status",
        "materialized",
        "verification_status",
    }

    def __init__(self, plan: dict[str, Any]):
        self.base_plan = copy.deepcopy(plan)
        ExperimentSelector(self.base_plan)

    def _selector(self, weights: dict[str, float]) -> ExperimentSelector:
        plan = copy.deepcopy(self.base_plan)
        plan["hypotheses"] = copy.deepcopy(weights)
        return ExperimentSelector(plan)

    @staticmethod
    def _prune(policies: list[Policy]) -> list[Policy]:
        """Keep the Pareto frontier in expected entropy / expected cost."""
        kept: list[Policy] = []
        for candidate in sorted(
            policies,
            key=lambda p: (
                p.expected_terminal_entropy_bits,
                p.expected_cost,
                p.first_experiment_id or "",
            ),
        ):
            dominated = False
            for other in policies:
                if other is candidate:
                    continue
                no_worse_entropy = (
                    other.expected_terminal_entropy_bits
                    <= candidate.expected_terminal_entropy_bits + 1e-12
                )
                no_worse_cost = other.expected_cost <= candidate.expected_cost + 1e-12
                strictly_better = (
                    other.expected_terminal_entropy_bits
                    < candidate.expected_terminal_entropy_bits - 1e-12
                    or other.expected_cost < candidate.expected_cost - 1e-12
                )
                if no_worse_entropy and no_worse_cost and strictly_better:
                    dominated = True
                    break
            if not dominated:
                kept.append(candidate)
        return kept

    def _frontier(
        self,
        weights: dict[str, float],
        completed: frozenset[str],
        *,
        budget: float,
        depth: int,
    ) -> list[Policy]:
        current_entropy = entropy_bits(weights)
        stop = Policy(
            first_experiment_id=None,
            expected_terminal_entropy_bits=current_entropy,
            expected_cost=0.0,
            max_path_cost=0.0,
        )
        if depth <= 0 or budget <= 0:
            return [stop]

        selector = self._selector(weights)
        candidates: list[Policy] = [stop]

        for experiment in selector.experiments:
            experiment_id = experiment["experiment_id"]
            cost = float(experiment["cost"])
            if experiment_id in completed or cost > budget + 1e-12:
                continue

            child_frontiers: list[tuple[dict[str, Any], list[Policy]]] = []
            for outcome in experiment["outcomes"]:
                posterior = selector.posterior_for_outcome(
                    experiment_id, outcome["name"]
                )
                children = self._frontier(
                    posterior,
                    completed | {experiment_id},
                    budget=budget - cost,
                    depth=depth - 1,
                )
                child_frontiers.append((outcome, children))

            for chosen_children in product(
                *(children for _, children in child_frontiers)
            ):
                expected_entropy = 0.0
                expected_future_cost = 0.0
                max_future_path_cost = 0.0
                branches: list[tuple[str, Policy | None]] = []

                for (outcome, _), child in zip(child_frontiers, chosen_children):
                    probability = float(outcome["planning_probability"])
                    expected_entropy += (
                        probability * child.expected_terminal_entropy_bits
                    )
                    expected_future_cost += probability * child.expected_cost
                    max_future_path_cost = max(
                        max_future_path_cost, child.max_path_cost
                    )
                    branches.append(
                        (outcome["name"], child if child.first_experiment_id else None)
                    )

                policy = Policy(
                    first_experiment_id=experiment_id,
                    expected_terminal_entropy_bits=expected_entropy,
                    expected_cost=cost + expected_future_cost,
                    max_path_cost=cost + max_future_path_cost,
                    branches=tuple(branches),
                )
                if policy.max_path_cost <= budget + 1e-12:
                    candidates.append(policy)

        return self._prune(candidates)

    def frontier(self, *, budget: float, max_depth: int = 2) -> list[Policy]:
        if not isinstance(budget, (int, float)) or budget <= 0:
            raise StrategyPlannerError("budget must be > 0")
        if not isinstance(max_depth, int) or not (1 <= max_depth <= 3):
            raise StrategyPlannerError("max_depth must be an integer in [1, 3]")
        return self._frontier(
            copy.deepcopy(self.base_plan["hypotheses"]),
            frozenset(),
            budget=float(budget),
            depth=max_depth,
        )

    def best(
        self,
        *,
        budget: float,
        max_depth: int = 2,
        objective: str = "information_gain",
    ) -> StrategyResult:
        prior_entropy = entropy_bits(self.base_plan["hypotheses"])
        policies = [p for p in self.frontier(budget=budget, max_depth=max_depth) if p.first_experiment_id]
        if not policies:
            raise StrategyPlannerError("budget cannot fund any experiment")

        if objective == "information_gain":
            best = min(
                policies,
                key=lambda p: (
                    p.expected_terminal_entropy_bits,
                    p.expected_cost,
                    p.first_experiment_id or "",
                ),
            )
        elif objective == "gain_per_cost":
            best = max(
                policies,
                key=lambda p: (
                    p.gain_per_cost(prior_entropy),
                    p.expected_information_gain(prior_entropy),
                    -(p.expected_cost),
                ),
            )
        else:
            raise StrategyPlannerError(f"unknown objective: {objective}")

        gain = best.expected_information_gain(prior_entropy)
        return StrategyResult(
            objective=objective,
            budget=float(budget),
            max_depth=max_depth,
            prior_entropy_bits=prior_entropy,
            first_experiment_id=best.first_experiment_id or "",
            expected_terminal_entropy_bits=best.expected_terminal_entropy_bits,
            expected_information_gain_bits=gain,
            expected_cost=best.expected_cost,
            max_path_cost=best.max_path_cost,
            gain_per_cost=best.gain_per_cost(prior_entropy),
            next_experiments_by_outcome=best.next_experiments(),
        )
