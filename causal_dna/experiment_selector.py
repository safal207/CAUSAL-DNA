#!/usr/bin/env python3
"""Information-gain experiment selector for CAUSAL-DNA.

Planner probabilities are *not biological evidence*. They are explicit planning
priors used only to rank candidate experiments by expected reduction of
hypothesis uncertainty. Material causal confidence changes only through the
append-only evidence/verification pipeline.

`posterior_weights` in an experiment outcome are treated as likelihood-like
planning multipliers. They are multiplied by the current planning weights, so
ranking can be recomputed after each hypothetical or observed planning outcome.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any


class ExperimentSelectorError(ValueError):
    pass


def _entropy(weights: dict[str, float]) -> float:
    total = sum(weights.values())
    if total <= 0:
        raise ExperimentSelectorError("hypothesis weights must sum to > 0")
    h = 0.0
    for value in weights.values():
        if value < 0:
            raise ExperimentSelectorError("hypothesis weights must be non-negative")
        if value == 0:
            continue
        p = value / total
        h -= p * math.log2(p)
    return h


def _apply_multipliers(
    prior: dict[str, float], multipliers: dict[str, float]
) -> dict[str, float]:
    if set(prior) != set(multipliers):
        raise ExperimentSelectorError("posterior_weights must cover all hypotheses")
    updated = {key: prior[key] * multipliers[key] for key in prior}
    if sum(updated.values()) <= 0:
        raise ExperimentSelectorError("outcome multipliers eliminate all hypotheses")
    return updated


@dataclass(frozen=True)
class ExperimentScore:
    experiment_id: str
    prior_entropy_bits: float
    expected_posterior_entropy_bits: float
    expected_information_gain_bits: float
    normalized_information_gain: float
    cost: float
    utility_per_cost: float


class ExperimentSelector:
    def __init__(self, plan: dict[str, Any]):
        self.plan = plan
        self.hypotheses = plan.get("hypotheses", {})
        self.experiments = plan.get("experiments", [])
        self.assert_valid()

    def assert_valid(self) -> None:
        if not self.hypotheses:
            raise ExperimentSelectorError("at least one hypothesis is required")
        _entropy(self.hypotheses)
        known = set(self.hypotheses)
        ids: set[str] = set()
        for exp in self.experiments:
            eid = exp.get("experiment_id")
            if not eid or eid in ids:
                raise ExperimentSelectorError("experiment_id values must be unique and non-empty")
            ids.add(eid)
            cost = exp.get("cost", 1.0)
            if not isinstance(cost, (int, float)) or cost <= 0:
                raise ExperimentSelectorError(f"{eid}: cost must be > 0")
            outcomes = exp.get("outcomes", [])
            if len(outcomes) < 2:
                raise ExperimentSelectorError(f"{eid}: at least two outcomes are required")
            outcome_probability_sum = 0.0
            for outcome in outcomes:
                probability = outcome.get("planning_probability")
                if not isinstance(probability, (int, float)) or not (0 <= probability <= 1):
                    raise ExperimentSelectorError(f"{eid}: invalid planning_probability")
                outcome_probability_sum += probability
                multipliers = outcome.get("posterior_weights", {})
                if set(multipliers) != known:
                    raise ExperimentSelectorError(f"{eid}: posterior_weights must cover all hypotheses")
                if any(not isinstance(v, (int, float)) or v < 0 for v in multipliers.values()):
                    raise ExperimentSelectorError(f"{eid}: posterior_weights must be non-negative numbers")
                _apply_multipliers(self.hypotheses, multipliers)
            if not math.isclose(outcome_probability_sum, 1.0, abs_tol=1e-9):
                raise ExperimentSelectorError(f"{eid}: planning probabilities must sum to 1")

    def experiment(self, experiment_id: str) -> dict[str, Any]:
        for experiment in self.experiments:
            if experiment["experiment_id"] == experiment_id:
                return experiment
        raise ExperimentSelectorError(f"unknown experiment: {experiment_id}")

    def posterior_for_outcome(
        self, experiment_id: str, outcome_name: str
    ) -> dict[str, float]:
        experiment = self.experiment(experiment_id)
        for outcome in experiment["outcomes"]:
            if outcome["name"] == outcome_name:
                return _apply_multipliers(self.hypotheses, outcome["posterior_weights"])
        raise ExperimentSelectorError(
            f"unknown outcome {outcome_name!r} for experiment {experiment_id}"
        )

    def score(self, experiment: dict[str, Any]) -> ExperimentScore:
        prior_h = _entropy(self.hypotheses)
        expected_h = 0.0
        for outcome in experiment["outcomes"]:
            updated = _apply_multipliers(self.hypotheses, outcome["posterior_weights"])
            expected_h += outcome["planning_probability"] * _entropy(updated)
        gain = max(0.0, prior_h - expected_h)
        normalized = gain / prior_h if prior_h else 0.0
        cost = float(experiment.get("cost", 1.0))
        return ExperimentScore(
            experiment_id=experiment["experiment_id"],
            prior_entropy_bits=prior_h,
            expected_posterior_entropy_bits=expected_h,
            expected_information_gain_bits=gain,
            normalized_information_gain=normalized,
            cost=cost,
            utility_per_cost=gain / cost,
        )

    def rank(
        self,
        *,
        objective: str = "information_gain",
        exclude: set[str] | None = None,
    ) -> list[ExperimentScore]:
        excluded = exclude or set()
        scores = [
            self.score(exp)
            for exp in self.experiments
            if exp["experiment_id"] not in excluded
        ]
        if objective == "information_gain":
            return sorted(scores, key=lambda s: (-s.expected_information_gain_bits, s.cost, s.experiment_id))
        if objective == "gain_per_cost":
            return sorted(scores, key=lambda s: (-s.utility_per_cost, -s.expected_information_gain_bits, s.experiment_id))
        raise ExperimentSelectorError(f"unknown objective: {objective}")

    def best(
        self,
        *,
        objective: str = "information_gain",
        exclude: set[str] | None = None,
    ) -> ExperimentScore:
        ranked = self.rank(objective=objective, exclude=exclude)
        if not ranked:
            raise ExperimentSelectorError("no experiments available")
        return ranked[0]
