#!/usr/bin/env python3
"""Adaptive experiment replanning for CAUSAL-DNA.

This layer updates *planning weights only*. An outcome may change which
experiment is most informative next, but it cannot promote a causal claim,
materialize a Bardo state, or set cause_found. Those changes remain reserved
for the append-only evidence processor plus independent verification.
"""
from __future__ import annotations

import copy
import math
from dataclasses import dataclass
from typing import Any

from .experiment_selector import ExperimentScore, ExperimentSelector, ExperimentSelectorError


class AdaptiveReplannerError(ValueError):
    pass


def entropy_bits(weights: dict[str, float]) -> float:
    total = sum(weights.values())
    if total <= 0:
        raise AdaptiveReplannerError("planning weights must sum to > 0")
    value = 0.0
    for weight in weights.values():
        if weight < 0:
            raise AdaptiveReplannerError("planning weights must be non-negative")
        if weight:
            p = weight / total
            value -= p * math.log2(p)
    return value


@dataclass(frozen=True)
class ReplanSnapshot:
    generation: int
    planning_weights: dict[str, float]
    entropy_bits: float
    completed_experiments: tuple[str, ...]
    best_information_gain: ExperimentScore | None
    best_gain_per_cost: ExperimentScore | None


class AdaptiveReplanner:
    """Replay outcome records and recompute the next experiment ranking."""

    FORBIDDEN_OUTCOME_FIELDS = {
        "cause_found",
        "causal_status",
        "edge_status",
        "materialized",
        "verification_status",
    }

    def __init__(self, plan: dict[str, Any], outcomes: list[dict[str, Any]] | None = None):
        self.base_plan = copy.deepcopy(plan)
        ExperimentSelector(self.base_plan)
        self.outcomes = copy.deepcopy(outcomes or [])
        self._validate_outcomes()

    def _validate_outcomes(self) -> None:
        last_generation = 1
        seen_experiments: set[str] = set()
        for record in self.outcomes:
            forbidden = self.FORBIDDEN_OUTCOME_FIELDS.intersection(record)
            if forbidden:
                raise AdaptiveReplannerError(
                    "planning outcome cannot mutate causal state: " + ", ".join(sorted(forbidden))
                )
            generation = record.get("generation")
            if not isinstance(generation, int) or generation <= last_generation:
                raise AdaptiveReplannerError("replanning generations must strictly increase after generation 1")
            last_generation = generation
            experiment_id = record.get("experiment_id")
            if experiment_id in seen_experiments:
                raise AdaptiveReplannerError("an experiment outcome may only be applied once in a replan history")
            seen_experiments.add(experiment_id)
            if not isinstance(record.get("simulation_only"), bool):
                raise AdaptiveReplannerError("simulation_only must be explicit boolean")

            selector = ExperimentSelector(self.base_plan)
            try:
                selector.posterior_for_outcome(experiment_id, record.get("outcome_name"))
            except ExperimentSelectorError as exc:
                raise AdaptiveReplannerError(str(exc)) from exc

    def replay(self, *, objective: str = "information_gain") -> ReplanSnapshot:
        plan = copy.deepcopy(self.base_plan)
        completed: set[str] = set()
        generation = 1

        for record in self.outcomes:
            selector = ExperimentSelector(plan)
            new_weights = selector.posterior_for_outcome(
                record["experiment_id"], record["outcome_name"]
            )
            plan["hypotheses"] = new_weights
            completed.add(record["experiment_id"])
            generation = record["generation"]

        selector = ExperimentSelector(plan)
        remaining = [
            exp for exp in selector.experiments
            if exp["experiment_id"] not in completed
        ]
        if remaining:
            best_ig = selector.best(objective="information_gain", exclude=completed)
            best_cost = selector.best(objective="gain_per_cost", exclude=completed)
        else:
            best_ig = None
            best_cost = None

        return ReplanSnapshot(
            generation=generation,
            planning_weights=copy.deepcopy(plan["hypotheses"]),
            entropy_bits=entropy_bits(plan["hypotheses"]),
            completed_experiments=tuple(sorted(completed)),
            best_information_gain=best_ig,
            best_gain_per_cost=best_cost,
        )

    def append_outcome(self, record: dict[str, Any]) -> ReplanSnapshot:
        candidate = self.outcomes + [copy.deepcopy(record)]
        validated = AdaptiveReplanner(self.base_plan, candidate)
        self.outcomes = candidate
        return validated.replay()
