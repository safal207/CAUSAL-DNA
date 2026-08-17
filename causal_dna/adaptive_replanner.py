#!/usr/bin/env python3
"""Adaptive experiment replanning for CAUSAL-DNA.

This layer updates *planning weights only*. An outcome may change which
experiment is most informative next, but it cannot promote a causal claim,
materialize a Bardo state, or set cause_found. Those changes remain reserved
for the append-only evidence processor plus independent verification.

Two modes are supported:
- simulation replay: planning-only scenarios, never appended to evidence history;
- real event mode: an evidence-backed `experiment_outcome` is appended as the
  next causal-processor generation, then the planning projection is recomputed.
"""
from __future__ import annotations

import copy
import math
from dataclasses import dataclass
from typing import Any

from .experiment_selector import ExperimentScore, ExperimentSelector, ExperimentSelectorError
from .processor import CausalProcessor, Projection


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


@dataclass(frozen=True)
class EventReplanResult:
    generation: int
    source_event_id: str
    experiment_id: str
    outcome_name: str
    prior_entropy_bits: float
    posterior_entropy_bits: float
    planning_weights: dict[str, float]
    best_information_gain: ExperimentScore | None
    best_gain_per_cost: ExperimentScore | None
    causal_projection_before: Projection
    causal_projection_after: Projection


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

    def _materialized_plan(self) -> tuple[dict[str, Any], set[str], int]:
        plan = copy.deepcopy(self.base_plan)
        completed: set[str] = set()
        generation = 1
        for record in self.outcomes:
            selector = ExperimentSelector(plan)
            plan["hypotheses"] = selector.posterior_for_outcome(
                record["experiment_id"], record["outcome_name"]
            )
            completed.add(record["experiment_id"])
            generation = record["generation"]
        return plan, completed, generation

    def replay(self, *, objective: str = "information_gain") -> ReplanSnapshot:
        plan, completed, generation = self._materialized_plan()
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
        """Append a planning record only; this does not touch causal evidence history."""
        candidate = self.outcomes + [copy.deepcopy(record)]
        validated = AdaptiveReplanner(self.base_plan, candidate)
        self.outcomes = candidate
        return validated.replay()

    def replan_from_event(
        self,
        processor: CausalProcessor,
        event: dict[str, Any],
    ) -> EventReplanResult:
        """Append one real outcome event and recompute the planning projection.

        The event must be evidence-backed and target exactly the next processor
        generation. The processor projection is checked before/after to prove
        that replanning changed no hypothesis status or verification state.
        """
        if event.get("event_type") != "experiment_outcome":
            raise AdaptiveReplannerError("real adaptive replanning requires experiment_outcome event")
        if event.get("observer") != "experiment":
            raise AdaptiveReplannerError("experiment_outcome must be asserted by experiment observer")
        if not event.get("evidence_refs"):
            raise AdaptiveReplannerError("experiment_outcome requires evidence provenance")
        if event.get("generation") != processor.generation + 1:
            raise AdaptiveReplannerError("experiment_outcome must target next processor generation")

        payload = event.get("payload") or {}
        experiment_id = payload.get("experiment_id")
        outcome_name = payload.get("outcome_name")
        if not experiment_id or not outcome_name:
            raise AdaptiveReplannerError(
                "experiment_outcome payload requires experiment_id and outcome_name"
            )

        plan, completed, planning_generation = self._materialized_plan()
        if planning_generation != processor.generation:
            raise AdaptiveReplannerError(
                "planning generation must match authoritative processor generation before event"
            )
        if experiment_id in completed:
            raise AdaptiveReplannerError("completed experiment cannot be applied twice")

        selector = ExperimentSelector(plan)
        prior_entropy = entropy_bits(plan["hypotheses"])
        try:
            updated_weights = selector.posterior_for_outcome(experiment_id, outcome_name)
        except ExperimentSelectorError as exc:
            raise AdaptiveReplannerError(str(exc)) from exc

        before = processor.project()
        processor.next_generation([copy.deepcopy(event)])
        after = processor.project()

        # `experiment_outcome` is evidence-bearing, but it is not verification.
        # It must not change open/rejected/superseded/verified scientific state.
        if (
            before.open_hypotheses != after.open_hypotheses
            or before.rejected_hypotheses != after.rejected_hypotheses
            or before.superseded_hypotheses != after.superseded_hypotheses
            or before.verified_subjects != after.verified_subjects
        ):
            raise AdaptiveReplannerError("experiment outcome mutated causal projection")

        record = {
            "generation": event["generation"],
            "experiment_id": str(experiment_id),
            "outcome_name": str(outcome_name),
            "simulation_only": False,
            "note": f"Derived from authoritative event {event['event_id']}",
        }
        self.outcomes.append(record)

        updated_plan = copy.deepcopy(plan)
        updated_plan["hypotheses"] = updated_weights
        updated_selector = ExperimentSelector(updated_plan)
        completed = set(completed) | {str(experiment_id)}
        remaining = [
            exp for exp in updated_selector.experiments
            if exp["experiment_id"] not in completed
        ]
        if remaining:
            best_ig = updated_selector.best(objective="information_gain", exclude=completed)
            best_cost = updated_selector.best(objective="gain_per_cost", exclude=completed)
        else:
            best_ig = None
            best_cost = None

        return EventReplanResult(
            generation=after.generation,
            source_event_id=str(event["event_id"]),
            experiment_id=str(experiment_id),
            outcome_name=str(outcome_name),
            prior_entropy_bits=prior_entropy,
            posterior_entropy_bits=entropy_bits(updated_weights),
            planning_weights=copy.deepcopy(updated_weights),
            best_information_gain=best_ig,
            best_gain_per_cost=best_cost,
            causal_projection_before=before,
            causal_projection_after=after,
        )
