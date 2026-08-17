#!/usr/bin/env python3
"""Phase-boundary mapping for CAUSAL-DNA experiment choice.

This layer maps where the preferred *planning action* changes as declared
assumptions vary. It is intentionally planning-only. A phase boundary is a
switch in experiment recommendation, not a biological transition and not new
evidence.

The map uses conservative one-step utility:

    utility = min_outcome_probability_shift(EIG) / effective_cost

The target experiment may vary in cost and assay reliability while outcome
probability uncertainty is swept globally. For every phase point we report the
winner under the baseline prior and the dominant winner across the same 1,140
planning-prior grid used by sensitivity analysis.
"""
from __future__ import annotations

import copy
from collections import Counter
from dataclasses import dataclass
from statistics import fmean
from typing import Any

from .distributional_robust import (
    DistributionalRobustError,
    _apply,
    _entropy,
    _probability_scenarios,
    _reliability_adjusted_multipliers,
)
from .experiment_selector import ExperimentSelector
from .sensitivity import StrategySensitivityAnalyzer


class PhaseMapError(ValueError):
    pass


@dataclass(frozen=True)
class PhasePoint:
    target_cost_multiplier: float
    target_reliability: float
    outcome_probability_max_transfer: float
    baseline_winner: str
    dominant_winner: str
    dominant_support: float
    winner_counts: dict[str, int]


@dataclass(frozen=True)
class BreakEvenPoint:
    target_reliability: float
    outcome_probability_max_transfer: float
    competitor_id: str
    competitor_utility: float
    target_conservative_gain_bits: float
    break_even_cost_multiplier: float | None
    break_even_absolute_cost: float | None


@dataclass(frozen=True)
class PhaseMapResult:
    analysis_id: str
    prior_points: int
    target_experiment_id: str
    phase_points: tuple[PhasePoint, ...]
    break_even_points: tuple[BreakEvenPoint, ...]


class PhaseBoundaryAnalyzer:
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
            raise PhaseMapError(
                "phase-map config cannot mutate causal state: "
                + ", ".join(sorted(forbidden))
            )
        if self.config.get("case_id") != self.plan.get("case_id"):
            raise PhaseMapError("case_id must match experiment plan")
        if self.config.get("sensitivity_analysis_ref") != self.sensitivity_config.get("analysis_id"):
            raise PhaseMapError("sensitivity_analysis_ref mismatch")
        if self.config.get("objective") != "conservative_information_gain_per_cost":
            raise PhaseMapError("unsupported phase-map objective")
        target = self.config.get("target_experiment_id")
        if target not in {e["experiment_id"] for e in self.plan.get("experiments", [])}:
            raise PhaseMapError("target_experiment_id is not in experiment plan")
        if self.config.get("include_stop") is not True:
            raise PhaseMapError("include_stop must be true")

        for field in (
            "target_cost_multipliers",
            "target_reliabilities",
            "outcome_probability_max_transfers",
        ):
            values = self.config.get(field)
            if not isinstance(values, list) or not values:
                raise PhaseMapError(f"{field} must be a non-empty list")
            if values != sorted(set(values)):
                raise PhaseMapError(f"{field} must be sorted and unique")

        if any(not isinstance(v, (int, float)) or v <= 0 for v in self.config["target_cost_multipliers"]):
            raise PhaseMapError("target_cost_multipliers must be > 0")
        if any(not isinstance(v, (int, float)) or not (0 <= v <= 1) for v in self.config["target_reliabilities"]):
            raise PhaseMapError("target_reliabilities must be in [0, 1]")
        if any(not isinstance(v, (int, float)) or not (0 <= v < 0.5) for v in self.config["outcome_probability_max_transfers"]):
            raise PhaseMapError("outcome_probability_max_transfers must be in [0, 0.5)")

        competitor_reliability = self.config.get("competitor_reliability")
        if not isinstance(competitor_reliability, (int, float)) or not (0 <= competitor_reliability <= 1):
            raise PhaseMapError("competitor_reliability must be in [0, 1]")

    def _experiment(self, experiment_id: str) -> dict[str, Any]:
        for experiment in self.plan["experiments"]:
            if experiment["experiment_id"] == experiment_id:
                return experiment
        raise PhaseMapError(f"unknown experiment: {experiment_id}")

    @staticmethod
    def conservative_gain(
        prior: dict[str, float],
        experiment: dict[str, Any],
        *,
        reliability: float,
        max_transfer: float,
    ) -> float:
        prior_entropy = _entropy(prior)
        base_probabilities = [float(o["planning_probability"]) for o in experiment["outcomes"]]
        gains: list[float] = []
        for probabilities in _probability_scenarios(base_probabilities, max_transfer):
            expected_entropy = 0.0
            for probability, outcome in zip(probabilities, experiment["outcomes"]):
                multipliers = _reliability_adjusted_multipliers(
                    outcome["posterior_weights"], reliability
                )
                posterior = _apply(prior, multipliers)
                expected_entropy += probability * _entropy(posterior)
            gains.append(max(0.0, prior_entropy - expected_entropy))
        return min(gains)

    def _utilities(
        self,
        prior: dict[str, float],
        *,
        target_cost_multiplier: float,
        target_reliability: float,
        max_transfer: float,
    ) -> dict[str, float]:
        target_id = self.config["target_experiment_id"]
        competitor_reliability = float(self.config["competitor_reliability"])
        utilities: dict[str, float] = {}
        for experiment in self.plan["experiments"]:
            is_target = experiment["experiment_id"] == target_id
            reliability = target_reliability if is_target else competitor_reliability
            gain = self.conservative_gain(
                prior,
                experiment,
                reliability=reliability,
                max_transfer=max_transfer,
            )
            cost_multiplier = target_cost_multiplier if is_target else 1.0
            effective_cost = float(experiment["cost"]) * cost_multiplier
            utilities[experiment["experiment_id"]] = gain / effective_cost
        utilities["STOP"] = 0.0
        return utilities

    @staticmethod
    def _winner(utilities: dict[str, float]) -> str:
        return min(
            utilities,
            key=lambda action_id: (-utilities[action_id], action_id),
        )

    def _break_even(
        self,
        baseline_prior: dict[str, float],
        *,
        target_reliability: float,
        max_transfer: float,
    ) -> BreakEvenPoint:
        target_id = self.config["target_experiment_id"]
        target = self._experiment(target_id)
        target_gain = self.conservative_gain(
            baseline_prior,
            target,
            reliability=target_reliability,
            max_transfer=max_transfer,
        )

        competitor_utilities: dict[str, float] = {"STOP": 0.0}
        for experiment in self.plan["experiments"]:
            if experiment["experiment_id"] == target_id:
                continue
            gain = self.conservative_gain(
                baseline_prior,
                experiment,
                reliability=float(self.config["competitor_reliability"]),
                max_transfer=max_transfer,
            )
            competitor_utilities[experiment["experiment_id"]] = gain / float(experiment["cost"])

        competitor_id = self._winner(competitor_utilities)
        competitor_utility = competitor_utilities[competitor_id]
        if competitor_utility <= 1e-15:
            multiplier = None
            absolute_cost = None
        else:
            absolute_cost = target_gain / competitor_utility
            multiplier = absolute_cost / float(target["cost"])

        return BreakEvenPoint(
            target_reliability=float(target_reliability),
            outcome_probability_max_transfer=float(max_transfer),
            competitor_id=competitor_id,
            competitor_utility=competitor_utility,
            target_conservative_gain_bits=target_gain,
            break_even_cost_multiplier=multiplier,
            break_even_absolute_cost=absolute_cost,
        )

    def analyze(self) -> PhaseMapResult:
        prior_analyzer = StrategySensitivityAnalyzer(self.plan, self.sensitivity_config)
        priors = list(prior_analyzer.prior_grid())
        if not priors:
            raise PhaseMapError("prior grid is empty")
        baseline_prior = copy.deepcopy(self.plan["hypotheses"])

        points: list[PhasePoint] = []
        for max_transfer in self.config["outcome_probability_max_transfers"]:
            for reliability in self.config["target_reliabilities"]:
                for cost_multiplier in self.config["target_cost_multipliers"]:
                    baseline_winner = self._winner(
                        self._utilities(
                            baseline_prior,
                            target_cost_multiplier=float(cost_multiplier),
                            target_reliability=float(reliability),
                            max_transfer=float(max_transfer),
                        )
                    )
                    counts: Counter[str] = Counter()
                    for prior in priors:
                        winner = self._winner(
                            self._utilities(
                                prior,
                                target_cost_multiplier=float(cost_multiplier),
                                target_reliability=float(reliability),
                                max_transfer=float(max_transfer),
                            )
                        )
                        counts[winner] += 1
                    dominant_winner, dominant_count = counts.most_common(1)[0]
                    points.append(
                        PhasePoint(
                            target_cost_multiplier=float(cost_multiplier),
                            target_reliability=float(reliability),
                            outcome_probability_max_transfer=float(max_transfer),
                            baseline_winner=baseline_winner,
                            dominant_winner=dominant_winner,
                            dominant_support=dominant_count / len(priors),
                            winner_counts=dict(sorted(counts.items())),
                        )
                    )

        break_evens = [
            self._break_even(
                baseline_prior,
                target_reliability=float(reliability),
                max_transfer=float(max_transfer),
            )
            for max_transfer in self.config["outcome_probability_max_transfers"]
            for reliability in self.config["target_reliabilities"]
        ]

        return PhaseMapResult(
            analysis_id=self.config["analysis_id"],
            prior_points=len(priors),
            target_experiment_id=self.config["target_experiment_id"],
            phase_points=tuple(points),
            break_even_points=tuple(break_evens),
        )
