#!/usr/bin/env python3
"""Adversarial robustness radius for CAUSAL-DNA planning decisions.

The radius measures how far declared planning assumptions must move from a
reference profile before the current minimax-regret first action changes.
It is a property of the decision model, not biological evidence or causal
confidence.
"""
from __future__ import annotations

import copy
from dataclasses import dataclass
from typing import Any

from .joint_assay_frontier import (
    JointAssayFrontierAnalyzer,
    JointAssayFrontierError,
    JointCell,
)


class AdversarialRadiusError(ValueError):
    pass


@dataclass(frozen=True)
class AxisThreshold:
    axis_id: str
    found: bool
    normalized_distance: float | None
    boundary_value: float | None
    winner_after_boundary: str | None
    regret_margin_after_boundary: float | None


@dataclass(frozen=True)
class AdversarialRadiusResult:
    analysis_id: str
    target_action_id: str
    baseline_winner_action_id: str
    baseline_regret_margin: float
    coordinated_switch_found: bool
    coordinated_radius: float | None
    coordinated_winner_after_boundary: str | None
    coordinated_regret_margin_after_boundary: float | None
    coordinated_settings: tuple[tuple[str, float, float], ...] | None
    axis_thresholds: tuple[AxisThreshold, ...]


class AdversarialRadiusAnalyzer:
    """Find first-action decision boundaries along declared adverse directions."""

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
        joint_config: dict[str, Any],
        radius_config: dict[str, Any],
    ):
        self.plan = copy.deepcopy(plan)
        self.joint_config = copy.deepcopy(joint_config)
        self.config = copy.deepcopy(radius_config)
        self.joint = JointAssayFrontierAnalyzer(self.plan, self.joint_config)
        self.target_action_id = self.config.get("target_action_id")
        self.axes = tuple(self.config.get("axes", ()))
        self.assert_valid()
        self.reference = self.joint._reference_settings()

    def assert_valid(self) -> None:
        forbidden = self.FORBIDDEN_CONFIG_FIELDS.intersection(self.config)
        if forbidden:
            raise AdversarialRadiusError(
                "radius config cannot mutate causal state: "
                + ", ".join(sorted(forbidden))
            )
        if self.config.get("case_id") != self.plan.get("case_id"):
            raise AdversarialRadiusError("case_id must match experiment plan")
        if self.config.get("distance_metric") != "weighted_linf":
            raise AdversarialRadiusError("distance_metric must be weighted_linf")
        if self.config.get("decision_rule") != "minimax_interval_regret":
            raise AdversarialRadiusError(
                "decision_rule must be minimax_interval_regret"
            )
        experiment_ids = {item["experiment_id"] for item in self.plan["experiments"]}
        if self.target_action_id not in experiment_ids:
            raise AdversarialRadiusError("target_action_id must be an experiment")
        if not isinstance(self.axes, tuple) or not self.axes:
            raise AdversarialRadiusError("axes must be a non-empty list")

        seen: set[str] = set()
        for axis in self.axes:
            if not isinstance(axis, dict):
                raise AdversarialRadiusError("every axis must be an object")
            axis_id = axis.get("axis_id")
            if not isinstance(axis_id, str) or not axis_id or axis_id in seen:
                raise AdversarialRadiusError("axis_id values must be unique strings")
            seen.add(axis_id)
            experiment_id = axis.get("experiment_id")
            if experiment_id not in self.joint.joint_ids:
                raise AdversarialRadiusError(
                    "radius axes must target independent joint experiment axes"
                )
            field = axis.get("field")
            if field not in {"cost_multiplier", "reliability"}:
                raise AdversarialRadiusError(
                    "axis field must be cost_multiplier or reliability"
                )
            direction = axis.get("direction")
            if direction not in {"increase", "decrease"}:
                raise AdversarialRadiusError(
                    "axis direction must be increase or decrease"
                )
            scale = axis.get("scale")
            limit = axis.get("limit")
            if not isinstance(scale, (int, float)) or scale <= 0:
                raise AdversarialRadiusError("axis scale must be > 0")
            if not isinstance(limit, (int, float)):
                raise AdversarialRadiusError("axis limit must be numeric")
            if field == "cost_multiplier" and limit <= 0:
                raise AdversarialRadiusError("cost limit must be > 0")
            if field == "reliability" and not (0 <= limit <= 1):
                raise AdversarialRadiusError("reliability limit must be in [0, 1]")

        tolerance = self.config.get("binary_tolerance")
        if not isinstance(tolerance, (int, float)) or not (0 < tolerance <= 0.05):
            raise AdversarialRadiusError("binary_tolerance must be in (0, 0.05]")

    def _reference_value(self, axis: dict[str, Any]) -> float:
        cost, reliability = self.reference[axis["experiment_id"]]
        return cost if axis["field"] == "cost_multiplier" else reliability

    def _axis_cap(self, axis: dict[str, Any]) -> float:
        reference = self._reference_value(axis)
        limit = float(axis["limit"])
        direction = axis["direction"]
        if direction == "increase" and limit < reference - 1e-12:
            raise AdversarialRadiusError(
                f"axis {axis['axis_id']} increase limit is below reference"
            )
        if direction == "decrease" and limit > reference + 1e-12:
            raise AdversarialRadiusError(
                f"axis {axis['axis_id']} decrease limit is above reference"
            )
        return abs(limit - reference) / float(axis["scale"])

    def _settings_at_radius(
        self,
        radius: float,
        active_axis_ids: set[str],
    ) -> dict[str, tuple[float, float]]:
        settings = copy.deepcopy(self.reference)
        for axis in self.axes:
            if axis["axis_id"] not in active_axis_ids:
                continue
            experiment_id = axis["experiment_id"]
            cost, reliability = settings[experiment_id]
            reference = self._reference_value(axis)
            step = min(
                radius * float(axis["scale"]),
                abs(float(axis["limit"]) - reference),
            )
            value = reference + step if axis["direction"] == "increase" else reference - step
            if axis["field"] == "cost_multiplier":
                cost = value
            else:
                reliability = value
            settings[experiment_id] = (float(cost), float(reliability))
        return settings

    def _cell_at_radius(
        self,
        radius: float,
        active_axis_ids: set[str],
    ) -> JointCell:
        return self.joint._cell(self._settings_at_radius(radius, active_axis_ids))

    def _first_switch(
        self,
        active_axis_ids: set[str],
    ) -> tuple[bool, float | None, JointCell | None]:
        active = [axis for axis in self.axes if axis["axis_id"] in active_axis_ids]
        if not active:
            raise AdversarialRadiusError("at least one active axis is required")
        max_radius = max(self._axis_cap(axis) for axis in active)
        at_limit = self._cell_at_radius(max_radius, active_axis_ids)
        if at_limit.winner_action_id == self.target_action_id:
            return False, None, None

        low = 0.0
        high = max_radius
        tolerance = float(self.config["binary_tolerance"])
        while high - low > tolerance:
            mid = (low + high) / 2.0
            cell = self._cell_at_radius(mid, active_axis_ids)
            if cell.winner_action_id == self.target_action_id:
                low = mid
            else:
                high = mid
        boundary = self._cell_at_radius(high, active_axis_ids)
        return True, high, boundary

    def _axis_threshold(self, axis: dict[str, Any]) -> AxisThreshold:
        found, radius, boundary = self._first_switch({axis["axis_id"]})
        if not found or radius is None or boundary is None:
            return AxisThreshold(
                axis_id=axis["axis_id"],
                found=False,
                normalized_distance=None,
                boundary_value=None,
                winner_after_boundary=None,
                regret_margin_after_boundary=None,
            )
        settings = dict(
            (experiment_id, (cost, reliability))
            for experiment_id, cost, reliability in boundary.settings
        )
        value_pair = settings[axis["experiment_id"]]
        boundary_value = value_pair[0] if axis["field"] == "cost_multiplier" else value_pair[1]
        return AxisThreshold(
            axis_id=axis["axis_id"],
            found=True,
            normalized_distance=radius,
            boundary_value=boundary_value,
            winner_after_boundary=boundary.winner_action_id,
            regret_margin_after_boundary=boundary.regret_margin,
        )

    def analyze(self) -> AdversarialRadiusResult:
        baseline = self.joint._cell(self.reference)
        if baseline.winner_action_id != self.target_action_id:
            raise AdversarialRadiusError(
                "target action must win at the declared reference settings"
            )

        active = {axis["axis_id"] for axis in self.axes}
        found, radius, boundary = self._first_switch(active)
        coordinated_settings = boundary.settings if boundary is not None else None
        return AdversarialRadiusResult(
            analysis_id=self.config["analysis_id"],
            target_action_id=self.target_action_id,
            baseline_winner_action_id=baseline.winner_action_id,
            baseline_regret_margin=baseline.regret_margin,
            coordinated_switch_found=found,
            coordinated_radius=radius,
            coordinated_winner_after_boundary=(
                boundary.winner_action_id if boundary is not None else None
            ),
            coordinated_regret_margin_after_boundary=(
                boundary.regret_margin if boundary is not None else None
            ),
            coordinated_settings=coordinated_settings,
            axis_thresholds=tuple(self._axis_threshold(axis) for axis in self.axes),
        )
