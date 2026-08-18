#!/usr/bin/env python3
"""Adversarial robustness radius for CAUSAL-DNA planning decisions.

The radius measures how far declared planning assumptions must move from a
reference profile before the current minimax-regret first action changes.
It is a property of the decision model, not biological evidence or causal
confidence.

Because minimax interval regret need not worsen monotonically when every
"adverse" coordinate is moved together, the analyzer searches a declared
family of normalized rays. Each ray is scanned for its *first* observed switch
before a local binary refinement, so a later return to the baseline winner does
not erase an earlier decision boundary.
"""
from __future__ import annotations

import copy
from dataclasses import dataclass
from itertools import product
from typing import Any

from .joint_assay_frontier import JointAssayFrontierAnalyzer, JointCell


class AdversarialRadiusError(ValueError):
    pass


@dataclass(frozen=True)
class AxisThreshold:
    axis_id: str
    found: bool
    normalized_distance: float | None
    boundary_value: float | None
    winner_after_boundary: str | None
    target_regret_gap_after_boundary: float | None


@dataclass(frozen=True)
class AdversarialRadiusResult:
    analysis_id: str
    target_action_id: str
    baseline_winner_action_id: str
    baseline_regret_margin: float
    scan_step: float
    coordinated_switch_found: bool
    coordinated_radius: float | None
    coordinated_winner_after_boundary: str | None
    coordinated_target_regret_gap_after_boundary: float | None
    coordinated_settings: tuple[tuple[str, float, float], ...] | None
    fractional_rays_tested: int
    minimum_ray_switch_found: bool
    minimum_ray_radius: float | None
    minimum_ray_fractions: tuple[tuple[str, float], ...] | None
    minimum_ray_winner_after_boundary: str | None
    minimum_ray_target_regret_gap_after_boundary: float | None
    minimum_ray_settings: tuple[tuple[str, float, float], ...] | None
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
        scan_step = self.config.get("coarse_scan_step")
        if not isinstance(scan_step, (int, float)) or not (0 < scan_step <= 0.25):
            raise AdversarialRadiusError("coarse_scan_step must be in (0, 0.25]")
        if scan_step < tolerance:
            raise AdversarialRadiusError(
                "coarse_scan_step must be >= binary_tolerance"
            )

        fractions = self.config.get("fraction_grid")
        if (
            not isinstance(fractions, list)
            or len(fractions) < 2
            or len(fractions) > 5
            or any(not isinstance(value, (int, float)) for value in fractions)
        ):
            raise AdversarialRadiusError(
                "fraction_grid must contain 2..5 numeric values"
            )
        numeric = [float(value) for value in fractions]
        if numeric != sorted(set(numeric)):
            raise AdversarialRadiusError(
                "fraction_grid must contain unique ascending values"
            )
        if numeric[0] != 0.0 or numeric[-1] != 1.0:
            raise AdversarialRadiusError("fraction_grid must span 0.0 to 1.0")
        if any(value < 0 or value > 1 for value in numeric):
            raise AdversarialRadiusError("fraction_grid values must be in [0, 1]")

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

    def _settings_at_profile(
        self,
        radius: float,
        fractions: dict[str, float],
    ) -> dict[str, tuple[float, float]]:
        settings = copy.deepcopy(self.reference)
        for axis in self.axes:
            fraction = float(fractions.get(axis["axis_id"], 0.0))
            if fraction <= 0:
                continue
            experiment_id = axis["experiment_id"]
            cost, reliability = settings[experiment_id]
            reference = self._reference_value(axis)
            step = min(
                radius * fraction * float(axis["scale"]),
                abs(float(axis["limit"]) - reference),
            )
            value = (
                reference + step
                if axis["direction"] == "increase"
                else reference - step
            )
            if axis["field"] == "cost_multiplier":
                cost = value
            else:
                reliability = value
            settings[experiment_id] = (float(cost), float(reliability))
        return settings

    def _cell_at_profile(
        self,
        radius: float,
        fractions: dict[str, float],
    ) -> JointCell:
        return self.joint._cell(self._settings_at_profile(radius, fractions))

    def _profile_max_radius(self, fractions: dict[str, float]) -> float:
        caps: list[float] = []
        for axis in self.axes:
            fraction = float(fractions.get(axis["axis_id"], 0.0))
            if fraction > 0:
                caps.append(self._axis_cap(axis) / fraction)
        if not caps:
            raise AdversarialRadiusError("profile must move at least one axis")
        return max(caps)

    def _first_switch_profile(
        self,
        fractions: dict[str, float],
        *,
        search_ceiling: float | None = None,
    ) -> tuple[bool, float | None, JointCell | None]:
        max_radius = self._profile_max_radius(fractions)
        if search_ceiling is not None:
            max_radius = min(max_radius, float(search_ceiling))
        if max_radius <= 0:
            return False, None, None

        step = float(self.config["coarse_scan_step"])
        previous_radius = 0.0
        current_radius = min(step, max_radius)
        first_non_target: JointCell | None = None
        while True:
            current = self._cell_at_profile(current_radius, fractions)
            if current.winner_action_id != self.target_action_id:
                first_non_target = current
                break
            if current_radius >= max_radius - 1e-12:
                return False, None, None
            previous_radius = current_radius
            current_radius = min(current_radius + step, max_radius)

        low = previous_radius
        high = current_radius
        tolerance = float(self.config["binary_tolerance"])
        while high - low > tolerance:
            mid = (low + high) / 2.0
            cell = self._cell_at_profile(mid, fractions)
            if cell.winner_action_id == self.target_action_id:
                low = mid
            else:
                high = mid
        boundary = self._cell_at_profile(high, fractions)
        if boundary.winner_action_id == self.target_action_id:
            boundary = first_non_target
            high = current_radius
        return True, high, boundary

    def _target_regret_gap(
        self,
        settings: dict[str, tuple[float, float]],
        winner_action_id: str,
    ) -> float:
        summaries = {
            item.action_id: item for item in self.joint._cell_summaries(settings)
        }
        return max(
            0.0,
            summaries[self.target_action_id].max_regret
            - summaries[winner_action_id].max_regret,
        )

    def _axis_threshold(self, axis: dict[str, Any]) -> AxisThreshold:
        fractions = {current["axis_id"]: 0.0 for current in self.axes}
        fractions[axis["axis_id"]] = 1.0
        found, radius, boundary = self._first_switch_profile(fractions)
        if not found or radius is None or boundary is None:
            return AxisThreshold(
                axis_id=axis["axis_id"],
                found=False,
                normalized_distance=None,
                boundary_value=None,
                winner_after_boundary=None,
                target_regret_gap_after_boundary=None,
            )
        settings = {
            experiment_id: (cost, reliability)
            for experiment_id, cost, reliability in boundary.settings
        }
        pair = settings[axis["experiment_id"]]
        boundary_value = pair[0] if axis["field"] == "cost_multiplier" else pair[1]
        return AxisThreshold(
            axis_id=axis["axis_id"],
            found=True,
            normalized_distance=radius,
            boundary_value=boundary_value,
            winner_after_boundary=boundary.winner_action_id,
            target_regret_gap_after_boundary=self._target_regret_gap(
                settings, boundary.winner_action_id
            ),
        )

    def _fractional_profiles(self) -> tuple[dict[str, float], ...]:
        values = [float(value) for value in self.config["fraction_grid"]]
        axis_ids = [axis["axis_id"] for axis in self.axes]
        profiles: list[dict[str, float]] = []
        for combination in product(values, repeat=len(axis_ids)):
            if max(combination) < 1.0 - 1e-12:
                continue
            profiles.append(dict(zip(axis_ids, combination)))
        return tuple(
            sorted(
                profiles,
                key=lambda profile: (
                    -sum(profile.values()),
                    tuple(-profile[axis_id] for axis_id in axis_ids),
                ),
            )
        )

    def analyze(self) -> AdversarialRadiusResult:
        baseline = self.joint._cell(self.reference)
        if baseline.winner_action_id != self.target_action_id:
            raise AdversarialRadiusError(
                "target action must win at the declared reference settings"
            )

        all_ones = {axis["axis_id"]: 1.0 for axis in self.axes}
        coordinated_found, coordinated_radius, coordinated_boundary = (
            self._first_switch_profile(all_ones)
        )
        coordinated_settings_dict = (
            {
                experiment_id: (cost, reliability)
                for experiment_id, cost, reliability in coordinated_boundary.settings
            }
            if coordinated_boundary is not None
            else None
        )

        profiles = self._fractional_profiles()
        axis_order = tuple(axis["axis_id"] for axis in self.axes)
        best: tuple[float, tuple[float, ...], dict[str, float], JointCell] | None = None
        for fractions in profiles:
            ceiling = best[0] if best is not None else None
            found, radius, boundary = self._first_switch_profile(
                fractions, search_ceiling=ceiling
            )
            if not found or radius is None or boundary is None:
                continue
            key = tuple(fractions[axis_id] for axis_id in axis_order)
            candidate = (radius, key, fractions, boundary)
            if best is None or (candidate[0], candidate[1]) < (best[0], best[1]):
                best = candidate

        if best is not None:
            minimum_radius, _key, minimum_fractions, minimum_boundary = best
            minimum_settings_dict = {
                experiment_id: (cost, reliability)
                for experiment_id, cost, reliability in minimum_boundary.settings
            }
            minimum_fraction_tuple = tuple(
                (axis_id, minimum_fractions[axis_id]) for axis_id in axis_order
            )
            minimum_target_gap = self._target_regret_gap(
                minimum_settings_dict, minimum_boundary.winner_action_id
            )
        else:
            minimum_radius = None
            minimum_boundary = None
            minimum_fraction_tuple = None
            minimum_target_gap = None

        return AdversarialRadiusResult(
            analysis_id=self.config["analysis_id"],
            target_action_id=self.target_action_id,
            baseline_winner_action_id=baseline.winner_action_id,
            baseline_regret_margin=baseline.regret_margin,
            scan_step=float(self.config["coarse_scan_step"]),
            coordinated_switch_found=coordinated_found,
            coordinated_radius=coordinated_radius,
            coordinated_winner_after_boundary=(
                coordinated_boundary.winner_action_id
                if coordinated_boundary is not None
                else None
            ),
            coordinated_target_regret_gap_after_boundary=(
                self._target_regret_gap(
                    coordinated_settings_dict,
                    coordinated_boundary.winner_action_id,
                )
                if coordinated_boundary is not None
                and coordinated_settings_dict is not None
                else None
            ),
            coordinated_settings=(
                coordinated_boundary.settings if coordinated_boundary is not None else None
            ),
            fractional_rays_tested=len(profiles),
            minimum_ray_switch_found=best is not None,
            minimum_ray_radius=minimum_radius,
            minimum_ray_fractions=minimum_fraction_tuple,
            minimum_ray_winner_after_boundary=(
                minimum_boundary.winner_action_id
                if minimum_boundary is not None
                else None
            ),
            minimum_ray_target_regret_gap_after_boundary=minimum_target_gap,
            minimum_ray_settings=(
                minimum_boundary.settings if minimum_boundary is not None else None
            ),
            axis_thresholds=tuple(self._axis_threshold(axis) for axis in self.axes),
        )
