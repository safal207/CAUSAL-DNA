#!/usr/bin/env python3
from __future__ import annotations

import json
import sys
from collections import Counter
from pathlib import Path

from jsonschema import Draft202012Validator

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from causal_dna.phase_map import PhaseBoundaryAnalyzer, PhaseMapError  # noqa: E402


def load(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def compress(points, field):
    ordered = sorted(points, key=lambda point: point.target_cost_multiplier)
    groups = []
    start = ordered[0].target_cost_multiplier
    end = start
    value = getattr(ordered[0], field)
    for point in ordered[1:]:
        current = getattr(point, field)
        if current == value:
            end = point.target_cost_multiplier
            continue
        groups.append((start, end, value))
        start = end = point.target_cost_multiplier
        value = current
    groups.append((start, end, value))
    return "; ".join(
        f"{start:.2f}-{end:.2f}:{value}" if start != end else f"{start:.2f}:{value}"
        for start, end, value in groups
    )


def main() -> int:
    plan = load(ROOT / "cases" / "CDNA-001.experiment-plan.json")
    config = load(ROOT / "cases" / "CDNA-001.phase-map.json")
    sensitivity = load(ROOT / "cases" / "CDNA-001.sensitivity-analysis.json")
    schema = load(ROOT / "schemas" / "phase-map.schema.json")

    errors = sorted(
        Draft202012Validator(schema).iter_errors(config),
        key=lambda error: list(error.path),
    )
    if errors:
        print("PHASE MAP VALIDATION: FAIL")
        for error in errors:
            location = ".".join(str(x) for x in error.path) or "<root>"
            print(f"- {location}: {error.message}")
        return 1

    try:
        result = PhaseBoundaryAnalyzer(plan, config, sensitivity).analyze()
    except (PhaseMapError, ValueError) as exc:
        print("PHASE MAP VALIDATION: FAIL")
        print(f"- {exc}")
        return 1

    baseline_counts = Counter(point.baseline_winner for point in result.phase_points)
    dominant_counts = Counter(point.dominant_winner for point in result.phase_points)

    print("PHASE MAP VALIDATION: PASS")
    print(f"- prior points: {result.prior_points}")
    print(f"- phase points: {len(result.phase_points)}")
    print(f"- target: {result.target_experiment_id}")
    print("- baseline phase occupancy: " + ", ".join(f"{k}={v}" for k, v in sorted(baseline_counts.items())))
    print("- dominant-prior phase occupancy: " + ", ".join(f"{k}={v}" for k, v in sorted(dominant_counts.items())))

    for transfer in (0.0, 0.15):
        selected = [
            point for point in result.break_even_points
            if abs(point.outcome_probability_max_transfer - transfer) <= 1e-12
        ]
        for point in selected:
            multiplier = (
                "none" if point.break_even_cost_multiplier is None
                else f"{point.break_even_cost_multiplier:.4f}"
            )
            cost = (
                "none" if point.break_even_absolute_cost is None
                else f"{point.break_even_absolute_cost:.4f}"
            )
            print(
                f"- break-even transfer={transfer:.2f}, reliability={point.target_reliability:.2f}: "
                f"competitor={point.competitor_id}, cost_multiplier={multiplier}, absolute_cost={cost}"
            )

    for transfer in (0.0, 0.15):
        for reliability in (0.6, 0.8, 1.0):
            slice_points = [
                point for point in result.phase_points
                if abs(point.outcome_probability_max_transfer - transfer) <= 1e-12
                and abs(point.target_reliability - reliability) <= 1e-12
            ]
            print(
                f"- slice transfer={transfer:.2f}, reliability={reliability:.2f}, baseline: "
                + compress(slice_points, "baseline_winner")
            )
            print(
                f"- slice transfer={transfer:.2f}, reliability={reliability:.2f}, dominant-prior: "
                + compress(slice_points, "dominant_winner")
            )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
