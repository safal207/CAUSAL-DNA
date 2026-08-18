#!/usr/bin/env python3
from __future__ import annotations

import json
import sys
from pathlib import Path

from jsonschema import Draft202012Validator

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from causal_dna.adversarial_radius import (  # noqa: E402
    AdversarialRadiusAnalyzer,
    AdversarialRadiusError,
)


def load(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> int:
    plan = load(ROOT / "cases" / "CDNA-001.experiment-plan.json")
    joint = load(ROOT / "cases" / "CDNA-001.joint-assay-frontier.json")
    config = load(ROOT / "cases" / "CDNA-001.adversarial-radius.json")
    schema = load(ROOT / "schemas" / "adversarial-radius.schema.json")

    errors = sorted(
        Draft202012Validator(schema).iter_errors(config),
        key=lambda error: list(error.path),
    )
    if errors:
        print("ADVERSARIAL RADIUS VALIDATION: FAIL")
        for error in errors:
            location = ".".join(str(x) for x in error.path) or "<root>"
            print(f"- {location}: {error.message}")
        return 1

    try:
        result = AdversarialRadiusAnalyzer(plan, joint, config).analyze()
    except AdversarialRadiusError as exc:
        print("ADVERSARIAL RADIUS VALIDATION: FAIL")
        print(f"- {exc}")
        return 1

    print("ADVERSARIAL RADIUS VALIDATION: PASS")
    print(f"- target: {result.target_action_id}")
    print(
        f"- baseline: winner={result.baseline_winner_action_id}, "
        f"regret_margin={result.baseline_regret_margin:.6f}"
    )
    print(f"- normalized fractional rays tested: {result.fractional_rays_tested}")

    if result.minimum_ray_switch_found:
        fractions = ", ".join(
            f"{axis_id}={fraction:.2f}"
            for axis_id, fraction in result.minimum_ray_fractions or ()
        )
        print(
            f"- minimum ray radius: {result.minimum_ray_radius:.6f} "
            f"-> {result.minimum_ray_winner_after_boundary}, "
            f"target_regret_gap={result.minimum_ray_target_regret_gap_after_boundary:.6f}"
        )
        print(f"- minimum ray fractions: {fractions}")
        print("- minimum ray boundary settings:")
        for experiment_id, cost, reliability in result.minimum_ray_settings or ():
            print(f"  {experiment_id}: cost={cost:.6f}x, reliability={reliability:.6f}")
    else:
        print("- minimum ray radius: no switch within declared ray limits")

    if result.coordinated_switch_found:
        print(
            f"- all-axis coordinated radius: {result.coordinated_radius:.6f} "
            f"-> {result.coordinated_winner_after_boundary}, "
            f"target_regret_gap={result.coordinated_target_regret_gap_after_boundary:.6f}"
        )
        print("- all-axis coordinated boundary settings:")
        for experiment_id, cost, reliability in result.coordinated_settings or ():
            print(f"  {experiment_id}: cost={cost:.6f}x, reliability={reliability:.6f}")
    else:
        print("- all-axis coordinated radius: not found within declared limits")

    print("- single-axis thresholds:")
    for threshold in result.axis_thresholds:
        if threshold.found:
            print(
                f"  {threshold.axis_id}: d={threshold.normalized_distance:.6f}, "
                f"boundary_value={threshold.boundary_value:.6f}, "
                f"winner={threshold.winner_after_boundary}, "
                f"target_regret_gap={threshold.target_regret_gap_after_boundary:.6f}"
            )
        else:
            print(f"  {threshold.axis_id}: no switch within declared limit")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
