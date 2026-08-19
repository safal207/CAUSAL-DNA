#!/usr/bin/env python3
from __future__ import annotations

import json
import sys
from collections import Counter
from pathlib import Path

from jsonschema import Draft202012Validator

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from causal_dna.phase_boundary import PhaseBoundaryError, RobustPhaseBoundaryAnalyzer  # noqa: E402


def load(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> int:
    plan = load(ROOT / "cases" / "CDNA-001.experiment-plan.json")
    config = load(ROOT / "cases" / "CDNA-001.phase-boundary.json")
    schema = load(ROOT / "schemas" / "phase-boundary.schema.json")

    errors = sorted(
        Draft202012Validator(schema).iter_errors(config),
        key=lambda error: list(error.path),
    )
    if errors:
        print("PHASE BOUNDARY VALIDATION: FAIL")
        for error in errors:
            location = ".".join(str(x) for x in error.path) or "<root>"
            print(f"- {location}: {error.message}")
        return 1

    try:
        result = RobustPhaseBoundaryAnalyzer(plan, config).analyze()
    except PhaseBoundaryError as exc:
        print("PHASE BOUNDARY VALIDATION: FAIL")
        print(f"- {exc}")
        return 1

    print("PHASE BOUNDARY VALIDATION: PASS")
    print(f"- cells: {len(result.cells)}")
    print(f"- prior points per cell: {result.prior_points_per_cell}")
    print(f"- prior-cell evaluations: {len(result.cells) * result.prior_points_per_cell}")
    print(f"- target: {result.target_experiment_id}")
    print(
        "- winner counts: "
        + ", ".join(
            f"{action_id}={count}" for action_id, count in result.winner_counts.items()
        )
    )

    transition_pairs = Counter(
        tuple(sorted((boundary.from_action_id, boundary.to_action_id)))
        for boundary in result.boundaries
    )
    print(f"- boundary edges: {len(result.boundaries)}")
    print(
        "- transition pairs: "
        + (
            ", ".join(
                f"{left}<->{right}={count}"
                for (left, right), count in sorted(transition_pairs.items())
            )
            if transition_pairs
            else "none"
        )
    )

    print("- target cost ceilings by transfer/reliability:")
    for transfer in config["outcome_probability_max_transfer_grid"]:
        entries = []
        for reliability in config["target_reliability_grid"]:
            ceiling = result.target_cost_ceiling(
                reliability=float(reliability),
                probability_transfer=float(transfer),
            )
            entries.append(
                f"r={float(reliability):.2f}:"
                + (f"{ceiling:.2f}x" if ceiling is not None else "none")
            )
        print(f"  transfer={float(transfer):.2f} -> " + ", ".join(entries))

    smallest_margins = sorted(
        result.cells,
        key=lambda item: (
            item.regret_margin,
            item.target_reliability,
            item.target_cost_multiplier,
            item.outcome_probability_max_transfer,
        ),
    )[:10]
    print("- nearest decision boundaries:")
    for cell in smallest_margins:
        print(
            "  "
            f"r={cell.target_reliability:.2f}, "
            f"cost={cell.target_cost_multiplier:.2f}x, "
            f"transfer={cell.outcome_probability_max_transfer:.2f}: "
            f"{cell.winner_action_id} vs {cell.runner_up_action_id}, "
            f"regret_margin={cell.regret_margin:.6f}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
