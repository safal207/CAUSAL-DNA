#!/usr/bin/env python3
from __future__ import annotations

import json
import sys
from collections import Counter
from pathlib import Path

from jsonschema import Draft202012Validator

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from causal_dna.joint_assay_frontier import (  # noqa: E402
    JointAssayFrontierAnalyzer,
    JointAssayFrontierError,
)


def load(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> int:
    plan = load(ROOT / "cases" / "CDNA-001.experiment-plan.json")
    config = load(ROOT / "cases" / "CDNA-001.joint-assay-frontier.json")
    schema = load(ROOT / "schemas" / "joint-assay-frontier.schema.json")

    errors = sorted(
        Draft202012Validator(schema).iter_errors(config),
        key=lambda error: list(error.path),
    )
    if errors:
        print("JOINT ASSAY FRONTIER VALIDATION: FAIL")
        for error in errors:
            location = ".".join(str(x) for x in error.path) or "<root>"
            print(f"- {location}: {error.message}")
        return 1

    try:
        result = JointAssayFrontierAnalyzer(plan, config).analyze()
    except JointAssayFrontierError as exc:
        print("JOINT ASSAY FRONTIER VALIDATION: FAIL")
        print(f"- {exc}")
        return 1

    print("JOINT ASSAY FRONTIER VALIDATION: PASS")
    print(f"- cells: {len(result.cells)}")
    print(f"- prior points per cell: {result.prior_points_per_cell}")
    print(
        f"- prior-cell evaluations: "
        f"{len(result.cells) * result.prior_points_per_cell}"
    )
    print(
        "- winner counts: "
        + ", ".join(f"{key}={value}" for key, value in result.winner_counts.items())
    )
    print(
        "- boundary pairs: "
        + (
            ", ".join(
                f"{key}={value}" for key, value in result.boundary_pair_counts.items()
            )
            or "none"
        )
    )
    print(
        "- boundary axes: "
        + (
            ", ".join(
                f"{key}={value}" for key, value in result.boundary_axis_counts.items()
            )
            or "none"
        )
    )

    smallest = sorted(result.cells, key=lambda cell: (cell.regret_margin, cell.winner_action_id))[:8]
    print("- nearest joint decision boundaries:")
    for cell in smallest:
        settings = ", ".join(
            f"{experiment_id}:cost={cost:.2f}x,r={reliability:.2f}"
            for experiment_id, cost, reliability in cell.settings
        )
        print(
            f"  {settings} -> {cell.winner_action_id} vs {cell.runner_up_action_id}, "
            f"margin={cell.regret_margin:.6f}"
        )

    print(f"- Pareto policies: {len(result.pareto_points)}")
    horizon_counts = Counter(point.horizon for point in result.pareto_points)
    first_counts = Counter(point.first_experiment_id for point in result.pareto_points)
    print(
        "- Pareto horizons: "
        + ", ".join(f"h{key}={value}" for key, value in sorted(horizon_counts.items()))
    )
    print(
        "- Pareto first actions: "
        + ", ".join(f"{key}={value}" for key, value in sorted(first_counts.items()))
    )
    for point in result.pareto_points[:12]:
        print(
            f"- Pareto h={point.horizon}, first={point.first_experiment_id}, "
            f"gain={point.expected_information_gain_bits:.4f} bits, "
            f"expected_cost={point.expected_cost:.4f}, "
            f"max_path_cost={point.max_path_cost:.4f}, "
            f"first_step_max_regret={point.first_step_max_regret:.4f}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
