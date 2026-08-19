#!/usr/bin/env python3
"""Validate and report CAUSAL-DNA robust first-step planning analysis."""
from __future__ import annotations

import json
import sys
from pathlib import Path

from jsonschema import Draft202012Validator

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from causal_dna.robust_decision import RobustDecisionError, RobustFirstStepAnalyzer  # noqa: E402


def load(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> int:
    plan = load(ROOT / "cases" / "CDNA-001.experiment-plan.json")
    config = load(ROOT / "cases" / "CDNA-001.robust-decision.json")
    sensitivity = load(ROOT / "cases" / "CDNA-001.sensitivity-analysis.json")
    schema = load(ROOT / "schemas" / "robust-decision.schema.json")

    errors = sorted(
        Draft202012Validator(schema).iter_errors(config),
        key=lambda error: list(error.path),
    )
    if errors:
        print("ROBUST DECISION VALIDATION: FAIL")
        for error in errors:
            location = ".".join(str(x) for x in error.path) or "<root>"
            print(f"- {location}: {error.message}")
        return 1

    try:
        result = RobustFirstStepAnalyzer(plan, config, sensitivity).analyze()
    except RobustDecisionError as exc:
        print("ROBUST DECISION VALIDATION: FAIL")
        print(f"- {exc}")
        return 1

    minimax = result.action(result.minimax_regret_action_id)
    maximin = result.action(result.maximin_gain_action_id)

    print("ROBUST DECISION VALIDATION: PASS")
    print(f"- grid points: {result.grid_points}")
    print(
        f"- minimax regret: {minimax.action_id} "
        f"(max_regret={minimax.max_regret_bits:.4f} bits, "
        f"mean_regret={minimax.mean_regret_bits:.4f} bits, "
        f"worst_gain={minimax.worst_case_gain_bits:.4f} bits)"
    )
    print(
        f"- maximin gain: {maximin.action_id} "
        f"(worst_gain={maximin.worst_case_gain_bits:.4f} bits, "
        f"max_regret={maximin.max_regret_bits:.4f} bits)"
    )
    for item in result.actions:
        print(
            f"- action {item.action_id}: worst_gain={item.worst_case_gain_bits:.4f}, "
            f"mean_gain={item.mean_gain_bits:.4f}, max_regret={item.max_regret_bits:.4f}, "
            f"mean_regret={item.mean_regret_bits:.4f}, optimal_points={item.optimal_count}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
