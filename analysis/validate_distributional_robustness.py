#!/usr/bin/env python3
from __future__ import annotations

import json
import sys
from pathlib import Path

import jsonschema

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from causal_dna.distributional_robust import DistributionallyRobustAnalyzer  # noqa: E402

PLAN_PATH = ROOT / "cases" / "CDNA-001.experiment-plan.json"
CONFIG_PATH = ROOT / "cases" / "CDNA-001.distributional-robustness.json"
SENSITIVITY_PATH = ROOT / "cases" / "CDNA-001.sensitivity-analysis.json"
SCHEMA_PATH = ROOT / "schemas" / "distributional-robustness.schema.json"


def load(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> int:
    plan = load(PLAN_PATH)
    config = load(CONFIG_PATH)
    sensitivity = load(SENSITIVITY_PATH)
    schema = load(SCHEMA_PATH)
    jsonschema.Draft202012Validator(schema).validate(config)

    result = DistributionallyRobustAnalyzer(plan, config, sensitivity).analyze()
    minimax = result.action(result.minimax_interval_regret_action_id)
    maximin = result.action(result.maximin_utility_action_id)

    print("DISTRIBUTIONAL ROBUSTNESS VALIDATION: PASS")
    print(f"- prior points: {result.prior_points}")
    print(f"- uncertainty model: {result.uncertainty_model}")
    print(
        "- minimax interval regret: "
        f"{minimax.action_id} "
        f"(max_regret={minimax.max_interval_regret:.4f} bits/cost-unit, "
        f"mean_regret={minimax.mean_interval_regret:.4f}, "
        f"worst_utility={minimax.worst_case_utility:.4f})"
    )
    print(
        "- maximin utility: "
        f"{maximin.action_id} "
        f"(worst_utility={maximin.worst_case_utility:.4f}, "
        f"max_regret={maximin.max_interval_regret:.4f})"
    )
    for action in result.actions:
        print(
            f"- action {action.action_id}: "
            f"worst_utility={action.worst_case_utility:.4f}, "
            f"best_utility={action.best_case_utility:.4f}, "
            f"mean_lower_utility={action.mean_lower_utility:.4f}, "
            f"max_regret={action.max_interval_regret:.4f}, "
            f"mean_regret={action.mean_interval_regret:.4f}, "
            f"robust_optimal_points={action.robust_optimal_count}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
