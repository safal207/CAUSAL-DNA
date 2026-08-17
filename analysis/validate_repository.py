#!/usr/bin/env python3
"""Validate CAUSAL-DNA machine-readable evidence and research-gap state.

This gate checks structure, cross-record identifiers, and a few research-
integrity invariants. It does not judge whether a biological claim is true.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

from jsonschema import Draft202012Validator, FormatChecker

ROOT = Path(__file__).resolve().parents[1]


def load(path: Path):
    with path.open(encoding="utf-8") as f:
        return json.load(f)


def validate(instance, schema, label: str) -> list[str]:
    validator = Draft202012Validator(schema, format_checker=FormatChecker())
    errors = sorted(validator.iter_errors(instance), key=lambda e: list(e.path))
    out = []
    for err in errors:
        loc = ".".join(str(x) for x in err.absolute_path) or "<root>"
        out.append(f"{label}:{loc}: {err.message}")
    return out


def main() -> int:
    errors: list[str] = []

    edge_schema = load(ROOT / "schemas" / "causal-edge.schema.json")
    gap_schema = load(ROOT / "schemas" / "gap-status.schema.json")
    edges = load(ROOT / "cases" / "CDNA-001-rs1421085.edges.json")
    gap = load(ROOT / "cases" / "CDNA-001-GAP-001.status.json")

    if not isinstance(edges, list) or not edges:
        errors.append("edge ledger must be a non-empty JSON array")
    else:
        seen = set()
        for i, edge in enumerate(edges):
            errors.extend(validate(edge, edge_schema, f"edges[{i}]"))
            edge_id = edge.get("edge_id")
            if edge_id in seen:
                errors.append(f"duplicate edge_id: {edge_id}")
            seen.add(edge_id)

    errors.extend(validate(gap, gap_schema, "GAP-001"))

    # Integrity invariant: an open missing edge cannot coexist with cause_found=true.
    if gap.get("missing_edge", {}).get("status") == "OPEN" and gap.get("cause_found") is True:
        errors.append("integrity: cause_found=true while missing_edge.status=OPEN")

    # Integrity invariant: every active hypothesis needs either supporting evidence
    # or an explicit next test; the schema guarantees next_test, this guards empty support
    # from being mistaken for established evidence in status wording.
    for hyp in gap.get("hypotheses", []):
        status = str(hyp.get("status", "")).upper()
        if "ESTABLISHED" in status and not hyp.get("support"):
            errors.append(f"integrity: {hyp.get('id')} says ESTABLISHED but has no support records")

    if errors:
        print("CAUSAL-DNA VALIDATION: FAIL")
        for err in errors:
            print(f"- {err}")
        return 1

    print("CAUSAL-DNA VALIDATION: PASS")
    print(f"- causal edges: {len(edges)}")
    print(f"- hypotheses: {len(gap.get('hypotheses', []))}")
    print(f"- GAP-001 status: {gap['missing_edge']['status']}")
    print(f"- cause_found: {str(gap['cause_found']).lower()}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
