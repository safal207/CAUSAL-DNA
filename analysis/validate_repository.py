#!/usr/bin/env python3
"""Validate CAUSAL-DNA machine-readable evidence and research-gap state.

This gate checks structure, cross-record identifiers, three-space graph
integrity, four-dimensional temporal-lattice integrity, append-only causal
processor history, and research-integrity invariants. It does not judge
whether a biological claim is true.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

from jsonschema import Draft202012Validator, FormatChecker

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from causal_dna.processor import CausalProcessor, CausalProcessorError  # noqa: E402
from causal_dna.space_graph import SpaceGraph  # noqa: E402
from causal_dna.temporal_lattice import TemporalLattice  # noqa: E402


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
    graph_schema = load(ROOT / "schemas" / "space-graph.schema.json")
    lattice_schema = load(ROOT / "schemas" / "temporal-lattice.schema.json")
    event_schema = load(ROOT / "schemas" / "causal-event.schema.json")

    edges = load(ROOT / "cases" / "CDNA-001-rs1421085.edges.json")
    gap = load(ROOT / "cases" / "CDNA-001-GAP-001.status.json")
    graph_doc = load(ROOT / "cases" / "CDNA-001-rs1421085.space-graph.json")
    lattice_doc = load(ROOT / "cases" / "CDNA-001-rs1421085.temporal-lattice.json")
    events = load(ROOT / "cases" / "CDNA-001.events.json")

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
    errors.extend(validate(graph_doc, graph_schema, "SPACE-GRAPH"))
    errors.extend(validate(lattice_doc, lattice_schema, "TEMPORAL-LATTICE"))
    for i, event in enumerate(events):
        errors.extend(validate(event, event_schema, f"EVENTS[{i}]"))

    graph = SpaceGraph(graph_doc)
    errors.extend(f"SPACE-GRAPH:{err}" for err in graph.errors())

    lattice = TemporalLattice(lattice_doc)
    errors.extend(f"TEMPORAL-LATTICE:{err}" for err in lattice.errors())

    try:
        processor = CausalProcessor(events)
    except CausalProcessorError as exc:
        errors.append(f"CAUSAL-PROCESSOR:{exc}")
        processor = None

    if gap.get("missing_edge", {}).get("status") == "OPEN" and gap.get("cause_found") is True:
        errors.append("integrity: cause_found=true while missing_edge.status=OPEN")

    for hyp in gap.get("hypotheses", []):
        status = str(hyp.get("status", "")).upper()
        if "ESTABLISHED" in status and not hyp.get("support"):
            errors.append(f"integrity: {hyp.get('id')} says ESTABLISHED but has no support records")

    projective_hypotheses = {
        ref
        for node in graph.nodes_in("projective")
        for ref in node.get("hypothesis_refs", [])
        if node.get("kind") == "hypothesis"
    }
    bardo_hypotheses = {
        ref
        for node in graph.nodes_in("bardo")
        for ref in node.get("hypothesis_refs", [])
    }
    missing_bardo = sorted(projective_hypotheses - bardo_hypotheses)
    if missing_bardo:
        errors.append(
            "three-space integrity: projective hypotheses missing Bardo representation: "
            + ", ".join(missing_bardo)
        )

    if gap.get("missing_edge", {}).get("status") == "OPEN":
        integration_gate = graph.node("B_INTEGRATION_GATE")
        if integration_gate.get("resolution", "open") != "open":
            errors.append("three-space integrity: GAP-001 is OPEN but B_INTEGRATION_GATE is resolved")

        future_material = lattice.states_at(space="material", time="future")
        if future_material:
            errors.append(
                "4D integrity: GAP-001 is OPEN but future material states already exist: "
                + ", ".join(s["id"] for s in future_material)
            )

        if processor is not None:
            projection = processor.project()
            if "GAP-001" in projection.verified_subjects:
                errors.append("processor integrity: open GAP-001 is already independently verified")

    contour_sources = [
        "P_H1_ARID5B", "P_H2_CUX1", "P_H3_ACCESS", "P_H4_CONTACT",
        "P_H6_COMPOSITION", "P_H7_SEX", "P_H8_GRAMMAR", "P_H9_TIME",
    ]
    for source in contour_sources:
        if not graph.contours(source, "M_IRX3_UP"):
            errors.append(
                f"three-space integrity: no projective→Bardo→material contour from {source} to M_IRX3_UP"
            )

    verifier_states = [
        s for s in lattice.states_at(time="future")
        if s.get("observer") == "independent_verifier"
    ]
    if not verifier_states:
        errors.append("4D integrity: no future independent verification gate is represented")

    if errors:
        print("CAUSAL-DNA VALIDATION: FAIL")
        for err in errors:
            print(f"- {err}")
        return 1

    summary = graph.summary()
    lattice_summary = lattice.summary()
    projection = processor.project() if processor is not None else None

    print("CAUSAL-DNA VALIDATION: PASS")
    print(f"- causal edges: {len(edges)}")
    print(f"- hypotheses: {len(gap.get('hypotheses', []))}")
    print(f"- GAP-001 status: {gap['missing_edge']['status']}")
    print(f"- cause_found: {str(gap['cause_found']).lower()}")
    print(
        "- three-space nodes: "
        f"projective={summary['nodes_by_space']['projective']}, "
        f"bardo={summary['nodes_by_space']['bardo']}, "
        f"material={summary['nodes_by_space']['material']}"
    )
    print(f"- open Bardo states: {summary['open_bardo']}")
    print("- materialization frontier: " + ", ".join(summary["materialization_frontier"]))
    print(f"- 4D lattice states: {lattice_summary['states']}")
    print(f"- 4D lattice transitions: {lattice_summary['transitions']}")
    print("- 4D future frontier: " + ", ".join(lattice_summary["future_frontier"]))
    print("- 4D materialized history: " + ", ".join(lattice_summary["materialized_history"]))
    if projection is not None:
        print(f"- causal processor events: {len(events)}")
        print(f"- causal processor generation: {projection.generation}")
        print("- processor open hypotheses: " + ", ".join(projection.open_hypotheses))
        print("- processor rejected hypotheses: " + ", ".join(projection.rejected_hypotheses))
    return 0


if __name__ == "__main__":
    sys.exit(main())
