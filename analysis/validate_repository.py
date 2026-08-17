#!/usr/bin/env python3
"""Validate CAUSAL-DNA machine-readable evidence and research-gap state.

This gate checks structure, cross-record identifiers, three-space graph
integrity, four-dimensional temporal-lattice integrity, append-only causal
processor history, active experiment planning, adaptive replanning, budgeted
multi-step strategy planning, and research-integrity invariants. It does not
judge whether a biological claim is true.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

from jsonschema import Draft202012Validator, FormatChecker

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from causal_dna.adaptive_replanner import AdaptiveReplanner, AdaptiveReplannerError  # noqa: E402
from causal_dna.experiment_selector import ExperimentSelector, ExperimentSelectorError  # noqa: E402
from causal_dna.processor import CausalProcessor, CausalProcessorError  # noqa: E402
from causal_dna.space_graph import SpaceGraph  # noqa: E402
from causal_dna.strategy_planner import StrategyPlanner, StrategyPlannerError  # noqa: E402
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
    plan_schema = load(ROOT / "schemas" / "experiment-plan.schema.json")
    replan_schema = load(ROOT / "schemas" / "adaptive-replan.schema.json")
    strategy_schema = load(ROOT / "schemas" / "strategy-config.schema.json")

    edges = load(ROOT / "cases" / "CDNA-001-rs1421085.edges.json")
    gap = load(ROOT / "cases" / "CDNA-001-GAP-001.status.json")
    graph_doc = load(ROOT / "cases" / "CDNA-001-rs1421085.space-graph.json")
    lattice_doc = load(ROOT / "cases" / "CDNA-001-rs1421085.temporal-lattice.json")
    events = load(ROOT / "cases" / "CDNA-001.events.json")
    experiment_plan = load(ROOT / "cases" / "CDNA-001.experiment-plan.json")
    replanning_scenarios = load(ROOT / "cases" / "CDNA-001.replanning-scenarios.json")
    strategy_configs = load(ROOT / "cases" / "CDNA-001.strategy-configs.json")

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
    errors.extend(validate(experiment_plan, plan_schema, "EXPERIMENT-PLAN"))
    for i, event in enumerate(events):
        errors.extend(validate(event, event_schema, f"EVENTS[{i}]"))
    for i, scenario in enumerate(replanning_scenarios):
        errors.extend(validate(scenario, replan_schema, f"REPLAN[{i}]"))
    for i, config in enumerate(strategy_configs):
        errors.extend(validate(config, strategy_schema, f"STRATEGY[{i}]"))

    graph = SpaceGraph(graph_doc)
    errors.extend(f"SPACE-GRAPH:{err}" for err in graph.errors())

    lattice = TemporalLattice(lattice_doc)
    errors.extend(f"TEMPORAL-LATTICE:{err}" for err in lattice.errors())

    try:
        processor = CausalProcessor(events)
    except CausalProcessorError as exc:
        errors.append(f"CAUSAL-PROCESSOR:{exc}")
        processor = None

    try:
        selector = ExperimentSelector(experiment_plan)
    except ExperimentSelectorError as exc:
        errors.append(f"EXPERIMENT-SELECTOR:{exc}")
        selector = None

    replans: dict[str, object] = {}
    for scenario in replanning_scenarios:
        try:
            replans[scenario["scenario_id"]] = AdaptiveReplanner(
                experiment_plan, scenario["outcomes"]
            ).replay()
        except AdaptiveReplannerError as exc:
            errors.append(f"ADAPTIVE-REPLANNER:{scenario.get('scenario_id')}:{exc}")

    strategies: dict[str, object] = {}
    try:
        strategy_planner = StrategyPlanner(experiment_plan)
    except Exception as exc:  # ExperimentSelector errors are surfaced by the planner.
        errors.append(f"STRATEGY-PLANNER:{exc}")
        strategy_planner = None

    if strategy_planner is not None:
        for config in strategy_configs:
            try:
                strategies[config["strategy_id"]] = strategy_planner.best(
                    budget=config["budget"],
                    max_depth=config["max_depth"],
                    objective=config["objective"],
                )
            except StrategyPlannerError as exc:
                errors.append(f"STRATEGY-PLANNER:{config.get('strategy_id')}:{exc}")

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

    # Active-discovery integrity: planner state is advisory only.
    if "biological evidence" not in experiment_plan.get("disclaimer", ""):
        errors.append("experiment-plan integrity: disclaimer must state planning priors are not biological evidence")

    if selector is not None and processor is not None:
        projection = processor.project()
        planned_hypotheses = set(experiment_plan.get("hypotheses", {}))
        open_processor_hypotheses = set(projection.open_hypotheses)
        missing_open = sorted(planned_hypotheses - open_processor_hypotheses)
        if missing_open:
            errors.append(
                "experiment-plan integrity: planner includes hypotheses not open in processor: "
                + ", ".join(missing_open)
            )

    # Adaptive-replanning integrity: simulations may alter planning weights and
    # rankings only. They must not touch the authoritative event stream or claim
    # causal closure.
    authoritative_event_count = len(events)
    for scenario in replanning_scenarios:
        if scenario.get("simulation_only") is not True:
            errors.append(f"adaptive-replan integrity: {scenario.get('scenario_id')} is not explicitly simulated")
        for outcome in scenario.get("outcomes", []):
            forbidden = {
                "cause_found", "causal_status", "edge_status", "materialized", "verification_status"
            }.intersection(outcome)
            if forbidden:
                errors.append(
                    f"adaptive-replan integrity: {scenario.get('scenario_id')} mutates causal state: "
                    + ", ".join(sorted(forbidden))
                )
    if len(events) != authoritative_event_count:
        errors.append("adaptive-replan integrity: authoritative event history was mutated")

    # Multi-step strategy planning is also advisory-only and must obey every
    # declared path budget.
    forbidden_strategy_fields = {
        "cause_found", "causal_status", "edge_status", "materialized", "verification_status"
    }
    for config in strategy_configs:
        if config.get("planning_only") is not True:
            errors.append(f"strategy integrity: {config.get('strategy_id')} is not planning_only")
        forbidden = forbidden_strategy_fields.intersection(config)
        if forbidden:
            errors.append(
                f"strategy integrity: {config.get('strategy_id')} mutates causal state: "
                + ", ".join(sorted(forbidden))
            )
        result = strategies.get(config.get("strategy_id"))
        if result is not None and result.max_path_cost > config["budget"] + 1e-9:
            errors.append(
                f"strategy integrity: {config.get('strategy_id')} exceeds path budget"
            )

    if errors:
        print("CAUSAL-DNA VALIDATION: FAIL")
        for err in errors:
            print(f"- {err}")
        return 1

    summary = graph.summary()
    lattice_summary = lattice.summary()
    projection = processor.project() if processor is not None else None
    ranked = selector.rank() if selector is not None else []

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
    if ranked:
        print(f"- experiment candidates: {len(ranked)}")
        print(
            "- max expected information gain: "
            f"{ranked[0].experiment_id} ({ranked[0].expected_information_gain_bits:.4f} bits)"
        )
        efficient = selector.best(objective="gain_per_cost")
        print(
            "- max information gain per cost: "
            f"{efficient.experiment_id} ({efficient.utility_per_cost:.4f} bits/cost-unit)"
        )
    for scenario_id, snapshot in sorted(replans.items()):
        best_ig = snapshot.best_information_gain
        best_cost = snapshot.best_gain_per_cost
        print(
            f"- adaptive replan {scenario_id}: generation={snapshot.generation}, "
            f"entropy={snapshot.entropy_bits:.4f} bits, "
            f"next_IG={best_ig.experiment_id if best_ig else 'NONE'}, "
            f"next_cost={best_cost.experiment_id if best_cost else 'NONE'}"
        )
    for strategy_id, result in sorted(strategies.items()):
        branches = ", ".join(
            f"{outcome}->{next_id or 'STOP'}"
            for outcome, next_id in sorted(result.next_experiments_by_outcome.items())
        )
        print(
            f"- strategy {strategy_id}: first={result.first_experiment_id}, "
            f"gain={result.expected_information_gain_bits:.4f} bits, "
            f"expected_cost={result.expected_cost:.4f}, "
            f"max_path_cost={result.max_path_cost:.4f}, "
            f"gain/cost={result.gain_per_cost:.4f}, branches=[{branches}]"
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
