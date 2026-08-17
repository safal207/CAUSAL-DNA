#!/usr/bin/env python3
"""Three-space causal graph engine for CAUSAL-DNA.

The model is epistemic, not metaphysical:

- projective: proposed structures, hypotheses, planned interventions;
- bardo: unresolved transition states where alternatives compete;
- material: observations/interventions backed by evidence artifacts.

The important property is provenance-preserving materialization: a proposed
transition cannot silently become a material fact. It must cross an explicit
Bardo state and carry evidence references when it reaches material space.
"""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Iterator

SPACES = ("projective", "bardo", "material")
BARDO_RESOLUTIONS = ("open", "materialized", "rejected", "superseded")


class SpaceGraphError(ValueError):
    """Raised when a three-space graph violates structural invariants."""


@dataclass(frozen=True)
class PathView:
    node_ids: tuple[str, ...]
    edge_ids: tuple[str, ...]
    spaces: tuple[str, ...]

    @property
    def compressed_spaces(self) -> tuple[str, ...]:
        out: list[str] = []
        for space in self.spaces:
            if not out or out[-1] != space:
                out.append(space)
        return tuple(out)


class SpaceGraph:
    """In-memory representation of a projective→Bardo→material causal graph."""

    def __init__(self, document: dict[str, Any]):
        self.document = document
        self.nodes = document.get("nodes", [])
        self.edges = document.get("edges", [])
        self._nodes = {n.get("id"): n for n in self.nodes}
        self._edges = {e.get("id"): e for e in self.edges}
        self._out: dict[str, list[dict[str, Any]]] = defaultdict(list)
        for edge in self.edges:
            self._out[edge.get("from")].append(edge)

    @classmethod
    def load(cls, path: str | Path) -> "SpaceGraph":
        with Path(path).open(encoding="utf-8") as handle:
            graph = cls(json.load(handle))
        graph.assert_valid()
        return graph

    def errors(self) -> list[str]:
        errors: list[str] = []

        node_ids = [n.get("id") for n in self.nodes]
        if any(not value for value in node_ids):
            errors.append("every node requires a non-empty id")
        if len(node_ids) != len(set(node_ids)):
            errors.append("node ids must be unique")

        edge_ids = [e.get("id") for e in self.edges]
        if any(not value for value in edge_ids):
            errors.append("every edge requires a non-empty id")
        if len(edge_ids) != len(set(edge_ids)):
            errors.append("edge ids must be unique")

        for node in self.nodes:
            node_id = node.get("id", "<unknown>")
            space = node.get("space")
            if space not in SPACES:
                errors.append(f"{node_id}: invalid space {space!r}")

            evidence = node.get("evidence_refs", [])
            if space == "material" and not evidence:
                errors.append(f"{node_id}: material node requires evidence_refs")

            if space == "bardo":
                resolution = node.get("resolution", "open")
                if resolution not in BARDO_RESOLUTIONS:
                    errors.append(f"{node_id}: invalid Bardo resolution {resolution!r}")
                if resolution == "materialized" and not node.get("resolved_by"):
                    errors.append(f"{node_id}: materialized Bardo node requires resolved_by")
                if resolution in {"rejected", "superseded"} and not node.get("resolution_reason"):
                    errors.append(
                        f"{node_id}: {resolution} Bardo node requires resolution_reason"
                    )

        known_nodes = set(node_ids)
        for edge in self.edges:
            edge_id = edge.get("id", "<unknown>")
            source = edge.get("from")
            target = edge.get("to")
            if source not in known_nodes:
                errors.append(f"{edge_id}: unknown source node {source!r}")
            if target not in known_nodes:
                errors.append(f"{edge_id}: unknown target node {target!r}")
            if source not in self._nodes or target not in self._nodes:
                continue

            source_space = self._nodes[source].get("space")
            target_space = self._nodes[target].get("space")
            evidence = edge.get("evidence_refs", [])

            # No hypothesis is allowed to jump straight into material fact.
            if source_space == "projective" and target_space == "material":
                errors.append(
                    f"{edge_id}: projective→material bypass forbidden; cross Bardo explicitly"
                )

            # Crossing into material space is a provenance boundary.
            if target_space == "material" and not evidence:
                errors.append(f"{edge_id}: edge entering material space requires evidence_refs")

            # An edge declared as materialization must actually cross Bardo→material.
            if edge.get("relation") == "materializes_as" and (
                source_space != "bardo" or target_space != "material"
            ):
                errors.append(
                    f"{edge_id}: materializes_as must connect Bardo→material"
                )

        return errors

    def assert_valid(self) -> None:
        errors = self.errors()
        if errors:
            raise SpaceGraphError("\n".join(errors))

    def node(self, node_id: str) -> dict[str, Any]:
        try:
            return self._nodes[node_id]
        except KeyError as exc:
            raise SpaceGraphError(f"unknown node: {node_id}") from exc

    def nodes_in(self, space: str) -> list[dict[str, Any]]:
        if space not in SPACES:
            raise SpaceGraphError(f"unknown space: {space}")
        return [node for node in self.nodes if node.get("space") == space]

    def outgoing(self, node_id: str) -> list[dict[str, Any]]:
        self.node(node_id)
        return list(self._out.get(node_id, ()))

    def paths(
        self,
        start: str,
        end: str,
        *,
        max_depth: int = 32,
    ) -> Iterator[PathView]:
        """Yield simple directed paths without revisiting a node."""
        self.node(start)
        self.node(end)

        stack: list[tuple[str, list[str], list[str]]] = [(start, [start], [])]
        while stack:
            current, node_path, edge_path = stack.pop()
            if len(edge_path) > max_depth:
                continue
            if current == end:
                yield PathView(
                    tuple(node_path),
                    tuple(edge_path),
                    tuple(self._nodes[n]["space"] for n in node_path),
                )
                continue
            for edge in reversed(self._out.get(current, [])):
                nxt = edge["to"]
                if nxt in node_path:
                    continue
                stack.append((nxt, node_path + [nxt], edge_path + [edge["id"]]))

    def contours(self, start: str, end: str) -> list[PathView]:
        """Return paths that explicitly traverse projective→Bardo→material.

        Extra same-space steps are allowed. Material observations may feed back
        into projective space elsewhere in the graph, but a contour itself is
        the forward materialization pass.
        """
        return [
            path
            for path in self.paths(start, end)
            if path.compressed_spaces == ("projective", "bardo", "material")
        ]

    def open_bardo(self) -> list[dict[str, Any]]:
        return [
            node
            for node in self.nodes_in("bardo")
            if node.get("resolution", "open") == "open"
        ]

    def materialization_frontier(self) -> list[dict[str, Any]]:
        """Open Bardo states with at least one planned materialization test."""
        frontier = []
        for node in self.open_bardo():
            relations = {edge.get("relation") for edge in self._out.get(node["id"], [])}
            if relations & {"tested_by", "bounded_by", "materializes_as"}:
                frontier.append(node)
        return frontier

    def summary(self) -> dict[str, Any]:
        by_space = {space: len(self.nodes_in(space)) for space in SPACES}
        return {
            "graph_id": self.document.get("graph_id"),
            "case_id": self.document.get("case_id"),
            "nodes": len(self.nodes),
            "edges": len(self.edges),
            "nodes_by_space": by_space,
            "open_bardo": len(self.open_bardo()),
            "materialization_frontier": [n["id"] for n in self.materialization_frontier()],
        }


def _format_summary(summary: dict[str, Any]) -> str:
    return json.dumps(summary, indent=2, ensure_ascii=False)


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("graph", type=Path)
    parser.add_argument("--start")
    parser.add_argument("--end")
    parser.add_argument("--contours", action="store_true")
    args = parser.parse_args(list(argv) if argv is not None else None)

    graph = SpaceGraph.load(args.graph)
    print(_format_summary(graph.summary()))

    if args.start or args.end:
        if not (args.start and args.end):
            parser.error("--start and --end must be supplied together")
        paths = graph.contours(args.start, args.end) if args.contours else list(
            graph.paths(args.start, args.end)
        )
        for path in paths:
            print(" -> ".join(path.node_ids))
            print("spaces:", " -> ".join(path.compressed_spaces))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
