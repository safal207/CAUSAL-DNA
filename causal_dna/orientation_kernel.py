#!/usr/bin/env python3
"""Causal Orientation Kernel (COK) detector for CAUSAL-DNA.

COK asks a narrower question than graph centrality: which smallest experimentally
supported set of intermediate states intercepts every currently validated
source -> target causal path?

Important boundary: observational/hypothesis edges can help plan experiments,
but they cannot establish an orientation kernel. A kernel is only promoted from
intervention-supported causal paths, and independent verification is tracked
separately.
"""
from __future__ import annotations

import argparse
import itertools
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

EDGE_RELATIONS = ("CAUSES", "ENABLES", "MODULATES", "FEEDBACK")
EDGE_STATUSES = (
    "HYPOTHESIS",
    "OBSERVATIONAL",
    "INTERVENTION_SUPPORTED",
    "RESCUE_SUPPORTED",
    "VERIFIED",
)
CAUSAL_STATUSES = {"INTERVENTION_SUPPORTED", "RESCUE_SUPPORTED", "VERIFIED"}


class OrientationKernelError(ValueError):
    """Raised when an exchange-process graph violates COK invariants."""


@dataclass(frozen=True)
class KernelCertificate:
    source: str
    target: str
    kernel: tuple[str, ...]
    causal_paths: int
    independently_verified: bool
    status: str


class OrientationKernelGraph:
    def __init__(self, document: dict[str, Any]):
        self.document = document
        self.nodes = {node["id"]: node for node in document.get("nodes", []) if node.get("id")}
        self.edges = document.get("edges", [])

    @classmethod
    def load(cls, path: str | Path) -> "OrientationKernelGraph":
        with Path(path).open(encoding="utf-8") as handle:
            graph = cls(json.load(handle))
        graph.assert_valid()
        return graph

    def errors(self) -> list[str]:
        errors: list[str] = []
        nodes = self.document.get("nodes", [])
        ids = [n.get("id") for n in nodes]
        if not nodes:
            errors.append("at least one node is required")
        if any(not node_id for node_id in ids):
            errors.append("every node requires a non-empty id")
        if len(ids) != len(set(ids)):
            errors.append("node ids must be unique")

        source = self.document.get("source")
        target = self.document.get("target")
        if source not in self.nodes:
            errors.append("source must reference a declared node")
        if target not in self.nodes:
            errors.append("target must reference a declared node")
        if source == target and source is not None:
            errors.append("source and target must differ")

        for index, edge in enumerate(self.edges):
            prefix = f"edge[{index}]"
            if edge.get("from") not in self.nodes:
                errors.append(f"{prefix}: unknown from node")
            if edge.get("to") not in self.nodes:
                errors.append(f"{prefix}: unknown to node")
            if edge.get("relation") not in EDGE_RELATIONS:
                errors.append(f"{prefix}: invalid relation {edge.get('relation')!r}")
            if edge.get("status") not in EDGE_STATUSES:
                errors.append(f"{prefix}: invalid status {edge.get('status')!r}")
            refs = edge.get("evidence_refs")
            if not isinstance(refs, list):
                errors.append(f"{prefix}: evidence_refs must be a list")
            elif edge.get("status") != "HYPOTHESIS" and not refs:
                errors.append(f"{prefix}: non-hypothesis edge requires evidence_refs")
            elif any(not isinstance(ref, str) or not ref.strip() for ref in refs):
                errors.append(f"{prefix}: evidence_refs must contain non-empty strings")

        declared = self.document.get("declared_kernel", [])
        if not isinstance(declared, list):
            errors.append("declared_kernel must be a list")
        elif any(node_id not in self.nodes for node_id in declared):
            errors.append("declared_kernel references unknown node")
        elif source in declared or target in declared:
            errors.append("declared_kernel cannot contain source or target")

        computed = [list(k) for k in self.minimal_kernels()]
        if declared and declared not in computed:
            errors.append("declared_kernel is not a minimal causal cut under current evidence")
        if self.document.get("kernel_found") is True and not declared:
            errors.append("kernel_found=true requires declared_kernel")
        if self.document.get("kernel_found") is True and not computed:
            errors.append("kernel_found=true requires at least one causal source-target path")
        if self.document.get("gap_status") == "OPEN" and self.document.get("kernel_found") is True:
            errors.append("gap_status=OPEN forbids kernel_found=true")
        return errors

    def assert_valid(self) -> None:
        errors = self.errors()
        if errors:
            raise OrientationKernelError("\n".join(errors))

    def _adjacency(self, causal_only: bool = True) -> dict[str, list[str]]:
        adjacency = {node_id: [] for node_id in self.nodes}
        for edge in self.edges:
            if edge.get("relation") == "FEEDBACK":
                continue
            if causal_only and edge.get("status") not in CAUSAL_STATUSES:
                continue
            start, end = edge.get("from"), edge.get("to")
            if start in adjacency and end in self.nodes:
                adjacency[start].append(end)
        return adjacency

    def causal_paths(self) -> list[tuple[str, ...]]:
        """Enumerate simple intervention-supported source -> target paths."""
        source, target = self.document.get("source"), self.document.get("target")
        if source not in self.nodes or target not in self.nodes:
            return []
        adjacency = self._adjacency(causal_only=True)
        paths: list[tuple[str, ...]] = []

        def visit(node: str, path: tuple[str, ...]) -> None:
            if node == target:
                paths.append(path)
                return
            for nxt in adjacency.get(node, []):
                if nxt not in path:
                    visit(nxt, path + (nxt,))

        visit(source, (source,))
        return paths

    def minimal_kernels(self) -> list[tuple[str, ...]]:
        """Return all smallest node cuts intercepting every validated causal path.

        Source and target are excluded. If no validated causal path exists, no
        kernel is returned: lack of evidence is not converted into a kernel.
        """
        paths = self.causal_paths()
        if not paths:
            return []
        source, target = self.document["source"], self.document["target"]
        interiors = [set(path[1:-1]) for path in paths]
        if any(not interior for interior in interiors):
            # A direct validated source -> target edge bypasses every mediator,
            # so no intermediate orientation kernel can yet be claimed.
            return []
        candidates = sorted(set().union(*interiors))
        for size in range(1, len(candidates) + 1):
            kernels = []
            for combo in itertools.combinations(candidates, size):
                cut = set(combo)
                if all(cut & interior for interior in interiors):
                    kernels.append(combo)
            if kernels:
                return kernels
        return []

    def _edge_verified(self, start: str, end: str) -> bool:
        return any(
            edge.get("from") == start
            and edge.get("to") == end
            and edge.get("status") == "VERIFIED"
            for edge in self.edges
        )

    def certificate(self, kernel: Iterable[str] | None = None) -> KernelCertificate | None:
        kernels = self.minimal_kernels()
        if kernel is None:
            if not kernels:
                return None
            selected = kernels[0]
        else:
            selected = tuple(kernel)
            if selected not in kernels:
                raise OrientationKernelError("requested kernel is not a minimal causal cut")

        paths = self.causal_paths()
        # Conservative verification rule: every edge on every current causal
        # path must be independently VERIFIED before the kernel certificate is.
        verified = all(
            all(self._edge_verified(a, b) for a, b in zip(path, path[1:]))
            for path in paths
        )
        return KernelCertificate(
            source=self.document["source"],
            target=self.document["target"],
            kernel=tuple(selected),
            causal_paths=len(paths),
            independently_verified=verified,
            status="VERIFIED_KERNEL" if verified else "INTERVENTION_SUPPORTED_KERNEL",
        )

    def summary(self) -> dict[str, Any]:
        paths = self.causal_paths()
        kernels = self.minimal_kernels()
        return {
            "graph_id": self.document.get("graph_id"),
            "case_id": self.document.get("case_id"),
            "gap_id": self.document.get("gap_id"),
            "source": self.document.get("source"),
            "target": self.document.get("target"),
            "validated_causal_paths": len(paths),
            "minimal_kernels": [list(k) for k in kernels],
            "kernel_found": self.document.get("kernel_found", False),
            "scientific_status": "OPEN" if not kernels else "CANDIDATE_CAUSAL_CUT",
        }


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("document", type=Path)
    parser.add_argument("--certificate", action="store_true")
    args = parser.parse_args(list(argv) if argv is not None else None)

    graph = OrientationKernelGraph.load(args.document)
    payload: Any = graph.summary()
    if args.certificate:
        cert = graph.certificate()
        payload = cert.__dict__ if cert else None
    print(json.dumps(payload, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
