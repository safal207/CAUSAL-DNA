#!/usr/bin/env python3
"""Causal Orientation Boundary (COB) gate for CAUSAL-DNA.

A COB is not a correlation threshold. It is a reproducible transition boundary
whose components precede the RNA transition, can be experimentally pushed
across the boundary, can be reversed/rescued, and are independently verified.
"""
from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

STATUSES = (
    "PROJECTIVE",
    "TEMPORALLY_SUPPORTED",
    "INTERVENTION_SUPPORTED",
    "REVERSIBLE",
    "VERIFIED",
)


class OrientationBoundaryError(ValueError):
    pass


@dataclass(frozen=True)
class BoundaryCertificate:
    boundary_id: str
    status: str
    components: tuple[str, ...]
    temporal_precedence: bool
    crossing_intervention: bool
    reversal_rescue: bool
    independent_verification: bool


class OrientationBoundaryGate:
    def __init__(self, document: dict[str, Any]):
        self.document = document
        self.boundaries = document.get("boundaries", [])

    @classmethod
    def load(cls, path: str | Path) -> "OrientationBoundaryGate":
        with Path(path).open(encoding="utf-8") as handle:
            gate = cls(json.load(handle))
        gate.assert_valid()
        return gate

    @staticmethod
    def _refs(boundary: dict[str, Any], key: str) -> tuple[str, ...]:
        evidence = boundary.get("evidence", {})
        refs = evidence.get(key, []) if isinstance(evidence, dict) else []
        if not isinstance(refs, list):
            return ()
        return tuple(ref for ref in refs if isinstance(ref, str) and ref.strip())

    def computed_status(self, boundary: dict[str, Any]) -> str:
        temporal = bool(self._refs(boundary, "temporal_precedence"))
        crossing = bool(self._refs(boundary, "boundary_crossing_intervention"))
        reversal = bool(self._refs(boundary, "reversal_rescue"))
        verification = bool(self._refs(boundary, "independent_verification"))
        if temporal and crossing and reversal and verification:
            return "VERIFIED"
        if temporal and crossing and reversal:
            return "REVERSIBLE"
        if temporal and crossing:
            return "INTERVENTION_SUPPORTED"
        if temporal:
            return "TEMPORALLY_SUPPORTED"
        return "PROJECTIVE"

    def certificate(self, boundary_id: str) -> BoundaryCertificate:
        boundary = next((b for b in self.boundaries if b.get("id") == boundary_id), None)
        if boundary is None:
            raise OrientationBoundaryError(f"unknown boundary: {boundary_id}")
        return BoundaryCertificate(
            boundary_id=boundary_id,
            status=self.computed_status(boundary),
            components=tuple(boundary.get("components", [])),
            temporal_precedence=bool(self._refs(boundary, "temporal_precedence")),
            crossing_intervention=bool(self._refs(boundary, "boundary_crossing_intervention")),
            reversal_rescue=bool(self._refs(boundary, "reversal_rescue")),
            independent_verification=bool(self._refs(boundary, "independent_verification")),
        )

    def errors(self) -> list[str]:
        errors: list[str] = []
        ids = [b.get("id") for b in self.boundaries]
        if not self.boundaries:
            errors.append("at least one candidate boundary is required")
        if any(not boundary_id for boundary_id in ids):
            errors.append("every boundary requires a non-empty id")
        if len(ids) != len(set(ids)):
            errors.append("boundary ids must be unique")

        for boundary in self.boundaries:
            boundary_id = boundary.get("id", "<unknown>")
            components = boundary.get("components")
            if not isinstance(components, list) or not components:
                errors.append(f"{boundary_id}: components must be a non-empty list")
            elif len(components) != len(set(components)):
                errors.append(f"{boundary_id}: components must be unique")

            evidence = boundary.get("evidence")
            if not isinstance(evidence, dict):
                errors.append(f"{boundary_id}: evidence must be an object")
                continue
            for key in (
                "temporal_precedence",
                "boundary_crossing_intervention",
                "reversal_rescue",
                "independent_verification",
            ):
                refs = evidence.get(key)
                if not isinstance(refs, list):
                    errors.append(f"{boundary_id}: evidence.{key} must be a list")
                elif any(not isinstance(ref, str) or not ref.strip() for ref in refs):
                    errors.append(f"{boundary_id}: evidence.{key} must contain non-empty strings")

            temporal = bool(self._refs(boundary, "temporal_precedence"))
            crossing = bool(self._refs(boundary, "boundary_crossing_intervention"))
            reversal = bool(self._refs(boundary, "reversal_rescue"))
            verification = bool(self._refs(boundary, "independent_verification"))
            if crossing and not temporal:
                errors.append(f"{boundary_id}: crossing intervention requires temporal precedence evidence")
            if reversal and not crossing:
                errors.append(f"{boundary_id}: reversal/rescue requires crossing intervention evidence")
            if verification and not reversal:
                errors.append(f"{boundary_id}: independent verification requires reversal/rescue evidence")

            computed = self.computed_status(boundary)
            declared = boundary.get("declared_status")
            if declared not in STATUSES:
                errors.append(f"{boundary_id}: invalid declared_status {declared!r}")
            elif declared != computed:
                errors.append(f"{boundary_id}: declared_status={declared} but evidence computes {computed}")

        verified = [b for b in self.boundaries if self.computed_status(b) == "VERIFIED"]
        if self.document.get("boundary_found") is True and not verified:
            errors.append("boundary_found=true requires at least one VERIFIED boundary")
        if self.document.get("gap_status") == "OPEN" and self.document.get("boundary_found") is True:
            errors.append("gap_status=OPEN forbids boundary_found=true")
        return errors

    def assert_valid(self) -> None:
        errors = self.errors()
        if errors:
            raise OrientationBoundaryError("\n".join(errors))

    def summary(self) -> dict[str, Any]:
        certs = [self.certificate(b["id"]) for b in self.boundaries]
        return {
            "gate_id": self.document.get("gate_id"),
            "case_id": self.document.get("case_id"),
            "gap_id": self.document.get("gap_id"),
            "gap_status": self.document.get("gap_status"),
            "boundary_found": self.document.get("boundary_found", False),
            "boundaries": {cert.boundary_id: cert.status for cert in certs},
        }


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("document", type=Path)
    args = parser.parse_args(list(argv) if argv is not None else None)
    gate = OrientationBoundaryGate.load(args.document)
    print(json.dumps(gate.summary(), indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
