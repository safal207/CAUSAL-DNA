#!/usr/bin/env python3
"""Causal Mediation Gate (CMG) for CAUSAL-DNA.

CMG prevents an associated molecular mediator from being promoted to a causal
materialized bridge. A mediator remains in Bardo until intervention and rescue
have both been observed in the declared context. Independent verification is
required for final materialization.

This module evaluates evidence structure; it does not decide whether a
biological result is true.
"""
from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

REQUIRED_EVIDENCE_KEYS = (
    "genotype_to_mediator",
    "genotype_to_outcome",
    "intervention_on_mediator",
    "rescue",
    "independent_verification",
)

STATUSES = (
    "PROJECTIVE",
    "BARDO_ASSOCIATED",
    "BARDO_INTERVENTION_SUPPORTED",
    "BARDO_RESCUE_SUPPORTED",
    "MATERIAL",
)

STATUS_TO_SPACE = {
    "PROJECTIVE": "projective",
    "BARDO_ASSOCIATED": "bardo",
    "BARDO_INTERVENTION_SUPPORTED": "bardo",
    "BARDO_RESCUE_SUPPORTED": "bardo",
    "MATERIAL": "material",
}


class MediationGateError(ValueError):
    """Raised when a mediation-gate document violates causal invariants."""


@dataclass(frozen=True)
class MediatorCertificate:
    mediator_id: str
    status: str
    space: str
    materialized: bool
    missing_requirements: tuple[str, ...]
    evidence_refs: tuple[str, ...]


class MediationGate:
    """Evaluate whether candidate mediators may cross Bardo into material space."""

    def __init__(self, document: dict[str, Any]):
        self.document = document
        self.mediators = document.get("mediators", [])

    @classmethod
    def load(cls, path: str | Path) -> "MediationGate":
        with Path(path).open(encoding="utf-8") as handle:
            gate = cls(json.load(handle))
        gate.assert_valid()
        return gate

    @staticmethod
    def _refs(candidate: dict[str, Any], key: str) -> tuple[str, ...]:
        evidence = candidate.get("evidence", {})
        values = evidence.get(key, [])
        if not isinstance(values, list):
            return ()
        return tuple(ref for ref in values if isinstance(ref, str) and ref.strip())

    def computed_status(self, candidate: dict[str, Any]) -> str:
        has_gm = bool(self._refs(candidate, "genotype_to_mediator"))
        has_gy = bool(self._refs(candidate, "genotype_to_outcome"))
        has_intervention = bool(self._refs(candidate, "intervention_on_mediator"))
        has_rescue = bool(self._refs(candidate, "rescue"))
        has_verification = bool(self._refs(candidate, "independent_verification"))

        # Materialization is deliberately conjunctive. No single evidence class
        # may silently substitute for another.
        if has_gm and has_gy and has_intervention and has_rescue and has_verification:
            return "MATERIAL"
        if has_gm and has_gy and has_intervention and has_rescue:
            return "BARDO_RESCUE_SUPPORTED"
        if has_gm and has_gy and has_intervention:
            return "BARDO_INTERVENTION_SUPPORTED"
        if has_gm and has_gy:
            return "BARDO_ASSOCIATED"
        return "PROJECTIVE"

    def certificate(self, mediator_id: str) -> MediatorCertificate:
        candidate = self.candidate(mediator_id)
        status = self.computed_status(candidate)
        refs: list[str] = []
        missing: list[str] = []
        for key in REQUIRED_EVIDENCE_KEYS:
            found = self._refs(candidate, key)
            refs.extend(found)
            if not found:
                missing.append(key)
        return MediatorCertificate(
            mediator_id=mediator_id,
            status=status,
            space=STATUS_TO_SPACE[status],
            materialized=status == "MATERIAL",
            missing_requirements=tuple(missing),
            evidence_refs=tuple(dict.fromkeys(refs)),
        )

    def candidate(self, mediator_id: str) -> dict[str, Any]:
        for candidate in self.mediators:
            if candidate.get("id") == mediator_id:
                return candidate
        raise MediationGateError(f"unknown mediator: {mediator_id}")

    def certificates(self) -> list[MediatorCertificate]:
        return [self.certificate(candidate["id"]) for candidate in self.mediators]

    def errors(self) -> list[str]:
        errors: list[str] = []
        ids = [candidate.get("id") for candidate in self.mediators]
        if not self.mediators:
            errors.append("at least one mediator candidate is required")
        if any(not value for value in ids):
            errors.append("every mediator requires a non-empty id")
        if len(ids) != len(set(ids)):
            errors.append("mediator ids must be unique")

        for candidate in self.mediators:
            mediator_id = candidate.get("id", "<unknown>")
            evidence = candidate.get("evidence")
            if not isinstance(evidence, dict):
                errors.append(f"{mediator_id}: evidence must be an object")
                continue

            for key in REQUIRED_EVIDENCE_KEYS:
                refs = evidence.get(key)
                if not isinstance(refs, list):
                    errors.append(f"{mediator_id}: evidence.{key} must be a list")
                elif any(not isinstance(ref, str) or not ref.strip() for ref in refs):
                    errors.append(f"{mediator_id}: evidence.{key} must contain non-empty strings")

            has_gm = bool(self._refs(candidate, "genotype_to_mediator"))
            has_gy = bool(self._refs(candidate, "genotype_to_outcome"))
            has_intervention = bool(self._refs(candidate, "intervention_on_mediator"))
            has_rescue = bool(self._refs(candidate, "rescue"))
            has_verification = bool(self._refs(candidate, "independent_verification"))

            if has_intervention and not (has_gm and has_gy):
                errors.append(
                    f"{mediator_id}: intervention cannot establish mediation before genotype→mediator and genotype→outcome evidence"
                )
            if has_rescue and not has_intervention:
                errors.append(f"{mediator_id}: rescue requires intervention evidence")
            if has_verification and not has_rescue:
                errors.append(f"{mediator_id}: independent verification requires rescue evidence")

            declared_status = candidate.get("declared_status")
            computed = self.computed_status(candidate)
            if declared_status not in STATUSES:
                errors.append(f"{mediator_id}: invalid declared_status {declared_status!r}")
            elif declared_status != computed:
                errors.append(
                    f"{mediator_id}: declared_status={declared_status} but evidence computes {computed}"
                )

            declared_space = candidate.get("declared_space")
            computed_space = STATUS_TO_SPACE[computed]
            if declared_space != computed_space:
                errors.append(
                    f"{mediator_id}: declared_space={declared_space!r} but evidence computes {computed_space!r}"
                )

        materialized = [c for c in self.certificates() if c.materialized] if self.mediators else []
        cause_found = self.document.get("cause_found")
        gap_status = self.document.get("gap_status")
        if cause_found is True and not materialized:
            errors.append("cause_found=true requires at least one MATERIAL mediator certificate")
        if gap_status == "OPEN" and cause_found is True:
            errors.append("gap_status=OPEN forbids cause_found=true")
        if gap_status == "OPEN" and materialized:
            errors.append("gap_status=OPEN forbids a MATERIAL mediator certificate")

        return errors

    def assert_valid(self) -> None:
        errors = self.errors()
        if errors:
            raise MediationGateError("\n".join(errors))

    def summary(self) -> dict[str, Any]:
        certs = self.certificates()
        counts = {status: 0 for status in STATUSES}
        for cert in certs:
            counts[cert.status] += 1
        return {
            "gate_id": self.document.get("gate_id"),
            "case_id": self.document.get("case_id"),
            "gap_id": self.document.get("gap_id"),
            "gap_status": self.document.get("gap_status"),
            "cause_found": self.document.get("cause_found"),
            "mediators": len(certs),
            "status_counts": counts,
            "materialized": [cert.mediator_id for cert in certs if cert.materialized],
        }


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("document", type=Path)
    parser.add_argument("--mediator")
    args = parser.parse_args(list(argv) if argv is not None else None)

    gate = MediationGate.load(args.document)
    if args.mediator:
        cert = gate.certificate(args.mediator)
        print(json.dumps(cert.__dict__, indent=2, ensure_ascii=False))
    else:
        print(json.dumps(gate.summary(), indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
