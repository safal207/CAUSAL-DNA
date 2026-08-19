#!/usr/bin/env python3
"""Necessity & Sufficiency Gate for CAUSAL-DNA.

This layer refines a candidate mediator after causal-mediation testing. It asks
whether the mediator is necessary, sufficient, partially mediating, redundant,
or non-load-bearing under an explicitly declared biological context.

The gate evaluates evidence structure and declared experimental outcomes. It
does not infer biological truth from observational association and never turns
planning assumptions into causal evidence.
"""
from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

NECESSITY_OUTCOMES = ("not_tested", "no_change", "attenuated", "abolished")
SUFFICIENCY_OUTCOMES = ("not_tested", "no_change", "partial", "recapitulated")
RESIDUAL_OUTCOMES = ("not_tested", "absent", "present")

ROLES = (
    "UNRESOLVED",
    "NECESSITY_SUPPORTED",
    "SUFFICIENCY_SUPPORTED",
    "PARTIAL_MEDIATOR",
    "NECESSARY_NOT_SUFFICIENT",
    "SUFFICIENT_NOT_NECESSARY",
    "NECESSARY_PARTIALLY_SUFFICIENT",
    "NECESSARY_AND_SUFFICIENT",
    "NON_LOAD_BEARING",
)


class NecessitySufficiencyError(ValueError):
    """Raised when a necessity/sufficiency document violates invariants."""


@dataclass(frozen=True)
class CausalRoleCertificate:
    mediator_id: str
    role: str
    independently_verified: bool
    verified_role: str | None
    necessity_outcome: str
    sufficiency_outcome: str
    residual_genotype_effect_after_block: str
    evidence_refs: tuple[str, ...]


class NecessitySufficiencyGate:
    def __init__(self, document: dict[str, Any]):
        self.document = document
        self.mediators = document.get("mediators", [])

    @classmethod
    def load(cls, path: str | Path) -> "NecessitySufficiencyGate":
        with Path(path).open(encoding="utf-8") as handle:
            gate = cls(json.load(handle))
        gate.assert_valid()
        return gate

    def candidate(self, mediator_id: str) -> dict[str, Any]:
        for candidate in self.mediators:
            if candidate.get("id") == mediator_id:
                return candidate
        raise NecessitySufficiencyError(f"unknown mediator: {mediator_id}")

    @staticmethod
    def _refs(candidate: dict[str, Any], section: str) -> tuple[str, ...]:
        block = candidate.get(section, {})
        values = block.get("evidence_refs", []) if isinstance(block, dict) else []
        if not isinstance(values, list):
            return ()
        return tuple(v for v in values if isinstance(v, str) and v.strip())

    def computed_role(self, candidate: dict[str, Any]) -> str:
        necessity = candidate.get("necessity_test", {}).get("outcome", "not_tested")
        sufficiency = candidate.get("sufficiency_test", {}).get("outcome", "not_tested")

        n_positive = necessity in {"attenuated", "abolished"}
        s_positive = sufficiency in {"partial", "recapitulated"}

        if necessity == "not_tested" and sufficiency == "not_tested":
            return "UNRESOLVED"

        if necessity == "no_change" and sufficiency == "no_change":
            return "NON_LOAD_BEARING"

        if necessity == "no_change" and s_positive:
            return "SUFFICIENT_NOT_NECESSARY"

        # Attenuation is definitionally partial mediation under this contract:
        # some genotype effect remains after mediator blockade. This precedence
        # must hold regardless of the sufficiency-test result.
        if necessity == "attenuated":
            return "PARTIAL_MEDIATOR"

        if n_positive and sufficiency == "no_change":
            return "NECESSARY_NOT_SUFFICIENT"

        if necessity == "abolished" and sufficiency == "partial":
            return "NECESSARY_PARTIALLY_SUFFICIENT"

        if necessity == "abolished" and sufficiency == "recapitulated":
            return "NECESSARY_AND_SUFFICIENT"

        if n_positive and sufficiency == "not_tested":
            return "NECESSITY_SUPPORTED"

        if necessity == "not_tested" and s_positive:
            return "SUFFICIENCY_SUPPORTED"

        return "UNRESOLVED"

    def certificate(self, mediator_id: str) -> CausalRoleCertificate:
        candidate = self.candidate(mediator_id)
        role = self.computed_role(candidate)
        verification_refs = self._refs(candidate, "independent_verification")
        refs = (
            list(self._refs(candidate, "necessity_test"))
            + list(self._refs(candidate, "sufficiency_test"))
            + list(verification_refs)
        )
        verified = bool(verification_refs) and role not in {"UNRESOLVED"}
        return CausalRoleCertificate(
            mediator_id=mediator_id,
            role=role,
            independently_verified=verified,
            verified_role=f"VERIFIED_{role}" if verified else None,
            necessity_outcome=candidate.get("necessity_test", {}).get("outcome", "not_tested"),
            sufficiency_outcome=candidate.get("sufficiency_test", {}).get("outcome", "not_tested"),
            residual_genotype_effect_after_block=candidate.get(
                "residual_genotype_effect_after_block", "not_tested"
            ),
            evidence_refs=tuple(dict.fromkeys(refs)),
        )

    def certificates(self) -> list[CausalRoleCertificate]:
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
            necessity = candidate.get("necessity_test")
            sufficiency = candidate.get("sufficiency_test")
            residual = candidate.get("residual_genotype_effect_after_block")

            if not isinstance(necessity, dict):
                errors.append(f"{mediator_id}: necessity_test must be an object")
                continue
            if not isinstance(sufficiency, dict):
                errors.append(f"{mediator_id}: sufficiency_test must be an object")
                continue

            n_outcome = necessity.get("outcome")
            s_outcome = sufficiency.get("outcome")
            if n_outcome not in NECESSITY_OUTCOMES:
                errors.append(f"{mediator_id}: invalid necessity outcome {n_outcome!r}")
            if s_outcome not in SUFFICIENCY_OUTCOMES:
                errors.append(f"{mediator_id}: invalid sufficiency outcome {s_outcome!r}")
            if residual not in RESIDUAL_OUTCOMES:
                errors.append(f"{mediator_id}: invalid residual genotype effect {residual!r}")

            for section, outcome in (("necessity_test", n_outcome), ("sufficiency_test", s_outcome)):
                refs = candidate.get(section, {}).get("evidence_refs")
                if not isinstance(refs, list):
                    errors.append(f"{mediator_id}: {section}.evidence_refs must be a list")
                elif any(not isinstance(ref, str) or not ref.strip() for ref in refs):
                    errors.append(f"{mediator_id}: {section}.evidence_refs must contain non-empty strings")
                elif outcome != "not_tested" and not refs:
                    errors.append(f"{mediator_id}: tested {section} requires evidence_refs")
                elif outcome == "not_tested" and refs:
                    errors.append(f"{mediator_id}: untested {section} cannot carry evidence_refs")

            verification = candidate.get("independent_verification")
            if not isinstance(verification, dict):
                errors.append(f"{mediator_id}: independent_verification must be an object")
            else:
                refs = verification.get("evidence_refs")
                if not isinstance(refs, list):
                    errors.append(f"{mediator_id}: independent_verification.evidence_refs must be a list")
                elif any(not isinstance(ref, str) or not ref.strip() for ref in refs):
                    errors.append(
                        f"{mediator_id}: independent_verification.evidence_refs must contain non-empty strings"
                    )

            if residual != "not_tested" and n_outcome == "not_tested":
                errors.append(
                    f"{mediator_id}: residual genotype effect requires a completed necessity test"
                )
            if n_outcome == "abolished" and residual == "present":
                errors.append(
                    f"{mediator_id}: necessity outcome abolished conflicts with residual genotype effect present"
                )
            if n_outcome == "attenuated" and residual == "absent":
                errors.append(
                    f"{mediator_id}: attenuated necessity outcome conflicts with residual genotype effect absent"
                )

            computed = self.computed_role(candidate)
            declared = candidate.get("declared_role")
            if declared not in ROLES:
                errors.append(f"{mediator_id}: invalid declared_role {declared!r}")
            elif declared != computed:
                errors.append(
                    f"{mediator_id}: declared_role={declared} but evidence computes {computed}"
                )

            verification_refs = self._refs(candidate, "independent_verification")
            if verification_refs and computed == "UNRESOLVED":
                errors.append(
                    f"{mediator_id}: independent verification cannot verify an unresolved causal role"
                )

        verified_load_bearing = [
            cert
            for cert in self.certificates()
            if cert.independently_verified and cert.role not in {"UNRESOLVED", "NON_LOAD_BEARING"}
        ] if self.mediators else []

        if self.document.get("causal_role_found") is True and not verified_load_bearing:
            errors.append(
                "causal_role_found=true requires an independently verified load-bearing mediator role"
            )
        if self.document.get("gap_status") == "OPEN" and self.document.get("causal_role_found") is True:
            errors.append("gap_status=OPEN forbids causal_role_found=true")

        return errors

    def assert_valid(self) -> None:
        errors = self.errors()
        if errors:
            raise NecessitySufficiencyError("\n".join(errors))

    def summary(self) -> dict[str, Any]:
        certs = self.certificates()
        counts = {role: 0 for role in ROLES}
        for cert in certs:
            counts[cert.role] += 1
        return {
            "gate_id": self.document.get("gate_id"),
            "case_id": self.document.get("case_id"),
            "gap_id": self.document.get("gap_id"),
            "gap_status": self.document.get("gap_status"),
            "causal_role_found": self.document.get("causal_role_found"),
            "mediators": len(certs),
            "role_counts": counts,
            "verified_roles": {
                cert.mediator_id: cert.verified_role
                for cert in certs
                if cert.verified_role is not None
            },
        }


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("document", type=Path)
    parser.add_argument("--mediator")
    args = parser.parse_args(list(argv) if argv is not None else None)

    gate = NecessitySufficiencyGate.load(args.document)
    if args.mediator:
        print(json.dumps(gate.certificate(args.mediator).__dict__, indent=2, ensure_ascii=False))
    else:
        print(json.dumps(gate.summary(), indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
