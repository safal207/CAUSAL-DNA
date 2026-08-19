#!/usr/bin/env python3
"""4D causal lattice for CAUSAL-DNA.

Axes:
- space: projective / bardo / material
- time: past / present / future
- evidence: none / supportive / perturbed / rescued / established
- observer: model / experiment / independent_verifier

This is an epistemic representation. It tracks how a proposed causal state is
projected, tested, materialized, and later used to update the model without
confusing future hypotheses with past or present observations.
"""

from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path
from typing import Any

TIME_ORDER = {"past": 0, "present": 1, "future": 2}
EVIDENCE_ORDER = {"none": 0, "supportive": 1, "perturbed": 2, "rescued": 3, "established": 4}


class TemporalLatticeError(ValueError):
    pass


class TemporalLattice:
    def __init__(self, document: dict[str, Any]):
        self.document = document
        self.states = document.get("states", [])
        self.transitions = document.get("transitions", [])
        self._states = {state.get("id"): state for state in self.states}
        self._out: dict[str, list[dict[str, Any]]] = defaultdict(list)
        for transition in self.transitions:
            self._out[transition.get("from")].append(transition)

    @classmethod
    def load(cls, path: str | Path) -> "TemporalLattice":
        with Path(path).open(encoding="utf-8") as handle:
            lattice = cls(json.load(handle))
        lattice.assert_valid()
        return lattice

    def errors(self) -> list[str]:
        errors: list[str] = []
        ids = [state.get("id") for state in self.states]
        if len(ids) != len(set(ids)):
            errors.append("state ids must be unique")

        known = set(ids)
        for state in self.states:
            sid = state.get("id", "<unknown>")
            evidence_refs = state.get("evidence_refs", [])
            if state.get("space") == "material" and not evidence_refs:
                errors.append(f"{sid}: material state requires evidence_refs")
            if state.get("time") == "future" and state.get("space") == "material":
                errors.append(f"{sid}: future state cannot be material fact")
            if state.get("evidence") in {"perturbed", "rescued", "established"} and not evidence_refs:
                errors.append(f"{sid}: strong evidence level requires provenance")
            if state.get("observer") == "independent_verifier" and not evidence_refs:
                errors.append(f"{sid}: independent verifier state requires provenance")

        for tr in self.transitions:
            tid = tr.get("id", "<unknown>")
            src = tr.get("from")
            dst = tr.get("to")
            if src not in known or dst not in known:
                errors.append(f"{tid}: transition references unknown state")
                continue
            a = self._states[src]
            b = self._states[dst]

            if tr.get("relation") == "materializes":
                if a.get("space") != "bardo" or b.get("space") != "material":
                    errors.append(f"{tid}: materializes must be Bardo→Material")
                if not tr.get("evidence_refs"):
                    errors.append(f"{tid}: materialization requires evidence_refs")

            if a.get("space") == "projective" and b.get("space") == "material":
                errors.append(f"{tid}: projective→material bypass forbidden")

            # Knowledge may point back from material evidence to a new future model,
            # but forward causal/testing transitions may not claim future evidence
            # happened before its source.
            if tr.get("relation") in {"projects", "enters_bardo", "tests", "materializes", "bounds"}:
                if TIME_ORDER[b.get("time")] < TIME_ORDER[a.get("time")]:
                    errors.append(f"{tid}: forward transition reverses temporal order")

        return errors

    def assert_valid(self) -> None:
        errors = self.errors()
        if errors:
            raise TemporalLatticeError("\n".join(errors))

    def states_at(self, *, space: str | None = None, time: str | None = None) -> list[dict[str, Any]]:
        out = self.states
        if space is not None:
            out = [s for s in out if s.get("space") == space]
        if time is not None:
            out = [s for s in out if s.get("time") == time]
        return out

    def future_frontier(self) -> list[dict[str, Any]]:
        return [
            s for s in self.states
            if s.get("time") == "future" and s.get("space") in {"projective", "bardo"}
        ]

    def materialized_history(self) -> list[dict[str, Any]]:
        return [
            s for s in self.states
            if s.get("space") == "material" and s.get("time") in {"past", "present"}
        ]

    def observer_views(self) -> dict[str, list[str]]:
        views: dict[str, list[str]] = defaultdict(list)
        for state in self.states:
            views[state.get("observer")].append(state.get("id"))
        return dict(views)

    def summary(self) -> dict[str, Any]:
        return {
            "lattice_id": self.document.get("lattice_id"),
            "states": len(self.states),
            "transitions": len(self.transitions),
            "future_frontier": [s["id"] for s in self.future_frontier()],
            "materialized_history": [s["id"] for s in self.materialized_history()],
            "observer_views": self.observer_views(),
        }
