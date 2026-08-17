#!/usr/bin/env python3
"""Append-only causal processor for CAUSAL-DNA.

Material evidence is immutable history. New evidence creates a new generation
of the model projection instead of rewriting prior generations.
"""
from __future__ import annotations

import copy
from dataclasses import dataclass
from typing import Any


class CausalProcessorError(ValueError):
    pass


@dataclass(frozen=True)
class Projection:
    generation: int
    open_hypotheses: tuple[str, ...]
    rejected_hypotheses: tuple[str, ...]
    superseded_hypotheses: tuple[str, ...]
    verified_subjects: tuple[str, ...]


class CausalProcessor:
    def __init__(self, events: list[dict[str, Any]] | None = None):
        self.events = copy.deepcopy(events or [])
        self.assert_valid()

    def assert_valid(self) -> None:
        ids = [e.get("event_id") for e in self.events]
        if len(ids) != len(set(ids)):
            raise CausalProcessorError("event_id must be unique")
        last_generation = -1
        for event in self.events:
            generation = event.get("generation")
            if not isinstance(generation, int) or generation < last_generation:
                raise CausalProcessorError("event generations must be monotonic")
            last_generation = generation
            if event.get("event_type") in {"observation","perturbation","verification"} and not event.get("evidence_refs"):
                raise CausalProcessorError("material/verification events require evidence_refs")
            if event.get("event_type") == "verification" and event.get("observer") != "independent_verifier":
                raise CausalProcessorError("verification must be asserted by independent_verifier")

    @property
    def generation(self) -> int:
        return max((e["generation"] for e in self.events), default=0)

    def append(self, event: dict[str, Any]) -> None:
        if any(e["event_id"] == event.get("event_id") for e in self.events):
            raise CausalProcessorError("cannot rewrite existing event_id")
        if event.get("generation", -1) < self.generation:
            raise CausalProcessorError("cannot append event to historical generation")
        candidate = self.events + [copy.deepcopy(event)]
        CausalProcessor(candidate)
        self.events.append(copy.deepcopy(event))

    def next_generation(self, events: list[dict[str, Any]]) -> int:
        gen = self.generation + 1
        for event in events:
            if event.get("generation") != gen:
                raise CausalProcessorError("all update events must target next generation")
        candidate = self.events + copy.deepcopy(events)
        CausalProcessor(candidate)
        self.events = candidate
        return gen

    def project(self) -> Projection:
        open_h, rejected, superseded, verified = set(), set(), set(), set()
        for event in self.events:
            subject = event["subject"]
            kind = event["event_type"]
            if kind == "hypothesis_opened":
                open_h.add(subject)
            elif kind == "hypothesis_rejected":
                open_h.discard(subject); rejected.add(subject)
            elif kind == "hypothesis_superseded":
                open_h.discard(subject); superseded.add(subject)
            elif kind == "verification":
                verified.add(subject)
        return Projection(
            self.generation,
            tuple(sorted(open_h)),
            tuple(sorted(rejected)),
            tuple(sorted(superseded)),
            tuple(sorted(verified)),
        )
