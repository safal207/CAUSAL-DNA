#!/usr/bin/env python3
"""Targeted MeCP2 ChIP/Input query around the mouse rs1421085 ortholog.

This parser reads processed WIG signal from GSE66868 and summarizes local
ChIP and matched Input signal around one coordinate. The result is a
feasibility/occupancy-context observation only. MeCP2 is broadly distributed
in mature neurons, so local ChIP/Input signal must not be promoted to causal
mediation or allele-specific binding.
"""
from __future__ import annotations

import argparse
import gzip
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, TextIO


@dataclass(frozen=True)
class WigPoint:
    chrom: str
    start0: int
    end0: int
    value: float


def _open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def iter_wig(path: Path) -> Iterator[WigPoint]:
    """Parse fixedStep, variableStep, or bedGraph-like WIG records."""
    mode = None
    chrom = None
    span = 1
    step = None
    next_start1 = None

    with _open_text(path) as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith(("track", "browser", "#")):
                continue
            if line.startswith("fixedStep"):
                fields = dict(x.split("=", 1) for x in line.split()[1:])
                mode = "fixed"
                chrom = fields["chrom"]
                next_start1 = int(fields["start"])
                step = int(fields.get("step", "1"))
                span = int(fields.get("span", "1"))
                continue
            if line.startswith("variableStep"):
                fields = dict(x.split("=", 1) for x in line.split()[1:])
                mode = "variable"
                chrom = fields["chrom"]
                span = int(fields.get("span", "1"))
                continue

            parts = line.split()
            if mode == "fixed":
                if len(parts) != 1 or chrom is None or next_start1 is None or step is None:
                    raise ValueError(f"invalid fixedStep record: {line}")
                start0 = next_start1 - 1
                value = float(parts[0])
                yield WigPoint(chrom, start0, start0 + span, value)
                next_start1 += step
            elif mode == "variable":
                if len(parts) != 2 or chrom is None:
                    raise ValueError(f"invalid variableStep record: {line}")
                start0 = int(parts[0]) - 1
                value = float(parts[1])
                yield WigPoint(chrom, start0, start0 + span, value)
            else:
                if len(parts) < 4:
                    raise ValueError(f"unsupported WIG/bedGraph record: {line}")
                yield WigPoint(parts[0], int(parts[1]), int(parts[2]), float(parts[3]))


def collect_window(path: Path, chrom: str, center0: int, radius: int) -> list[WigPoint]:
    low, high = center0 - radius, center0 + radius + 1
    out: list[WigPoint] = []
    for p in iter_wig(path):
        if p.chrom != chrom:
            continue
        if p.end0 <= low:
            continue
        if p.start0 >= high:
            # WIGs from GEO are coordinate sorted; safe early stop after the target window.
            if out:
                break
            continue
        out.append(p)
    return out


def weighted_mean(points: Iterable[WigPoint], low: int, high: int) -> float | None:
    weighted = 0.0
    covered = 0
    for p in points:
        overlap = max(0, min(p.end0, high) - max(p.start0, low))
        if overlap:
            weighted += p.value * overlap
            covered += overlap
    return weighted / covered if covered else None


def value_at(points: Iterable[WigPoint], pos0: int) -> float | None:
    vals = [p.value for p in points if p.start0 <= pos0 < p.end0]
    if not vals:
        return None
    return sum(vals) / len(vals)


def percentile_rank(values: list[float], value: float | None) -> float | None:
    if value is None or not values:
        return None
    return sum(v <= value for v in values) / len(values)


def summarize(chip: list[WigPoint], inp: list[WigPoint], center0: int) -> dict:
    radii = [250, 500, 2000, 10000, 50000]
    result = {"windows": {}}
    for r in radii:
        low, high = center0 - r, center0 + r + 1
        c = weighted_mean(chip, low, high)
        i = weighted_mean(inp, low, high)
        ratio = None
        log2_ratio = None
        if c is not None and i is not None and i > 0:
            ratio = c / i
            if ratio > 0:
                log2_ratio = math.log2(ratio)
        result["windows"][str(r)] = {
            "radius_bp": r,
            "chip_mean": c,
            "input_mean": i,
            "chip_over_input": ratio,
            "log2_chip_over_input": log2_ratio,
        }

    c0 = value_at(chip, center0)
    i0 = value_at(inp, center0)
    result["exact_base_signal"] = {
        "chip": c0,
        "input": i0,
        "chip_over_input": (c0 / i0 if c0 is not None and i0 not in (None, 0) else None),
    }

    # Contextualize the exact-bin ChIP signal against the +/-50 kb local neighborhood.
    local_chip_values = [p.value for p in chip if p.end0 > center0 - 50000 and p.start0 < center0 + 50001]
    local_input_values = [p.value for p in inp if p.end0 > center0 - 50000 and p.start0 < center0 + 50001]
    result["local_percentile"] = {
        "chip": percentile_rank(local_chip_values, c0),
        "input": percentile_rank(local_input_values, i0),
    }
    return result


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--chip", type=Path, required=True)
    ap.add_argument("--input", type=Path, required=True)
    ap.add_argument("--chrom", default="chr8")
    ap.add_argument("--position0", type=int, required=True)
    ap.add_argument("--assembly", default="mm9")
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    chip = collect_window(args.chip, args.chrom, args.position0, 50000)
    inp = collect_window(args.input, args.chrom, args.position0, 50000)
    if not chip or not inp:
        raise RuntimeError("target window absent from ChIP or Input WIG")

    payload = {
        "dataset": "GSE66868",
        "assay": "MeCP2 ChIP-seq versus matched Input",
        "tissue": "adult mouse hypothalamus",
        "assembly": args.assembly,
        "query": {"chrom": args.chrom, "position_0_based": args.position0, "position_1_based": args.position0 + 1},
        "summary": summarize(chip, inp, args.position0),
        "interpretation_guard": (
            "This is bulk hypothalamus, non-allele-specific ChIP context. MeCP2 is broadly distributed in mature neurons; "
            "local ChIP/Input enrichment is occupancy-compatible evidence only and cannot establish rs1421085 mediation, "
            "cell-type specificity, enhancer target, or causality."
        ),
        "discovery_level_ceiling": "D1_OCCUPANCY_CONTEXT",
        "cause_found": False,
    }
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
