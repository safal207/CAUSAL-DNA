#!/usr/bin/env python3
"""Scan the compact GSE246791 enhancer/cCRE archive at one exact locus.

The scanner treats genomic intervals as BED-style 0-based half-open coordinates.
It searches every text-like member of the tar archive for intervals overlapping
the exact mouse rs1421085 ortholog and reports both direct overlaps and nearest
intervals when no overlap is present.

A cCRE/accessibility overlap is descriptive evidence only. It does not establish
allele dependence, target gene, TF occupancy, or causal mediation.
"""
from __future__ import annotations

import argparse
import gzip
import heapq
import io
import json
import tarfile
from pathlib import Path
from typing import Any, BinaryIO


def interval_distance(position0: int, start: int, end: int) -> int:
    if start <= position0 < end:
        return 0
    if position0 < start:
        return start - position0
    return position0 - end + 1


def parse_interval(fields: list[str], offset: int = 0) -> tuple[str, int, int] | None:
    if len(fields) < offset + 3:
        return None
    chrom = fields[offset]
    try:
        start = int(fields[offset + 1])
        end = int(fields[offset + 2])
    except ValueError:
        return None
    if start < 0 or end <= start:
        return None
    return chrom, start, end


def scan_member(
    name: str,
    handle: BinaryIO,
    *,
    chrom: str,
    position0: int,
    nearest_k: int,
) -> tuple[list[dict[str, Any]], list[tuple[int, str, dict[str, Any]]], int]:
    overlaps: list[dict[str, Any]] = []
    nearest: list[tuple[int, str, dict[str, Any]]] = []
    parsed_intervals = 0
    for line_number, raw in enumerate(handle, start=1):
        if not raw or raw.startswith(b"#"):
            continue
        try:
            text = raw.decode("utf-8").strip()
        except UnicodeDecodeError:
            continue
        if not text:
            continue
        fields = text.split("\t")
        # Scan first BED triplet, and the second BEDPE triplet when present.
        for side, offset in (("left", 0), ("right", 3)):
            interval = parse_interval(fields, offset)
            if interval is None:
                continue
            current_chrom, start, end = interval
            if current_chrom != chrom:
                continue
            parsed_intervals += 1
            distance = interval_distance(position0, start, end)
            record = {
                "member": name,
                "line_number": line_number,
                "side": side,
                "chrom": current_chrom,
                "start_0_based": start,
                "end_0_based_exclusive": end,
                "distance_bp": distance,
                "raw_record": text[:2000],
            }
            if distance == 0:
                overlaps.append(record)
            else:
                key = f"{name}:{line_number}:{side}"
                item = (-distance, key, record)
                if len(nearest) < nearest_k:
                    heapq.heappush(nearest, item)
                elif distance < -nearest[0][0]:
                    heapq.heapreplace(nearest, item)
    return overlaps, nearest, parsed_intervals


def member_stream(name: str, extracted: BinaryIO) -> BinaryIO:
    if name.endswith(".gz"):
        return gzip.GzipFile(fileobj=extracted, mode="rb")
    return extracted


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--archive", required=True, type=Path)
    parser.add_argument("--chrom", default="chr8")
    parser.add_argument("--position0", type=int, default=91374371)
    parser.add_argument("--nearest-k", type=int, default=20)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    all_overlaps: list[dict[str, Any]] = []
    global_nearest: list[tuple[int, str, dict[str, Any]]] = []
    scanned_members = 0
    parsed_intervals = 0

    with tarfile.open(args.archive, "r:*") as archive:
        for member in archive.getmembers():
            if not member.isfile():
                continue
            extracted = archive.extractfile(member)
            if extracted is None:
                continue
            scanned_members += 1
            stream = member_stream(member.name, extracted)
            try:
                overlaps, nearest, count = scan_member(
                    member.name,
                    stream,
                    chrom=args.chrom,
                    position0=args.position0,
                    nearest_k=args.nearest_k,
                )
            except (OSError, EOFError):
                # Ignore non-text members that happen to use a .gz suffix or are corrupt;
                # the result reports the number of parsed chromosome intervals.
                continue
            parsed_intervals += count
            all_overlaps.extend(overlaps)
            for neg_distance, key, record in nearest:
                distance = -neg_distance
                item = (-distance, key, record)
                if len(global_nearest) < args.nearest_k:
                    heapq.heappush(global_nearest, item)
                elif distance < -global_nearest[0][0]:
                    heapq.heapreplace(global_nearest, item)

    nearest_sorted = [
        item[2]
        for item in sorted(global_nearest, key=lambda x: (-x[0], x[1]))
    ]
    all_overlaps.sort(
        key=lambda item: (item["member"], item["start_0_based"], item["line_number"])
    )
    result = {
        "dataset": "GSE246791",
        "assembly": "mm10",
        "query": {
            "chrom": args.chrom,
            "position_0_based": args.position0,
            "position_1_based": args.position0 + 1,
            "variant": "mouse ortholog of human rs1421085 T>C",
        },
        "archive_members_scanned": scanned_members,
        "chromosome_intervals_parsed": parsed_intervals,
        "overlap_count": len(all_overlaps),
        "overlaps": all_overlaps,
        "nearest_intervals": nearest_sorted,
        "interpretation_guard": (
            "An overlap means the coordinate falls inside a published accessible/cCRE "
            "interval represented in this archive. It does not show allele-specific "
            "accessibility, TF occupancy, enhancer-to-Irx3 targeting, or causality."
        ),
        "discovery_level_ceiling": "D1_DESCRIPTIVE",
    }
    text = json.dumps(result, indent=2, sort_keys=True)
    print(text)
    if args.output:
        args.output.write_text(text + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
