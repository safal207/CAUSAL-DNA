#!/usr/bin/env python3
"""Cross-check the exact mouse rs1421085 ortholog across mm10 and GRCm39.

This is a coordinate/provenance validator only. It does not promote any causal
hypothesis. The mm10 source coordinate is the unique published-guide match
resolved by gap001_map_mouse_ortholog.py; UCSC liftOver supplies the mm39 map.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

MM10_SITE = ("chr8", 91374371, 91374372)
MM10_TILE = ("chr8", 91374000, 91374500)
# MGI:8245447 / Rr695574 / Fto-Irx::hibE1, GRCm39.
HIBE1_MM39 = ("chr8", 92114859, 92119641)
# NCBI Gene 26383, GRCm39 current genomic context.
FTO_MM39 = ("chr8", 92039995, 92395061)


def read_bed(path: Path) -> tuple[str, int, int]:
    rows = [line.strip() for line in path.read_text().splitlines() if line.strip() and not line.startswith("#")]
    if len(rows) != 1:
        raise RuntimeError(f"expected exactly one BED row in {path}, found {len(rows)}")
    fields = rows[0].split("\t")
    if len(fields) < 3:
        raise RuntimeError(f"invalid BED row: {rows[0]}")
    chrom, start, end = fields[0], int(fields[1]), int(fields[2])
    if start < 0 or end <= start:
        raise RuntimeError(f"invalid interval: {(chrom, start, end)}")
    return chrom, start, end


def overlaps(a: tuple[str, int, int], b: tuple[str, int, int]) -> bool:
    return a[0] == b[0] and max(a[1], b[1]) < min(a[2], b[2])


def interval_distance(a: tuple[str, int, int], b: tuple[str, int, int]) -> int | None:
    if a[0] != b[0]:
        return None
    if overlaps(a, b):
        return 0
    if a[2] <= b[1]:
        return b[1] - a[2]
    return a[1] - b[2]


def inside(inner: tuple[str, int, int], outer: tuple[str, int, int]) -> bool:
    return inner[0] == outer[0] and outer[1] <= inner[1] and inner[2] <= outer[2]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--lifted-site", required=True, type=Path)
    ap.add_argument("--lifted-tile", required=True, type=Path)
    ap.add_argument("--output", required=True, type=Path)
    args = ap.parse_args()

    site39 = read_bed(args.lifted_site)
    tile39 = read_bed(args.lifted_tile)

    if site39[2] - site39[1] != 1:
        raise RuntimeError(f"1-bp site did not remain 1 bp after liftOver: {site39}")
    if not inside(site39, tile39):
        raise RuntimeError(f"lifted exact site is not inside lifted target tile: {site39} vs {tile39}")
    if not inside(site39, FTO_MM39):
        raise RuntimeError(f"lifted exact site is not inside current GRCm39 Fto interval: {site39}")

    result = {
        "human_variant": "rs1421085 T>C",
        "claim_scope": "assembly and coordinate provenance only; no causal promotion",
        "source": {
            "assembly": "GRCm38/mm10",
            "exact_site_0_based_half_open": list(MM10_SITE),
            "target_tile_0_based_half_open": list(MM10_TILE),
            "basis": "unique match of published WT CRISPR guide in UCSC mm10 chr8 FASTA",
        },
        "lifted": {
            "assembly": "GRCm39/mm39",
            "exact_site_0_based_half_open": list(site39),
            "exact_site_1_based": site39[1] + 1,
            "target_tile_0_based_half_open": list(tile39),
        },
        "cross_checks": {
            "inside_current_Fto_GRCm39": inside(site39, FTO_MM39),
            "Fto_GRCm39_0_based_half_open": list(FTO_MM39),
            "hibE1_GRCm39_0_based_half_open": list(HIBE1_MM39),
            "exact_site_overlaps_hibE1": overlaps(site39, HIBE1_MM39),
            "target_tile_overlaps_hibE1": overlaps(tile39, HIBE1_MM39),
            "exact_site_distance_to_hibE1_bp": interval_distance(site39, HIBE1_MM39),
            "target_tile_distance_to_hibE1_bp": interval_distance(tile39, HIBE1_MM39),
        },
        "provenance": [
            "UCSC mm10ToMm39.over.chain.gz",
            "NCBI Gene 26383 Fto: GRCm39 chr8:92039995-92395061",
            "MGI:8245447 Rr695574/Fto-Irx::hibE1: GRCm39 chr8:92114859-92119641",
            "Laber et al. 2021 PMID:34290091 published WT guide",
        ],
    }
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
