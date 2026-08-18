#!/usr/bin/env python3
"""Aggregate per-sample GSE246791 locus queries without pooling away replicates."""
from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("input_dir", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    files = sorted(args.input_dir.rglob("*-rs1421085-locus.json"))
    if not files:
        raise RuntimeError("no locus result JSON files found")

    rows = []
    by_region = defaultdict(list)
    for path in files:
        p = json.loads(path.read_text(encoding="utf-8"))
        sample = p["sample_id"]
        region = sample.rsplit("_", 1)[-1]
        target = p["target_500bp_tile"]
        insertion = p["raw_insertion_window"]
        neighbor_fracs = [float(x["fraction_nuclei_nonzero"]) for x in p["neighbor_tiles"]]
        sorted_fracs = sorted(neighbor_fracs)
        median = sorted_fracs[len(sorted_fracs) // 2]
        row = {
            "region": region,
            "sample_id": sample,
            "geo_accession": p["geo_accession"],
            "nuclei": int(p["nuclei"]),
            "target_tile_nuclei": int(target["nuclei_nonzero"]),
            "target_tile_fraction": float(target["fraction_nuclei_nonzero"]),
            "target_tile_counts": int(target["sum_tile_counts"]),
            "target_tile_selected": bool(target["selected"]),
            "local_21bin_median_fraction": median,
            "target_over_local_median": (
                float(target["fraction_nuclei_nonzero"]) / median if median > 0 else None
            ),
            "insertion_window_nuclei": int(insertion["nuclei_with_any_insertion"]),
            "insertion_window_fraction": float(insertion["fraction_nuclei_with_any_insertion"]),
            "insertion_window_counts": int(insertion["total_insertion_counts"]),
            "insertion_positions_with_signal": int(insertion["positions_with_signal"]),
            "exact_variant_base_nuclei": int(insertion["exact_variant_base_nuclei"]),
        }
        rows.append(row)
        by_region[region].append(row)

    region_summary = []
    for region, reps in sorted(by_region.items()):
        region_summary.append(
            {
                "region": region,
                "replicate_count": len(reps),
                "replicates_with_target_tile_signal": sum(r["target_tile_nuclei"] > 0 for r in reps),
                "replicates_with_insertion_window_signal": sum(r["insertion_window_nuclei"] > 0 for r in reps),
                "target_tile_fraction_range": [
                    min(r["target_tile_fraction"] for r in reps),
                    max(r["target_tile_fraction"] for r in reps),
                ],
                "insertion_window_fraction_range": [
                    min(r["insertion_window_fraction"] for r in reps),
                    max(r["insertion_window_fraction"] for r in reps),
                ],
            }
        )

    result = {
        "dataset": "GSE246791",
        "assembly": "mm10",
        "query": {
            "chrom": "chr8",
            "position_0_based": 91374371,
            "position_1_based": 91374372,
            "variant": "mouse ortholog of human rs1421085 T>C",
        },
        "sample_count": len(rows),
        "samples_with_target_tile_signal": sum(r["target_tile_nuclei"] > 0 for r in rows),
        "samples_with_insertion_window_signal": sum(r["insertion_window_nuclei"] > 0 for r in rows),
        "samples": rows,
        "regions": region_summary,
        "interpretation_guard": (
            "Replicate consistency within one wild-type atlas strengthens a descriptive D1 signal "
            "but is not independent-dataset replication (D2), allele specificity, cell-type "
            "localization, or causal mediation."
        ),
        "cause_found": False,
    }
    text = json.dumps(result, indent=2, sort_keys=True)
    print(text)
    if args.output:
        args.output.write_text(text + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
