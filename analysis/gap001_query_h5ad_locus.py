#!/usr/bin/env python3
"""Query exact-locus accessibility in one GSE246791 SnapATAC2 h5ad.

Reads only persisted sparse-matrix indices plus small targeted data slices.
Reports the exact 500-bp tile and raw Tn5 insertion activity around the mapped
mouse ortholog of human rs1421085, preserving the exact H5AD nucleus barcodes
for downstream joins to the authors' official cell annotations.

This is descriptive wild-type atlas evidence. It cannot establish allele
specificity, enhancer target, TF occupancy, or causal mediation.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import h5py
import numpy as np


def decode(value):
    if isinstance(value, bytes):
        return value.decode("utf-8")
    return str(value)


def interval_layout(var_index: h5py.Dataset, x_cols: int):
    names = [decode(v) for v in var_index.attrs["names"]]
    intervals = np.asarray(var_index.attrs["intervals"], dtype=np.int64)
    if len(names) != len(intervals):
        raise RuntimeError("interval names/metadata length mismatch")
    counts = []
    for start, end, step, width in intervals:
        if step <= 0 or width <= 0 or end <= start:
            raise RuntimeError("invalid interval index metadata")
        counts.append(math.ceil((int(end) - int(start)) / int(step)))
    if sum(counts) != x_cols:
        raise RuntimeError(
            f"derived tile count {sum(counts)} does not match X columns {x_cols}"
        )
    return names, intervals, counts


def tile_column(names, intervals, counts, chrom: str, position0: int):
    try:
        chrom_i = names.index(chrom)
    except ValueError as exc:
        raise RuntimeError(f"chromosome not found in tile index: {chrom}") from exc
    start, end, step, width = [int(v) for v in intervals[chrom_i]]
    if not (start <= position0 < end):
        raise RuntimeError("position lies outside chromosome interval metadata")
    local_i = (position0 - start) // step
    col = sum(counts[:chrom_i]) + local_i
    tile_start = start + local_i * step
    tile_end = tile_start + width
    return col, tile_start, tile_end


def csr_hits(group: h5py.Group, low_col: int, high_col_exclusive: int):
    indices = group["indices"][:]
    mask = (indices >= low_col) & (indices < high_col_exclusive)
    positions = np.flatnonzero(mask)
    cols = indices[positions]
    indptr = group["indptr"][:]
    rows = np.searchsorted(indptr, positions, side="right") - 1
    values = group["data"][positions] if len(positions) else np.array([], dtype=np.uint8)
    return rows.astype(np.int64), cols.astype(np.int64), np.asarray(values)


def cell_records(rows: np.ndarray, obs_index: list[str]) -> list[dict]:
    records = []
    for row in sorted(int(v) for v in np.unique(rows)):
        if row < 0 or row >= len(obs_index):
            raise RuntimeError(f"sparse row {row} is outside obs/index")
        records.append({"row_index": row, "barcode": obs_index[row]})
    return records


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("h5ad", type=Path)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--geo-accession", required=True)
    parser.add_argument("--chrom", default="chr8")
    parser.add_argument("--position0", type=int, default=91374371)
    parser.add_argument("--neighbor-bins", type=int, default=10)
    parser.add_argument("--insertion-radius", type=int, default=250)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    with h5py.File(args.h5ad, "r") as f:
        n_cells, n_tiles = [int(v) for v in f["X"].attrs["shape"]]
        obs_index = [decode(v) for v in f["obs/index"][:]]
        if len(obs_index) != n_cells:
            raise RuntimeError(
                f"obs/index has {len(obs_index)} rows but X has {n_cells} nuclei"
            )
        if len(set(obs_index)) != len(obs_index):
            raise RuntimeError("obs/index nucleus barcodes are not unique within sample")

        names, intervals, counts = interval_layout(f["var/index"], n_tiles)
        target_col, tile_start, tile_end = tile_column(
            names, intervals, counts, args.chrom, args.position0
        )

        lo_col = max(0, target_col - args.neighbor_bins)
        hi_col = min(n_tiles, target_col + args.neighbor_bins + 1)
        rows, cols, values = csr_hits(f["X"], lo_col, hi_col)

        bin_summaries = []
        target_rows = np.array([], dtype=np.int64)
        for col in range(lo_col, hi_col):
            m = cols == col
            hit_rows = rows[m]
            hit_values = values[m]
            if col == target_col:
                target_rows = np.unique(hit_rows)
            offset = col - target_col
            start = tile_start + offset * 500
            end = start + 500
            bin_summaries.append(
                {
                    "column": col,
                    "offset_bins": offset,
                    "start_0_based": start,
                    "end_0_based_exclusive": end,
                    "nuclei_nonzero": int(len(np.unique(hit_rows))),
                    "fraction_nuclei_nonzero": (
                        float(len(np.unique(hit_rows)) / n_cells) if n_cells else 0.0
                    ),
                    "sum_tile_counts": int(hit_values.sum()) if len(hit_values) else 0,
                }
            )

        target_summary = next(x for x in bin_summaries if x["offset_bins"] == 0)
        target_summary["var_count"] = float(f["var/count"][target_col])
        target_summary["selected"] = bool(f["var/selected"][target_col])
        try:
            target_summary["var_index_label"] = decode(f["var/index"][target_col])
        except Exception:
            target_summary["var_index_label"] = None

        ref_names = [decode(v) for v in f["uns/reference_sequences/reference_seq_name"][:]]
        ref_lengths = [int(v) for v in f["uns/reference_sequences/reference_seq_length"][:]]
        insertion_shape = [int(v) for v in f["obsm/insertion"].attrs["shape"]]
        if sum(ref_lengths) != insertion_shape[1]:
            raise RuntimeError(
                "reference sequence lengths do not match insertion matrix columns"
            )
        try:
            ref_i = ref_names.index(args.chrom)
        except ValueError as exc:
            raise RuntimeError(f"chromosome not found in insertion reference: {args.chrom}") from exc
        chrom_offset = sum(ref_lengths[:ref_i])
        low_pos = max(0, args.position0 - args.insertion_radius)
        high_pos = min(ref_lengths[ref_i], args.position0 + args.insertion_radius + 1)
        low_flat = chrom_offset + low_pos
        high_flat = chrom_offset + high_pos
        i_rows, i_cols, i_values = csr_hits(f["obsm/insertion"], low_flat, high_flat)
        genomic_positions = i_cols - chrom_offset
        insertion_rows = np.unique(i_rows)

        insertion_by_position = []
        if len(genomic_positions):
            unique_positions = np.unique(genomic_positions)
            for pos in unique_positions:
                m = genomic_positions == pos
                insertion_by_position.append(
                    {
                        "position_0_based": int(pos),
                        "relative_to_variant_bp": int(pos - args.position0),
                        "nuclei_with_insertion": int(len(np.unique(i_rows[m]))),
                        "sum_insertion_counts": int(i_values[m].sum()),
                    }
                )

        insertion_window = {
            "radius_bp": args.insertion_radius,
            "start_0_based": low_pos,
            "end_0_based_exclusive": high_pos,
            "nuclei_with_any_insertion": int(len(insertion_rows)),
            "fraction_nuclei_with_any_insertion": (
                float(len(insertion_rows) / n_cells) if n_cells else 0.0
            ),
            "total_insertion_counts": int(i_values.sum()) if len(i_values) else 0,
            "positions_with_signal": int(len(np.unique(genomic_positions))),
            "exact_variant_base_nuclei": int(
                len(np.unique(i_rows[genomic_positions == args.position0]))
            ),
            "by_position": insertion_by_position,
        }

        both_rows = np.intersect1d(target_rows, insertion_rows)
        union_rows = np.union1d(target_rows, insertion_rows)
        signal_cells = {
            "barcode_source": "H5AD obs/index",
            "target_tile": cell_records(target_rows, obs_index),
            "insertion_window": cell_records(insertion_rows, obs_index),
            "both_target_tile_and_insertion_window": cell_records(both_rows, obs_index),
            "either_signal": cell_records(union_rows, obs_index),
        }

        result = {
            "dataset": "GSE246791",
            "sample_id": args.sample_id,
            "geo_accession": args.geo_accession,
            "assembly": "mm10",
            "query": {
                "chrom": args.chrom,
                "position_0_based": args.position0,
                "position_1_based": args.position0 + 1,
                "variant": "mouse ortholog of human rs1421085 T>C",
            },
            "nuclei": n_cells,
            "tile_matrix_columns": n_tiles,
            "target_500bp_tile": target_summary,
            "neighbor_tiles": bin_summaries,
            "raw_insertion_window": insertion_window,
            "signal_carrying_nuclei": signal_cells,
            "interpretation_guard": (
                "Wild-type adult male hypothalamus atlas signal is descriptive. "
                "Presence does not establish allele dependence or enhancer-to-Irx3 causality; "
                "absence does not exclude rare-cell or T>C-induced accessibility. "
                "Barcode identities are preserved only for an official annotation join."
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
