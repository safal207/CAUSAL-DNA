#!/usr/bin/env python3
"""GAP-001: join exact-locus ATAC fragments to RNA in the same GSE226277 nuclei.

The dataset is wild-type and not an rs1421085 genotype experiment. This pass
therefore tests only whether the exact mouse-ortholog locus and Irx3/candidate
transcripts can be observed in the same processed Multiome nucleus. It cannot
establish allele effects, TF occupancy, enhancer direction, or causality.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
from pathlib import Path
import urllib.request

import h5py
import numpy as np
from scipy.sparse import csc_matrix

GENES = {
    "IRX3": "Irx3",
    "ARID5B": "Arid5b",
    "CUX1": "Cux1",
    "TET1": "Tet1",
    "KDM2B": "Kdm2b",
    "DNMT1": "Dnmt1",
    "MECP2": "Mecp2",
    "CXXC1": "Cxxc1",
}
FORBIDDEN_CAUSAL_FIELDS = {
    "cause_found", "causal_status", "materialized", "verified", "mediator_found"
}


def _decode(values) -> np.ndarray:
    return np.asarray([
        x.decode("utf-8") if isinstance(x, (bytes, np.bytes_)) else str(x)
        for x in values
    ])


def load_rna(path: Path) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    with h5py.File(path, "r") as h5:
        if "matrix" not in h5:
            raise ValueError(f"{path}: missing /matrix")
        g = h5["matrix"]
        names = _decode(g["features"]["name"][:])
        rows: dict[str, int] = {}
        for label, gene in GENES.items():
            hits = np.flatnonzero(names == gene)
            if len(hits) != 1:
                raise ValueError(f"{path}: expected exactly one {gene}, got {len(hits)}")
            rows[label] = int(hits[0])
        shape = tuple(int(x) for x in g["shape"][:])
        matrix = csc_matrix((g["data"][:], g["indices"][:], g["indptr"][:]), shape=shape)
        barcodes = _decode(g["barcodes"][:])
        if len(barcodes) != shape[1]:
            raise ValueError(f"{path}: barcode count does not match cell count")
        expr = {
            label: np.asarray(matrix[row, :].todense()).ravel()
            for label, row in rows.items()
        }
    return barcodes, expr


def stream_locus_fragments(url: str, chrom: str, start: int, end: int) -> tuple[dict[str, dict], dict]:
    if start < 0 or end <= start:
        raise ValueError("invalid locus interval")
    request = urllib.request.Request(url, headers={"User-Agent": "CAUSAL-DNA/0.1 locus-scan"})
    overlap: dict[str, dict] = {}
    lines_read = 0
    target_lines = 0
    seen_target = False
    stopped_after_target = False
    with urllib.request.urlopen(request, timeout=120) as response:
        with gzip.GzipFile(fileobj=response, mode="rb") as gz:
            for raw in gz:
                lines_read += 1
                if not raw or raw.startswith(b"#"):
                    continue
                parts = raw.rstrip(b"\n").split(b"\t")
                if len(parts) < 5:
                    raise ValueError(f"malformed fragment row at line {lines_read}")
                row_chrom = parts[0].decode("utf-8")
                if row_chrom != chrom:
                    if seen_target:
                        stopped_after_target = True
                        break
                    continue
                seen_target = True
                f_start = int(parts[1])
                f_end = int(parts[2])
                if f_start >= end:
                    stopped_after_target = True
                    break
                target_lines += 1
                if f_end <= start or f_start >= end:
                    continue
                barcode = parts[3].decode("utf-8")
                support = int(parts[4])
                state = overlap.setdefault(barcode, {"fragment_records": 0, "read_support": 0})
                state["fragment_records"] += 1
                state["read_support"] += support
    return overlap, {
        "lines_read_until_stop": lines_read,
        "target_chromosome_rows_examined": target_lines,
        "target_chromosome_seen": seen_target,
        "early_stop_after_locus": stopped_after_target,
    }


def build_join(barcodes: np.ndarray, expr: dict[str, np.ndarray], locus: dict[str, dict]) -> tuple[list[dict], dict]:
    index = {barcode: i for i, barcode in enumerate(barcodes.tolist())}
    joined = []
    for barcode in sorted(set(locus) & set(index)):
        i = index[barcode]
        row = {
            "barcode": barcode,
            "fragment_records": locus[barcode]["fragment_records"],
            "read_support": locus[barcode]["read_support"],
        }
        for label in GENES:
            row[f"{label.lower()}_umi"] = int(expr[label][i])
        joined.append(row)

    irx_all = expr["IRX3"] > 0
    summary = {
        "rna_filtered_nuclei": int(len(barcodes)),
        "rna_irx3_positive": int(irx_all.sum()),
        "raw_locus_fragment_barcodes": int(len(locus)),
        "locus_barcodes_in_filtered_rna": int(len(joined)),
        "locus_irx3_positive": int(sum(r["irx3_umi"] > 0 for r in joined)),
        "candidate_same_nucleus": {},
    }
    for label in GENES:
        if label == "IRX3":
            continue
        ckey = f"{label.lower()}_umi"
        candidate_locus = sum(r[ckey] > 0 for r in joined)
        candidate_irx3_locus = sum(r[ckey] > 0 and r["irx3_umi"] > 0 for r in joined)
        summary["candidate_same_nucleus"][label] = {
            "candidate_positive_among_locus_barcodes": int(candidate_locus),
            "irx3_candidate_positive_among_locus_barcodes": int(candidate_irx3_locus),
        }
    return joined, summary


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--pair-id", required=True)
    p.add_argument("--rna-h5", type=Path, required=True)
    p.add_argument("--atac-url", required=True)
    p.add_argument("--chrom", default="chr8")
    p.add_argument("--start", type=int, default=91374000)
    p.add_argument("--end", type=int, default=91374500)
    p.add_argument("--output-dir", type=Path, required=True)
    args = p.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    barcodes, expr = load_rna(args.rna_h5)
    locus, stream = stream_locus_fragments(args.atac_url, args.chrom, args.start, args.end)
    joined, summary = build_join(barcodes, expr, locus)
    summary.update({
        "dataset": "GSE226277",
        "pair_id": args.pair_id,
        "assembly": "mm10",
        "locus": {"chrom": args.chrom, "start": args.start, "end": args.end},
        "stream": stream,
        "claim_scope": (
            "same-nucleus wild-type Multiome feasibility only; not allele-specific accessibility, "
            "TF occupancy, enhancer direction, target-gene causality, or mediator identification"
        ),
    })
    if FORBIDDEN_CAUSAL_FIELDS & set(summary):
        raise RuntimeError("summary leaked causal-promotion fields")

    fields = ["barcode", "fragment_records", "read_support"] + [f"{x.lower()}_umi" for x in GENES]
    with (args.output_dir / "locus_barcodes.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(joined)
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
