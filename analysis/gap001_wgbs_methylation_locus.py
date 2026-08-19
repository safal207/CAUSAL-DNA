#!/usr/bin/env python3
"""Summarize adult hypothalamus WGBS around the exact rs1421085 mouse ortholog.

The input is a Bismark coverage file: chrom, start, end, methylation %, methylated,
unmethylated. Sequence context is classified from an mm9 reference window so
both forward- and reverse-strand cytosines can be recognized as CG, CAC, or
other contexts.

This is bulk wild-type adult-male hypothalamus context. It can show whether
methylated CG/CAC substrate exists near the ortholog, but cannot establish
cell-specific methylation, allele dependence, TF occupancy, or causal mediation.
"""
from __future__ import annotations

import argparse
import gzip
import json
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Call:
    pos0: int
    meth: int
    unmeth: int
    context: str

    @property
    def cov(self) -> int:
        return self.meth + self.unmeth

    @property
    def fraction(self) -> float:
        return self.meth / self.cov if self.cov else 0.0


def revcomp(seq: str) -> str:
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(table)[::-1]


def read_fasta(path: Path) -> tuple[str, str]:
    name = None
    chunks: list[str] = []
    for raw in path.read_text().splitlines():
        if raw.startswith(">"):
            if name is not None:
                raise ValueError("reference FASTA must contain exactly one record")
            name = raw[1:].split()[0]
        elif raw.strip():
            chunks.append(raw.strip().upper())
    if name is None or not chunks:
        raise ValueError("empty reference FASTA")
    return name, "".join(chunks)


def classify_context(seq: str, seq_start0: int, pos0: int) -> str:
    i = pos0 - seq_start0
    if i < 0 or i >= len(seq):
        return "OUTSIDE_REFERENCE_WINDOW"
    base = seq[i]
    if base == "C":
        tri = seq[i:i+3]
        if len(tri) >= 2 and tri[:2] == "CG":
            return "CG"
        if len(tri) == 3 and tri == "CAC":
            return "CAC"
        return "CH_OTHER"
    if base == "G":
        tri = seq[max(0, i-2):i+1]
        if len(tri) != 3:
            return "REVERSE_C_OTHER"
        rc = revcomp(tri)
        if rc[:2] == "CG":
            return "CG"
        if rc == "CAC":
            return "CAC"
        return "CH_OTHER"
    return "NON_CYTOSINE_REFERENCE"


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8", errors="replace") if path.suffix == ".gz" else path.open("r")


def scan_cov(path: Path, chrom: str, center0: int, radius: int, seq: str, seq_start0: int, min_cov: int) -> list[Call]:
    low1 = center0 - radius + 1
    high1 = center0 + radius + 1
    calls: list[Call] = []
    with open_text(path) as fh:
        for raw in fh:
            if not raw.strip() or raw.startswith("#"):
                continue
            f = raw.split()
            if len(f) < 6 or f[0] != chrom:
                continue
            try:
                start1 = int(f[1])
                meth = int(float(f[4]))
                unmeth = int(float(f[5]))
            except ValueError:
                continue
            if start1 < low1 or start1 > high1:
                continue
            pos0 = start1 - 1
            cov = meth + unmeth
            if cov < min_cov:
                continue
            calls.append(Call(pos0, meth, unmeth, classify_context(seq, seq_start0, pos0)))
    return calls


def aggregate(calls: list[Call], center0: int, radius: int, context: str | None = None) -> dict:
    subset = [c for c in calls if abs(c.pos0-center0) <= radius and (context is None or c.context == context)]
    meth = sum(c.meth for c in subset)
    total = sum(c.cov for c in subset)
    nearest = None
    if subset:
        c = min(subset, key=lambda x: (abs(x.pos0-center0), x.pos0))
        nearest = {
            "position_0_based": c.pos0,
            "distance_bp": c.pos0-center0,
            "coverage": c.cov,
            "methylation_fraction": c.fraction,
            "context": c.context,
        }
    return {
        "sites": len(subset),
        "methylated_reads": meth,
        "total_reads": total,
        "coverage_weighted_methylation": (meth/total if total else None),
        "nearest_site": nearest,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cov", type=Path, required=True)
    ap.add_argument("--fasta", type=Path, required=True)
    ap.add_argument("--fasta-start0", type=int, required=True)
    ap.add_argument("--chrom", default="chr8")
    ap.add_argument("--position0", type=int, required=True)
    ap.add_argument("--radius", type=int, default=5000)
    ap.add_argument("--min-cov", type=int, default=5)
    ap.add_argument("--sample-id", required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    _, seq = read_fasta(args.fasta)
    calls = scan_cov(args.cov, args.chrom, args.position0, args.radius, seq, args.fasta_start0, args.min_cov)
    if not calls:
        raise RuntimeError("no qualifying WGBS calls found in target window")

    contexts: dict[str, int] = {}
    for c in calls:
        contexts[c.context] = contexts.get(c.context, 0) + 1

    windows = {}
    for r in (250, 500, 2000, 5000):
        windows[str(r)] = {
            "all_calls": aggregate(calls, args.position0, r),
            "CG": aggregate(calls, args.position0, r, "CG"),
            "CAC": aggregate(calls, args.position0, r, "CAC"),
        }

    payload = {
        "dataset": "GSE84533",
        "sample_id": args.sample_id,
        "assay": "adult male mouse hypothalamus WGBS",
        "assembly": "mm9",
        "query": {"chrom": args.chrom, "position_0_based": args.position0, "position_1_based": args.position0+1},
        "minimum_coverage": args.min_cov,
        "window_radius_bp": args.radius,
        "qualifying_calls": len(calls),
        "context_counts": contexts,
        "windows": windows,
        "exact_variant_base_note": "The wild-type ortholog base is T, so it is not itself a methylatable cytosine; methylation context is evaluated at nearby cytosines.",
        "interpretation_guard": (
            "Bulk wild-type adult male hypothalamus WGBS can establish local methylated-CG/CAC substrate only. "
            "It cannot establish cell-specific methylation, allele dependence, MeCP2/TET1/DNMT1 occupancy, enhancer activity, target-gene assignment, or mediation."
        ),
        "discovery_level_ceiling": "D1_METHYLATION_CONTEXT",
        "cause_found": False,
    }
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
