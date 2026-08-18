#!/usr/bin/env python3
"""Map the published mouse ortholog sequence of human rs1421085.

The exact coordinate is resolved by scanning a full chromosome FASTA for the
published WT CRISPR guide (or its reverse complement). A separately lifted
published fosmid interval is used only as a provenance cross-check, never as the
source of truth for the SNP coordinate.

This utility maps sequence coordinates only. It makes no causal claim.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

WT_GUIDE = "TAATCAATACGATGCCTT"
WT_GUIDE_RC = "AAGGCATCGTATTGATTA"
TARGET_OFFSET = {WT_GUIDE: 8, WT_GUIDE_RC: 9}
ALLELES = {WT_GUIDE: ("A", "G"), WT_GUIDE_RC: ("T", "C")}


def read_single_bed(path: Path) -> tuple[str, int, int]:
    rows = [
        line.strip()
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.startswith("#")
    ]
    if len(rows) != 1:
        raise ValueError(f"expected exactly one mapped BED row, got {len(rows)}")
    fields = rows[0].split("\t")
    if len(fields) < 3:
        raise ValueError("BED row must contain chrom, start, end")
    chrom, start, end = fields[0], int(fields[1]), int(fields[2])
    if start < 0 or end <= start:
        raise ValueError("invalid BED interval")
    return chrom, start, end


def read_fasta(path: Path) -> tuple[str, str]:
    header: str | None = None
    sequence_parts: list[str] = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                raise ValueError("expected a single-record chromosome FASTA")
            header = line[1:].split()[0]
        else:
            sequence_parts.append(line.upper())
    if header is None or not sequence_parts:
        raise ValueError("invalid FASTA")
    return header, "".join(sequence_parts)


def find_all(sequence: str, needle: str) -> list[int]:
    hits: list[int] = []
    offset = 0
    while True:
        index = sequence.find(needle, offset)
        if index < 0:
            return hits
        hits.append(index)
        offset = index + 1


def map_site(
    genome: str,
    chrom: str,
    sequence: str,
    lifted_interval: tuple[str, int, int] | None,
) -> dict:
    candidates: list[tuple[str, int]] = []
    for anchor in (WT_GUIDE, WT_GUIDE_RC):
        for index in find_all(sequence, anchor):
            candidates.append((anchor, index))
    if len(candidates) != 1:
        raise RuntimeError(
            "published WT guide must map exactly once on the chromosome; "
            f"found {len(candidates)} matches: {candidates[:20]}"
        )

    anchor, anchor_start0 = candidates[0]
    offset = TARGET_OFFSET[anchor]
    ref, alt = ALLELES[anchor]
    base0 = anchor_start0 + offset
    observed = sequence[base0]
    if observed != ref:
        raise RuntimeError(f"anchor matched but target base is {observed}, expected {ref}")

    crosscheck = None
    if lifted_interval is not None:
        lifted_chrom, lifted_start, lifted_end = lifted_interval
        crosscheck = {
            "lifted_fosmid_chrom": lifted_chrom,
            "lifted_fosmid_interval_0_based_half_open": [lifted_start, lifted_end],
            "same_chromosome": lifted_chrom == chrom,
            "exact_site_inside_lifted_fosmid": lifted_chrom == chrom and lifted_start <= base0 < lifted_end,
        }
        if lifted_chrom == chrom:
            if base0 < lifted_start:
                crosscheck["distance_to_lifted_fosmid_bp"] = lifted_start - base0
            elif base0 >= lifted_end:
                crosscheck["distance_to_lifted_fosmid_bp"] = base0 - lifted_end + 1
            else:
                crosscheck["distance_to_lifted_fosmid_bp"] = 0

    return {
        "genome": genome,
        "chrom": chrom,
        "position_0_based": base0,
        "position_1_based": base0 + 1,
        "reference_on_assembly_strand": ref,
        "alternate_on_assembly_strand": alt,
        "human_variant": "rs1421085 T>C",
        "anchor_start_0_based": anchor_start0,
        "anchor_orientation": (
            "published_guide"
            if anchor == WT_GUIDE
            else "reverse_complement_of_published_guide"
        ),
        "published_wt_guide": WT_GUIDE,
        "matched_sequence": anchor,
        "fosmid_crosscheck": crosscheck,
        "provenance": [
            "Laber et al. Science Advances 2021 PMID:34290091",
            "published guide: TAATCAATACGATGCCTT; PAM AGG",
            "published rs1421085-enhancer fosmid: mm9 chr8:93911541-93950923",
            "UCSC mm10 chromosome FASTA",
        ],
        "claim_scope": "coordinate mapping only; not biological evidence",
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", required=True, type=Path)
    parser.add_argument("--mapped-bed", type=Path)
    parser.add_argument("--genome", default="mm10")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    chrom, sequence = read_fasta(args.fasta)
    lifted = read_single_bed(args.mapped_bed) if args.mapped_bed else None
    result = map_site(args.genome, chrom, sequence, lifted)
    text = json.dumps(result, indent=2, sort_keys=True)
    print(text)
    if args.output:
        args.output.write_text(text + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
