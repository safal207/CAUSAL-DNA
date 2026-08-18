#!/usr/bin/env python3
"""Map the published mouse ortholog sequence of human rs1421085.

Input is a BED interval already lifted to the requested assembly. The script
fetches that interval sequence from the UCSC Genome Browser API and requires a
unique exact match to the published WT CRISPR guide sequence (or its reverse
complement). It then emits the exact genomic base corresponding to the human
T>C site.

This is a coordinate/provenance utility only. It makes no causal claim.
"""
from __future__ import annotations

import argparse
import json
import urllib.parse
import urllib.request
from pathlib import Path

WT_GUIDE = "TAATCAATACGATGCCTT"
WT_GUIDE_RC = "AAGGCATCGTATTGATTA"
# The homologous T>C site is represented as A>G on WT_GUIDE orientation and
# T>C on its reverse-complement orientation.
TARGET_OFFSET = {WT_GUIDE: 8, WT_GUIDE_RC: 9}
ALLELES = {WT_GUIDE: ("A", "G"), WT_GUIDE_RC: ("T", "C")}


def read_single_bed(path: Path) -> tuple[str, int, int]:
    rows = [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip() and not line.startswith("#")]
    if len(rows) != 1:
        raise ValueError(f"expected exactly one mapped BED row, got {len(rows)}")
    fields = rows[0].split("\t")
    if len(fields) < 3:
        raise ValueError("BED row must contain chrom, start, end")
    chrom, start, end = fields[0], int(fields[1]), int(fields[2])
    if start < 0 or end <= start:
        raise ValueError("invalid BED interval")
    return chrom, start, end


def fetch_sequence(genome: str, chrom: str, start: int, end: int) -> str:
    params = urllib.parse.urlencode({"genome": genome, "chrom": chrom, "start": start, "end": end})
    url = f"https://api.genome.ucsc.edu/getData/sequence?{params}"
    with urllib.request.urlopen(url, timeout=60) as response:
        payload = json.load(response)
    dna = payload.get("dna")
    if not isinstance(dna, str) or not dna:
        raise RuntimeError(f"UCSC sequence API returned no DNA for {chrom}:{start}-{end}")
    return dna.upper()


def find_all(sequence: str, needle: str) -> list[int]:
    hits: list[int] = []
    offset = 0
    while True:
        index = sequence.find(needle, offset)
        if index < 0:
            return hits
        hits.append(index)
        offset = index + 1


def map_site(genome: str, chrom: str, start: int, end: int, sequence: str) -> dict:
    candidates: list[tuple[str, int]] = []
    for anchor in (WT_GUIDE, WT_GUIDE_RC):
        for index in find_all(sequence, anchor):
            candidates.append((anchor, index))
    if len(candidates) != 1:
        raise RuntimeError(
            "published WT guide must map exactly once inside lifted interval; "
            f"found {len(candidates)} matches: {candidates}"
        )
    anchor, local_start = candidates[0]
    offset = TARGET_OFFSET[anchor]
    ref, alt = ALLELES[anchor]
    base0 = start + local_start + offset
    observed = sequence[local_start + offset]
    if observed != ref:
        raise RuntimeError(f"anchor matched but target base is {observed}, expected {ref}")
    return {
        "genome": genome,
        "chrom": chrom,
        "position_0_based": base0,
        "position_1_based": base0 + 1,
        "reference_on_assembly_strand": ref,
        "alternate_on_assembly_strand": alt,
        "human_variant": "rs1421085 T>C",
        "anchor_orientation": "published_guide" if anchor == WT_GUIDE else "reverse_complement_of_published_guide",
        "published_wt_guide": WT_GUIDE,
        "matched_sequence": anchor,
        "lifted_search_interval_0_based_half_open": [start, end],
        "provenance": [
            "Laber et al. Science Advances 2021 PMID:34290091",
            "published guide: TAATCAATACGATGCCTT; PAM AGG",
            "published rs1421085-enhancer fosmid: mm9 chr8:93911541-93950923",
        ],
        "claim_scope": "coordinate mapping only; not biological evidence",
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mapped-bed", required=True, type=Path)
    parser.add_argument("--genome", default="mm10")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    chrom, start, end = read_single_bed(args.mapped_bed)
    sequence = fetch_sequence(args.genome, chrom, start, end)
    result = map_site(args.genome, chrom, start, end, sequence)
    text = json.dumps(result, indent=2, sort_keys=True)
    print(text)
    if args.output:
        args.output.write_text(text + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
