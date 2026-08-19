#!/usr/bin/env python3
"""Build a provenance-first paired RNA/ATAC manifest for GSE226277 hypothalamus Multiome.

This step intentionally does not analyze biology. It verifies that RNA and ATAC
samples can be paired by tissue, condition, replicate and published sample
metadata before any same-nucleus join is attempted.
"""

from __future__ import annotations

import argparse
import gzip
import json
import re
from pathlib import Path
from urllib.parse import urlparse

TITLE_RE = re.compile(
    r"^Hypothalamus\s+(?P<condition>RH|AH),\s+replicate\s+(?P<replicate>[12]),\s+Multiome\s+\((?P<modality>RNA|ATAC)\)$"
)
EXPECTED = {(condition, replicate) for condition in ("RH", "AH") for replicate in (1, 2)}


def _https_ftp(url: str) -> str:
    if url.startswith("ftp://ftp.ncbi.nlm.nih.gov/"):
        return "https://ftp.ncbi.nlm.nih.gov/" + url.split("ftp://ftp.ncbi.nlm.nih.gov/", 1)[1]
    return url


def parse_soft(text: str) -> list[dict]:
    samples: list[dict] = []
    current: dict | None = None
    for raw in text.splitlines():
        line = raw.rstrip("\n")
        if line.startswith("^SAMPLE = "):
            if current is not None:
                samples.append(current)
            current = {
                "accession": line.split("=", 1)[1].strip(),
                "characteristics": {},
                "supplementary_files": [],
            }
            continue
        if current is None:
            continue
        if line.startswith("!Sample_title = "):
            current["title"] = line.split("=", 1)[1].strip()
        elif line.startswith("!Sample_characteristics_ch1 = "):
            value = line.split("=", 1)[1].strip()
            if ":" in value:
                key, val = value.split(":", 1)
                current["characteristics"][key.strip().lower()] = val.strip()
        elif line.startswith("!Sample_supplementary_file"):
            url = line.split("=", 1)[1].strip()
            if url and url.upper() != "NONE":
                current["supplementary_files"].append(_https_ftp(url))
    if current is not None:
        samples.append(current)
    return samples


def select_hypothalamus_multiome(samples: list[dict]) -> list[dict]:
    selected = []
    for sample in samples:
        title = sample.get("title", "")
        match = TITLE_RE.match(title)
        if not match:
            continue
        ch = sample.get("characteristics", {})
        tissue = ch.get("tissue", "")
        sex = ch.get("sex", "")
        genotype = ch.get("genotype", "")
        if tissue.lower() != "hypothalamus":
            raise ValueError(f"{sample['accession']}: expected hypothalamus, got {tissue!r}")
        if sex.lower() != "male":
            raise ValueError(f"{sample['accession']}: expected male, got {sex!r}")
        if genotype.upper() != "WT":
            raise ValueError(f"{sample['accession']}: expected WT, got {genotype!r}")
        selected.append(
            {
                "accession": sample["accession"],
                "title": title,
                "condition": match.group("condition"),
                "replicate": int(match.group("replicate")),
                "modality": match.group("modality"),
                "sex": sex,
                "genotype": genotype,
                "tissue": tissue,
                "supplementary_files": sample.get("supplementary_files", []),
            }
        )
    return selected


def pair_samples(selected: list[dict]) -> list[dict]:
    by_key: dict[tuple[str, int], dict[str, dict]] = {}
    for row in selected:
        key = (row["condition"], row["replicate"])
        by_key.setdefault(key, {})[row["modality"]] = row

    if set(by_key) != EXPECTED:
        raise ValueError(f"expected pair keys {sorted(EXPECTED)}, got {sorted(by_key)}")

    pairs = []
    for condition, replicate in sorted(EXPECTED):
        modalities = by_key[(condition, replicate)]
        if set(modalities) != {"RNA", "ATAC"}:
            raise ValueError(f"{condition} replicate {replicate}: missing RNA/ATAC pair")
        rna, atac = modalities["RNA"], modalities["ATAC"]
        if (rna["sex"], rna["genotype"], rna["tissue"]) != (atac["sex"], atac["genotype"], atac["tissue"]):
            raise ValueError(f"{condition} replicate {replicate}: metadata mismatch")
        pairs.append(
            {
                "pair_id": f"HYPO-{condition}-R{replicate}",
                "condition": condition,
                "replicate": replicate,
                "sex": rna["sex"],
                "genotype": rna["genotype"],
                "rna_accession": rna["accession"],
                "atac_accession": atac["accession"],
                "rna_files": rna["supplementary_files"],
                "atac_files": atac["supplementary_files"],
            }
        )
    return pairs


def read_soft(path: Path) -> str:
    if path.suffix == ".gz":
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
            return handle.read()
    return path.read_text(encoding="utf-8", errors="replace")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--soft", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    samples = parse_soft(read_soft(args.soft))
    selected = select_hypothalamus_multiome(samples)
    pairs = pair_samples(selected)
    payload = {
        "dataset": "GSE226277",
        "assembly": "mm10",
        "pair_count": len(pairs),
        "pairs": pairs,
        "claim_scope": "provenance/pairing only; no biological or causal claim",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
