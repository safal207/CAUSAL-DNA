#!/usr/bin/env python3
"""Normalize Nature Codebook Supplementary Data 1 into one Lx4 PWM per TF.

Supported text encodings are intentionally explicit: simple Lx4 matrices,
MEME motif blocks, JASPAR-style A/C/G/T row matrices, and TRANSFAC-like blocks.
The conversion is fail-closed for the published full atlas: the number of
unique parsed TFs must equal --expected-tfs (1,421 in CDNA-001 Pass B).

Matrix values are preserved. Downstream allele scanning normalizes each motif
against its own theoretical row-wise min/max, so this conversion does not claim
that scores are comparable between different PWM construction methods.
"""
from __future__ import annotations

import argparse
import csv
import io
import json
from pathlib import Path
import re
import zipfile

BASES = ("A", "C", "G", "T")
FLOAT_RE = re.compile(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?")


def floats(line: str) -> list[float]:
    return [float(x) for x in FLOAT_RE.findall(line)]


def sanitize_tf(name: str) -> str:
    name = name.strip().lstrip(">").strip()
    # Prefer the final gene-like token if an accession precedes it.
    tokens = [t for t in re.split(r"[\s|]+", name) if t]
    if len(tokens) > 1 and re.fullmatch(r"(?:MA|M)\d+(?:\.\d+)?", tokens[0], re.I):
        name = tokens[-1]
    else:
        name = tokens[0] if tokens else ""
    name = re.sub(r"[^A-Za-z0-9_.-]", "_", name)
    if not name:
        raise ValueError("empty TF name")
    return name


def parse_simple_lx4(text: str, fallback: str) -> list[tuple[str, list[list[float]]]]:
    lines = [x.strip() for x in text.splitlines() if x.strip() and not x.lstrip().startswith("#")]
    if not lines:
        return []
    name = fallback
    start = 0
    if lines[0].startswith(">"):
        name = sanitize_tf(lines[0])
        start = 1
    rows: list[list[float]] = []
    for line in lines[start:]:
        vals = floats(line)
        if len(vals) != 4:
            return []
        rows.append(vals)
    return [(name, rows)] if len(rows) >= 2 else []


def parse_jaspar_rows(text: str, fallback: str) -> list[tuple[str, list[list[float]]]]:
    lines = [x.strip() for x in text.splitlines() if x.strip()]
    if not lines:
        return []
    name = sanitize_tf(lines[0]) if lines[0].startswith(">") else fallback
    by_base: dict[str, list[float]] = {}
    for line in lines[1:] if lines[0].startswith(">") else lines:
        m = re.match(r"^([ACGT])\s*(?:\[)?(.*?)(?:\])?\s*$", line, re.I)
        if not m:
            continue
        by_base[m.group(1).upper()] = floats(m.group(2))
    if set(by_base) != set(BASES):
        return []
    lengths = {len(by_base[b]) for b in BASES}
    if len(lengths) != 1 or next(iter(lengths)) < 2:
        return []
    length = next(iter(lengths))
    return [(name, [[by_base[b][i] for b in BASES] for i in range(length)])]


def parse_meme(text: str) -> list[tuple[str, list[list[float]]]]:
    lines = text.splitlines()
    out: list[tuple[str, list[list[float]]]] = []
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if not line.startswith("MOTIF "):
            i += 1
            continue
        parts = line.split(maxsplit=2)
        name = sanitize_tf(parts[2] if len(parts) > 2 else parts[1])
        i += 1
        while i < len(lines) and "letter-probability matrix" not in lines[i].lower():
            if lines[i].strip().startswith("MOTIF "):
                break
            i += 1
        if i >= len(lines) or "letter-probability matrix" not in lines[i].lower():
            continue
        i += 1
        rows: list[list[float]] = []
        while i < len(lines):
            vals = floats(lines[i])
            if len(vals) != 4:
                break
            rows.append(vals)
            i += 1
        if len(rows) >= 2:
            out.append((name, rows))
    return out


def parse_transfac(text: str, fallback: str) -> list[tuple[str, list[list[float]]]]:
    lines = text.splitlines()
    out: list[tuple[str, list[list[float]]]] = []
    name = fallback
    rows: list[list[float]] = []
    in_matrix = False
    for raw in lines + ["//"]:
        line = raw.strip()
        if line.startswith("ID") or line.startswith("NA"):
            parts = line.split(maxsplit=1)
            if len(parts) == 2:
                name = sanitize_tf(parts[1])
        elif line.startswith("P0") or line.startswith("PO"):
            in_matrix = True
        elif line == "//":
            if len(rows) >= 2:
                out.append((name, rows))
            rows = []
            in_matrix = False
            name = fallback
        elif in_matrix:
            vals = floats(line)
            # Numbered rows usually have position first, then A,C,G,T.
            if len(vals) >= 5:
                rows.append(vals[1:5])
            elif len(vals) == 4:
                rows.append(vals)
    return out


def parse_member(name: str, data: bytes) -> list[tuple[str, list[list[float]]]]:
    if len(data) > 20_000_000:
        return []
    try:
        text = data.decode("utf-8-sig")
    except UnicodeDecodeError:
        try:
            text = data.decode("latin-1")
        except UnicodeDecodeError:
            return []
    fallback = sanitize_tf(Path(name).stem)
    parsers = (
        lambda: parse_meme(text),
        lambda: parse_jaspar_rows(text, fallback),
        lambda: parse_transfac(text, fallback),
        lambda: parse_simple_lx4(text, fallback),
    )
    for parser in parsers:
        motifs = parser()
        if motifs:
            return motifs
    return []


def matrix_signature(rows: list[list[float]]) -> tuple[tuple[float, ...], ...]:
    return tuple(tuple(round(v, 12) for v in row) for row in rows)


def write_pwm(path: Path, tf: str, rows: list[list[float]]) -> None:
    path.write_text(
        ">" + tf + "\n" + "\n".join("\t".join(f"{v:.12g}" for v in row) for row in rows) + "\n",
        encoding="utf-8",
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-zip", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--expected-tfs", type=int, default=1421)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    parsed: dict[str, tuple[list[list[float]], str]] = {}
    duplicates: list[dict[str, str]] = []
    unparsed: list[str] = []
    source_member_count = 0

    with zipfile.ZipFile(args.input_zip) as zf:
        for info in zf.infolist():
            if info.is_dir():
                continue
            source_member_count += 1
            data = zf.read(info)
            motifs = parse_member(info.filename, data)
            if not motifs:
                unparsed.append(info.filename)
                continue
            for tf, rows in motifs:
                if len(rows) < 2 or any(len(row) != 4 for row in rows):
                    continue
                if tf in parsed:
                    same = matrix_signature(parsed[tf][0]) == matrix_signature(rows)
                    duplicates.append({"tf": tf, "first": parsed[tf][1], "second": info.filename, "identical": str(same).lower()})
                    if not same:
                        continue
                else:
                    parsed[tf] = (rows, info.filename)

    # Non-identical duplicate representatives make the one-TF/one-PWM contract ambiguous.
    conflicting = [x for x in duplicates if x["identical"] == "false"]
    for tf, (rows, _) in parsed.items():
        write_pwm(args.output_dir / f"{tf}.pwm", tf, rows)

    manifest = {
        "expected_unique_tfs": args.expected_tfs,
        "parsed_unique_tfs": len(parsed),
        "source_member_count": source_member_count,
        "unparsed_member_count": len(unparsed),
        "unparsed_members": unparsed,
        "duplicate_count": len(duplicates),
        "conflicting_duplicate_count": len(conflicting),
        "duplicates": duplicates,
        "tf_sources": {tf: src for tf, (_, src) in sorted(parsed.items())},
        "full_atlas_contract_pass": len(parsed) == args.expected_tfs and not conflicting,
        "claim_scope": "motif normalization only; not TF occupancy or causality",
    }
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    args.manifest.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({k: v for k, v in manifest.items() if k not in {"tf_sources", "duplicates", "unparsed_members"}}, indent=2))

    if conflicting:
        raise RuntimeError(f"ambiguous non-identical duplicate TF representatives: {len(conflicting)}")
    if len(parsed) != args.expected_tfs:
        raise RuntimeError(f"full-atlas contract failed: expected {args.expected_tfs} unique TFs, parsed {len(parsed)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
