#!/usr/bin/env python3
"""Normalize Nature Codebook Supplementary Data 1 into Lx4 PWMs.

The published archive covers 1,421 unique human TFs but can contain more than
one representative motif record for a TF.  This normalizer therefore preserves
*all* published motif representatives and makes the fail-closed atlas contract
about unique TF coverage, not a false one-TF/one-PWM assumption.

Supported encodings are intentionally explicit: the official Codebook TSV
format (TF/Motif/Pos+A+C+G+T), simple Lx4 matrices, MEME motif blocks,
JASPAR-style A/C/G/T row matrices, and TRANSFAC-like blocks.

Matrix values are preserved. Downstream allele scanning normalizes each motif
against its own theoretical row-wise min/max. Multiple representatives for the
same TF are scanned independently and only then collapsed at the TF decision
layer; normalization itself makes no occupancy or causality claim.
"""
from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path
import re
import zipfile

BASES = ("A", "C", "G", "T")
FLOAT_RE = re.compile(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?")


def floats(line: str) -> list[float]:
    return [float(x) for x in FLOAT_RE.findall(line)]


def sanitize_tf(name: str) -> str:
    name = name.strip().lstrip(">").strip()
    tokens = [t for t in re.split(r"[\s|]+", name) if t]
    if len(tokens) > 1 and re.fullmatch(r"(?:MA|M)\d+(?:\.\d+)?", tokens[0], re.I):
        name = tokens[-1]
    else:
        name = tokens[0] if tokens else ""
    name = re.sub(r"[^A-Za-z0-9_.-]", "_", name)
    if not name:
        raise ValueError("empty TF name")
    return name


def parse_codebook_tsv(text: str) -> list[tuple[str, list[list[float]]]]:
    """Parse official SupplementaryData1 files.

    Example:
      TF\tNOBOX
      Motif\tM01232_1.94d
      Pos\tA\tC\tG\tT
      1\t0.21\t0.30\t0.26\t0.22
    """
    lines = [line.rstrip("\r\n") for line in text.splitlines() if line.strip()]
    if len(lines) < 5:
        return []
    first = lines[0].split("\t")
    second = lines[1].split("\t")
    header = lines[2].split("\t")
    if len(first) < 2 or first[0].strip().lower() != "tf":
        return []
    if len(second) < 2 or second[0].strip().lower() != "motif":
        return []
    if [x.strip().upper() for x in header[:5]] != ["POS", "A", "C", "G", "T"]:
        return []
    tf = sanitize_tf(first[1])
    rows: list[list[float]] = []
    expected_pos = 1
    for raw in lines[3:]:
        parts = raw.split("\t")
        if len(parts) < 5:
            return []
        try:
            pos = int(float(parts[0]))
            vals = [float(parts[i]) for i in range(1, 5)]
        except ValueError:
            return []
        if pos != expected_pos:
            return []
        rows.append(vals)
        expected_pos += 1
    return [(tf, rows)] if len(rows) >= 2 else []


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
        lambda: parse_codebook_tsv(text),
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
    records: list[tuple[str, list[list[float]], str]] = []
    unparsed: list[str] = []
    source_member_count = 0

    with zipfile.ZipFile(args.input_zip) as zf:
        for info in zf.infolist():
            if info.is_dir():
                continue
            source_member_count += 1
            motifs = parse_member(info.filename, zf.read(info))
            if not motifs:
                unparsed.append(info.filename)
                continue
            for tf, rows in motifs:
                if len(rows) >= 2 and all(len(row) == 4 for row in rows):
                    records.append((tf, rows, info.filename))

    by_tf: dict[str, list[tuple[list[list[float]], str]]] = defaultdict(list)
    for tf, rows, source in records:
        by_tf[tf].append((rows, source))

    output_records: list[dict[str, object]] = []
    for tf in sorted(by_tf):
        representatives = by_tf[tf]
        for index, (rows, source) in enumerate(representatives, start=1):
            if len(representatives) == 1:
                filename = f"{tf}.pwm"
            else:
                filename = f"{tf}__rep{index:02d}.pwm"
            write_pwm(args.output_dir / filename, tf, rows)
            output_records.append({
                "tf": tf,
                "representative_index": index,
                "representative_count_for_tf": len(representatives),
                "source_member": source,
                "output_pwm": filename,
                "motif_length": len(rows),
            })

    unique_tf_count = len(by_tf)
    motif_record_count = len(records)
    multi_tf = {tf: len(reps) for tf, reps in by_tf.items() if len(reps) > 1}
    manifest = {
        "expected_unique_tfs": args.expected_tfs,
        "parsed_unique_tfs": unique_tf_count,
        "parsed_motif_records": motif_record_count,
        "source_member_count": source_member_count,
        "unparsed_member_count": len(unparsed),
        "unparsed_members": unparsed,
        "multiple_representative_tf_count": len(multi_tf),
        "max_representatives_per_tf": max((len(x) for x in by_tf.values()), default=0),
        "multiple_representative_counts": dict(sorted(multi_tf.items())),
        "motif_records": output_records,
        "full_atlas_contract_pass": unique_tf_count == args.expected_tfs,
        "contract_semantics": "1,421 unique TFs; all published representative motif records preserved",
        "claim_scope": "motif normalization only; not TF occupancy or causality",
    }
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    args.manifest.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({k: v for k, v in manifest.items() if k not in {"motif_records", "unparsed_members", "multiple_representative_counts"}}, indent=2))

    if unique_tf_count != args.expected_tfs:
        raise RuntimeError(
            f"full-atlas contract failed: expected {args.expected_tfs} unique TFs, parsed {unique_tf_count} "
            f"from {motif_record_count} motif records"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
