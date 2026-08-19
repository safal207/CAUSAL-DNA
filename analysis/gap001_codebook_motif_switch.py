#!/usr/bin/env python3
"""Allele-aware Codebook motif scan for GAP-001 / rs1421085.

Only motif placements that physically overlap the variant are compared, at the
same placement and strand for T and C. Multiple published motif representatives
for the same TF are scanned independently and retain a stable motif_record_id.
PWM changes are sequence-model signals, not evidence of in-vivo occupancy,
target-gene regulation, or causality. VPH expression is a feasibility filter.
"""
from __future__ import annotations

import argparse
import csv
import json
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Iterable

import h5py
import numpy as np
from scipy.sparse import csc_matrix

BASES = "ACGT"
BASE_INDEX = {b: i for i, b in enumerate(BASES)}
RC = str.maketrans("ACGT", "TGCA")
SAMPLES = {"AJ18003": "male", "AJ18004": "female", "AJ19001": "female", "AJ19002": "male"}


@dataclass(frozen=True)
class MotifSwitch:
    tf: str
    motif_length: int
    strand: str
    start_offset: int
    ref_score: float
    alt_score: float
    ref_normalized: float
    alt_normalized: float
    delta_normalized: float
    abs_delta_normalized: float
    max_site_strength: float
    direction: str
    threshold_075_crossing: bool
    ref_window: str
    alt_window: str
    motif_record_id: str = ""


def reverse_complement(seq: str) -> str:
    return seq.translate(RC)[::-1]


def load_pwm(path: Path) -> tuple[str, np.ndarray]:
    lines = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    if not lines or not lines[0].startswith(">"):
        raise ValueError(f"invalid PWM header: {path}")
    tf = lines[0][1:].strip() or path.stem
    pwm = np.asarray([[float(x) for x in line.split()] for line in lines[1:]], dtype=float)
    if pwm.ndim != 2 or pwm.shape[1] != 4 or pwm.shape[0] < 2:
        raise ValueError(f"expected Lx4 PWM: {path}, got {pwm.shape}")
    return tf, pwm


def score_window(window: str, pwm: np.ndarray) -> float:
    if len(window) != len(pwm):
        raise ValueError("window/PWM length mismatch")
    return float(sum(pwm[i, BASE_INDEX[b]] for i, b in enumerate(window)))


def normalize_score(score: float, pwm: np.ndarray) -> float:
    lo = float(np.min(pwm, axis=1).sum())
    hi = float(np.max(pwm, axis=1).sum())
    return 0.0 if hi <= lo else (score - lo) / (hi - lo)


def best_variant_overlapping_switch(
    tf: str,
    pwm: np.ndarray,
    ref_seq: str,
    alt_seq: str,
    variant_index: int,
    motif_record_id: str = "",
) -> MotifSwitch:
    length = len(pwm)
    candidates: list[MotifSwitch] = []
    min_start = max(0, variant_index - length + 1)
    max_start = min(variant_index, len(ref_seq) - length)
    for start in range(min_start, max_start + 1):
        ref_forward = ref_seq[start:start + length]
        alt_forward = alt_seq[start:start + length]
        for strand, rw, aw in (("+", ref_forward, alt_forward), ("-", reverse_complement(ref_forward), reverse_complement(alt_forward))):
            rs, ac = score_window(rw, pwm), score_window(aw, pwm)
            rn, an = normalize_score(rs, pwm), normalize_score(ac, pwm)
            delta = an - rn
            direction = "gain_on_C" if delta > 0 else "loss_on_C" if delta < 0 else "no_change"
            candidates.append(MotifSwitch(
                tf=tf, motif_length=length, strand=strand, start_offset=start - variant_index,
                ref_score=rs, alt_score=ac, ref_normalized=rn, alt_normalized=an,
                delta_normalized=delta, abs_delta_normalized=abs(delta),
                max_site_strength=max(rn, an), direction=direction,
                threshold_075_crossing=(rn >= 0.75) != (an >= 0.75),
                ref_window=rw, alt_window=aw, motif_record_id=motif_record_id,
            ))
    if not candidates:
        raise ValueError(f"motif {tf} cannot overlap variant in supplied sequence")
    return max(candidates, key=lambda x: (x.abs_delta_normalized, x.max_site_strength))


def scan_motifs(motif_dir: Path, ref_seq: str, alt_seq: str, variant_index: int) -> list[MotifSwitch]:
    results = []
    for path in sorted(motif_dir.glob("*.pwm")):
        tf, pwm = load_pwm(path)
        results.append(best_variant_overlapping_switch(tf, pwm, ref_seq, alt_seq, variant_index, path.stem))
    if not results:
        raise ValueError(f"no .pwm files found under {motif_dir}")
    return sorted(results, key=lambda x: (x.abs_delta_normalized, x.max_site_strength), reverse=True)


def mouse_symbol_candidates(human_tf: str) -> list[str]:
    if not human_tf:
        return []
    out = [human_tf[0].upper() + human_tf[1:].lower()]
    if human_tf.startswith("ZNF") and human_tf[3:].isdigit():
        out.append("Zfp" + human_tf[3:])
    return list(dict.fromkeys(out))


def decode(values: Iterable[object]) -> np.ndarray:
    return np.asarray([x.decode() if isinstance(x, (bytes, np.bytes_)) else str(x) for x in values])


def load_candidate_presence(path: Path, human_tfs: list[str]) -> tuple[dict[str, np.ndarray], dict[str, str | None], int]:
    with h5py.File(path, "r") as h5:
        g = h5["matrix"]
        names = decode(g["features"]["name"][:])
        shape = tuple(int(x) for x in g["shape"][:])
        matrix = csc_matrix((g["data"][:], g["indices"][:], g["indptr"][:]), shape=shape)
        vectors, resolved = {}, {}
        for tf in human_tfs:
            match = None
            for symbol in mouse_symbol_candidates(tf):
                hits = np.flatnonzero(names == symbol)
                if len(hits) == 1:
                    match = (symbol, int(hits[0]))
                    break
            if match is None:
                vectors[tf] = np.zeros(shape[1], dtype=bool)
                resolved[tf] = None
            else:
                symbol, idx = match
                vectors[tf] = np.asarray(matrix[idx, :].todense()).ravel() > 0
                resolved[tf] = symbol
        return vectors, resolved, shape[1]


def expression_summary(input_root: Path, tfs: list[str]) -> list[dict[str, object]]:
    tfs = list(dict.fromkeys(tfs))
    aggregate = {tf: {"cells": 0, "positive": 0, "male_cells": 0, "male_positive": 0, "female_cells": 0, "female_positive": 0, "symbols": set()} for tf in tfs}
    for sample, sex in SAMPLES.items():
        vectors, resolved, n_cells = load_candidate_presence(input_root / sample / "filtered_feature_bc_matrix.h5", tfs)
        for tf in tfs:
            rec, pos = aggregate[tf], int(vectors[tf].sum())
            rec["cells"] += n_cells
            rec["positive"] += pos
            rec[f"{sex}_cells"] += n_cells
            rec[f"{sex}_positive"] += pos
            if resolved[tf]:
                rec["symbols"].add(resolved[tf])
    rows = []
    for tf in tfs:
        rec = aggregate[tf]
        rows.append({
            "tf": tf,
            "mouse_symbols": ";".join(sorted(rec["symbols"])),
            "vph_cells": rec["cells"], "vph_positive": rec["positive"],
            "vph_positive_pct": 100.0 * rec["positive"] / rec["cells"] if rec["cells"] else 0.0,
            "male_positive_pct": 100.0 * rec["male_positive"] / rec["male_cells"] if rec["male_cells"] else 0.0,
            "female_positive_pct": 100.0 * rec["female_positive"] / rec["female_cells"] if rec["female_cells"] else 0.0,
        })
    return rows


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if rows:
        with path.open("w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--motif-dir", type=Path, required=True)
    p.add_argument("--sequence-json", type=Path, required=True)
    p.add_argument("--variant-index", type=int, required=True)
    p.add_argument("--ref", default="T")
    p.add_argument("--alt", default="C")
    p.add_argument("--expression-root", type=Path)
    p.add_argument("--top-expression", type=int, default=30)
    p.add_argument("--output-dir", type=Path, required=True)
    args = p.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    ref_seq = str(json.loads(args.sequence_json.read_text()).get("dna", "")).upper()
    if not ref_seq or any(b not in BASES for b in ref_seq):
        raise ValueError("UCSC sequence JSON missing clean A/C/G/T dna")
    if ref_seq[args.variant_index] != args.ref:
        raise ValueError(f"reference base mismatch at variant: expected {args.ref}, got {ref_seq[args.variant_index]}")
    alt_seq = ref_seq[:args.variant_index] + args.alt + ref_seq[args.variant_index + 1:]

    switches = scan_motifs(args.motif_dir, ref_seq, alt_seq, args.variant_index)
    write_csv(args.output_dir / "motif_switches.csv", [asdict(x) for x in switches])
    strong = [x for x in switches if x.max_site_strength >= 0.75]
    crossings = [x for x in switches if x.threshold_075_crossing]
    unique_tfs = len({x.tf for x in switches})

    selected_tfs = list(dict.fromkeys(
        [x.tf for x in switches[:args.top_expression]] + [x.tf for x in crossings] + [x.tf for x in strong]
    ))
    expr_rows, expr_by_tf = [], {}
    if args.expression_root:
        expr_rows = expression_summary(args.expression_root, selected_tfs)
        write_csv(args.output_dir / "vph_expression_feasibility.csv", expr_rows)
        expr_by_tf = {str(x["tf"]): x for x in expr_rows}

    rank_by_record = {(x.tf, x.motif_record_id): i + 1 for i, x in enumerate(switches)}
    integrated = []
    for sw in switches[:args.top_expression]:
        row = {"rank": rank_by_record[(sw.tf, sw.motif_record_id)], **asdict(sw)}
        row.update({f"expr_{k}": v for k, v in expr_by_tf.get(sw.tf, {}).items() if k != "tf"})
        integrated.append(row)
    write_csv(args.output_dir / "top_integrated_candidates.csv", integrated)

    crossing_rows = []
    for sw in crossings:
        row = {"rank": rank_by_record[(sw.tf, sw.motif_record_id)], **asdict(sw)}
        row.update({f"expr_{k}": v for k, v in expr_by_tf.get(sw.tf, {}).items() if k != "tf"})
        crossing_rows.append(row)
    write_csv(args.output_dir / "threshold_crossings.csv", crossing_rows)

    report = [
        "# GAP-001 Codebook allele motif-switch scan", "",
        "**Status: sequence-model discovery screen; not occupancy or causality.**", "",
        f"- Motif records scanned: **{len(switches)}**",
        f"- Unique TFs covered: **{unique_tfs}**",
        f"- Exact edit: **{args.ref}>{args.alt}** at supplied variant index {args.variant_index}",
        f"- Motif records with >=0.75 normalized site strength on either allele: **{len(strong)}**",
        f"- Heuristic 0.75 threshold-crossing motif records: **{len(crossings)}**", "",
        "## Heuristic 0.75 threshold crossings", "",
        "|rank|TF|motif record|direction|ref|alt|abs Δnorm|VPH +%|male +%|female +%|",
        "|---:|---|---|---|---:|---:|---:|---:|---:|---:|",
    ]
    for row in crossing_rows:
        report.append(
            f"|{row['rank']}|{row['tf']}|{row['motif_record_id']}|{row['direction']}|"
            f"{row['ref_normalized']:.4f}|{row['alt_normalized']:.4f}|{row['abs_delta_normalized']:.4f}|"
            f"{float(row.get('expr_vph_positive_pct',0)):.3f}|{float(row.get('expr_male_positive_pct',0)):.3f}|"
            f"{float(row.get('expr_female_positive_pct',0)):.3f}|"
        )
    report += [
        "", "## Top allele-sensitive motif records", "",
        "|rank|TF|motif record|direction|abs Δnorm|max strength|strand|offset|VPH +%|mouse symbol|",
        "|---:|---|---|---|---:|---:|:---:|---:|---:|---|",
    ]
    for row in integrated[:20]:
        report.append(
            f"|{row['rank']}|{row['tf']}|{row['motif_record_id']}|{row['direction']}|"
            f"{row['abs_delta_normalized']:.4f}|{row['max_site_strength']:.4f}|{row['strand']}|{row['start_offset']}|"
            f"{float(row.get('expr_vph_positive_pct',0)):.3f}|{row.get('expr_mouse_symbols','')}|"
        )
    report += [
        "", "## Interpretation boundary", "",
        "- PWM score change means the DNA substitution changes a sequence model for TF recognition.",
        "- Multiple representatives of one TF are independent published motif records, not independent biological replications.",
        "- VPH transcript detection is only a feasibility prerequisite; absence can demote a direct adult-VPH interpreter, presence cannot establish binding.",
        "- The 0.75 crossing is a declared heuristic for triage, not a biochemical affinity threshold.",
        "- No candidate is promoted to causal status without allele-specific occupancy and perturbation in the relevant PH cell state.",
    ]
    text = "\n".join(report) + "\n"
    (args.output_dir / "report.md").write_text(text)
    print(text)


if __name__ == "__main__":
    main()
