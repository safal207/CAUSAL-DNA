#!/usr/bin/env python3
"""Cross-species rs1421085 motif-grammar stability for GAP-001.

This gate compares the exact human rs1421085 T>C sequence-model effect with
the exact mouse-ortholog T>C effect used in the 2025 knock-in study.

The same published motif_record_id is compared across species. For TF-level
prioritization, the most allele-sensitive HUMAN representative is frozen first
and that SAME record is evaluated in mouse. Strict grammar portability further
requires the best variant-overlapping placement to preserve strand and offset.

Important: for an additive PWM and a one-base edit, delta magnitude at a fixed
motif column is largely determined by that column's allele weights. Identical
human/mouse delta therefore is not independent replication by itself. Flanking
sequence is tested through placement and full-site strength.

Outputs are sequence-model robustness signals only. They do not establish
in-vivo occupancy, Irx3 regulation, or causality.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


def truthy(value: object) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes"}


def f(row: dict[str, str], key: str) -> float:
    try:
        return float(row.get(key, "") or 0.0)
    except ValueError:
        return 0.0


def i(row: dict[str, str], key: str) -> int:
    return int(str(row.get(key, "0") or 0))


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        raise ValueError(f"empty motif-switch CSV: {path}")
    required = {
        "tf", "motif_record_id", "direction", "delta_normalized",
        "abs_delta_normalized", "max_site_strength", "threshold_075_crossing",
        "strand", "start_offset",
    }
    missing = required.difference(rows[0])
    if missing:
        raise ValueError(f"missing columns in {path}: {sorted(missing)}")
    return rows


def by_record(rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    out: dict[str, dict[str, str]] = {}
    for row in rows:
        rid = row["motif_record_id"]
        if not rid:
            raise ValueError("empty motif_record_id")
        if rid in out:
            raise ValueError(f"duplicate motif_record_id: {rid}")
        out[rid] = row
    return out


def human_selected_by_tf(rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    grouped: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        grouped.setdefault(row["tf"], []).append(row)
    return {
        tf: max(
            reps,
            key=lambda r: (
                f(r, "abs_delta_normalized"),
                f(r, "max_site_strength"),
                r["motif_record_id"],
            ),
        )
        for tf, reps in grouped.items()
    }


def same_nonzero_direction(a: float, b: float) -> bool:
    return a != 0.0 and b != 0.0 and ((a > 0) == (b > 0))


def pearson(xs: list[float], ys: list[float]) -> float | None:
    if len(xs) < 2 or len(xs) != len(ys):
        return None
    mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
    dx, dy = [x - mx for x in xs], [y - my for y in ys]
    denom = math.sqrt(sum(x*x for x in dx) * sum(y*y for y in dy))
    return None if denom == 0 else sum(x*y for x, y in zip(dx, dy)) / denom


def compare_species(
    human_rows: list[dict[str, str]],
    mouse_rows: list[dict[str, str]],
    min_abs_delta: float = 0.10,
) -> tuple[list[dict[str, object]], dict[str, object]]:
    human_records = by_record(human_rows)
    mouse_records = by_record(mouse_rows)
    if set(human_records) != set(mouse_records):
        missing_mouse = sorted(set(human_records) - set(mouse_records))
        missing_human = sorted(set(mouse_records) - set(human_records))
        raise ValueError(
            "human/mouse motif_record_id sets differ: "
            f"missing_mouse={missing_mouse[:10]} missing_human={missing_human[:10]}"
        )

    selected = human_selected_by_tf(human_rows)
    tf_rows: list[dict[str, object]] = []
    signed_h, signed_m, abs_h, abs_m = [], [], [], []
    all_record_direction_ok = 0

    for rid, hr in human_records.items():
        mr = mouse_records[rid]
        if hr["tf"] != mr["tf"]:
            raise ValueError(f"TF mismatch for {rid}: {hr['tf']} vs {mr['tf']}")
        hd, md = f(hr, "delta_normalized"), f(mr, "delta_normalized")
        signed_h.append(hd); signed_m.append(md)
        abs_h.append(abs(hd)); abs_m.append(abs(md))
        if same_nonzero_direction(hd, md):
            all_record_direction_ok += 1

    for tf, hr in selected.items():
        rid = hr["motif_record_id"]
        mr = mouse_records[rid]
        hd, md = f(hr, "delta_normalized"), f(mr, "delta_normalized")
        direction_ok = same_nonzero_direction(hd, md)
        h_abs, m_abs = abs(hd), abs(md)
        h_site, m_site = f(hr, "max_site_strength"), f(mr, "max_site_strength")
        same_placement = hr["strand"] == mr["strand"] and i(hr, "start_offset") == i(mr, "start_offset")
        robust_delta = direction_ok and h_abs >= min_abs_delta and m_abs >= min_abs_delta
        strict_grammar = robust_delta and same_placement
        row = {
            "tf": tf,
            "motif_record_id": rid,
            "human_direction": hr["direction"],
            "mouse_direction": mr["direction"],
            "direction_conserved": direction_ok,
            "human_delta": hd,
            "mouse_delta": md,
            "human_abs_delta": h_abs,
            "mouse_abs_delta": m_abs,
            "min_abs_delta": min(h_abs, m_abs),
            "human_max_site": h_site,
            "mouse_max_site": m_site,
            "min_site_strength": min(h_site, m_site),
            "mouse_to_human_site_ratio": (m_site / h_site) if h_site > 0 else None,
            "human_strand": hr["strand"],
            "mouse_strand": mr["strand"],
            "human_start_offset": i(hr, "start_offset"),
            "mouse_start_offset": i(mr, "start_offset"),
            "placement_conserved": same_placement,
            "human_crossing": truthy(hr["threshold_075_crossing"]),
            "mouse_crossing": truthy(mr["threshold_075_crossing"]),
            "crossing_conserved": bool(
                truthy(hr["threshold_075_crossing"])
                and truthy(mr["threshold_075_crossing"])
                and direction_ok
                and same_placement
            ),
            "robust_cross_species_delta": robust_delta,
            "strict_grammar_portability": strict_grammar,
        }
        tf_rows.append(row)

    tf_rows.sort(key=lambda r: (
        not bool(r["strict_grammar_portability"]),
        -float(r["min_abs_delta"]),
        -float(r["min_site_strength"]),
        str(r["tf"]),
    ))
    for rank, row in enumerate(tf_rows, 1):
        row["cross_species_rank"] = rank

    robust_delta_rows = [r for r in tf_rows if r["robust_cross_species_delta"]]
    strict_rows = [r for r in tf_rows if r["strict_grammar_portability"]]
    direction_rows = [r for r in tf_rows if r["direction_conserved"]]
    placement_rows = [r for r in tf_rows if r["placement_conserved"]]
    crossing_rows = [r for r in tf_rows if r["crossing_conserved"]]
    summary = {
        "motif_records": len(human_records),
        "unique_tfs": len(tf_rows),
        "record_direction_conserved": all_record_direction_ok,
        "record_direction_conserved_fraction": all_record_direction_ok / len(human_records),
        "record_signed_delta_pearson": pearson(signed_h, signed_m),
        "record_abs_delta_pearson": pearson(abs_h, abs_m),
        "delta_magnitude_caveat": "for an additive PWM and one-base edit, fixed-column delta is determined by allele weights at that motif column; identical delta is not independent biological replication",
        "tf_direction_conserved": len(direction_rows),
        "tf_direction_conserved_fraction": len(direction_rows) / len(tf_rows),
        "tf_placement_conserved": len(placement_rows),
        "tf_placement_conserved_fraction": len(placement_rows) / len(tf_rows),
        "min_abs_delta_threshold": min_abs_delta,
        "robust_cross_species_delta_tfs": len(robust_delta_rows),
        "strict_grammar_portable_tfs": len(strict_rows),
        "conserved_crossing_tfs": len(crossing_rows),
        "strict_grammar_candidates": [str(r["tf"]) for r in strict_rows],
        "claim_scope": "human-mouse motif grammar portability only; not occupancy, Irx3 regulation, or causality",
    }
    return tf_rows, summary


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--human-csv", type=Path, required=True)
    p.add_argument("--mouse-csv", type=Path, required=True)
    p.add_argument("--min-abs-delta", type=float, default=0.10)
    p.add_argument("--output-dir", type=Path, required=True)
    args = p.parse_args()
    if not 0 <= args.min_abs_delta <= 1:
        raise ValueError("--min-abs-delta must be in [0,1]")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    rows, summary = compare_species(read_rows(args.human_csv), read_rows(args.mouse_csv), args.min_abs_delta)
    write_csv(args.output_dir / "cross_species_motif_stability.csv", rows)
    (args.output_dir / "cross_species_motif_stability_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    strict = [r for r in rows if r["strict_grammar_portability"]]
    report = [
        "# GAP-001 human-mouse motif grammar stability", "",
        "**Status: sequence-model cross-species robustness analysis; not occupancy or causality.**", "",
        f"- Motif records compared: **{summary['motif_records']}**",
        f"- Unique TFs: **{summary['unique_tfs']}**",
        f"- Record-level direction conserved: **{summary['record_direction_conserved']}** ({summary['record_direction_conserved_fraction']:.3f})",
        f"- Record-level signed Δ Pearson: **{summary['record_signed_delta_pearson']}**",
        f"- Record-level |Δ| Pearson: **{summary['record_abs_delta_pearson']}**",
        f"- Human-selected TF representatives with direction conserved: **{summary['tf_direction_conserved']}** ({summary['tf_direction_conserved_fraction']:.3f})",
        f"- Human-selected TF representatives with placement conserved: **{summary['tf_placement_conserved']}** ({summary['tf_placement_conserved_fraction']:.3f})",
        f"- Same-record direction+|Δ| >= {args.min_abs_delta}: **{summary['robust_cross_species_delta_tfs']}**",
        f"- Strict same-record + same-placement grammar candidates: **{summary['strict_grammar_portable_tfs']}**",
        f"- Conserved 0.75 crossings under strict placement: **{summary['conserved_crossing_tfs']}**", "",
        "## Top strict human→mouse grammar candidates", "",
        "|rank|TF|record|dir|min |Δ||human site|mouse site|strand|offset|",
        "|---:|---|---|---|---:|---:|---:|:---:|---:|",
    ]
    for row in strict[:30]:
        report.append(
            f"|{row['cross_species_rank']}|{row['tf']}|{row['motif_record_id']}|{row['human_direction']}|"
            f"{float(row['min_abs_delta']):.4f}|{float(row['human_max_site']):.4f}|{float(row['mouse_max_site']):.4f}|"
            f"{row['human_strand']}|{row['human_start_offset']}|"
        )

    by_tf = {str(r["tf"]): r for r in rows}
    report += ["", "## Named GAP-001 candidates", ""]
    for tf in ("CXXC4", "ARID5B", "CUX1", "CXXC1", "TET1", "KDM2B", "DNMT1"):
        row = by_tf.get(tf)
        if row is None:
            report.append(f"- `{tf}`: not present in the compared atlas.")
            continue
        report.append(
            f"- `{tf}`: {row['human_direction']} human Δ={float(row['human_delta']):.4f}, mouse Δ={float(row['mouse_delta']):.4f}; "
            f"placement human={row['human_strand']}/{row['human_start_offset']} mouse={row['mouse_strand']}/{row['mouse_start_offset']}; "
            f"site human={float(row['human_max_site']):.4f}, mouse={float(row['mouse_max_site']):.4f}; "
            f"strict_grammar_portability={row['strict_grammar_portability']}."
        )
    report += [
        "", "## Interpretation boundary", "",
        "- Identical Δ magnitude is not treated as independent replication: under an additive PWM, a one-base edit at the same motif column can mechanically preserve Δ across different flanks.",
        "- Strict portability therefore requires the same published motif record, same allele-effect direction, |Δ| threshold, and the same best variant-overlapping strand/offset.",
        "- Full-site strength is reported separately and is not given an invented biochemical cutoff.",
        "- A placement-changing candidate such as CXXC4 is demoted only as a literal same-grammar explanation of the mouse knock-in; it is not biologically ruled out.",
        "- None of these signals establishes adult-PH occupancy, enhancer→Irx3 regulation, methylation dependence, or mediation.",
    ]
    text = "\n".join(report) + "\n"
    (args.output_dir / "cross_species_motif_stability_report.md").write_text(text, encoding="utf-8")
    print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
