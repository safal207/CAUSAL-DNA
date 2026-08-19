#!/usr/bin/env python3
"""Compare rs1421085 allele-switch calls across motif sources.

This gate exists because the 235-motif manual Codebook subset and the official
full 1,421-TF representative atlas can use different PWM representatives for the
same TF. Candidate priority must not depend silently on motif-source choice.

The analysis collapses each source to the most allele-sensitive representative
per TF, then compares direction, effect magnitude, site strength, and the
explicit 0.75 crossing heuristic. It never converts agreement into occupancy or
causality.
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


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        raise ValueError(f"empty motif-switch CSV: {path}")
    required = {"tf", "direction", "abs_delta_normalized", "max_site_strength", "threshold_075_crossing"}
    missing = required.difference(rows[0])
    if missing:
        raise ValueError(f"missing columns in {path}: {sorted(missing)}")
    return rows


def best_by_tf(rows: list[dict[str, str]]) -> dict[str, dict[str, object]]:
    grouped: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        grouped.setdefault(row["tf"], []).append(row)
    out: dict[str, dict[str, object]] = {}
    for tf, reps in grouped.items():
        best = max(reps, key=lambda r: (f(r, "abs_delta_normalized"), f(r, "max_site_strength"), r.get("motif_record_id", "")))
        out[tf] = {
            "tf": tf,
            "representative_count": len(reps),
            "motif_record_id": best.get("motif_record_id", "") or tf,
            "direction": best["direction"],
            "delta": f(best, "delta_normalized"),
            "abs_delta": f(best, "abs_delta_normalized"),
            "max_site": f(best, "max_site_strength"),
            "crossing": truthy(best["threshold_075_crossing"]),
        }
    return out


def same_nonzero_direction(a: dict[str, object], b: dict[str, object]) -> bool:
    da, db = float(a["delta"]), float(b["delta"])
    return da != 0.0 and db != 0.0 and ((da > 0) == (db > 0))


def pearson(xs: list[float], ys: list[float]) -> float | None:
    if len(xs) < 2 or len(xs) != len(ys):
        return None
    mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
    dx = [x - mx for x in xs]
    dy = [y - my for y in ys]
    denom = math.sqrt(sum(x*x for x in dx) * sum(y*y for y in dy))
    return None if denom == 0 else sum(x*y for x, y in zip(dx, dy)) / denom


def compare_sources(
    a_rows: list[dict[str, str]],
    b_rows: list[dict[str, str]],
    source_a: str,
    source_b: str,
    min_abs_delta: float = 0.10,
) -> tuple[list[dict[str, object]], dict[str, object]]:
    a = best_by_tf(a_rows)
    b = best_by_tf(b_rows)
    overlap = sorted(set(a) & set(b))
    comparisons: list[dict[str, object]] = []
    signed_a: list[float] = []
    signed_b: list[float] = []
    abs_a: list[float] = []
    abs_b: list[float] = []

    for tf in overlap:
        ra, rb = a[tf], b[tf]
        sign_ok = same_nonzero_direction(ra, rb)
        crossing_both = bool(ra["crossing"] and rb["crossing"] and sign_ok)
        robust_delta = sign_ok and float(ra["abs_delta"]) >= min_abs_delta and float(rb["abs_delta"]) >= min_abs_delta
        row = {
            "tf": tf,
            f"{source_a}_record": ra["motif_record_id"],
            f"{source_b}_record": rb["motif_record_id"],
            f"{source_a}_direction": ra["direction"],
            f"{source_b}_direction": rb["direction"],
            "direction_consistent": sign_ok,
            f"{source_a}_abs_delta": ra["abs_delta"],
            f"{source_b}_abs_delta": rb["abs_delta"],
            "min_abs_delta": min(float(ra["abs_delta"]), float(rb["abs_delta"])),
            f"{source_a}_max_site": ra["max_site"],
            f"{source_b}_max_site": rb["max_site"],
            f"{source_a}_crossing": ra["crossing"],
            f"{source_b}_crossing": rb["crossing"],
            "crossing_replicated": crossing_both,
            "robust_delta_direction": robust_delta,
            "source_sensitive_crossing": bool(ra["crossing"] != rb["crossing"]),
            f"{source_a}_representatives": ra["representative_count"],
            f"{source_b}_representatives": rb["representative_count"],
        }
        comparisons.append(row)
        signed_a.append(float(ra["delta"])); signed_b.append(float(rb["delta"]))
        abs_a.append(float(ra["abs_delta"])); abs_b.append(float(rb["abs_delta"]))

    comparisons.sort(key=lambda r: (
        not bool(r["crossing_replicated"]),
        not bool(r["robust_delta_direction"]),
        -float(r["min_abs_delta"]),
        str(r["tf"]),
    ))
    for i, row in enumerate(comparisons, 1):
        row["stability_rank"] = i

    a_cross = {tf for tf, r in a.items() if r["crossing"]}
    b_cross = {tf for tf, r in b.items() if r["crossing"]}
    overlap_a_cross = a_cross & set(overlap)
    overlap_b_cross = b_cross & set(overlap)
    crossing_repl = {r["tf"] for r in comparisons if r["crossing_replicated"]}
    robust = {r["tf"] for r in comparisons if r["robust_delta_direction"]}
    direction_ok = {r["tf"] for r in comparisons if r["direction_consistent"]}

    summary = {
        "source_a": source_a,
        "source_b": source_b,
        "source_a_unique_tfs": len(a),
        "source_b_unique_tfs": len(b),
        "overlap_tfs": len(overlap),
        "direction_consistent_tfs": len(direction_ok),
        "direction_consistency_fraction": len(direction_ok) / len(overlap) if overlap else 0.0,
        "signed_delta_pearson": pearson(signed_a, signed_b),
        "abs_delta_pearson": pearson(abs_a, abs_b),
        "min_abs_delta_threshold": min_abs_delta,
        "robust_delta_direction_tfs": len(robust),
        "source_a_crossings_in_overlap": len(overlap_a_cross),
        "source_b_crossings_in_overlap": len(overlap_b_cross),
        "crossing_replicated_tfs": len(crossing_repl),
        "source_a_only_crossings": sorted(overlap_a_cross - overlap_b_cross),
        "source_b_only_crossings": sorted(overlap_b_cross - overlap_a_cross),
        "replicated_crossings": sorted(crossing_repl),
        "robust_delta_direction_candidates": sorted(robust),
        "claim_scope": "motif-source robustness only; not occupancy, target-gene regulation, or causality",
    }
    return comparisons, summary


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fields = list(rows[0])
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader(); writer.writerows(rows)


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--source-a-csv", type=Path, required=True)
    p.add_argument("--source-b-csv", type=Path, required=True)
    p.add_argument("--source-a-name", default="manual235")
    p.add_argument("--source-b-name", default="full1421")
    p.add_argument("--min-abs-delta", type=float, default=0.10)
    p.add_argument("--output-dir", type=Path, required=True)
    args = p.parse_args()
    if not 0 <= args.min_abs_delta <= 1:
        raise ValueError("--min-abs-delta must be in [0,1]")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    comparisons, summary = compare_sources(
        read_rows(args.source_a_csv), read_rows(args.source_b_csv),
        args.source_a_name, args.source_b_name, args.min_abs_delta,
    )
    write_csv(args.output_dir / "motif_source_stability.csv", comparisons)
    (args.output_dir / "motif_source_stability_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")

    stable = [r for r in comparisons if r["robust_delta_direction"]]
    report = [
        "# GAP-001 motif-source stability", "",
        "**Status: sequence-model robustness analysis; not occupancy or causality.**", "",
        f"- TF overlap: **{summary['overlap_tfs']}**",
        f"- Direction-consistent TFs: **{summary['direction_consistent_tfs']}** ({summary['direction_consistency_fraction']:.3f})",
        f"- Signed Δ Pearson: **{summary['signed_delta_pearson']}**",
        f"- |Δ| Pearson: **{summary['abs_delta_pearson']}**",
        f"- Robust same-direction TFs with |Δ| >= {args.min_abs_delta} in both sources: **{summary['robust_delta_direction_tfs']}**",
        f"- Replicated 0.75 crossings: **{summary['crossing_replicated_tfs']}**",
        "",
        "## Top cross-source stable candidates", "",
        "|rank|TF|min |Δ||crossing replicated|source-sensitive crossing|",
        "|---:|---|---:|:---:|:---:|",
    ]
    for row in stable[:30]:
        report.append(f"|{row['stability_rank']}|{row['tf']}|{float(row['min_abs_delta']):.4f}|{'yes' if row['crossing_replicated'] else 'no'}|{'yes' if row['source_sensitive_crossing'] else 'no'}|")
    report += [
        "", "## Crossing sensitivity", "",
        f"- {args.source_a_name}-only crossings in overlap: `{', '.join(summary['source_a_only_crossings']) or 'none'}`",
        f"- {args.source_b_name}-only crossings in overlap: `{', '.join(summary['source_b_only_crossings']) or 'none'}`",
        f"- replicated crossings: `{', '.join(summary['replicated_crossings']) or 'none'}`",
        "", "## Boundary", "",
        "- Different PWM representatives can change threshold-crossing status without changing the underlying genome.",
        "- Direction agreement is more stable than a single arbitrary threshold and is reported separately.",
        "- Cross-source agreement cannot establish TF occupancy or mediation of Irx3.",
    ]
    text = "\n".join(report) + "\n"
    (args.output_dir / "motif_source_stability_report.md").write_text(text)
    print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
