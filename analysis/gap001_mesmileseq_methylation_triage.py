#!/usr/bin/env python3
"""Join published meSMiLE methylation classes to rs1421085 motif switches.

Supplementary Data 5 contains a worksheet titled ``Classification by
methylation`` with columns TF, DBD, classification, misc. This script reads that
published classification verbatim, normalizes it into a small machine enum, and
joins it to an allele-switch CSV produced by gap001_codebook_motif_switch.py.

Important boundaries:
- absence from the published classification sheet => ``unknown``, not negative;
- methylation class is an annotation axis, never multiplied into a causal score;
- motif score + VPH transcript presence + methylation class remain feasibility
  evidence only, not in-vivo occupancy or a causal edge.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import zipfile

from gap001_probe_mesmileseq_xlsx import read_shared_strings, workbook_sheet_targets, read_rows

CLASS_SHEET = "Classification by methylation"
NORMAL_CLASSES = {"methyl_plus", "methyl_minus", "little_effect_no_cg", "unknown"}


def normalize_tf(value: object) -> str:
    return str(value).strip().upper()


def normalize_class(raw: object) -> tuple[str, str]:
    value = " ".join(str(raw).strip().lower().split())
    if value == "methyl plus":
        return "methyl_plus", "strong"
    if value == "weak methyl plus":
        return "methyl_plus", "weak"
    if value == "methyl minus":
        return "methyl_minus", "strong"
    if value == "weak methyl minus":
        return "methyl_minus", "weak"
    if value == "little effect":
        return "little_effect_no_cg", "little_effect"
    if value == "no cg":
        return "little_effect_no_cg", "no_cg"
    return "unknown", value or "unclassified"


def published_classification_rows(xlsx: Path) -> list[dict[str, object]]:
    with zipfile.ZipFile(xlsx) as zf:
        shared = read_shared_strings(zf)
        target = None
        for title, sheet_path in workbook_sheet_targets(zf):
            if title == CLASS_SHEET:
                target = sheet_path
                break
        if target is None:
            raise ValueError(f"missing published worksheet: {CLASS_SHEET}")
        rows = read_rows(zf, target, shared, max_rows=10000)
    if not rows:
        raise ValueError("classification worksheet is empty")
    header = [str(x).strip() for x in rows[0]]
    required = ["TF", "DBD", "classification", "misc"]
    missing = [x for x in required if x not in header]
    if missing:
        raise ValueError(f"missing classification columns: {missing}; header={header}")
    index = {name: header.index(name) for name in required}
    out: list[dict[str, object]] = []
    seen: set[str] = set()
    for raw_row in rows[1:]:
        def get(name: str) -> object:
            idx = index[name]
            return raw_row[idx] if idx < len(raw_row) else ""
        tf = normalize_tf(get("TF"))
        if not tf:
            continue
        if tf in seen:
            raise ValueError(f"duplicate TF in published methylation classification: {tf}")
        seen.add(tf)
        raw_class = str(get("classification")).strip()
        norm, strength = normalize_class(raw_class)
        if norm not in NORMAL_CLASSES:
            raise AssertionError(norm)
        out.append({
            "tf": tf,
            "dbd": str(get("DBD")).strip(),
            "published_classification": raw_class,
            "mesmile_class": norm,
            "class_strength": strength,
            "misc": str(get("misc")).strip(),
            "mesmile_profiled": True,
            "source_sheet": CLASS_SHEET,
        })
    if len(out) < 40:
        raise ValueError(f"classification sheet unexpectedly small: {len(out)} TFs")
    return out


def read_csv(path: Path | None) -> list[dict[str, str]]:
    if path is None:
        return []
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def truthy(value: object) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


def as_float(row: dict[str, object], key: str, default: float = 0.0) -> float:
    try:
        return float(row.get(key, default) or default)
    except (TypeError, ValueError):
        return default


def join_triage(
    motif_rows: list[dict[str, str]],
    classes: list[dict[str, object]],
    expression_rows: list[dict[str, str]] | None = None,
) -> list[dict[str, object]]:
    by_class = {normalize_tf(x["tf"]): x for x in classes}
    by_expr = {normalize_tf(x.get("tf", "")): x for x in (expression_rows or [])}
    out: list[dict[str, object]] = []
    for motif in motif_rows:
        tf = normalize_tf(motif.get("tf", ""))
        cls = by_class.get(tf)
        expr = by_expr.get(tf, {})
        row: dict[str, object] = dict(motif)
        row.update({
            "mesmile_profiled": bool(cls),
            "mesmile_class": cls["mesmile_class"] if cls else "unknown",
            "mesmile_class_strength": cls["class_strength"] if cls else "unprofiled",
            "mesmile_published_classification": cls["published_classification"] if cls else "",
            "mesmile_dbd": cls["dbd"] if cls else "",
            "mesmile_misc": cls["misc"] if cls else "",
            "methylation_relevant_class": bool(cls and cls["mesmile_class"] in {"methyl_plus", "methyl_minus"}),
            "expr_vph_positive_pct": expr.get("vph_positive_pct", ""),
            "expr_mouse_symbols": expr.get("mouse_symbols", ""),
        })
        out.append(row)

    # Transparent lexicographic triage, not a composite biological score.
    out.sort(key=lambda r: (
        not truthy(r.get("threshold_075_crossing", False)),
        -as_float(r, "max_site_strength"),
        not bool(r.get("methylation_relevant_class")),
        -(as_float(r, "expr_vph_positive_pct", -1.0) if str(r.get("expr_vph_positive_pct", "")) != "" else -1.0),
        -as_float(r, "abs_delta_normalized"),
        str(r.get("tf", "")),
    ))
    for idx, row in enumerate(out, 1):
        row["triage_rank"] = idx
    return out


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--xlsx", type=Path, required=True)
    parser.add_argument("--motif-switches", type=Path, required=True)
    parser.add_argument("--expression-csv", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--shortlist-size", type=int, default=30)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    classes = published_classification_rows(args.xlsx)
    motifs = read_csv(args.motif_switches)
    if not motifs:
        raise ValueError("motif-switch CSV is empty")
    expression = read_csv(args.expression_csv) if args.expression_csv else []
    joined = join_triage(motifs, classes, expression)

    write_csv(args.output_dir / "mesmile_published_classes.csv", classes)
    write_csv(args.output_dir / "motif_methylation_triage.csv", joined)
    shortlist = joined[: args.shortlist_size]
    write_csv(args.output_dir / "methylation_shortlist.csv", shortlist)

    counts: dict[str, int] = {}
    for row in classes:
        key = str(row["mesmile_class"])
        counts[key] = counts.get(key, 0) + 1
    overlap = sum(bool(row["mesmile_profiled"]) for row in joined)
    relevant = sum(bool(row["methylation_relevant_class"]) for row in joined)
    crossings = [r for r in joined if truthy(r.get("threshold_075_crossing", False))]
    crossing_profiled = [r for r in crossings if r["mesmile_profiled"]]

    summary = {
        "published_classification_tfs": len(classes),
        "published_class_counts": counts,
        "motifs_joined": len(joined),
        "motifs_with_mesmileseq_classification": overlap,
        "motifs_with_methyl_plus_or_minus_class": relevant,
        "threshold_crossings": len(crossings),
        "threshold_crossings_profiled_by_mesmileseq": len(crossing_profiled),
        "shortlist_size": len(shortlist),
        "claim_scope": "triage only; no TF occupancy, Irx3 regulation, or causality inferred",
    }
    (args.output_dir / "methylation_triage_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    report = [
        "# GAP-001 methylation-aware motif triage",
        "",
        "**Status: published methylation-class annotation + sequence/VPH feasibility; not occupancy or causality.**",
        "",
        f"- Published meSMiLE classifications parsed: **{len(classes)} TFs**",
        f"- Motif candidates joined: **{len(joined)}**",
        f"- Candidates with meSMiLE classification: **{overlap}**",
        f"- Candidates with methyl-plus/minus class: **{relevant}**",
        f"- 0.75 heuristic crossings: **{len(crossings)}**; profiled by meSMiLE: **{len(crossing_profiled)}**",
        "",
        "## Top transparent triage rows",
        "",
        "|rank|TF|allele direction|Δnorm|max site|crossing|VPH +%|meSMiLE class|strength|",
        "|---:|---|---|---:|---:|:---:|---:|---|---|",
    ]
    for row in shortlist[:20]:
        report.append(
            f"|{row['triage_rank']}|{row.get('tf','')}|{row.get('direction','')}|"
            f"{as_float(row,'abs_delta_normalized'):.4f}|{as_float(row,'max_site_strength'):.4f}|"
            f"{'yes' if truthy(row.get('threshold_075_crossing')) else 'no'}|"
            f"{str(row.get('expr_vph_positive_pct',''))}|{row.get('mesmile_class','unknown')}|"
            f"{row.get('mesmile_class_strength','')}|"
        )
    report += [
        "",
        "## Interpretation boundary",
        "",
        "- `unknown` means unprofiled/unclassified here, never methylation-insensitive.",
        "- `weak methyl plus/minus` is preserved through the `class_strength=weak` field.",
        "- `little effect` and `no CG` are grouped only in the machine enum; the published raw label is retained.",
        "- The sort is lexicographic triage, not a fitted score or biological probability.",
        "- Decisive evidence still requires allele-specific occupancy and perturbation in the relevant adult PH cell state.",
    ]
    text = "\n".join(report) + "\n"
    (args.output_dir / "methylation_triage_report.md").write_text(text, encoding="utf-8")
    print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
