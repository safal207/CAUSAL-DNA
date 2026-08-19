#!/usr/bin/env python3
"""Replicate-aware subclass analysis for the exact rs1421085 target tile.

Inputs are a frozen target-barcode evidence pack and the authors' official
GSE246791 Supplementary Table 2. The biological sample is the replication
stratum. All subclasses are tested; Benjamini-Hochberg correction is therefore
applied across the full subclass family.

This is descriptive/statistical prioritization in wild-type atlas data. It does
not establish allele dependence, enhancer function, target gene or causality.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from statsmodels.stats.contingency_tables import StratifiedTable
from statsmodels.stats.multitest import multipletests


def clean_subclass(value: str) -> str:
    return re.sub(r"^\d+\s+", "", value.strip())


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--barcodes", required=True, type=Path)
    ap.add_argument("--metadata", required=True, type=Path)
    ap.add_argument("--output", required=True, type=Path)
    args = ap.parse_args()

    pack = json.loads(args.barcodes.read_text())
    if pack.get("dataset") != "GSE246791":
        raise RuntimeError("unexpected barcode evidence dataset")
    if pack.get("totals", {}).get("samples") != 8:
        raise RuntimeError("analysis requires the frozen 8-sample evidence pack")

    sample_barcodes = {
        row["sample"]: set(row["target_tile_barcodes"])
        for row in pack["samples"]
    }
    if len(sample_barcodes) != 8:
        raise RuntimeError("sample IDs must be unique")

    sample_total = Counter()
    subclass_total = Counter()
    subclass_target = Counter()
    ss_total = Counter()
    ss_target = Counter()
    matched = defaultdict(set)

    with args.metadata.open("r", encoding="utf-8", errors="replace", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        needed = {"Sample", "Barcode", "Subclass"}
        if not needed.issubset(reader.fieldnames or []):
            raise RuntimeError(f"metadata missing columns: {needed - set(reader.fieldnames or [])}")
        for row in reader:
            sample = row["Sample"]
            if sample not in sample_barcodes:
                continue
            subclass = clean_subclass(row["Subclass"])
            barcode = row["Barcode"]
            sample_total[sample] += 1
            subclass_total[subclass] += 1
            ss_total[(sample, subclass)] += 1
            if barcode in sample_barcodes[sample]:
                subclass_target[subclass] += 1
                ss_target[(sample, subclass)] += 1
                matched[sample].add(barcode)

    for sample, expected in sample_barcodes.items():
        if matched[sample] != expected:
            missing = sorted(expected - matched[sample])
            extra = sorted(matched[sample] - expected)
            raise RuntimeError(
                f"barcode join mismatch for {sample}: missing={missing[:5]} extra={extra[:5]}"
            )

    metadata_n = sum(sample_total.values())
    target_n = sum(len(v) for v in sample_barcodes.values())
    if metadata_n != pack["totals"]["nuclei"]:
        raise RuntimeError(f"metadata nuclei {metadata_n} != frozen total {pack['totals']['nuclei']}")
    if target_n != pack["totals"]["target_tile_nuclei"]:
        raise RuntimeError("target barcode total mismatch")

    overall_rate = target_n / metadata_n
    rows = []
    subclasses = sorted(subclass_total)
    raw_p = []

    for subclass in subclasses:
        tables = []
        sample_hits = 0
        for sample in sorted(sample_barcodes):
            a = ss_target[(sample, subclass)]
            sub_n = ss_total[(sample, subclass)]
            target = len(sample_barcodes[sample])
            total = sample_total[sample]
            b = sub_n - a
            c = target - a
            d = total - sub_n - c
            if min(a, b, c, d) < 0:
                raise RuntimeError(f"invalid 2x2 table for {sample}/{subclass}")
            tables.append(np.array([[a, b], [c, d]], dtype=float))
            sample_hits += int(a > 0)

        arr = np.stack(tables, axis=2)
        st = StratifiedTable(arr)
        null_test = st.test_null_odds()
        try:
            heterogeneity = st.test_equal_odds()
            heterogeneity_p = float(heterogeneity.pvalue)
        except Exception:
            heterogeneity_p = None
        ci_low, ci_high = st.oddsratio_pooled_confint()
        target = subclass_target[subclass]
        total = subclass_total[subclass]
        pooled_enrichment = (target / total) / overall_rate if total else 0.0
        p = float(null_test.pvalue) if np.isfinite(null_test.pvalue) else 1.0
        raw_p.append(p)
        rows.append({
            "subclass": subclass,
            "target_nuclei": int(target),
            "total_nuclei": int(total),
            "samples_with_target": int(sample_hits),
            "samples_total": 8,
            "pooled_enrichment": float(pooled_enrichment),
            "mantel_haenszel_or": float(st.oddsratio_pooled) if np.isfinite(st.oddsratio_pooled) else None,
            "mh_or_95ci": [
                float(ci_low) if np.isfinite(ci_low) else None,
                float(ci_high) if np.isfinite(ci_high) else None,
            ],
            "mh_null_p": p,
            "breslow_day_equal_odds_p": heterogeneity_p,
        })

    reject, qvals, _, _ = multipletests(raw_p, alpha=0.05, method="fdr_bh")
    for row, reject_i, q in zip(rows, reject, qvals):
        row["bh_q"] = float(q)
        row["bh_fdr_0_05"] = bool(reject_i)

    rows.sort(key=lambda r: (r["bh_q"], -(r["samples_with_target"]), -(r["target_nuclei"])))
    result = {
        "analysis": "GAP-001 GSE246791 all8 replicate-aware subclass enrichment",
        "claim_scope": "wild-type atlas prioritization only; no allele-specific or causal promotion",
        "method": {
            "replication_unit": "GSE246791 biological sample",
            "stratified_test": "Mantel-Haenszel common odds ratio/test via statsmodels StratifiedTable",
            "heterogeneity": "Breslow-Day equal-odds test where estimable",
            "multiple_testing": "Benjamini-Hochberg across every official subclass represented in the 8 samples",
        },
        "totals": {
            "samples": 8,
            "nuclei": metadata_n,
            "target_tile_nuclei": target_n,
            "overall_target_rate": overall_rate,
            "official_subclasses_tested": len(subclasses),
        },
        "priority_rule": (
            "Prefer replicated subclasses over one-sample extreme enrichments; statistical significance alone "
            "does not identify the causal cell of action."
        ),
        "subclasses": rows,
    }
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")

    print(json.dumps(result["totals"], indent=2))
    print("TOP REPLICATED")
    for row in sorted(rows, key=lambda r: (-r["samples_with_target"], r["bh_q"]))[:12]:
        print(
            row["subclass"],
            f"target={row['target_nuclei']}/{row['total_nuclei']}",
            f"samples={row['samples_with_target']}/8",
            f"enrichment={row['pooled_enrichment']:.3f}",
            f"MH_OR={row['mantel_haenszel_or']}",
            f"q={row['bh_q']:.6g}",
            f"BD_p={row['breslow_day_equal_odds_p']}",
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
