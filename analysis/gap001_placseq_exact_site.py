#!/usr/bin/env python3
"""Query processed adult hypothalamus PLAC-seq contacts at exact rs1421085 site.

This analysis keeps three authorities separate:

1. GSE295850 10-kb significant BEDPE calls: primary evidence that a contact was
   called at the authors' published threshold.
2. GSE295850 whitespace interaction-model table: confirms the same 10-kb pair
   and exposes 1D H3K27ac peak flags plus model statistics.
3. GSE295853 5-kb Fed/Fasted EdgeR table: asks whether the exact target/TSS bins
   are represented across metabolic states and preserves replicate-level calls.

All evidence is wild-type adult hypothalamus. A positive contact cannot establish
T>C allele dependence, exact-base enhancer activity, or causal mediation.
"""
from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path
from typing import Iterable

TARGET_CHROM = "chr8"
TARGET_SITE0 = 91374371
TARGET_TILE = (91374000, 91374500)
TARGET_BIN_10K = (91370000, 91380000)
TARGET_BIN_5K = (91370000, 91375000)

# NCBI Gene 16373, GRCm38.p6: 91798511..91801654, complement (1-based).
IRX3_GENE = (91798510, 91801654)  # 0-based half-open
IRX3_TSS0 = 91801653
IRX3_PROMOTER_BIN_10K = (91800000, 91810000)
IRX3_PROMOTER_BIN_5K = (91800000, 91805000)


def overlap(a0: int, a1: int, b0: int, b1: int) -> bool:
    return max(a0, b0) < min(a1, b1)


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def normalize_float(value: str):
    try:
        return float(value)
    except (TypeError, ValueError):
        return value


def normalize_int(value: str):
    try:
        return int(value)
    except (TypeError, ValueError):
        return value


def pair_matches(
    chrom1: str,
    start1: int,
    end1: int,
    chrom2: str,
    start2: int,
    end2: int,
    left: tuple[str, int, int],
    right: tuple[str, int, int],
) -> bool:
    lchrom, l0, l1 = left
    rchrom, r0, r1 = right
    forward = (
        chrom1 == lchrom
        and overlap(start1, end1, l0, l1)
        and chrom2 == rchrom
        and overlap(start2, end2, r0, r1)
    )
    reverse = (
        chrom1 == rchrom
        and overlap(start1, end1, r0, r1)
        and chrom2 == lchrom
        and overlap(start2, end2, l0, l1)
    )
    return forward or reverse


def parse_significant_bedpe(path: Path):
    header = None
    preview = []
    parsed = 0
    target_contacts = []
    promoter_contacts = []
    exact_pair = []

    with open_text(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(preview) < 5:
                preview.append(fields)
            if header is None:
                try:
                    int(fields[1])
                except (ValueError, IndexError):
                    header = fields
                    continue
            if len(fields) < 6:
                continue
            try:
                s1, e1, s2, e2 = map(int, [fields[1], fields[2], fields[4], fields[5]])
            except ValueError:
                continue
            parsed += 1
            row = {
                "anchor1": [fields[0], s1, e1],
                "anchor2": [fields[3], s2, e2],
            }
            if header and len(fields) == len(header):
                for key, value in zip(header[6:], fields[6:]):
                    if key in {"count", "ClusterSize", "ClusterSummit"}:
                        row[key] = normalize_int(value)
                    elif key in {"expected", "fdr", "ClusterNegLog10P"}:
                        row[key] = normalize_float(value)
                    else:
                        row[key] = value
            else:
                row["extra"] = fields[6:]

            a1_target = fields[0] == TARGET_CHROM and overlap(s1, e1, *TARGET_BIN_10K)
            a2_target = fields[3] == TARGET_CHROM and overlap(s2, e2, *TARGET_BIN_10K)
            a1_prom = fields[0] == TARGET_CHROM and overlap(s1, e1, *IRX3_PROMOTER_BIN_10K)
            a2_prom = fields[3] == TARGET_CHROM and overlap(s2, e2, *IRX3_PROMOTER_BIN_10K)
            if a1_target or a2_target:
                target_contacts.append(row)
            if a1_prom or a2_prom:
                promoter_contacts.append(row)
            if (a1_target and a2_prom) or (a2_target and a1_prom):
                exact_pair.append(row)

    if parsed == 0:
        raise RuntimeError(f"no significant BEDPE contacts parsed from {path}")
    return {
        "header": header,
        "preview": preview,
        "parsed_contacts": parsed,
        "contacts_with_target_bin_anchor": len(target_contacts),
        "contacts_with_irx3_promoter_bin_anchor": len(promoter_contacts),
        "target_bin_to_irx3_promoter_bin_contacts": len(exact_pair),
        "target_bin_contacts": target_contacts,
        "target_to_irx3_contacts": exact_pair,
    }


def parse_interaction_model(path: Path):
    """Parse the authors' whitespace-delimited 10-kb interaction model table."""
    header = None
    preview = []
    parsed = 0
    target_rows = []
    exact_pair = []
    target_mid = (TARGET_BIN_10K[0] + TARGET_BIN_10K[1]) // 2
    promoter_mid = (IRX3_PROMOTER_BIN_10K[0] + IRX3_PROMOTER_BIN_10K[1]) // 2

    with open_text(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            if len(preview) < 5:
                preview.append(fields)
            if header is None:
                if fields[:3] == ["bin1_mid", "bin2_mid", "count"]:
                    header = fields
                    continue
                raise RuntimeError("interaction-model header not recognized")
            if len(fields) != len(header):
                raise RuntimeError(
                    f"interaction-model row width {len(fields)} != header width {len(header)}"
                )
            row = dict(zip(header, fields))
            try:
                b1 = int(row["bin1_mid"])
                b2 = int(row["bin2_mid"])
            except ValueError as exc:
                raise RuntimeError(f"invalid interaction midpoints: {row}") from exc
            parsed += 1
            if row["chr"] != TARGET_CHROM:
                continue
            if b1 == target_mid or b2 == target_mid:
                target_rows.append({
                    "bin1_mid": b1,
                    "bin2_mid": b2,
                    "count": normalize_int(row["count"]),
                    "X1D_peak_bin1": normalize_int(row["X1D_peak_bin1"]),
                    "X1D_peak_bin2": normalize_int(row["X1D_peak_bin2"]),
                    "expected": normalize_float(row["expected"]),
                    "expected2": normalize_float(row["expected2"]),
                    "ratio2": normalize_float(row["ratio2"]),
                    "p_val_reg2": normalize_float(row["p_val_reg2"]),
                    "fdr": normalize_float(row["fdr"]),
                    "ClusterType": row["ClusterType"],
                    "lab": row["lab"],
                    "ClusterSize": normalize_int(row["ClusterSize"]),
                    "NegLog10P": normalize_float(row["NegLog10P"]),
                })
            if {b1, b2} == {target_mid, promoter_mid}:
                exact_pair.append({
                    "bin1_mid": b1,
                    "bin2_mid": b2,
                    "count": normalize_int(row["count"]),
                    "X1D_peak_bin1": normalize_int(row["X1D_peak_bin1"]),
                    "X1D_peak_bin2": normalize_int(row["X1D_peak_bin2"]),
                    "expected": normalize_float(row["expected"]),
                    "expected2": normalize_float(row["expected2"]),
                    "ratio2": normalize_float(row["ratio2"]),
                    "p_val": normalize_float(row["p_val"]),
                    "p_val_reg2": normalize_float(row["p_val_reg2"]),
                    "p_bonferroni": normalize_float(row["p_bonferroni"]),
                    "fdr": normalize_float(row["fdr"]),
                    "CountNei": normalize_int(row["CountNei"]),
                    "label": row["label"],
                    "lab": row["lab"],
                    "NegLog10P": normalize_float(row["NegLog10P"]),
                    "ClusterSize": normalize_int(row["ClusterSize"]),
                    "ClusterType": row["ClusterType"],
                })

    if parsed == 0:
        raise RuntimeError(f"no interaction-model rows parsed from {path}")
    return {
        "header": header,
        "preview": preview,
        "parsed_rows": parsed,
        "target_midpoint": target_mid,
        "irx3_promoter_midpoint": promoter_mid,
        "rows_with_target_midpoint": len(target_rows),
        "target_midpoint_rows": target_rows,
        "target_to_irx3_rows": exact_pair,
    }


def parse_fed_fasted(path: Path):
    """Parse GSE295853 5-kb Fed/Fasted EdgeR loop table at the exact pair."""
    header = None
    preview = []
    parsed = 0
    target_anchor_rows = 0
    promoter_anchor_rows = 0
    exact_pair = []

    with open_text(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(preview) < 5:
                preview.append(fields)
            if header is None:
                if fields[:6] == ["chr1", "start1", "end1", "chr2", "start2", "end2"]:
                    header = fields
                    continue
                raise RuntimeError("Fed/Fasted loop header not recognized")
            if len(fields) != len(header):
                raise RuntimeError(f"Fed/Fasted row width {len(fields)} != header width {len(header)}")
            row = dict(zip(header, fields))
            try:
                s1, e1 = int(row["start1"]), int(row["end1"])
                s2, e2 = int(row["start2"]), int(row["end2"])
            except ValueError as exc:
                raise RuntimeError(f"invalid Fed/Fasted coordinates: {row}") from exc
            parsed += 1

            a1_target = row["chr1"] == TARGET_CHROM and overlap(s1, e1, *TARGET_BIN_5K)
            a2_target = row["chr2"] == TARGET_CHROM and overlap(s2, e2, *TARGET_BIN_5K)
            a1_prom = row["chr1"] == TARGET_CHROM and overlap(s1, e1, *IRX3_PROMOTER_BIN_5K)
            a2_prom = row["chr2"] == TARGET_CHROM and overlap(s2, e2, *IRX3_PROMOTER_BIN_5K)
            if a1_target or a2_target:
                target_anchor_rows += 1
            if a1_prom or a2_prom:
                promoter_anchor_rows += 1
            if not ((a1_target and a2_prom) or (a2_target and a1_prom)):
                continue

            keep = {
                "anchor1": [row["chr1"], s1, e1],
                "anchor2": [row["chr2"], s2, e2],
                "Fed_ChIPCov1": normalize_float(row["Fed_ChIPCov1"]),
                "Fed_ChIPCov2": normalize_float(row["Fed_ChIPCov2"]),
                "Fasted_ChIPCov1": normalize_float(row["Fasted_ChIPCov1"]),
                "Fasted_ChIPCov2": normalize_float(row["Fasted_ChIPCov2"]),
                "Bin1_Label": row["Bin1_Label"],
                "Bin2_Label": row["Bin2_Label"],
                "logFC": normalize_float(row["logFC"]),
                "logCPM": normalize_float(row["logCPM"]),
                "PValue": normalize_float(row["PValue"]),
                "FDR": normalize_float(row["FDR"]),
                "Fed_SigRepl": normalize_int(row["Fed_SigRepl"]),
                "Fasted_SigRepl": normalize_int(row["Fasted_SigRepl"]),
                "fed_replicates": [],
                "fasted_replicates": [],
            }
            for i in range(1, 8):
                keep["fed_replicates"].append({
                    "replicate": i,
                    "RawCC": normalize_int(row[f"Fed_Fed{i}_RawCC"]),
                    "QVal": normalize_float(row[f"Fed_Fed{i}_QVal"]),
                })
                keep["fasted_replicates"].append({
                    "replicate": i,
                    "RawCC": normalize_int(row[f"Fasted_Fasted{i}_RawCC"]),
                    "QVal": normalize_float(row[f"Fasted_Fasted{i}_QVal"]),
                })
            exact_pair.append(keep)

    if parsed == 0:
        raise RuntimeError(f"no Fed/Fasted rows parsed from {path}")
    return {
        "header": header,
        "preview": preview,
        "parsed_rows": parsed,
        "target_5kb_bin": [TARGET_CHROM, *TARGET_BIN_5K],
        "irx3_promoter_5kb_bin": [TARGET_CHROM, *IRX3_PROMOTER_BIN_5K],
        "rows_with_target_anchor": target_anchor_rows,
        "rows_with_irx3_promoter_anchor": promoter_anchor_rows,
        "target_to_irx3_rows": exact_pair,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--contacts", type=Path, required=True)
    ap.add_argument("--model", "--peaks", dest="model", type=Path, required=True)
    ap.add_argument("--fasted-loops", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    contact = parse_significant_bedpe(args.contacts)
    model = parse_interaction_model(args.model)
    fed_fasted = parse_fed_fasted(args.fasted_loops)

    result = {
        "analysis": "GAP-001 exact rs1421085-bin adult hypothalamus PLAC-seq query",
        "claim_scope": "wild-type processed contact-map evidence only; no allele-specific or causal promotion",
        "authority": {
            "significant_calls": {
                "series": "GSE295850",
                "file": args.contacts.name,
                "assay": "H3K27ac+ PLAC-seq",
                "published_threshold": "FDR <1%, total read count >=12, observed/expected >2; n=10",
                "assembly": "GRCm38/mm10",
            },
            "interaction_model": {
                "series": "GSE295850",
                "file": args.model.name,
                "assembly": "GRCm38/mm10",
            },
            "fed_fasted": {
                "series": "GSE295853",
                "file": args.fasted_loops.name,
                "resolution": "5 kb",
                "assembly": "GRCm38/mm10",
                "role": "state/replicate context, not primary significant-call authority",
            },
        },
        "coordinates": {
            "exact_site_0_based": [TARGET_CHROM, TARGET_SITE0, TARGET_SITE0 + 1],
            "target_tile_0_based": [TARGET_CHROM, *TARGET_TILE],
            "target_contact_bin_10kb": [TARGET_CHROM, *TARGET_BIN_10K],
            "target_contact_bin_5kb": [TARGET_CHROM, *TARGET_BIN_5K],
            "irx3_gene_0_based": [TARGET_CHROM, *IRX3_GENE],
            "irx3_tss_0_based": IRX3_TSS0,
            "irx3_promoter_contact_bin_10kb": [TARGET_CHROM, *IRX3_PROMOTER_BIN_10K],
            "irx3_promoter_contact_bin_5kb": [TARGET_CHROM, *IRX3_PROMOTER_BIN_5K],
            "target_to_irx3_tss_bp": IRX3_TSS0 - TARGET_SITE0,
        },
        "significant_contact_callset": contact,
        "interaction_model_table": model,
        "fed_fasted_state_table": fed_fasted,
        "interpretation": {
            "positive_rule": (
                "A published significant GSE295850 call contains the exact-site 10-kb bin on one anchor "
                "and the Irx3 promoter 10-kb bin on the other, and the independent interaction-model table "
                "contains the same midpoint pair."
            ),
            "resolution_boundary": (
                "The PLAC-seq result localizes a 10-kb anchor containing rs1421085; it does not identify "
                "the single base or 500-bp tile as the physical contact endpoint."
            ),
            "allele_boundary": (
                "All PLAC-seq evidence is wild-type. A positive contact cannot show that T>C creates, "
                "strengthens, weakens, or mediates the contact."
            ),
            "cell_boundary": (
                "The PLAC-seq tissue is adult hypothalamus and is not resolved to the prioritized "
                "LHA-AHN-PVH Otp Trh Glut subclass or to the male posterior-hypothalamic IRX3+ circuit."
            ),
        },
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")

    print("SIGNIFICANT_CONTACTS", contact["parsed_contacts"])
    print("TARGET_BIN_CONTACTS", contact["contacts_with_target_bin_anchor"])
    print("IRX3_PROMOTER_CONTACTS", contact["contacts_with_irx3_promoter_bin_anchor"])
    print("TARGET_TO_IRX3_SIGNIFICANT", contact["target_bin_to_irx3_promoter_bin_contacts"])
    print("MODEL_TARGET_TO_IRX3", len(model["target_to_irx3_rows"]))
    print("FED_FASTED_TARGET_TO_IRX3", len(fed_fasted["target_to_irx3_rows"]))
    for row in contact["target_to_irx3_contacts"]:
        print("SIGNIFICANT_ROW", json.dumps(row, sort_keys=True))
    for row in model["target_to_irx3_rows"]:
        print("MODEL_ROW", json.dumps(row, sort_keys=True))
    for row in fed_fasted["target_to_irx3_rows"]:
        print("FED_FASTED_ROW", json.dumps(row, sort_keys=True))

    if contact["target_bin_to_irx3_promoter_bin_contacts"] and model["target_to_irx3_rows"]:
        print("H4_WILDTYPE_3D_SUPPORT: POSITIVE_AT_10KB_RESOLUTION")
    else:
        print("H4_WILDTYPE_3D_SUPPORT: NOT_CONFIRMED")
    print("GAP-001 remains OPEN; no allele dependence or mediator is established.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
