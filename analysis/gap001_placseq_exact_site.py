#!/usr/bin/env python3
"""Query adult hypothalamus PLAC-seq at the exact rs1421085-containing bin.

Authorities are kept distinct:
- GSE295850 significant 10-kb BEDPE callset = primary called-contact evidence.
- GSE295850 whitespace model table = semantic/statistical cross-check of that call.
- GSE295853 5-kb Fed/Fasted table = metabolic-state / replicate context.

All inputs are wild-type adult hypothalamus. A positive 10-kb contact does not
identify the single rs1421085 base as the contacting nucleotide and cannot show
T>C allele dependence or causal mediation.
"""
from __future__ import annotations

import argparse
import gzip
import json
import math
from pathlib import Path

TARGET_CHROM = "chr8"
TARGET_SITE0 = 91374371
TARGET_TILE = (91374000, 91374500)
TARGET_BIN_10K = (91370000, 91380000)
TARGET_BIN_5K = (91370000, 91375000)

# NCBI Gene 16373, GRCm38.p6: 91798511..91801654, complement (1-based).
IRX3_GENE = (91798510, 91801654)
IRX3_TSS0 = 91801653
IRX3_BIN_10K = (91800000, 91810000)
IRX3_BIN_5K = (91800000, 91805000)


def overlap(a0: int, a1: int, b0: int, b1: int) -> bool:
    return max(a0, b0) < min(a1, b1)


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def number(value: str):
    try:
        return float(value)
    except ValueError:
        return value


def integral_number(value: str) -> int:
    """Accept ordinary or scientific-notation integer coordinates."""
    x = float(value)
    r = round(x)
    if not math.isclose(x, r, rel_tol=0.0, abs_tol=1e-6):
        raise ValueError(f"expected integral numeric value, got {value}")
    return int(r)


def pair_overlap(fields, left: tuple[int, int], right: tuple[int, int]) -> bool:
    try:
        c1, s1, e1, c2, s2, e2 = (
            fields[0], int(fields[1]), int(fields[2]),
            fields[3], int(fields[4]), int(fields[5]),
        )
    except (ValueError, IndexError):
        return False
    forward = c1 == TARGET_CHROM and c2 == TARGET_CHROM and overlap(s1, e1, *left) and overlap(s2, e2, *right)
    reverse = c1 == TARGET_CHROM and c2 == TARGET_CHROM and overlap(s1, e1, *right) and overlap(s2, e2, *left)
    return forward or reverse


def parse_significant(path: Path):
    header = None
    parsed = 0
    target_anchor_n = 0
    irx3_anchor_n = 0
    exact = []
    target_preview = []

    with open_text(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if header is None:
                try:
                    int(f[1])
                except (ValueError, IndexError):
                    header = f
                    continue
            if len(f) < 6:
                continue
            try:
                s1, e1, s2, e2 = int(f[1]), int(f[2]), int(f[4]), int(f[5])
            except ValueError:
                continue
            parsed += 1
            a1t = f[0] == TARGET_CHROM and overlap(s1, e1, *TARGET_BIN_10K)
            a2t = f[3] == TARGET_CHROM and overlap(s2, e2, *TARGET_BIN_10K)
            a1i = f[0] == TARGET_CHROM and overlap(s1, e1, *IRX3_BIN_10K)
            a2i = f[3] == TARGET_CHROM and overlap(s2, e2, *IRX3_BIN_10K)
            if a1t or a2t:
                target_anchor_n += 1
            if a1i or a2i:
                irx3_anchor_n += 1
            if not ((a1t and a2i) or (a2t and a1i)):
                continue
            row = {"anchor1": [f[0], s1, e1], "anchor2": [f[3], s2, e2]}
            if header and len(f) == len(header):
                for key, value in zip(header[6:], f[6:]):
                    row[key] = number(value)
            else:
                row["extra"] = f[6:]
            exact.append(row)

    if not parsed:
        raise RuntimeError("no significant contacts parsed")
    return {
        "header": header,
        "parsed_contacts": parsed,
        "contacts_with_target_bin_anchor": target_anchor_n,
        "contacts_with_irx3_promoter_bin_anchor": irx3_anchor_n,
        "target_bin_to_irx3_promoter_bin_contacts": len(exact),
        "target_to_irx3_contacts": exact,
    }


def parse_model(path: Path, primary: dict):
    header = None
    parsed = 0
    target_coordinate_candidates = []
    stat_matches = []

    p_count = float(primary["count"])
    p_expected = float(primary["expected"])
    p_fdr = float(primary["fdr"])

    with open_text(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.split()
            if header is None:
                if f[:3] != ["bin1_mid", "bin2_mid", "count"]:
                    raise RuntimeError("unrecognized interaction-model header")
                header = f
                continue
            if len(f) != len(header):
                raise RuntimeError(f"model row width {len(f)} != header width {len(header)}")
            row = dict(zip(header, f))
            parsed += 1
            if row["chr"] != TARGET_CHROM:
                continue
            b1, b2 = integral_number(row["bin1_mid"]), integral_number(row["bin2_mid"])
            # Preserve nearby coordinate candidates without assuming whether the
            # author's '*_mid' is encoded as a BED start, center, or rounded bin coordinate.
            target_near = min(abs(b1 - TARGET_SITE0), abs(b2 - TARGET_SITE0)) <= 10000
            irx3_near = min(abs(b1 - IRX3_TSS0), abs(b2 - IRX3_TSS0)) <= 10000
            compact = {
                "bin1_mid": b1,
                "bin2_mid": b2,
                "count": integral_number(row["count"]),
                "X1D_peak_bin1": integral_number(row["X1D_peak_bin1"]),
                "X1D_peak_bin2": integral_number(row["X1D_peak_bin2"]),
                "expected": float(row["expected"]),
                "expected2": float(row["expected2"]),
                "ratio2": float(row["ratio2"]),
                "p_val": float(row["p_val"]),
                "p_val_reg2": float(row["p_val_reg2"]),
                "p_bonferroni": float(row["p_bonferroni"]),
                "fdr": float(row["fdr"]),
                "CountNei": integral_number(row["CountNei"]),
                "label": row["label"],
                "NegLog10P": float(row["NegLog10P"]),
                "ClusterSize": integral_number(row["ClusterSize"]),
                "lab": row["lab"],
                "ClusterType": row["ClusterType"],
            }
            if target_near and irx3_near:
                target_coordinate_candidates.append(compact)
            if (
                math.isclose(compact["count"], p_count, rel_tol=0.0, abs_tol=0.0)
                and math.isclose(compact["expected"], p_expected, rel_tol=1e-12, abs_tol=1e-12)
                and math.isclose(compact["fdr"], p_fdr, rel_tol=1e-10, abs_tol=1e-20)
            ):
                stat_matches.append(compact)

    if not parsed:
        raise RuntimeError("no interaction-model rows parsed")
    return {
        "header": header,
        "parsed_rows": parsed,
        "coordinate_candidates_near_target_and_irx3": target_coordinate_candidates,
        "rows_matching_primary_count_expected_fdr": stat_matches,
        "note": (
            "The model table is a semantic/statistical cross-check of the GSE295850 callset, "
            "not an independent biological experiment."
        ),
    }


def parse_fed_fasted(path: Path):
    header = None
    parsed = 0
    target_anchor_n = 0
    irx3_anchor_n = 0
    exact = []

    with open_text(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if header is None:
                if f[:6] != ["chr1", "start1", "end1", "chr2", "start2", "end2"]:
                    raise RuntimeError("unrecognized Fed/Fasted header")
                header = f
                continue
            if len(f) != len(header):
                raise RuntimeError(f"Fed/Fasted row width {len(f)} != header width {len(header)}")
            row = dict(zip(header, f))
            parsed += 1
            s1, e1, s2, e2 = int(row["start1"]), int(row["end1"]), int(row["start2"]), int(row["end2"])
            a1t = row["chr1"] == TARGET_CHROM and overlap(s1, e1, *TARGET_BIN_5K)
            a2t = row["chr2"] == TARGET_CHROM and overlap(s2, e2, *TARGET_BIN_5K)
            a1i = row["chr1"] == TARGET_CHROM and overlap(s1, e1, *IRX3_BIN_5K)
            a2i = row["chr2"] == TARGET_CHROM and overlap(s2, e2, *IRX3_BIN_5K)
            if a1t or a2t:
                target_anchor_n += 1
            if a1i or a2i:
                irx3_anchor_n += 1
            if not ((a1t and a2i) or (a2t and a1i)):
                continue

            out = {
                "anchor1": [row["chr1"], s1, e1],
                "anchor2": [row["chr2"], s2, e2],
                "Fed_ChIPCov1": number(row["Fed_ChIPCov1"]),
                "Fed_ChIPCov2": number(row["Fed_ChIPCov2"]),
                "Fasted_ChIPCov1": number(row["Fasted_ChIPCov1"]),
                "Fasted_ChIPCov2": number(row["Fasted_ChIPCov2"]),
                "Bin1_Label": row["Bin1_Label"],
                "Bin2_Label": row["Bin2_Label"],
                "logFC": number(row["logFC"]),
                "logCPM": number(row["logCPM"]),
                "PValue": number(row["PValue"]),
                "FDR": number(row["FDR"]),
                "Fed_SigRepl": integral_number(row["Fed_SigRepl"]),
                "Fasted_SigRepl": integral_number(row["Fasted_SigRepl"]),
                "fed_replicates": [],
                "fasted_replicates": [],
            }
            for i in range(1, 8):
                out["fed_replicates"].append({
                    "replicate": i,
                    "RawCC": integral_number(row[f"Fed_Fed{i}_RawCC"]),
                    "QVal": number(row[f"Fed_Fed{i}_QVal"]),
                })
                out["fasted_replicates"].append({
                    "replicate": i,
                    "RawCC": integral_number(row[f"Fasted_Fasted{i}_RawCC"]),
                    "QVal": number(row[f"Fasted_Fasted{i}_QVal"]),
                })
            exact.append(out)

    if not parsed:
        raise RuntimeError("no Fed/Fasted rows parsed")
    return {
        "header": header,
        "parsed_rows": parsed,
        "target_5kb_bin": [TARGET_CHROM, *TARGET_BIN_5K],
        "irx3_promoter_5kb_bin": [TARGET_CHROM, *IRX3_BIN_5K],
        "rows_with_target_anchor": target_anchor_n,
        "rows_with_irx3_promoter_anchor": irx3_anchor_n,
        "target_to_irx3_rows": exact,
        "role": "metabolic-state and replicate context; not primary significant-call authority",
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--contacts", type=Path, required=True)
    ap.add_argument("--model", "--peaks", dest="model", type=Path, required=True)
    ap.add_argument("--fasted-loops", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    contact = parse_significant(args.contacts)
    if len(contact["target_to_irx3_contacts"]) != 1:
        raise RuntimeError(
            "expected one published 10-kb target→Irx3 significant contact; "
            f"found {len(contact['target_to_irx3_contacts'])}"
        )
    primary = contact["target_to_irx3_contacts"][0]
    model = parse_model(args.model, primary)
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
            "interaction_model": {"series": "GSE295850", "file": args.model.name, "assembly": "GRCm38/mm10"},
            "fed_fasted": {"series": "GSE295853", "file": args.fasted_loops.name, "assembly": "GRCm38/mm10"},
        },
        "coordinates": {
            "exact_site_0_based": [TARGET_CHROM, TARGET_SITE0, TARGET_SITE0 + 1],
            "target_tile_0_based": [TARGET_CHROM, *TARGET_TILE],
            "target_contact_bin_10kb": [TARGET_CHROM, *TARGET_BIN_10K],
            "target_contact_bin_5kb": [TARGET_CHROM, *TARGET_BIN_5K],
            "irx3_gene_0_based": [TARGET_CHROM, *IRX3_GENE],
            "irx3_tss_0_based": IRX3_TSS0,
            "irx3_promoter_contact_bin_10kb": [TARGET_CHROM, *IRX3_BIN_10K],
            "irx3_promoter_contact_bin_5kb": [TARGET_CHROM, *IRX3_BIN_5K],
            "target_to_irx3_tss_bp": IRX3_TSS0 - TARGET_SITE0,
        },
        "significant_contact_callset": contact,
        "interaction_model_table": model,
        "fed_fasted_state_table": fed_fasted,
        "interpretation": {
            "supported": (
                "The published adult-hypothalamus GSE295850 callset contains a significant 10-kb interaction "
                "between the bin containing rs1421085 and the 10-kb bin containing the Irx3 promoter/TSS."
            ),
            "resolution_boundary": (
                "The contact is localized to a 10-kb anchor that contains the SNP. It does not identify the "
                "single base or the 500-bp snATAC tile as the physical contact endpoint."
            ),
            "allele_boundary": (
                "The PLAC-seq data are wild-type. They cannot show that T>C creates, strengthens, weakens, "
                "or causally mediates the interaction."
            ),
            "cell_boundary": (
                "The assay is adult hypothalamus, not the prioritized LHA-AHN-PVH Otp Trh Glut subclass "
                "or the exact male posterior-hypothalamic IRX3+ responding cells."
            ),
        },
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")

    print("SIGNIFICANT_CONTACTS", contact["parsed_contacts"])
    print("TARGET_BIN_CONTACTS", contact["contacts_with_target_bin_anchor"])
    print("IRX3_PROMOTER_CONTACTS", contact["contacts_with_irx3_promoter_bin_anchor"])
    print("TARGET_TO_IRX3_SIGNIFICANT", contact["target_bin_to_irx3_promoter_bin_contacts"])
    print("MODEL_STAT_MATCHES", len(model["rows_matching_primary_count_expected_fdr"]))
    print("MODEL_COORDINATE_CANDIDATES", len(model["coordinate_candidates_near_target_and_irx3"]))
    print("FED_FASTED_TARGET_TO_IRX3", len(fed_fasted["target_to_irx3_rows"]))
    print("SIGNIFICANT_ROW", json.dumps(primary, sort_keys=True))
    for row in model["rows_matching_primary_count_expected_fdr"]:
        print("MODEL_STAT_MATCH", json.dumps(row, sort_keys=True))
    for row in fed_fasted["target_to_irx3_rows"]:
        print("FED_FASTED_ROW", json.dumps(row, sort_keys=True))
    print("H4_WILDTYPE_3D_SUPPORT: POSITIVE_AT_10KB_RESOLUTION")
    print("GAP-001 remains OPEN; no allele dependence or mediator is established.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
