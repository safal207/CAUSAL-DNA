#!/usr/bin/env python3
"""GAP-001 broad cell-context screen in GSE146692 VPH raw counts.

Goal: localize Irx3-containing candidate regulatory states to broad cell
contexts without pretending to reconstruct the curated Mickelsen taxonomy.

Published VPH markers used here include:
- astrocyte: Aqp4, Agt
- tanycyte: Rax
- ependymal: Ccdc153

General neuronal markers Snap25/Syt1/Rbfox3 are used as a neuronal-presence
signature. Marker detection is a coarse feasibility/localization screen only.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy.sparse import csc_matrix

SAMPLES = {
    "AJ18003": "male",
    "AJ18004": "female",
    "AJ19001": "female",
    "AJ19002": "male",
}
GENES = (
    "Irx3", "Arid5b", "Cux1",
    "Snap25", "Syt1", "Rbfox3",
    "Aqp4", "Agt",
    "Rax", "Ccdc153",
    "Tmem119", "Pecam1", "Pdgfra", "Mbp",
)


def decode(values):
    return np.asarray([x.decode() if isinstance(x, (bytes, np.bytes_)) else str(x) for x in values])


def load(path: Path):
    with h5py.File(path, "r") as h5:
        g = h5["matrix"]
        names = decode(g["features"]["name"][:])
        shape = tuple(int(x) for x in g["shape"][:])
        matrix = csc_matrix((g["data"][:], g["indices"][:], g["indptr"][:]), shape=shape)
        out = {}
        for gene in GENES:
            hits = np.flatnonzero(names == gene)
            out[gene] = (
                np.asarray(matrix[int(hits[0]), :].todense()).ravel()
                if len(hits) == 1 else np.zeros(shape[1], dtype=np.int64)
            )
        return out


def any_pos(g, names):
    return np.logical_or.reduce([g[x] > 0 for x in names])


def pct(n, d):
    return 100.0 * n / d if d else 0.0


def make_frame(sample, sex, g):
    f = pd.DataFrame({"sample": sample, "sex": sex})
    n = len(g["Irx3"])
    f = pd.DataFrame({"sample": [sample] * n, "sex": [sex] * n})
    f["irx3"] = g["Irx3"] > 0
    f["arid5b"] = g["Arid5b"] > 0
    f["cux1"] = g["Cux1"] > 0
    f["neuronal_signature"] = any_pos(g, ("Snap25", "Syt1", "Rbfox3"))
    f["astrocyte_signature"] = any_pos(g, ("Aqp4", "Agt"))
    f["tanycyte_signature"] = g["Rax"] > 0
    f["ependymal_signature"] = g["Ccdc153"] > 0
    f["microglia_signature"] = g["Tmem119"] > 0
    f["endothelial_signature"] = g["Pecam1"] > 0
    f["opc_signature"] = g["Pdgfra"] > 0
    f["oligodendrocyte_signature"] = g["Mbp"] > 0
    return f


def summarize_subset(cells, subset_name, mask):
    f = cells.loc[mask]
    rows = []
    signatures = [c for c in cells.columns if c.endswith("_signature")]
    for sig in signatures:
        count = int(f[sig].sum())
        rows.append({
            "subset": subset_name,
            "signature": sig,
            "n_subset": len(f),
            "positive": count,
            "pct_subset": pct(count, len(f)),
        })
    return rows


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input-root", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    args = p.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    frames = []
    for sample, sex in SAMPLES.items():
        frames.append(make_frame(sample, sex, load(args.input_root / sample / "filtered_feature_bc_matrix.h5")))
    cells = pd.concat(frames, ignore_index=True)

    subsets = {
        "all_cells": np.ones(len(cells), dtype=bool),
        "irx3_positive": cells.irx3,
        "irx3_arid5b": cells.irx3 & cells.arid5b,
        "irx3_cux1": cells.irx3 & cells.cux1,
        "irx3_arid5b_cux1": cells.irx3 & cells.arid5b & cells.cux1,
    }

    rows = []
    for name, mask in subsets.items():
        rows.extend(summarize_subset(cells, name, mask.to_numpy(bool) if hasattr(mask, "to_numpy") else mask))
    pooled = pd.DataFrame(rows)
    pooled.to_csv(args.output_dir / "pooled_context_summary.csv", index=False)

    sex_rows = []
    for sex, sf in cells.groupby("sex"):
        sex_subsets = {
            "irx3_positive": sf.irx3,
            "irx3_arid5b": sf.irx3 & sf.arid5b,
            "irx3_cux1": sf.irx3 & sf.cux1,
            "irx3_arid5b_cux1": sf.irx3 & sf.arid5b & sf.cux1,
        }
        for name, mask in sex_subsets.items():
            for row in summarize_subset(sf, name, mask.to_numpy(bool)):
                row["sex"] = sex
                sex_rows.append(row)
    sex_df = pd.DataFrame(sex_rows)
    sex_df.to_csv(args.output_dir / "sex_context_summary.csv", index=False)

    focus = pooled[
        (pooled.subset != "all_cells") &
        pooled.signature.isin(["neuronal_signature", "astrocyte_signature", "tanycyte_signature", "ependymal_signature"])
    ].copy()

    report = "# GAP-001 broad cell-context screen\n\n"
    report += "This screen uses marker **detection signatures**, not curated cell-type labels. Cells can match more than one signature.\n\n"
    report += "## Candidate subsets\n\n" + focus.to_markdown(index=False) + "\n\n"
    report += "## Guardrails\n\n"
    report += "- Marker dropout can undercount a true cell type.\n"
    report += "- OR-of-markers increases sensitivity but is not a formal classifier.\n"
    report += "- This does not replace the curated JAX annotation mapping.\n"
    report += "- The purpose is to prioritize which broad cellular compartment deserves exact-allele chromatin/occupancy experiments.\n"
    (args.output_dir / "report.md").write_text(report, encoding="utf-8")
    print(report)


if __name__ == "__main__":
    main()
