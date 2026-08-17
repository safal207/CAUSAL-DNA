#!/usr/bin/env python3
"""GAP-001 motif-candidate presence screen in GSE146692 VPH cells.

The original rs1421085 motif analysis reported disruption of motifs for
ARID-family factors, LHX6 and NKX6-3. This script asks a deliberately narrow
question: are transcripts for these candidate interpreters present in the same
VPH cells as Irx3?

Presence/co-detection is only a feasibility filter. It does not establish DNA
binding, allele-specific occupancy, or causality.
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
TARGET = "Irx3"
CANDIDATES = ("Arid5b", "Cux1", "Lhx6", "Nkx6-3")


def decode(values):
    return np.asarray([x.decode() if isinstance(x, (bytes, np.bytes_)) else str(x) for x in values])


def load_genes(path: Path):
    with h5py.File(path, "r") as h5:
        g = h5["matrix"]
        names = decode(g["features"]["name"][:])
        shape = tuple(int(x) for x in g["shape"][:])
        matrix = csc_matrix((g["data"][:], g["indices"][:], g["indptr"][:]), shape=shape)
        out = {}
        availability = {}
        for gene in (TARGET,) + CANDIDATES:
            hits = np.flatnonzero(names == gene)
            availability[gene] = int(len(hits))
            if len(hits) == 1:
                out[gene] = np.asarray(matrix[int(hits[0]), :].todense()).ravel()
            else:
                out[gene] = np.zeros(shape[1], dtype=np.int64)
        return out, availability


def pct(num: int, den: int) -> float:
    return 100.0 * num / den if den else 0.0


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--input-root", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    args = p.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    cell_frames = []
    availability_rows = []

    for sample, sex in SAMPLES.items():
        genes, availability = load_genes(args.input_root / sample / "filtered_feature_bc_matrix.h5")
        availability_rows.append({"sample": sample, **availability})
        frame = pd.DataFrame({"sample": sample, "sex": sex, TARGET: genes[TARGET] > 0})
        for candidate in CANDIDATES:
            frame[candidate] = genes[candidate] > 0
        cell_frames.append(frame)

        irx = frame[TARGET].to_numpy(bool)
        for candidate in CANDIDATES:
            cand = frame[candidate].to_numpy(bool)
            both = irx & cand
            rows.append({
                "sample": sample,
                "sex": sex,
                "candidate": candidate,
                "n_cells": len(frame),
                "irx3_positive": int(irx.sum()),
                "candidate_positive": int(cand.sum()),
                "double_positive": int(both.sum()),
                "candidate_pct_all": pct(int(cand.sum()), len(frame)),
                "candidate_given_irx3_pct": pct(int(both.sum()), int(irx.sum())),
            })

    cells = pd.concat(cell_frames, ignore_index=True)
    per_sample = pd.DataFrame(rows)
    per_sample.to_csv(args.output_dir / "sample_candidate_summary.csv", index=False)
    pd.DataFrame(availability_rows).to_csv(args.output_dir / "feature_availability.csv", index=False)

    pooled_rows = []
    irx = cells[TARGET].to_numpy(bool)
    for candidate in CANDIDATES:
        cand = cells[candidate].to_numpy(bool)
        both = irx & cand
        pooled_rows.append({
            "candidate": candidate,
            "n_cells": len(cells),
            "irx3_positive": int(irx.sum()),
            "candidate_positive": int(cand.sum()),
            "double_positive": int(both.sum()),
            "candidate_pct_all": pct(int(cand.sum()), len(cells)),
            "candidate_given_irx3_pct": pct(int(both.sum()), int(irx.sum())),
        })
    pooled = pd.DataFrame(pooled_rows)
    pooled.to_csv(args.output_dir / "pooled_candidate_summary.csv", index=False)

    sex_rows = []
    for sex, sf in cells.groupby("sex"):
        irx_s = sf[TARGET].to_numpy(bool)
        for candidate in CANDIDATES:
            cand = sf[candidate].to_numpy(bool)
            both = irx_s & cand
            sex_rows.append({
                "sex": sex,
                "candidate": candidate,
                "n_cells": len(sf),
                "irx3_positive": int(irx_s.sum()),
                "candidate_positive": int(cand.sum()),
                "double_positive": int(both.sum()),
                "candidate_pct_all": pct(int(cand.sum()), len(sf)),
                "candidate_given_irx3_pct": pct(int(both.sum()), int(irx_s.sum())),
            })
    sex_df = pd.DataFrame(sex_rows)
    sex_df.to_csv(args.output_dir / "sex_candidate_summary.csv", index=False)

    report = "# GAP-001 motif-candidate presence screen\n\n"
    report += "This is a **feasibility screen only**. Co-detection cannot establish TF occupancy or a causal edge.\n\n"
    report += "## Pooled VPH cells\n\n" + pooled.to_markdown(index=False) + "\n\n"
    report += "## By sex\n\n" + sex_df.to_markdown(index=False) + "\n\n"
    report += "## Feature availability\n\n" + pd.DataFrame(availability_rows).to_markdown(index=False) + "\n\n"
    report += "## Interpretation rules\n\n"
    report += "- Candidate absent from the feature matrix or nearly absent from Irx3+ cells -> demote as a direct adult-VPH interpreter.\n"
    report += "- Candidate co-detected with Irx3 -> prerequisite survives, but no causal upgrade.\n"
    report += "- Sex imbalance is descriptive because sex is confounded with batch/10x chemistry in this dataset.\n"
    report += "- Decisive evidence requires allele-specific occupancy/perturbation in the relevant PH cell type.\n"
    (args.output_dir / "report.md").write_text(report, encoding="utf-8")
    print(report)


if __name__ == "__main__":
    main()
