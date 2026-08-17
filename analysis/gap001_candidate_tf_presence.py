#!/usr/bin/env python3
"""GAP-001 P1: candidate-TF presence in GSE146692 raw VPH cells.

Tests whether Cux1 is detected in the same VPH cells as Irx3 and in the
Irx3+/Arid5b+ subset. This is a prerequisite test only; 3' scRNA-seq cannot
resolve functional CUX1 protein forms or prove site occupancy.
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
GENES = ("Irx3", "Arid5b", "Cux1")


def decode(values):
    return np.asarray([x.decode() if isinstance(x, (bytes, np.bytes_)) else str(x) for x in values])


def load(path: Path) -> dict[str, np.ndarray]:
    with h5py.File(path, "r") as h5:
        g = h5["matrix"]
        names = decode(g["features"]["name"][:])
        idx = {}
        for gene in GENES:
            hits = np.flatnonzero(names == gene)
            if len(hits) != 1:
                raise ValueError(f"{path}: expected one {gene}, got {len(hits)}")
            idx[gene] = int(hits[0])
        shape = tuple(int(x) for x in g["shape"][:])
        matrix = csc_matrix((g["data"][:], g["indices"][:], g["indptr"][:]), shape=shape)
        return {gene: np.asarray(matrix[row, :].todense()).ravel() for gene, row in idx.items()}


def pct(a: int, b: int) -> float:
    return 100.0 * a / b if b else 0.0


def summarize(sample: str, sex: str, x: dict[str, np.ndarray]) -> dict:
    irx = x["Irx3"] > 0
    arid = x["Arid5b"] > 0
    cux = x["Cux1"] > 0
    double = irx & arid
    triple = double & cux
    irx_cux = irx & cux
    return {
        "sample": sample,
        "sex": sex,
        "n_cells": len(irx),
        "irx3_positive": int(irx.sum()),
        "arid5b_positive": int(arid.sum()),
        "cux1_positive": int(cux.sum()),
        "irx3_arid5b_double": int(double.sum()),
        "irx3_cux1_double": int(irx_cux.sum()),
        "irx3_arid5b_cux1_triple": int(triple.sum()),
        "cux1_pct_all": pct(int(cux.sum()), len(cux)),
        "cux1_given_irx3_pct": pct(int(irx_cux.sum()), int(irx.sum())),
        "cux1_given_irx3_arid5b_pct": pct(int(triple.sum()), int(double.sum())),
    }


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--input-root", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    args = p.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    cell_frames = []
    for sample, sex in SAMPLES.items():
        x = load(args.input_root / sample / "filtered_feature_bc_matrix.h5")
        rows.append(summarize(sample, sex, x))
        cell_frames.append(pd.DataFrame({
            "sample": sample,
            "sex": sex,
            "irx3": x["Irx3"] > 0,
            "arid5b": x["Arid5b"] > 0,
            "cux1": x["Cux1"] > 0,
        }))

    sample_df = pd.DataFrame(rows)
    sample_df.to_csv(args.output_dir / "sample_summary.csv", index=False)

    cells = pd.concat(cell_frames, ignore_index=True)
    sex_rows = []
    for sex, f in cells.groupby("sex"):
        irx = f.irx3.to_numpy(bool)
        arid = f.arid5b.to_numpy(bool)
        cux = f.cux1.to_numpy(bool)
        da = irx & arid
        sex_rows.append({
            "sex": sex,
            "n_cells": len(f),
            "irx3_positive": int(irx.sum()),
            "cux1_positive": int(cux.sum()),
            "irx3_cux1_double": int((irx & cux).sum()),
            "irx3_arid5b_double": int(da.sum()),
            "triple_positive": int((da & cux).sum()),
            "cux1_pct_all": pct(int(cux.sum()), len(f)),
            "cux1_given_irx3_pct": pct(int((irx & cux).sum()), int(irx.sum())),
            "cux1_given_irx3_arid5b_pct": pct(int((da & cux).sum()), int(da.sum())),
        })
    sex_df = pd.DataFrame(sex_rows)
    sex_df.to_csv(args.output_dir / "sex_summary.csv", index=False)

    all_irx = cells.irx3.to_numpy(bool)
    all_arid = cells.arid5b.to_numpy(bool)
    all_cux = cells.cux1.to_numpy(bool)
    all_da = all_irx & all_arid
    pooled = {
        "n_cells": len(cells),
        "irx3_positive": int(all_irx.sum()),
        "arid5b_positive": int(all_arid.sum()),
        "cux1_positive": int(all_cux.sum()),
        "irx3_cux1_double": int((all_irx & all_cux).sum()),
        "irx3_arid5b_double": int(all_da.sum()),
        "triple_positive": int((all_da & all_cux).sum()),
        "cux1_pct_all": pct(int(all_cux.sum()), len(cells)),
        "cux1_given_irx3_pct": pct(int((all_irx & all_cux).sum()), int(all_irx.sum())),
        "cux1_given_irx3_arid5b_pct": pct(int((all_da & all_cux).sum()), int(all_da.sum())),
    }

    report = "# GAP-001 P1 — Cux1 candidate-TF presence\n\n"
    report += "This is a **presence prerequisite**, not evidence of rs1421085 occupancy or causality.\n\n"
    report += "## Pooled\n\n"
    for k, v in pooled.items():
        report += f"- {k}: **{v:.3f}**\n" if isinstance(v, float) else f"- {k}: **{v}**\n"
    report += "\n## Per sample\n\n" + sample_df.to_markdown(index=False) + "\n"
    report += "\n## By sex\n\n" + sex_df.to_markdown(index=False) + "\n"
    report += "\n## Guardrails\n\n"
    report += "- Cux1 transcript detection cannot distinguish P110/P200 functional protein forms.\n"
    report += "- Dropout can produce false-negative transcript detection.\n"
    report += "- A positive result only preserves H2a; it does not establish binding at rs1421085.\n"
    report += "- A decisive test requires exact-allele, PH-cell-specific CUX1/ARID5B occupancy and perturbation.\n"
    (args.output_dir / "report.md").write_text(report, encoding="utf-8")
    print(report)


if __name__ == "__main__":
    main()
