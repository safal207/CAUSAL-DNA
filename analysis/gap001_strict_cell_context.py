#!/usr/bin/env python3
"""GAP-001 strict broad cell-context screen for GSE146692 VPH cells.

This is a refinement of the high-sensitivity OR-marker screen. It requires
multiple concordant markers for a broad cellular program and explicitly keeps
an `ambiguous` category instead of forcing every cell into one identity.

The result is still a prioritization aid, not a substitute for the published
curated cell taxonomy or a causal rs1421085 experiment.
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

PROGRAMS = {
    # Require >=2 detected genes for each program.
    "neuronal": ("Snap25", "Syt1", "Rbfox3"),
    "astrocyte": ("Aqp4", "Agt", "Slc1a2", "Gja1", "Aldoc"),
    "oligodendrocyte": ("Mbp", "Mog", "Plp1", "Cnp"),
    "opc": ("Pdgfra", "Cspg4", "Olig1", "Olig2"),
    "microglia": ("Tmem119", "P2ry12", "Cx3cr1", "Aif1"),
    "endothelial": ("Pecam1", "Cldn5", "Kdr", "Emcn"),
}
TARGETS = ("Irx3", "Arid5b", "Cux1")
ALL_GENES = tuple(dict.fromkeys(TARGETS + tuple(g for genes in PROGRAMS.values() for g in genes)))


def decode(values):
    return np.asarray([x.decode() if isinstance(x, (bytes, np.bytes_)) else str(x) for x in values])


def load(path: Path):
    with h5py.File(path, "r") as h5:
        g = h5["matrix"]
        names = decode(g["features"]["name"][:])
        shape = tuple(int(x) for x in g["shape"][:])
        matrix = csc_matrix((g["data"][:], g["indices"][:], g["indptr"][:]), shape=shape)
        out = {}
        availability = {}
        for gene in ALL_GENES:
            hits = np.flatnonzero(names == gene)
            availability[gene] = int(len(hits))
            out[gene] = (
                np.asarray(matrix[int(hits[0]), :].todense()).ravel()
                if len(hits) == 1 else np.zeros(shape[1], dtype=np.int64)
            )
        return out, availability


def program_call(g, genes, min_detected=2):
    detected = np.vstack([g[x] > 0 for x in genes])
    return detected.sum(axis=0) >= min_detected


def pct(n, d):
    return 100.0 * n / d if d else 0.0


def build_cells(sample, sex, g):
    n = len(g["Irx3"])
    f = pd.DataFrame({"sample": [sample] * n, "sex": [sex] * n})
    for gene in TARGETS:
        f[gene.lower()] = g[gene] > 0
    for name, genes in PROGRAMS.items():
        f[name] = program_call(g, genes)

    broad = []
    for row in f[list(PROGRAMS)].itertuples(index=False, name=None):
        on = [name for name, value in zip(PROGRAMS, row) if value]
        broad.append(on[0] if len(on) == 1 else ("ambiguous" if len(on) > 1 else "unresolved"))
    f["broad_context"] = broad
    return f


def subset_rows(cells, subset_name, mask):
    f = cells.loc[mask]
    rows = []
    for context in list(PROGRAMS) + ["ambiguous", "unresolved"]:
        n = int((f.broad_context == context).sum())
        rows.append({
            "subset": subset_name,
            "broad_context": context,
            "n_subset": len(f),
            "n_cells": n,
            "pct_subset": pct(n, len(f)),
        })
    # Also retain non-exclusive program-call rates, useful for diagnosing ambiguity.
    for name in PROGRAMS:
        n = int(f[name].sum())
        rows.append({
            "subset": subset_name,
            "broad_context": f"{name}_program_any",
            "n_subset": len(f),
            "n_cells": n,
            "pct_subset": pct(n, len(f)),
        })
    return rows


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input-root", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    args = p.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    frames = []
    availability_rows = []
    for sample, sex in SAMPLES.items():
        g, availability = load(args.input_root / sample / "filtered_feature_bc_matrix.h5")
        availability_rows.append({"sample": sample, **availability})
        frames.append(build_cells(sample, sex, g))
    cells = pd.concat(frames, ignore_index=True)
    pd.DataFrame(availability_rows).to_csv(args.output_dir / "feature_availability.csv", index=False)

    subsets = {
        "irx3_positive": cells.irx3,
        "irx3_arid5b": cells.irx3 & cells.arid5b,
        "irx3_cux1": cells.irx3 & cells.cux1,
        "irx3_arid5b_cux1": cells.irx3 & cells.arid5b & cells.cux1,
    }
    pooled_rows = []
    for name, mask in subsets.items():
        pooled_rows.extend(subset_rows(cells, name, mask.to_numpy(bool)))
    pooled = pd.DataFrame(pooled_rows)
    pooled.to_csv(args.output_dir / "pooled_strict_context.csv", index=False)

    sex_rows = []
    for sex, sf in cells.groupby("sex"):
        sex_subsets = {
            "irx3_positive": sf.irx3,
            "irx3_arid5b": sf.irx3 & sf.arid5b,
            "irx3_cux1": sf.irx3 & sf.cux1,
            "irx3_arid5b_cux1": sf.irx3 & sf.arid5b & sf.cux1,
        }
        for name, mask in sex_subsets.items():
            for row in subset_rows(sf, name, mask.to_numpy(bool)):
                row["sex"] = sex
                sex_rows.append(row)
    sex_df = pd.DataFrame(sex_rows)
    sex_df.to_csv(args.output_dir / "sex_strict_context.csv", index=False)

    key_contexts = ["neuronal", "astrocyte", "ambiguous", "unresolved",
                    "neuronal_program_any", "astrocyte_program_any"]
    focus = pooled[pooled.broad_context.isin(key_contexts)].copy()

    report = "# GAP-001 strict broad cell-context screen\n\n"
    report += "Each broad program requires **at least 2 concordant markers**. Exclusive labels are assigned only when exactly one program passes; overlaps are reported as `ambiguous`.\n\n"
    report += "## Key result\n\n" + focus.to_markdown(index=False) + "\n\n"
    report += "## Marker programs\n\n"
    for name, genes in PROGRAMS.items():
        report += f"- **{name}**: {', '.join(genes)} (>=2 detected)\n"
    report += "\n## Guardrails\n\n"
    report += "- This is a deliberately conservative broad-context screen, not a trained cell classifier.\n"
    report += "- scRNA dropout can move true cells into `unresolved`; ambient RNA can still create false positives.\n"
    report += "- Broad programs do not identify the exact PH neuronal subtype implicated by Sullivan et al.\n"
    report += "- No result here establishes the cell in which rs1421085 exerts its native causal effect.\n"
    report += "- The purpose is to prioritize exact-allele, cell-resolved ATAC/occupancy/contact experiments.\n"
    (args.output_dir / "report.md").write_text(report, encoding="utf-8")
    print(report)


if __name__ == "__main__":
    main()
