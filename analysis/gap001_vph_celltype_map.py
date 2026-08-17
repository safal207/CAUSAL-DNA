#!/usr/bin/env python3
"""GAP-001 T1b: map Irx3/Arid5b co-detection onto curated VPH annotations.

Input is the curated JAX VPH h5ad (`vph-2020_global-clean.h5ad`). The script
uses published cell annotations only to localize the co-expression prerequisite.
It makes no causal claim about rs1421085 or ARID5B occupancy.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

TARGETS = ("Irx3", "Arid5b")
HINTS = ("cluster", "type", "class", "annot", "sub", "region", "sex", "sample", "batch", "library")


def gene_vector(adata: ad.AnnData, gene: str) -> np.ndarray:
    sources = []
    if adata.raw is not None:
        sources.append(("raw", adata.raw))
    sources.append(("X", adata))

    for _, source in sources:
        names = pd.Index(source.var_names.astype(str))
        hits = np.flatnonzero(names.to_numpy() == gene)
        if len(hits) == 1:
            x = source[:, int(hits[0])].X
            if sparse.issparse(x):
                x = x.toarray()
            return np.asarray(x).ravel()
    raise ValueError(f"Gene {gene!r} not found uniquely in adata.raw.var_names or adata.var_names")


def candidate_columns(obs: pd.DataFrame) -> list[str]:
    selected = []
    for col in obs.columns:
        nunique = int(obs[col].nunique(dropna=True))
        if nunique < 2 or nunique > 100:
            continue
        lower = str(col).lower()
        dtype_ok = isinstance(obs[col].dtype, pd.CategoricalDtype) or obs[col].dtype == object
        if dtype_ok or any(h in lower for h in HINTS):
            selected.append(str(col))
    return selected


def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_") or "group"


def summarize_group(frame: pd.DataFrame, column: str) -> pd.DataFrame:
    rows = []
    for label, part in frame.groupby(column, observed=True, dropna=False):
        irx = part["irx3_positive"].to_numpy(dtype=bool)
        arid = part["arid5b_positive"].to_numpy(dtype=bool)
        both = irx & arid
        rows.append(
            {
                "annotation_column": column,
                "group": str(label),
                "n_cells": len(part),
                "irx3_positive": int(irx.sum()),
                "arid5b_positive": int(arid.sum()),
                "double_positive": int(both.sum()),
                "irx3_pct": 100.0 * float(irx.mean()),
                "arid5b_pct": 100.0 * float(arid.mean()),
                "double_positive_pct": 100.0 * float(both.mean()),
                "arid5b_given_irx3_pct": (100.0 * int(both.sum()) / int(irx.sum())) if irx.any() else 0.0,
            }
        )
    return pd.DataFrame(rows).sort_values(
        ["double_positive", "arid5b_given_irx3_pct", "n_cells"], ascending=[False, False, False]
    )


def markdown_table(df: pd.DataFrame, n: int = 25) -> str:
    if df.empty:
        return "_No rows._"
    cols = [
        "annotation_column", "group", "n_cells", "irx3_positive", "arid5b_positive",
        "double_positive", "irx3_pct", "arid5b_pct", "double_positive_pct",
        "arid5b_given_irx3_pct",
    ]
    x = df[cols].head(n).copy()
    for col in ["irx3_pct", "arid5b_pct", "double_positive_pct", "arid5b_given_irx3_pct"]:
        x[col] = x[col].map(lambda v: f"{v:.3f}")
    return x.to_markdown(index=False)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--h5ad", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    args = p.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    adata = ad.read_h5ad(args.h5ad)
    irx3 = gene_vector(adata, "Irx3")
    arid5b = gene_vector(adata, "Arid5b")

    obs = adata.obs.copy()
    obs["irx3_positive"] = irx3 > 0
    obs["arid5b_positive"] = arid5b > 0
    obs["double_positive"] = obs["irx3_positive"] & obs["arid5b_positive"]

    candidates = candidate_columns(obs)
    schema = {
        "n_obs": int(adata.n_obs),
        "n_vars": int(adata.n_vars),
        "obs_columns": [str(c) for c in adata.obs.columns],
        "candidate_annotation_columns": candidates,
        "candidate_values": {
            c: [str(x) for x in obs[c].dropna().unique()[:100]] for c in candidates
        },
    }
    (args.output_dir / "obs_schema.json").write_text(json.dumps(schema, indent=2) + "\n", encoding="utf-8")

    all_tables = []
    for col in candidates:
        table = summarize_group(obs, col)
        table.to_csv(args.output_dir / f"by_{safe_name(col)}.csv", index=False)
        all_tables.append(table)

    combined = pd.concat(all_tables, ignore_index=True) if all_tables else pd.DataFrame()
    if not combined.empty:
        combined = combined.sort_values(
            ["double_positive", "arid5b_given_irx3_pct", "n_cells"], ascending=[False, False, False]
        )
        combined.to_csv(args.output_dir / "all_annotation_groups.csv", index=False)

    global_summary = {
        "n_cells": int(len(obs)),
        "irx3_positive": int(obs["irx3_positive"].sum()),
        "arid5b_positive": int(obs["arid5b_positive"].sum()),
        "double_positive": int(obs["double_positive"].sum()),
        "arid5b_given_irx3_pct": (
            100.0 * int(obs["double_positive"].sum()) / int(obs["irx3_positive"].sum())
            if obs["irx3_positive"].any() else 0.0
        ),
    }
    (args.output_dir / "global_summary.json").write_text(
        json.dumps(global_summary, indent=2) + "\n", encoding="utf-8"
    )

    report = f"""# GAP-001 T1b — curated VPH cell-type map

This report maps `Irx3`/`Arid5b` detection onto annotation columns already present in the published JAX VPH h5ad. It is a localization test, **not a causal test**.

## Dataset

- observations: **{adata.n_obs:,}**
- variables: **{adata.n_vars:,}**
- candidate annotation columns: **{len(candidates)}**
- Irx3+ cells: **{global_summary['irx3_positive']:,}**
- Arid5b+ cells: **{global_summary['arid5b_positive']:,}**
- double-positive cells: **{global_summary['double_positive']:,}**
- Arid5b+ among Irx3+: **{global_summary['arid5b_given_irx3_pct']:.3f}%**

## Highest double-positive annotation groups

{markdown_table(combined, 40)}

## Guardrails

- h5ad annotations are used as published; this script does not re-cluster cells.
- Detection threshold is expression value > 0; normalization/storage details may differ from raw UMI matrices.
- A cell-type overlap can keep H1 feasible but cannot demonstrate ARID5B binding or rs1421085 causality.
- Exact-SNP, sex-aware regulatory assays remain required for a causal edge.
"""
    (args.output_dir / "report.md").write_text(report, encoding="utf-8")
    print(report)


if __name__ == "__main__":
    main()
