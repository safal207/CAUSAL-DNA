#!/usr/bin/env python3
"""Inspect the official JAX VPH neuronal RData without requiring an R runtime.

The JAX portal documents that its RData exports contain normalized expression,
gene information, and a dataframe with dimensional-reduction coordinates and
cluster labels. This script inventories objects readable by pyreadr and exports
small cell-metadata-like tables for downstream barcode/cluster mapping.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import pandas as pd
import pyreadr


def safe(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", name) or "unnamed"


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--rdata", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    args = p.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    result = pyreadr.read_r(str(args.rdata))
    inventory = []
    candidates = []

    for name, obj in result.items():
        row = {
            "object": str(name),
            "python_type": type(obj).__name__,
            "nrow": int(obj.shape[0]) if hasattr(obj, "shape") and len(obj.shape) >= 1 else None,
            "ncol": int(obj.shape[1]) if hasattr(obj, "shape") and len(obj.shape) >= 2 else None,
            "columns": [str(c) for c in getattr(obj, "columns", [])[:100]],
            "index_head": [str(x) for x in getattr(obj, "index", [])[:10]],
        }
        inventory.append(row)

        if isinstance(obj, pd.DataFrame):
            cols_lower = [str(c).lower() for c in obj.columns]
            keyword = any(
                any(k in c for k in ("cluster", "class", "sex", "sample", "batch", "umap", "chem", "cell"))
                for c in cols_lower
            )
            cellshape = 1_000 <= len(obj) <= 100_000 and obj.shape[1] <= 100
            if (keyword or cellshape) and obj.shape[1] <= 100:
                out = obj.copy()
                out.insert(0, ".rowname", out.index.astype(str))
                path = args.output_dir / f"candidate_{safe(str(name))}.csv"
                out.to_csv(path, index=False)
                candidates.append({"object": str(name), "path": path.name, "rows": len(out), "cols": out.shape[1]})

    (args.output_dir / "object_inventory.json").write_text(
        json.dumps(inventory, indent=2) + "\n", encoding="utf-8"
    )
    pd.DataFrame(inventory).drop(columns=["columns", "index_head"]).to_csv(
        args.output_dir / "object_inventory.csv", index=False
    )
    (args.output_dir / "metadata_candidates.json").write_text(
        json.dumps(candidates, indent=2) + "\n", encoding="utf-8"
    )

    lines = ["# VPH neuronal RData inspection", ""]
    for row in inventory:
        lines += [
            f"## `{row['object']}`",
            f"- python type: `{row['python_type']}`",
            f"- shape: `{row['nrow']} x {row['ncol']}`",
            f"- columns/head: `{', '.join(row['columns'][:30])}`",
            f"- index/head: `{', '.join(row['index_head'])}`",
            "",
        ]
    lines += ["## Exported metadata candidates", "", "```json", json.dumps(candidates, indent=2), "```", ""]
    (args.output_dir / "report.md").write_text("\n".join(lines), encoding="utf-8")

    print("\n".join(lines))


if __name__ == "__main__":
    main()
