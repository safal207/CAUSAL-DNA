#!/usr/bin/env python3
"""Inspect one GSE246791 SnapATAC2 h5ad without loading the full matrix.

The goal is to identify the persisted 500-bp tile matrix and coordinate encoding
so a later locus query can be implemented with h5py against only the required
rows/columns. This script prints metadata only and makes no biological claim.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import h5py


def simple(value):
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    if hasattr(value, "tolist"):
        value = value.tolist()
    if isinstance(value, (list, tuple)):
        return [simple(v) for v in value[:20]]
    return value


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("h5ad", type=Path)
    parser.add_argument("--max-depth", type=int, default=4)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    records = []
    with h5py.File(args.h5ad, "r") as f:
        def visit(name, obj):
            depth = name.count("/") + 1 if name else 0
            if depth > args.max_depth:
                return
            record = {
                "path": "/" + name,
                "kind": "dataset" if isinstance(obj, h5py.Dataset) else "group",
                "attrs": {str(k): simple(v) for k, v in obj.attrs.items()},
            }
            if isinstance(obj, h5py.Dataset):
                record["shape"] = list(obj.shape)
                record["dtype"] = str(obj.dtype)
                if obj.ndim == 1 and obj.shape[0] and obj.shape[0] <= 20:
                    record["preview"] = simple(obj[:])
            records.append(record)
        f.visititems(visit)

    result = {
        "file": args.h5ad.name,
        "inspection_only": True,
        "records": records,
    }
    text = json.dumps(result, indent=2, sort_keys=True)
    print(text)
    if args.output:
        args.output.write_text(text + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
