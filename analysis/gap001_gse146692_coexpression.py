#!/usr/bin/env python3
"""GAP-001 T1: test Irx3/Arid5b co-expression in GSE146692 VPH scRNA-seq.

This script intentionally answers a narrow prerequisite question:

    Are Irx3 and Arid5b detected in the same ventral-posterior-hypothalamus
    cells, and is that pattern different between male and female samples?

It does NOT establish rs1421085 causality, ARID5B occupancy, or an enhancer
mechanism. A positive result keeps H1 feasible; a negative/rare result weakens
H1 and motivates alternative-TF/chromatin hypotheses.

Expected input: Cell Ranger filtered_feature_bc_matrix.h5 files from
TheJacksonLaboratory/ventral-posterior-hypothalamus-scrnaseq for samples:
AJ18003, AJ18004, AJ19001, AJ19002.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

import h5py
import numpy as np
import pandas as pd
from scipy.sparse import csc_matrix
from scipy.stats import fisher_exact


SAMPLES = {
    "AJ18003": "male",
    "AJ18004": "female",
    "AJ19001": "female",
    "AJ19002": "male",
}

TARGET_GENES = ("Irx3", "Arid5b")


@dataclass
class SampleResult:
    sample: str
    sex: str
    n_cells: int
    irx3_positive: int
    arid5b_positive: int
    double_positive: int
    irx3_pct: float
    arid5b_pct: float
    double_positive_pct: float
    arid5b_given_irx3_pct: float
    irx3_given_arid5b_pct: float
    fisher_odds_ratio: float
    fisher_p_value: float


def _decode(values: Iterable[bytes]) -> np.ndarray:
    return np.asarray(
        [v.decode("utf-8") if isinstance(v, (bytes, np.bytes_)) else str(v) for v in values]
    )


def load_target_expression(path: Path) -> tuple[np.ndarray, np.ndarray, int]:
    """Return per-cell UMI counts for Irx3 and Arid5b from a Cell Ranger H5."""
    with h5py.File(path, "r") as h5:
        if "matrix" not in h5:
            raise ValueError(f"{path}: missing Cell Ranger /matrix group")

        group = h5["matrix"]
        features = group["features"]
        names = _decode(features["name"][:])

        hits = {gene: np.flatnonzero(names == gene) for gene in TARGET_GENES}
        missing = [gene for gene, idx in hits.items() if len(idx) == 0]
        duplicate = [gene for gene, idx in hits.items() if len(idx) > 1]
        if missing:
            raise ValueError(f"{path}: target gene(s) absent: {', '.join(missing)}")
        if duplicate:
            raise ValueError(f"{path}: target gene(s) duplicated: {', '.join(duplicate)}")

        shape = tuple(int(x) for x in group["shape"][:])
        matrix = csc_matrix(
            (
                group["data"][:],
                group["indices"][:],
                group["indptr"][:],
            ),
            shape=shape,
        )
        n_cells = shape[1]

        irx3 = np.asarray(matrix[int(hits["Irx3"][0]), :].todense()).ravel()
        arid5b = np.asarray(matrix[int(hits["Arid5b"][0]), :].todense()).ravel()

    return irx3, arid5b, n_cells


def summarize(sample: str, sex: str, irx3: np.ndarray, arid5b: np.ndarray) -> SampleResult:
    irx3_pos = irx3 > 0
    arid5b_pos = arid5b > 0
    both = irx3_pos & arid5b_pos

    a = int(both.sum())
    b = int((irx3_pos & ~arid5b_pos).sum())
    c = int((~irx3_pos & arid5b_pos).sum())
    d = int((~irx3_pos & ~arid5b_pos).sum())
    n = len(irx3_pos)

    odds_ratio, p_value = fisher_exact([[a, b], [c, d]], alternative="two-sided")

    return SampleResult(
        sample=sample,
        sex=sex,
        n_cells=n,
        irx3_positive=int(irx3_pos.sum()),
        arid5b_positive=int(arid5b_pos.sum()),
        double_positive=a,
        irx3_pct=100.0 * irx3_pos.mean(),
        arid5b_pct=100.0 * arid5b_pos.mean(),
        double_positive_pct=100.0 * both.mean(),
        arid5b_given_irx3_pct=(100.0 * a / int(irx3_pos.sum())) if irx3_pos.any() else 0.0,
        irx3_given_arid5b_pct=(100.0 * a / int(arid5b_pos.sum())) if arid5b_pos.any() else 0.0,
        fisher_odds_ratio=float(odds_ratio),
        fisher_p_value=float(p_value),
    )


def aggregate_by_sex(cell_tables: list[pd.DataFrame]) -> pd.DataFrame:
    cells = pd.concat(cell_tables, ignore_index=True)
    rows = []
    for sex, frame in cells.groupby("sex", sort=True):
        irx3 = frame["irx3_positive"].to_numpy(dtype=bool)
        arid5b = frame["arid5b_positive"].to_numpy(dtype=bool)
        both = irx3 & arid5b
        a = int(both.sum())
        b = int((irx3 & ~arid5b).sum())
        c = int((~irx3 & arid5b).sum())
        d = int((~irx3 & ~arid5b).sum())
        odds, p = fisher_exact([[a, b], [c, d]], alternative="two-sided")
        rows.append(
            {
                "sex": sex,
                "n_cells": len(frame),
                "irx3_positive": int(irx3.sum()),
                "arid5b_positive": int(arid5b.sum()),
                "double_positive": a,
                "irx3_pct": 100.0 * irx3.mean(),
                "arid5b_pct": 100.0 * arid5b.mean(),
                "double_positive_pct": 100.0 * both.mean(),
                "arid5b_given_irx3_pct": 100.0 * a / int(irx3.sum()) if irx3.any() else 0.0,
                "fisher_odds_ratio": float(odds),
                "fisher_p_value": float(p),
            }
        )
    return pd.DataFrame(rows)


def render_markdown(sample_df: pd.DataFrame, sex_df: pd.DataFrame) -> str:
    def table(df: pd.DataFrame, columns: list[str]) -> str:
        pretty = df[columns].copy()
        for col in pretty.columns:
            if col.endswith("_pct") or col in {"fisher_odds_ratio", "fisher_p_value"}:
                pretty[col] = pretty[col].map(lambda x: f"{x:.6g}")
        return pretty.to_markdown(index=False)

    sample_cols = [
        "sample", "sex", "n_cells", "irx3_positive", "arid5b_positive",
        "double_positive", "irx3_pct", "arid5b_pct", "double_positive_pct",
        "arid5b_given_irx3_pct", "fisher_odds_ratio", "fisher_p_value",
    ]
    sex_cols = [
        "sex", "n_cells", "irx3_positive", "arid5b_positive", "double_positive",
        "irx3_pct", "arid5b_pct", "double_positive_pct",
        "arid5b_given_irx3_pct", "fisher_odds_ratio", "fisher_p_value",
    ]

    return f"""# GAP-001 T1 — GSE146692 Irx3/Arid5b co-expression\n\n"
"This report tests **cell-level feasibility**, not causality. Detection means ≥1 UMI in a filtered Cell Ranger cell.\n\n"
"## Per-sample\n\n{table(sample_df, sample_cols)}\n\n"
"## By sex (pooled cells)\n\n{table(sex_df, sex_cols)}\n\n"
"## Interpretation guardrails\n\n"
"- Co-expression supports only the prerequisite that ARID5B and IRX3 can coexist in sampled VPH cells.\n"
"- Dropout can make true co-expression look rarer than it is.\n"
"- Pooled-cell Fisher tests are descriptive and do not replace biological-replicate inference.\n"
"- This wild-type dataset cannot test rs1421085 genotype effects or ARID5B binding at the variant.\n"
"- A causal H1 test still requires allele/context-specific occupancy or perturbation evidence.\n"
"""


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input-root",
        type=Path,
        required=True,
        help="Directory containing <sample>/filtered_feature_bc_matrix.h5",
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    sample_results: list[SampleResult] = []
    cell_tables: list[pd.DataFrame] = []

    for sample, sex in SAMPLES.items():
        path = args.input_root / sample / "filtered_feature_bc_matrix.h5"
        if not path.exists():
            raise FileNotFoundError(path)
        irx3, arid5b, n_cells = load_target_expression(path)
        if len(irx3) != n_cells or len(arid5b) != n_cells:
            raise RuntimeError(f"{sample}: expression length does not match matrix cell count")

        result = summarize(sample, sex, irx3, arid5b)
        sample_results.append(result)
        cell_tables.append(
            pd.DataFrame(
                {
                    "sample": sample,
                    "sex": sex,
                    "irx3_positive": irx3 > 0,
                    "arid5b_positive": arid5b > 0,
                }
            )
        )

    sample_df = pd.DataFrame([asdict(x) for x in sample_results])
    sex_df = aggregate_by_sex(cell_tables)

    sample_df.to_csv(args.output_dir / "sample_summary.csv", index=False)
    sex_df.to_csv(args.output_dir / "sex_summary.csv", index=False)
    (args.output_dir / "report.md").write_text(render_markdown(sample_df, sex_df), encoding="utf-8")
    (args.output_dir / "provenance.json").write_text(
        json.dumps(
            {
                "dataset": "GSE146692",
                "test": "GAP-001/T1",
                "genes": list(TARGET_GENES),
                "samples": SAMPLES,
                "detection_rule": "UMI count > 0 in Cell Ranger filtered cells",
                "causal_claim": False,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )

    print(sample_df.to_string(index=False))
    print("\nBy sex:\n", sex_df.to_string(index=False))


if __name__ == "__main__":
    main()
