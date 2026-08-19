#!/usr/bin/env python3
"""GAP-001 decisive cell-context screen in adult VPH GSE146692.

For a fixed candidate set emerging from prior motif/grammar passes, measure
whether each mouse ortholog transcript is detected in the same filtered adult
VPH cells as Irx3, separately by biological sample and sex.

The screen is deliberately a feasibility gate. Transcript co-detection does not
establish protein abundance, TF occupancy at rs1421085, enhancer->Irx3 control,
or mediation of the phenotype. Pooled-cell statistics are descriptive only.
"""
from __future__ import annotations

import argparse
import json
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

# Human motif labels -> mouse transcript symbols in GSE146692.
CANDIDATES = {
    "ARID5B": "Arid5b",
    "CUX1": "Cux1",
    "TET1": "Tet1",
    "KDM2B": "Kdm2b",
    "DNMT1": "Dnmt1",
    "MECP2": "Mecp2",
    "CXXC1": "Cxxc1",
}

FORBIDDEN_CAUSAL_FIELDS = {
    "cause_found", "causal_status", "materialized", "verified", "mediator_found"
}


def decode(values: Iterable[object]) -> np.ndarray:
    return np.asarray([
        x.decode("utf-8") if isinstance(x, (bytes, np.bytes_)) else str(x)
        for x in values
    ])


def load_vectors(path: Path, genes: list[str]) -> tuple[dict[str, np.ndarray], int]:
    """Load exact feature-name vectors; never silently convert a missing gene to zero."""
    with h5py.File(path, "r") as h5:
        if "matrix" not in h5:
            raise ValueError(f"{path}: missing /matrix")
        g = h5["matrix"]
        names = decode(g["features"]["name"][:])
        rows: dict[str, int] = {}
        for gene in genes:
            hits = np.flatnonzero(names == gene)
            if len(hits) != 1:
                raise ValueError(f"{path}: expected exactly one {gene}, got {len(hits)}")
            rows[gene] = int(hits[0])
        shape = tuple(int(x) for x in g["shape"][:])
        matrix = csc_matrix(
            (g["data"][:], g["indices"][:], g["indptr"][:]), shape=shape
        )
        return {
            gene: np.asarray(matrix[row, :].todense()).ravel()
            for gene, row in rows.items()
        }, shape[1]


def pct(num: int, den: int) -> float:
    return 100.0 * num / den if den else 0.0


def safe_or(a: int, b: int, c: int, d: int) -> tuple[float, float]:
    odds, p_value = fisher_exact([[a, b], [c, d]], alternative="two-sided")
    return float(odds), float(p_value)


def summarize_candidate(
    sample: str,
    sex: str,
    candidate_tf: str,
    mouse_gene: str,
    irx_counts: np.ndarray,
    candidate_counts: np.ndarray,
) -> dict[str, object]:
    if len(irx_counts) != len(candidate_counts):
        raise ValueError("Irx3 and candidate vectors differ in length")
    irx = irx_counts > 0
    cand = candidate_counts > 0
    non_irx = ~irx
    both = irx & cand

    a = int(both.sum())
    b = int((irx & ~cand).sum())
    c = int((non_irx & cand).sum())
    d = int((non_irx & ~cand).sum())
    odds, p_value = safe_or(a, b, c, d)
    irx_n = int(irx.sum())
    non_irx_n = int(non_irx.sum())
    p_irx = pct(a, irx_n)
    p_non = pct(c, non_irx_n)

    return {
        "sample": sample,
        "sex": sex,
        "candidate_tf": candidate_tf,
        "mouse_gene": mouse_gene,
        "n_cells": len(irx),
        "irx3_positive": irx_n,
        "candidate_positive": int(cand.sum()),
        "irx3_candidate_double": a,
        "candidate_pct_all": pct(int(cand.sum()), len(cand)),
        "candidate_given_irx3_pct": p_irx,
        "candidate_given_non_irx3_pct": p_non,
        "irx3_context_delta_pp": p_irx - p_non,
        "descriptive_fisher_odds_ratio": odds,
        "descriptive_fisher_p_value": p_value,
    }


def aggregate_group(frame: pd.DataFrame, group_name: str, group_value: str) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for tf, mouse_gene in CANDIDATES.items():
        sub = frame[frame["candidate_tf"] == tf]
        # The frame contains one row per cell per candidate.
        irx = sub["irx3_positive"].to_numpy(dtype=bool)
        cand = sub["candidate_positive"].to_numpy(dtype=bool)
        non_irx = ~irx
        a = int((irx & cand).sum())
        b = int((irx & ~cand).sum())
        c = int((non_irx & cand).sum())
        d = int((non_irx & ~cand).sum())
        odds, p_value = safe_or(a, b, c, d)
        p_irx = pct(a, int(irx.sum()))
        p_non = pct(c, int(non_irx.sum()))
        rows.append({
            "group": group_name,
            "value": group_value,
            "candidate_tf": tf,
            "mouse_gene": mouse_gene,
            "n_cells": len(sub),
            "irx3_positive": int(irx.sum()),
            "candidate_positive": int(cand.sum()),
            "irx3_candidate_double": a,
            "candidate_pct_all": pct(int(cand.sum()), len(sub)),
            "candidate_given_irx3_pct": p_irx,
            "candidate_given_non_irx3_pct": p_non,
            "irx3_context_delta_pp": p_irx - p_non,
            "descriptive_fisher_odds_ratio": odds,
            "descriptive_fisher_p_value": p_value,
        })
    return rows


def build_cell_frame(sample: str, sex: str, vectors: dict[str, np.ndarray]) -> pd.DataFrame:
    irx = vectors["Irx3"] > 0
    pieces = []
    for tf, gene in CANDIDATES.items():
        cand = vectors[gene] > 0
        pieces.append(pd.DataFrame({
            "sample": sample,
            "sex": sex,
            "candidate_tf": tf,
            "mouse_gene": gene,
            "irx3_positive": irx,
            "candidate_positive": cand,
        }))
    return pd.concat(pieces, ignore_index=True)


def replicate_summary(sample_df: pd.DataFrame, pooled_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    pooled = pooled_df.set_index("candidate_tf")
    for tf, mouse_gene in CANDIDATES.items():
        sub = sample_df[sample_df.candidate_tf == tf].copy()
        male = sub[sub.sex == "male"]["candidate_given_irx3_pct"]
        female = sub[sub.sex == "female"]["candidate_given_irx3_pct"]
        values = sub["candidate_given_irx3_pct"]
        pooled_row = pooled.loc[tf]
        rows.append({
            "candidate_tf": tf,
            "mouse_gene": mouse_gene,
            "replicates_with_irx3_codetection": int((sub.irx3_candidate_double > 0).sum()),
            "replicate_mean_candidate_given_irx3_pct": float(values.mean()),
            "replicate_min_candidate_given_irx3_pct": float(values.min()),
            "replicate_max_candidate_given_irx3_pct": float(values.max()),
            "male_replicate_mean_given_irx3_pct": float(male.mean()),
            "female_replicate_mean_given_irx3_pct": float(female.mean()),
            "male_minus_female_pp": float(male.mean() - female.mean()),
            "pooled_candidate_given_irx3_pct": float(pooled_row["candidate_given_irx3_pct"]),
            "pooled_candidate_given_non_irx3_pct": float(pooled_row["candidate_given_non_irx3_pct"]),
            "pooled_irx3_context_delta_pp": float(pooled_row["irx3_context_delta_pp"]),
        })
    out = pd.DataFrame(rows)
    return out.sort_values(
        ["replicates_with_irx3_codetection", "replicate_mean_candidate_given_irx3_pct"],
        ascending=[False, False],
    ).reset_index(drop=True)


def render_report(
    sample_df: pd.DataFrame,
    sex_df: pd.DataFrame,
    pooled_df: pd.DataFrame,
    repl_df: pd.DataFrame,
) -> str:
    show = repl_df.copy()
    for col in show.columns:
        if col.endswith("_pct") or col.endswith("_pp"):
            show[col] = show[col].map(lambda x: f"{float(x):.3f}")

    sex_show = sex_df[[
        "value", "candidate_tf", "irx3_positive", "irx3_candidate_double",
        "candidate_given_irx3_pct", "candidate_given_non_irx3_pct", "irx3_context_delta_pp"
    ]].copy()
    for col in ["candidate_given_irx3_pct", "candidate_given_non_irx3_pct", "irx3_context_delta_pp"]:
        sex_show[col] = sex_show[col].map(lambda x: f"{float(x):.3f}")

    sample_show = sample_df[[
        "sample", "sex", "candidate_tf", "irx3_positive", "irx3_candidate_double",
        "candidate_given_irx3_pct", "candidate_given_non_irx3_pct", "irx3_context_delta_pp"
    ]].copy()
    for col in ["candidate_given_irx3_pct", "candidate_given_non_irx3_pct", "irx3_context_delta_pp"]:
        sample_show[col] = sample_show[col].map(lambda x: f"{float(x):.3f}")

    return (
        "# GAP-001 decisive Irx3+ candidate context — adult VPH\n\n"
        "**Status: transcript co-detection feasibility screen; not TF occupancy or causality.**\n\n"
        "Candidates: `ARID5B`, `CUX1`, `TET1`, `KDM2B`, `DNMT1`, `MECP2`, `CXXC1`.\n\n"
        "Detection rule: >0 UMI in a filtered GSE146692 Cell Ranger cell. No expression cutoff is invented.\n\n"
        "## Replicate-aware shortlist\n\n"
        + show.to_markdown(index=False)
        + "\n\n## By sex (pooled cells; descriptive only)\n\n"
        + sex_show.to_markdown(index=False)
        + "\n\n## Per biological sample\n\n"
        + sample_show.to_markdown(index=False)
        + "\n\n## Interpretation guardrails\n\n"
        "- Co-detection only tests whether a candidate transcript is available in sampled adult VPH Irx3+ cells.\n"
        "- 3' scRNA-seq dropout can create false negatives; transcript presence does not prove protein abundance or isoform state.\n"
        "- Pooled-cell Fisher statistics are descriptive because cells within a biological sample are not independent biological replicates.\n"
        "- A sex difference in this screen is a prioritization signal, not a genotype-by-sex causal effect.\n"
        "- This wild-type dataset cannot establish allele-specific occupancy at rs1421085 or enhancer->Irx3 regulation.\n"
        "- Causal promotion still requires exact-allele, adult-PH-cell-specific occupancy/perturbation and target-gene evidence.\n"
    )


def run(input_root: Path, output_dir: Path) -> dict[str, object]:
    output_dir.mkdir(parents=True, exist_ok=True)
    genes = ["Irx3", *CANDIDATES.values()]
    sample_rows: list[dict[str, object]] = []
    cell_frames: list[pd.DataFrame] = []

    for sample, sex in SAMPLES.items():
        path = input_root / sample / "filtered_feature_bc_matrix.h5"
        if not path.exists():
            raise FileNotFoundError(path)
        vectors, n_cells = load_vectors(path, genes)
        if any(len(v) != n_cells for v in vectors.values()):
            raise RuntimeError(f"{sample}: vector length mismatch")
        for tf, gene in CANDIDATES.items():
            sample_rows.append(summarize_candidate(sample, sex, tf, gene, vectors["Irx3"], vectors[gene]))
        cell_frames.append(build_cell_frame(sample, sex, vectors))

    sample_df = pd.DataFrame(sample_rows)
    cells = pd.concat(cell_frames, ignore_index=True)

    pooled_rows = aggregate_group(cells, "pooled", "all")
    pooled_df = pd.DataFrame(pooled_rows)
    sex_rows: list[dict[str, object]] = []
    for sex, frame in cells.groupby("sex", sort=True):
        sex_rows.extend(aggregate_group(frame, "sex", str(sex)))
    sex_df = pd.DataFrame(sex_rows)
    repl_df = replicate_summary(sample_df, pooled_df)

    sample_df.to_csv(output_dir / "sample_candidate_context.csv", index=False)
    sex_df.to_csv(output_dir / "sex_candidate_context.csv", index=False)
    pooled_df.to_csv(output_dir / "pooled_candidate_context.csv", index=False)
    repl_df.to_csv(output_dir / "replicate_shortlist.csv", index=False)
    (output_dir / "report.md").write_text(
        render_report(sample_df, sex_df, pooled_df, repl_df), encoding="utf-8"
    )

    summary = {
        "dataset": "GSE146692",
        "samples": SAMPLES,
        "candidate_count": len(CANDIDATES),
        "candidates": CANDIDATES,
        "detection_rule": "UMI count > 0 in Cell Ranger filtered cells",
        "ranking_rule": "replicate coverage, then mean candidate detection among Irx3+ cells",
        "top_candidate_order": repl_df["candidate_tf"].tolist(),
        "claim_scope": "adult VPH Irx3+ transcript feasibility only; not occupancy, target-gene regulation, or causality",
    }
    if FORBIDDEN_CAUSAL_FIELDS & set(summary):
        raise RuntimeError("summary unexpectedly contains causal promotion fields")
    (output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return summary


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--input-root", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    args = p.parse_args()
    summary = run(args.input_root, args.output_dir)
    print((args.output_dir / "report.md").read_text(encoding="utf-8"))
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
