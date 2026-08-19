# GAP-001 broad cell-context screen

Status: **PASS AS HIGH-SENSITIVITY SCREEN — identity not resolved**

GitHub Actions run: `31988538020`
Artifact: `gap001-broad-cell-context`
Artifact digest: `sha256:2a928e1c67cfee9a3c309539fb5d1940135e35bb7090fae5405d1f8c80304469`
Head SHA: `a8bbffeef4074ef81149a8570b57ef00631d18b2`

This screen uses marker **detection signatures**, not curated cell-type labels. Cells can match more than one signature.

## Candidate subsets

| subset | signature | n_subset | positive | pct_subset |
|---|---|---:|---:|---:|
| Irx3+ | neuronal OR signature | 1,791 | 1,396 | **77.945%** |
| Irx3+ | astrocyte OR signature | 1,791 | 896 | **50.028%** |
| Irx3+/Arid5b+ | neuronal OR signature | 509 | 470 | **92.338%** |
| Irx3+/Arid5b+ | astrocyte OR signature | 509 | 200 | **39.293%** |
| Irx3+/Cux1+ | neuronal OR signature | 632 | 552 | **87.342%** |
| Irx3+/Cux1+ | astrocyte OR signature | 632 | 291 | **46.044%** |
| Irx3+/Arid5b+/Cux1+ | neuronal OR signature | 277 | 263 | **94.946%** |
| Irx3+/Arid5b+/Cux1+ | astrocyte OR signature | 277 | 107 | **38.628%** |

Tanycyte and ependymal signatures were rare in all Irx3-containing candidate subsets (~1–2%).

## Interpretation

The candidate state becomes progressively more compatible with a neuronal compartment as `Arid5b` and `Cux1` are added to `Irx3`: the high-sensitivity neuronal signature rises from ~78% in all `Irx3+` cells to ~95% in the triple-positive subset.

This does **not** establish that the native rs1421085 effect occurs in neurons. The OR-based signatures overlap heavily, and ambient RNA/dropout can distort broad marker calls. In particular, single-marker oligodendrocyte signals were implausibly broad and are therefore not used for biological inference from this screen.

A stricter multi-marker follow-up is required before cell-context prioritization.

## Guardrails

- Marker dropout can undercount a true cell type.
- OR-of-markers increases sensitivity but is not a formal classifier.
- Cells can simultaneously satisfy neuronal and astrocyte signatures.
- This does not replace the curated Mickelsen/JAX annotation mapping.
- This is a **prioritization screen**, not a causal or cell-identity proof.
