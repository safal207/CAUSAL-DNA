# GAP-001 T1 — GSE146692 Irx3/Arid5b co-expression

Status: **PASS — prerequisite co-expression is present; causality remains unresolved**

Generated from public Cell Ranger filtered matrices for GSE146692 using `analysis/gap001_gse146692_coexpression.py`.

GitHub Actions run: `31987734506`
Artifact: `gap001-gse146692-coexpression`
Artifact digest: `sha256:0fa9620a6c669a71ecc946d0bbf561cce2175543bb0d8bc03ec5cd4c2d6af957`

## Question

Does a basic cell-level prerequisite for H1 hold?

> Are `Irx3` and `Arid5b` detectable in the same ventral-posterior-hypothalamus cells?

Detection rule: **UMI count > 0** in a Cell Ranger filtered cell.

## Per-sample results

| sample | sex | cells | Irx3+ | Arid5b+ | double+ | Irx3+ % | Arid5b+ % | double+ % | Arid5b+ among Irx3+ | Fisher OR | Fisher p |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| AJ18003 | male | 11,552 | 431 | 1,520 | 109 | 3.731% | 13.158% | 0.944% | 25.290% | 2.330 | 4.04e-12 |
| AJ18004 | female | 9,389 | 260 | 1,047 | 37 | 2.769% | 11.151% | 0.394% | 14.231% | 1.334 | 0.110 |
| AJ19001 | female | 10,623 | 561 | 2,187 | 193 | 5.281% | 20.587% | 1.817% | 34.403% | 2.122 | 5.28e-15 |
| AJ19002 | male | 12,585 | 539 | 2,378 | 170 | 4.283% | 18.896% | 1.351% | 31.540% | 2.053 | 7.58e-13 |

## Pooled descriptive result

Across all four samples:

- cells: **44,149**
- `Irx3+`: **1,791 (4.057%)**
- `Arid5b+`: **7,132 (16.154%)**
- double-positive: **509 (1.153%)**
- `P(Arid5b+ | Irx3+)`: **28.420%**

By sex, pooled cells give:

| sex | cells | Irx3+ | Arid5b+ | double+ | Arid5b+ among Irx3+ | Fisher OR | Fisher p |
|---|---:|---:|---:|---:|---:|---:|---:|
| female | 20,012 | 821 | 3,234 | 230 | 28.015% | 2.097 | 1.82e-18 |
| male | 24,137 | 970 | 3,898 | 279 | 28.763% | 2.181 | 7.15e-24 |

## Interpretation

### H1 — ARID5B reuse

**Update: `FEASIBLE_BUT_UNPROVEN`**.

The prerequisite is satisfied: `Arid5b` and `Irx3` are co-detected in hundreds of sampled VPH cells, and co-detection is enriched relative to their marginal detection rates in the pooled cell tables.

This does **not** demonstrate:

- ARID5B binding at the rs1421085-containing element;
- allele-specific ARID5B occupancy;
- enhancer-to-IRX3 causality;
- that the relevant double-positive cells are the exact PH IRX3+ neuronal population driving feeding;
- that the human rs1421085 mechanism is conserved in these wild-type mouse cells.

The 2021 DEL82 negative hypothalamus result therefore remains load-bearing counter-evidence to a simple tissue-independent ARID5B model.

### H7 — sex/context interaction

**Update: simple expression-availability explanation weakened.**

The coarse pooled value `P(Arid5b+ | Irx3+)` is very similar in male and female cells (28.76% vs 28.01%). Therefore the male-specific PH `Irx3` response reported for the exact SNP is not readily explained by a gross absence of `Arid5b`/`Irx3` co-expression in females.

This does not eliminate sex dependence. The causal difference may instead involve cell subtype composition, hormone/cofactor state, chromatin accessibility, TF occupancy, enhancer contact, developmental timing, or another context variable.

## Important confounders

- This is a wild-type dataset: there is no rs1421085 genotype contrast.
- Sex is partially confounded with batch/harvest and Cell Ranger chemistry (v2/v3).
- Single-cell dropout makes zero counts imperfect evidence of biological absence.
- Pooled-cell Fisher tests are descriptive; cells are not biological replicates.
- The VPH dissection spans multiple cell classes and anatomical subregions.

## Next discriminating test

Move from **co-expression** to **cell identity + regulatory state**:

1. map `Irx3+/Arid5b+` cells onto the published neuronal/non-neuronal subclusters;
2. identify whether double-positive cells overlap the PH IRX3+ neuronal classes implicated in feeding;
3. test candidate TF/cofactor expression within those cells;
4. prioritize exact-SNP WT-vs-C assays for ATAC/CUT&RUN/Capture-C in the surviving cell population.

The target causal edge remains unresolved:

```text
rs1421085 T>C
  -> ???
  -> Irx3 ↑ in male posterior hypothalamus
```
