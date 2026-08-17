# GAP-001 P1 — CUX1 candidate-TF presence

Status: **PASS — Cux1 is present and enriched in Irx3-containing VPH cells; site occupancy remains untested**

GitHub Actions run: `31987988505`
Artifact: `gap001-cux1-presence`
Artifact digest: `sha256:df3bcdea409ba1220b6142211857df123e2ab00e724fdfcad874a13299d910de`

## Question

Does the basic transcript-presence prerequisite hold for the exact-site CUX1 switch hypothesis H2a?

> Is `Cux1` detected in the same wild-type VPH cells as `Irx3`, especially in the `Irx3+/Arid5b+` subset?

## Pooled result

Across 44,149 filtered GSE146692 VPH cells:

- `Irx3+`: **1,791**
- `Arid5b+`: **7,132**
- `Cux1+`: **9,279 (21.017%)**
- `Irx3+ Cux1+`: **632**
- `Irx3+ Arid5b+`: **509**
- `Irx3+ Arid5b+ Cux1+`: **277**
- `P(Cux1+ | Irx3+)`: **35.288%**
- `P(Cux1+ | Irx3+ & Arid5b+)`: **54.420%**

Descriptive cell-level Fisher enrichment:

- CUX1 detection in `Irx3+` vs `Irx3-`: **OR ≈ 2.13**, p ≈ 7.65e-46
- CUX1 detection in `Irx3+/Arid5b+` vs other cells: **OR ≈ 4.59**, p ≈ 5.40e-62

These p-values are descriptive only: individual cells are not independent biological replicates.

## Per-sample result

| sample | sex | cells | Irx3+ | Arid5b+ | Cux1+ | Irx3/Arid5b double+ | Irx3/Cux1 double+ | triple+ | Cux1+ among Irx3+ | Cux1+ among Irx3+/Arid5b+ |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| AJ18003 | male | 11,552 | 431 | 1,520 | 2,251 | 109 | 131 | 48 | 30.394% | 44.037% |
| AJ18004 | female | 9,389 | 260 | 1,047 | 1,309 | 37 | 71 | 11 | 27.308% | 29.730% |
| AJ19001 | female | 10,623 | 561 | 2,187 | 2,683 | 193 | 224 | 123 | 39.929% | 63.731% |
| AJ19002 | male | 12,585 | 539 | 2,378 | 3,036 | 170 | 206 | 95 | 38.219% | 55.882% |

## By sex (pooled descriptive cells)

| sex | cells | Irx3+ | Cux1+ | Irx3/Cux1+ | Irx3/Arid5b+ | triple+ | Cux1+ among Irx3+ | Cux1+ among Irx3+/Arid5b+ |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| female | 20,012 | 821 | 3,992 | 295 | 230 | 134 | 35.932% | 58.261% |
| male | 24,137 | 970 | 5,287 | 337 | 279 | 143 | 34.742% | 51.254% |

## Interpretation

### H2a — CUX1 factor-state switch

**Update: `CELLULARLY_FEASIBLE_AND_ENRICHED`, still non-causal.**

The candidate is not excluded by cellular compartmentalization. `Cux1` is detected in more than one third of `Irx3+` cells and over half of the `Irx3+/Arid5b+` subset in this wild-type VPH dataset.

Together with the 2016 exact-site biochemical evidence that rs1421085 alleles alter CUX1 P110/P200 binding preferences, this makes CUX1 a concrete competitor/cofactor candidate for the missing CNS regulatory edge.

However, this analysis does **not** establish:

- CUX1 protein abundance or functional P110/P200 state;
- CUX1 occupancy at rs1421085 in posterior-hypothalamic cells;
- competition or cooperation between CUX1 and ARID5B;
- allele-specific enhancer activity;
- a `CUX1 -> Irx3` causal edge.

## Sex implication

The coarse transcript-presence result again does not explain the male-specific exact-SNP PH phenotype by simple absence in females. `Cux1` co-detection with `Irx3` is similar in pooled male/female cells, and CUX1 among `Irx3+/Arid5b+` is not lower in pooled female cells.

This shifts attention toward **regulatory state**, not mere gene presence: protein processing, hormone-dependent cofactors, occupancy, accessibility, 3D contact, developmental timing, or cell-subtype composition.

## Next decisive discriminator

The highest-value wet-lab / perturbational test is now increasingly specific:

```text
WT-T PH cells vs exact-C PH cells, male/female
  -> ARID5B occupancy at rs1421085
  -> CUX1 occupancy / functional state at rs1421085
  -> local accessibility
  -> enhancer-to-Irx3 contact
  -> Irx3 expression
```

Then perturb `CUX1` and `ARID5B` separately and jointly while preserving the SNP-containing element. This would distinguish loss-of-repressor, factor-switch, and cooperative models.
