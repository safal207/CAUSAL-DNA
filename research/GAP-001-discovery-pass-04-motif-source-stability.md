# GAP-001 discovery pass 04 — motif-source stability

## Question

How stable are exact `rs1421085 T>C` motif-switch calls when the PWM source / representative changes?

This pass compares:

- `manual235`: the 235-motif manual Codebook subset used in discovery pass A;
- `full1421`: the official Nature Codebook atlas selected from exact-DOI Springer `MOESM4`, normalized as 1,685 motif records covering 1,421 unique TFs.

The analysis intentionally compares sequence-model behavior only. VPH expression, methylation classes, accessibility, occupancy and target-gene evidence are not used to define motif-source stability.

## Provenance

- Exact variant: human `rs1421085 T>C`.
- Human reference context: hg38, exact SNP-centered sequence fetched from UCSC.
- Full Codebook archive SHA256: `76e2d9b5f3dbe9f399bc70a89e68c58ebc082c9b8c872532ef6423c396e13588`.
- Stability workflow: GitHub Actions run `32232128194`, job `96004044716` — **SUCCESS**.
- Stability unit tests: **7/7 PASS**.
- Motif scanner unit tests: **7/7 PASS**.
- Main causal/evidence pipeline on the same head: run `32232128080`, job `96004044710` — **SUCCESS**.

## Comparison contract

For each TF within each source:

1. scan only motif placements that overlap the SNP;
2. compare T and C at the same placement and strand;
3. retain all motif records in the raw scan;
4. for source-to-source comparison, select the most allele-sensitive representative of each TF within that source;
5. compare sign/direction of `Δ`, absolute `|Δ|`, maximum normalized site score and the declared `0.75` crossing heuristic;
6. call a TF `robust_delta_direction` only when both sources give the same non-zero direction and `|Δ| >= 0.10` in both.

`0.10` is a declared robustness threshold for this pass, not a biochemical affinity threshold.

## Main result

The two sources share **202 TFs** by the same TF label.

- direction-consistent TFs: **165 / 202 = 81.7%**;
- signed `Δ` Pearson correlation: **0.6733**;
- `|Δ|` Pearson correlation: **0.7979**;
- same-direction TFs with `|Δ| >= 0.10` in both sources: **77**;
- replicated `0.75` crossings: **0**.

### Crossing instability

`manual235` crossings within the shared TF set:

- `RARA`
- `SOX2`
- `SRY`
- `ZNF516`

`full1421` crossing within the shared TF set:

- `RORB`

Replicated crossings:

- **none**.

This means that the binary `0.75` crossing criterion is highly source / representative sensitive in this comparison. It should not be used as the primary candidate-selection rule for GAP-001.

## Top cross-source stable allele-sensitive TFs

Ranked by the minimum absolute allele effect observed across the two motif sources, requiring the same direction and `|Δ| >= 0.10` in both:

| Rank | TF | minimum `|Δ|` across sources |
|---:|---|---:|
| 1 | **CXXC4** | **0.2888** |
| 2 | SP100 | 0.2634 |
| 3 | SP140 | 0.2405 |
| 4 | GLYR1 | 0.2250 |
| 5 | KDM2A | 0.2184 |
| 6 | ZBTB8A | 0.2016 |
| 7 | TCF20 | 0.1977 |
| 8 | MBNL2 | 0.1934 |
| 9 | MKX | 0.1933 |
| 10 | AHCTF1 | 0.1925 |
| 11 | FBXL19 | 0.1917 |
| 12 | DNTTIP1 | 0.1895 |
| 13 | ZBED2 | 0.1893 |
| 14 | ZGLP1 | 0.1869 |
| 15 | MBD1 | 0.1817 |
| 16 | MBD3 | 0.1805 |
| 17 | SP140L | 0.1765 |
| 18 | ZNF703 | 0.1742 |
| 19 | ZBED4 | 0.1735 |
| 20 | ZNF395 | 0.1678 |

## CXXC4 convergence

`CXXC4` is currently the strongest cross-source stable candidate in this specific sequence-model comparison:

- `gain_on_C` in the manual motif source;
- `gain_on_C` in the full Nature Codebook source;
- minimum cross-source `|Δ| = 0.2888`;
- adult pooled VPH transcript feasibility from the earlier GSE146692 pass: approximately **22.6%** of cells detected `Cxxc4`;
- published meSMiLE-seq class: **methyl minus**;
- the rs1421085 mouse-ortholog locus has reproducible low-frequency raw accessibility across **8/8** adult hypothalamus snATAC samples in GSE246791.

These observations make `CXXC4` a **high-priority candidate for the next locus/cell-state-specific test**.

They do **not** establish:

- CXXC4 occupancy at rs1421085;
- allele-specific CXXC4 binding in adult PH;
- CXXC4 regulation of `Irx3`;
- methylation dependence at this locus;
- mediation of the mouse or human phenotype.

## Methodological implication

For this variant, continuous direction / magnitude of motif change is substantially more reproducible across motif sources than a single binary score threshold.

That is a useful robustness observation for CAUSAL-DNA, but it is **not claimed as a novel general result** without a broader variant benchmark and literature comparison.

The practical rule for GAP-001 is therefore:

> prioritize candidates by cross-source continuous effect stability, then require biological context filters; do not promote candidates because one PWM crosses an arbitrary threshold.

## Next decisive computational gate

The highest-information next analysis is **cross-species motif stability** across the exact human SNP and the exact mouse ortholog used in the 2025 T>C knock-in study.

For high-priority cross-source candidates, especially `CXXC4`, ask whether the same T>C allele direction survives the mouse flanking sequence:

`human T/C motif switch`
→ `same candidate / mouse ortholog T/C motif switch`
→ `adult hypothalamus cell-state feasibility`
→ eventual allele-specific occupancy / perturbation.

If a candidate's allele switch disappears or reverses in the mouse sequence context, it becomes a poorer explanation for the mouse `Irx3` phenotype even if it remains a human-specific regulatory candidate.

## Causal status

`GAP-001 = OPEN`

`cause_found = false`

No causal edge is promoted by this stability pass.
