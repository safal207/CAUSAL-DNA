# GAP-001 discovery pass 03 — full Codebook allele motif switch

## Question

Does the exact human `rs1421085 T>C` substitution alter sequence-recognition models for transcription factors beyond the previously prioritized ARID5B/CUX1 candidates, and which of those candidates remain feasible in adult ventral/posterior hypothalamus (VPH)?

This pass also asks whether the published meSMiLE-seq methylation classifications directly annotate the strongest allele-switch candidates.

## Provenance

- Nature Codebook article: DOI `10.1038/s41586-026-10798-9`.
- Exact-DOI Springer media object selected by parser contract: `MOESM4`.
- Selected archive SHA256: `76e2d9b5f3dbe9f399bc70a89e68c58ebc082c9b8c872532ef6423c396e13588`.
- Source members: 1,687.
- Parsed motif records: **1,685**.
- Unique human TFs: **1,421**.
- TFs with multiple published representatives: **264**; maximum two representatives per TF.
- Unparsed members: archive bookkeeping (`tmp.tmp`) and a macOS metadata file; neither is a motif record.
- Adult VPH expression feasibility: GSE146692, four samples (`AJ18003`, `AJ18004`, `AJ19001`, `AJ19002`).
- meSMiLE-seq published classification: Supplementary Data 5, worksheet `Classification by methylation`, DOI `10.1038/s41467-026-71387-y`.
- Executed workflow: GitHub Actions run `32231711305`, job `96002757014` — **SUCCESS**.

## Method contract

For each published motif representative:

1. use the exact hg38 rs1421085 sequence context;
2. compare `T` and `C` only at motif placements that physically overlap the SNP;
3. compare the same placement and strand for both alleles;
4. scan forward and reverse-complement strands;
5. normalize each motif against its own row-wise theoretical min/max;
6. retain a stable `motif_record_id` so multiple representatives of one TF are not collapsed;
7. use VPH transcript detection only as a feasibility filter;
8. join published meSMiLE classes as an annotation axis, not as a causal or composite score.

The `0.75` normalized-site crossing remains a declared heuristic for triage. It is not a biochemical affinity threshold.

## Full-atlas result

The full scan covered **1,685 motif records / 1,421 unique TFs**.

- motif records with maximum normalized site strength `>=0.75`: **17**;
- heuristic `0.75` threshold-crossing motif records: **17**.

### Threshold crossings

| TF | Representative | Direction | T | C | abs Δnorm | VPH +% |
|---|---|---|---:|---:|---:|---:|
| TBP | TBP | loss on C | 0.8711 | 0.6642 | 0.2069 | 6.895 |
| TBPL2 | TBPL2 | loss on C | 0.8711 | 0.6642 | 0.2069 | 0.002 |
| MEIS1 | MEIS1 | gain on C | 0.5665 | 0.7651 | 0.1986 | 2.507 |
| GATA2 | GATA2 | loss on C | 0.8347 | 0.6582 | 0.1765 | 1.268 |
| SIX1 | SIX1__rep02 | loss on C | 0.8092 | 0.6351 | 0.1741 | 0.059 |
| NR2C1 | NR2C1__rep01 | gain on C | 0.6364 | 0.8098 | 0.1733 | 4.048 |
| ZNF649 | ZNF649 | loss on C | 0.8333 | 0.6667 | 0.1667 | 0.000 |
| GATA5 | GATA5 | loss on C | 0.7752 | 0.6110 | 0.1642 | 0.000 |
| NR5A1 | NR5A1 | gain on C | 0.6963 | 0.8390 | 0.1427 | 0.034 |
| ESRRG | ESRRG | gain on C | 0.7116 | 0.8498 | 0.1383 | 2.945 |
| FOSL1 | FOSL1 | gain on C | 0.6270 | 0.7621 | 0.1352 | 0.007 |
| NR4A1 | NR4A1__rep02 | gain on C | 0.6657 | 0.7991 | 0.1333 | 4.902 |
| ESRRA | ESRRA | gain on C | 0.7273 | 0.8593 | 0.1320 | 11.321 |
| RORB | RORB | gain on C | 0.6289 | 0.7573 | 0.1284 | 20.540 |
| NR5A2 | NR5A2 | gain on C | 0.7394 | 0.8591 | 0.1198 | 0.503 |
| GATA6 | GATA6 | loss on C | 0.7872 | 0.6711 | 0.1161 | 0.025 |
| ZNF236 | ZNF236__rep02 | loss on C | 0.7773 | 0.6688 | 0.1084 | 6.476 |

Several crossing candidates are poor adult-VPH direct-interpreter candidates by transcript feasibility alone (`TBPL2`, `ZNF649`, `GATA5`, `NR5A1`, `FOSL1`, `GATA6`). This is demotion only; transcript dropout and context dependence prevent an absolute exclusion claim.

The most immediately interesting crossing candidates combining a threshold crossing with non-trivial VPH detection are therefore currently:

- `TBP` — loss on C; VPH ~6.9%;
- `MEIS1` — gain on C; VPH ~2.5%;
- `NR2C1` — gain on C; VPH ~4.0%;
- `NR4A1` — gain on C; VPH ~4.9%, with a descriptive male>female difference in this confounded dataset;
- `ESRRA` — gain on C; VPH ~11.3%;
- `RORB` — gain on C; VPH ~20.5%;
- `ZNF236` — loss on C; VPH ~6.5%.

These are **sequence-model candidates only**.

## Strongest allele-sensitive motifs regardless of the 0.75 crossing heuristic

The first full-atlas scan unexpectedly puts several chromatin / epigenetic regulators near the top of absolute allele sensitivity:

1. `CXXC1` — gain on C, abs Δnorm 0.3470, VPH ~10.99%;
2. `TET1` — gain on C, 0.3054, VPH ~8.32%;
3. `KDM2B` — gain on C, 0.2925, VPH ~9.05%;
4. `DNMT1` — gain on C, 0.2911, VPH ~11.55%;
5. `CXXC4` — gain on C, 0.2888, VPH ~22.64%;
6. `DBX2` — loss on C, 0.2823, VPH ~9.87%;
7. `GMEB1` — gain on C, 0.2713, VPH ~5.09%;
8. `SP100` — gain on C, 0.2700, VPH ~1.43%;
9. `MECP2` — loss on C, 0.2662, VPH ~24.88%;
10. `KMT2A` — gain on C, 0.2637, VPH ~27.57%;
11. `ARID5B` — **loss on C**, 0.2577, VPH ~16.15%.

`ARID5B loss_on_C` is directionally concordant with the known adipocyte motif-disruption literature, which is a useful positive-control behavior of this sequence-model scan. It does **not** establish ARID5B occupancy or mediation in adult PH.

The apparent concentration of chromatin / methylation-associated proteins near the top is a **pattern to test**, not an enrichment claim. No enrichment statistic or matched background analysis has yet been performed.

## Methylation-aware result

Published meSMiLE classifications parsed: **48 TFs**.

Across the 1,685 Codebook motif records:

- **51 motif records** correspond to meSMiLE-profiled TFs;
- **29 motif records** correspond to published `methyl plus` or `methyl minus` classes;
- **0 / 17 threshold-crossing motif records** are profiled in the meSMiLE classification sheet.

Therefore the new methylation study **does not directly explain the 17 strongest threshold-crossing candidates** in this pass. `unknown` means unprofiled/unclassified by that experiment, not methylation-insensitive.

Methylation remains relevant to lower-threshold candidates and to the broader chromatin-state hypothesis, but it cannot currently be used to promote one of the 17 crossings.

## Important source-sensitivity observation

The earlier 235-motif manual subset produced a different set of 0.75 crossings than the full 1,421-TF representative atlas. This means motif-source / representative choice is a load-bearing variable in candidate triage.

The full atlas is preferable for coverage, but a candidate should not be promoted merely because it ranks highly in one motif representation. The next computational gate should therefore measure **cross-source / cross-representative stability** before prioritizing occupancy experiments.

## What changed scientifically

Before this pass, adult-PH candidate space was dominated by manually carried hypotheses (`ARID5B`, `CUX1`, accessibility, 3D contact).

After this pass, the sequence-level candidate space is broader:

`rs1421085 T>C`
→ allele-sensitive TF recognition models including `TBP/MEIS1/NR2C1/NR4A1/ESRRA/RORB/...`
→ plus a non-crossing high-Δ group containing `CXXC1/TET1/KDM2B/DNMT1/CXXC4/MECP2/KMT2A/ARID5B`
→ **unknown in-vivo occupancy / cell-state specificity**
→ `Irx3 ↑` remains unresolved.

This is a real narrowing/expansion of the mechanistic search space, but **not a discovered mediator**.

## Required next evidence

1. cross-source / cross-representative stability analysis for allele-switch candidates;
2. identify which candidates are expressed in the exact adult-PH signal-carrying cell class rather than pooled VPH;
3. test motif/footprint support in the low-frequency accessible rs1421085-ortholog nuclei;
4. obtain allele-specific occupancy or chromatin evidence in matched `T/C × sex × PH cell type` context;
5. connect candidate occupancy to `Irx3` regulation and perturb/rescue the mediator.

## Causal status

`GAP-001 = OPEN`

`cause_found = false`

No causal edge was promoted by this discovery pass.
