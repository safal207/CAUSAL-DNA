# GAP-001 discovery pass 02 — official adult hypothalamus cell-state localization

## Question

The previous discovery pass established reproducible low-frequency raw chromatin accessibility at the exact mouse rs1421085 ortholog in all eight queried adult hypothalamus snATAC samples, despite the locus being absent from the thresholded union cCRE catalogue.

This pass asks the next load-bearing question:

> Which officially annotated adult-brain cell states carry that raw locus signal?

This is a **D1 descriptive localization analysis**. It is not allele-specific occupancy, Irx3 expression, enhancer-to-gene contact, or causal mediation.

## Provenance

- Dataset: `GSE246791` / adult mouse whole-brain snATAC atlas.
- Variant coordinate: mm10 `chr8:91,374,372` (1-based), exact mouse ortholog of human `rs1421085 T>C`.
- Signal source: workflow run `32227940175`, eight successful per-sample jobs plus successful aggregate.
- Signal definition:
  - target 500-bp tile `chr8:91374000-91374500`;
  - raw Tn5 insertion window ±250 bp around the exact ortholog;
  - `both` = nucleus is positive for both measurements.
- Official annotation source: Nature whole-brain supplementary metadata artifact from workflow run `32228624673`, artifact `WMB-supplementary-cellmeta-candidates` (`9356434107`, SHA256 `b61b67b5a1b7a7fac17529f946b2beb5a657998fc2d6b6b9de92dc808d483e6e`).
- Official files used:
  - SI Table 2: metadata for all 2.3 million snATAC nuclei (`CellID`, `Sample`, `Barcode`, `L1-L4`, neurotransmitter, subclass);
  - SI Table 4: L4-to-subclass and major-region integration annotation.

The join key was exactly `Sample + Barcode`; no marker-based re-annotation was performed.

## Join integrity

Across the eight hypothalamus samples:

- background nuclei represented in official SI Table 2: **83,321**;
- unique nuclei carrying at least one locus signal: **166**;
- official metadata matches: **166/166**;
- unmatched signal barcodes: **0**.

Signal-set sizes:

- target-tile positive: **154** nuclei;
- ±250-bp insertion-window positive: **142** nuclei;
- positive for both: **130** nuclei;
- positive for either: **166** nuclei.

## Main localization result

The strongest reproducible enrichment in the stricter `both` set is:

### `131 LHA-AHN-PVH Otp Trh Glut`

- signal nuclei: **11/130 = 8.46%**;
- same-sample background frequency: **2.21%**;
- one-sided Fisher exact odds ratio: **~4.11**;
- nominal `p ≈ 1.58e-4`;
- Benjamini-Hochberg `q ≈ 0.00584` across observed subclasses;
- observed in **6 of 8 samples**.

The same subclass remains the leading signal under other definitions:

- target tile: 11/154, OR ~3.42, BH q ~0.027;
- insertion window: 11/142, OR ~3.73, BH q ~0.0129;
- either signal: 11/166, OR ~3.15, BH q ~0.051.

This cross-definition persistence is useful because tile and insertion sets are not identical.

## L4 concentration inside subclass 131

All 11 signal nuclei assigned to subclass 131 are also positive for both tile and insertion-window signal.

Their official L4 assignments are:

- `2-8-4-1`: **8/11**;
- `2-8-10-1`: **2/11**;
- `2-8-14-1`: **1/11**.

SI Table 4 independently maps all three L4 clusters to `LHA-AHN-PVH Otp Trh Glut` and reports major-region weights strongly concentrated in hypothalamus:

- `2-8-4-1`: ATAC L4 major region `HY:0.95; PAL:0.04; TH:0.01`;
- `2-8-10-1`: `HY:0.86; PAL:0.07; AMY:0.04`;
- `2-8-14-1`: `HY:0.95; PAL:0.04; AMY:0.01`.

Thus the leading signal is not produced by importing an annotation from another brain region.

## Secondary signals worth retaining

For the strict `both` set, additional enrichments include:

- `127 DMH-LHA Vgll2 Glut`: 4/130, background 0.46%, OR ~6.95, BH q ~0.045; present in 3 samples;
- `116 AVPV-MEPO-SFO Tbr1 Glut`: 3/130, background 0.24%, OR ~10.16, BH q ~0.045, but all three observations come from one sample and therefore this is not replicate-stable;
- `090 BST-MPN Six3 Nrgn Gaba`: 6/130, background 1.44%, OR ~3.32, BH q ~0.108;
- `100 AHN Onecut3 Gaba`: broad recurrence across samples but weaker multiple-testing support in the strict set.

`VMH Nr5a1 Glut` is frequent among signal nuclei (9/130 strict; 12/166 either) and recurs across five samples, but its enrichment is weaker after multiple-testing correction than subclass 131.

## Negative / balancing observations

The signal is **not exclusive to one cell class**. It also occurs in astrocytes, oligodendrocytes, tanycytes, and multiple neuronal subclasses. In particular, common non-neuronal classes are not enriched in the same direction:

- `Astro-NT NN`: 11/130 strict versus ~16.26% background;
- `Oligo NN`: 7/130 strict versus ~15.24% background.

This argues against interpreting the locus as a universally accessible housekeeping element and is compatible with a sparse cell-state-biased regulatory element, but it does not prove that interpretation.

## Updated causal boundary

The previous state was:

`exact locus -> reproducible low-frequency adult-HY accessibility -> ?`

This pass supports a narrower descriptive frontier:

`exact locus -> sparse accessibility enriched in LHA-AHN-PVH Otp Trh Glut (especially L4 2-8-4-1) -> ?`

It **does not** establish:

- that the risk C allele causes this accessibility state;
- that the signal nuclei express `Irx3`;
- that subclass 131 is the native cell of action for the 2025 phenotype;
- that any Codebook TF binds the locus in these cells;
- that the locus contacts the `Irx3` promoter;
- that the enrichment is independent-dataset replication.

The eight samples are biological replicates within one atlas and therefore remain D1, not D2.

## Next decisive joins

1. **Irx3 alignment:** test `Irx3` expression in the RNA counterpart of subclass 131 / L4-linked cell states using an independent adult hypothalamus or whole-brain RNA atlas.
2. **TF feasibility:** intersect the exact T>C Codebook motif-switch candidates with TF expression in subclass 131 rather than pooled VPH.
3. **Independent chromatin replication:** test the orthologous locus in a separate adult hypothalamus ATAC/multiome dataset.
4. **Mechanism:** allele-specific accessibility/TF occupancy plus enhancer-to-`Irx3` contact and perturbation in the same cell state.

## Status

- `GAP-001`: **OPEN**
- `cause_found`: **false**
- discovery ceiling: **D1 cell-state localization**

The result changes prioritization, not causal status.
