# GAP-001 Discovery Pass 07 — structural and methylation convergence at the exact rs1421085 ortholog

Date: 2026-08-19

## Question

What adult posterior-hypothalamic molecular regulatory state connects the exact rs1421085 T>C substitution to increased Irx3 expression, particularly in males?

The unresolved edge remains:

`rs1421085 T>C -> M? -> Irx3 up in adult male posterior hypothalamus`

This pass asks whether public adult-hypothalamus data support three prerequisites of a local regulatory mechanism around the exact mouse ortholog:

1. a reproducible cell-resolved accessible state at the SNP-containing element;
2. structural proximity from the SNP-containing contact bin to the Irx3 promoter;
3. a local methylated-C substrate compatible with methylation-sensitive regulatory readers.

None of these prerequisites alone establishes allele dependence or a molecular mediator.

## 1. Exact coordinate and element

The sequence-verified mouse ortholog is:

- mm10 / GRCm38: `chr8:91,374,372` (1-based), T>C;
- exact 0-based site: `chr8:91374371`;
- 500-bp snATAC tile: `chr8:91374000-91374500`.

An independent assembly cross-check maps the exact base to mm39 `chr8:92,101,000` (1-based). The exact ortholog is not the previously functionally deleted `Fto-Irx::hibE1` element; after assembly reconciliation the exact base is about 13.9 kb from hibE1. Therefore the 2025 hibE deletion experiments do not constitute a perturbation of the exact rs1421085 ortholog.

## 2. Adult hypothalamic accessibility is rare but reproducible

GSE246791 adult male mouse snATAC analysis preserves biological replicates across four hypothalamic dissections (6H, 7J, 8K, 9L; two replicates each).

After provenance-safe repair of the eighth sample using the official GEO SOFT record:

- intended samples: 8/8 available;
- total nuclei across the all-8 subclass analysis: 83,321;
- exact 500-bp tile-positive nuclei: 154;
- overall tile-positive rate: ~0.1848%.

The exact locus is therefore not a strong constitutive atlas peak, but it is not closed: low-frequency raw accessibility is reproducible across the intended hypothalamic sample set.

### Replicate-aware cell-state prioritization

The strongest reproducible primary neuronal subclass is:

`LHA-AHN-PVH Otp Trh Glut`

- 11 / 1,842 nuclei tile-positive;
- signal present in 6 / 8 biological samples;
- pooled enrichment ~3.231x;
- Mantel-Haenszel OR ~3.640;
- 95% CI ~1.952-6.790;
- BH-adjusted q ~0.00137;
- Breslow-Day p ~0.032, indicating cross-sample heterogeneity that must be retained.

`DMH-LHA Vgll2 Glut` shows a larger point enrichment (~5.67x) but only 3/8-sample support, so it remains a secondary, less replicate-stable candidate state.

This prioritizes a neuronal cell context. It does not identify the native causal cell of action.

## 3. Independent same-nucleus feasibility

In independent paired RNA+ATAC hypothalamus Multiome data (GSE226277), exact-locus ATAC and Irx3 RNA occur in the same nucleus in 2/4 biological pairs.

Only two such nuclei were observed, so the result is sparse. Both co-detect ARID5B, CUX1, TET1 and KDM2B RNA; MECP2 is detected in one of the two. DNMT1 and CXXC1 are not detected in these two nuclei.

Interpretation ceiling: orthogonal same-nucleus feasibility, not mediator identification.

## 4. Exact-site-containing 10-kb bin contacts the Irx3 promoter bin

Reanalysis of the published adult hypothalamus H3K27ac+ PLAC-seq callset (GSE295850) finds one significant interaction between:

- target bin containing the exact ortholog: `chr8:91370000-91380000`;
- Irx3 promoter/TSS bin: `chr8:91800000-91810000`.

Primary published callset statistics:

- observed count: 49;
- expected: 12.6272643735288;
- FDR: 2.8483024481983e-14.

This is strong wild-type adult-hypothalamus structural support for the target-bin -> Irx3-promoter-bin relationship.

### Required conflict preservation

The secondary interaction-model table contains coordinate candidates near the same regions but does **not** reproduce the exact primary count/expected/FDR tuple.

A separate Fed/Fasted 5-kb comparison contains the corresponding target/Irx3 pair but does not show a state-dependent difference (FDR = 1.0 in the processed table).

Therefore the correct claim is:

> the primary published 10-kb callset supports a wild-type adult-hypothalamus contact between the bin containing the exact rs1421085 ortholog and the Irx3 promoter bin; secondary statistical/state layers do not independently reproduce the same significance.

This cannot resolve whether the exact 500-bp element or exact base is the contacting anchor, whether the contact is allele-dependent, or which cell subclass carries it.

## 5. Local methylated-C substrate is reproducible in 3/3 adult-male hypothalamus WGBS replicates

GSE84533 WGBS was queried after dynamic mm10 -> mm9 liftOver of the exact ortholog. The exact wild-type ortholog base is T and therefore is not itself a methylatable cytosine; the analysis asks whether nearby cytosines provide methylation-sensitive regulatory substrate.

Exact mm9 ortholog used in this dataset:

- `chr8:93,898,271` (1-based);
- 0-based `chr8:93898270`.

### +/-250 bp mCG

Coverage-weighted mCG methylation:

- replicate 1 (GSM2241593): 95.918%;
- replicate 2 (GSM2241594): 97.222%;
- replicate 3 (GSM2241595): 93.103%.

The nearest qualifying CG lies **1 bp upstream** of the ortholog in all three replicates:

- replicate 1: coverage 11, methylation 100%;
- replicate 2: coverage 6, methylation 83.33%;
- replicate 3: coverage 6, methylation 83.33%.

### +/-250 bp mCAC

Coverage-weighted mCAC is lower and more variable:

- replicate 1: 12.903%;
- replicate 2: 25.974%;
- replicate 3: 13.043%.

The nearest qualifying CAC is +19 bp from the ortholog; methylation at that individual site is low/absent in these bulk replicates.

Interpretation:

- a highly methylated CG substrate immediately adjacent to the ortholog is reproducible in 3/3 adult male hypothalamus WGBS replicates;
- this makes methylation-sensitive reader/epigenetic-state hypotheses biophysically feasible;
- it does **not** establish MECP2, TET1, DNMT1, KDM2B or another factor as the mediator;
- it is bulk wild-type data, not allele-resolved or cell-subclass-resolved.

## 6. What this changes

Before this pass, the adult PH gap could still be framed as an almost unconstrained missing mediator. The public-data evidence now supports a much narrower regulatory surface:

`exact rs1421085 ortholog`

`-> rare adult hypothalamic accessible state`

`-> replicate-prioritized LHA/AHN/PVH neuronal context`

`-> structural compatibility with Irx3 promoter contact`

`+ strongly methylated adjacent CG substrate`

`-> M still unknown`

`-> Irx3 up after exact T>C in vivo`

This is **convergence of prerequisites**, not a completed causal chain.

## 7. Current candidate interpretation

### ARID5B / CUX1

Remain high-priority exact-site candidates because:

- their allele-sensitive motif grammar survives the strict human->mouse motif-slot portability test;
- both are transcriptionally feasible in adult Irx3+ VPH cells;
- both appear in the two sparse same-nucleus locus+/Irx3+ Multiome nuclei.

A targeted search has not yet identified a direct public adult-hypothalamus ARID5B or CUX1 occupancy dataset at this site. Direct allele-resolved occupancy remains a decisive missing measurement.

### Methylation-sensitive / epigenetic-state branch

The 3/3 WGBS result materially increases the plausibility of a methylation-sensitive local state, because a highly methylated CG sits immediately adjacent to the exact ortholog.

MECP2 remains a candidate reader because its strict motif score, adult Irx3+ expression feasibility and local methylated-C substrate are compatible. However the direct adult-hypothalamus MeCP2 ChIP/Input exact-locus gate is still pending at the time of this note and must decide whether local occupancy is actually enriched above matched input.

TET1/KDM2B/DNMT1-related state hypotheses also remain open; bulk methylation alone cannot distinguish them.

## 8. Novelty boundary

Do **not** claim:

- first Fto-region -> Irx3 contact;
- first adult hypothalamic regulation within the Fto-Irx TAD;
- rs1421085-specific 3D contact causality;
- methylation-mediated rs1421085 mechanism;
- a discovered TF mediator.

Prior work already establishes broader Fto-Irx regulatory contacts and functional CREs in hypothalamic contexts.

The defensible computational contribution is narrower:

> exact rs1421085-ortholog-focused integration of reproducible 500-bp adult hypothalamic accessibility, replicate-aware neuronal subclass enrichment, sparse independent same-nucleus locus/Irx3 feasibility, a published significant 10-kb target-bin/Irx3-promoter-bin contact, and reproducible local methylated-C substrate.

## 9. Highest-information next evidence

The decisive experiment remains an exact-allele, sex-stratified, cell-resolved adult posterior-hypothalamus assay measuring in matched biological material:

1. Irx3 transcription;
2. exact-locus accessibility;
3. ARID5B, CUX1 and top methylation-sensitive reader occupancy;
4. local methylation state;
5. exact-element -> Irx3 contact;
6. followed by mediator perturbation/rescue.

The first public-data priority before wet-lab work is to finish the adult hypothalamus MeCP2 ChIP/Input exact-locus check and to determine whether existing H3K27ac/ATAC processed peaks can narrow the 10-kb PLAC anchor to the exact 500-bp element.

## Status

`GAP-001 = OPEN`

`cause_found = false`

No causal edge is promoted by this pass.
