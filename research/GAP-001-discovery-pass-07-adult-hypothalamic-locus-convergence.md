# GAP-001 discovery pass 07 — adult hypothalamic exact-locus convergence

## Status

**Level-1 computational discovery candidate: PASS.**

**Mechanistic / causal discovery: NOT PASSED.**

`GAP-001 = OPEN` and `cause_found = false` remain unchanged.

## Narrow claim that passed the Level-1 gate

Open-data reanalysis supports a previously unreported adult-mouse hypothalamic pattern at the 500-bp interval containing the exact mouse ortholog of human `rs1421085` (`mm10 chr8:91,374,372`, T>C):

1. the interval shows rare raw chromatin accessibility below the atlas-wide candidate-enhancer/cCRE calling threshold;
2. accessibility is non-uniform across officially annotated hypothalamic cell subclasses, with a reproducible neuronal component;
3. an independent paired hypothalamus Multiome dataset contains nuclei in which an ATAC fragment overlapping the same interval and `Irx3` RNA are observed in the same filtered nucleus.

This claim does **not** say that rs1421085 causes the accessibility pattern, that the interval is an active enhancer in those nuclei, or that any candidate TF mediates `Irx3` regulation.

## Exact coordinate provenance

The mouse ortholog was resolved from the published exact-edit CRISPR sequence rather than guessed from a broad Fto interval:

- assembly: mm10
- exact base: `chr8:91,374,372` (1-based)
- zero-based base: `chr8:91,374,371`
- queried snATAC tile: `chr8:91,374,000-91,374,500`

## Evidence line 1 — adult whole-brain snATAC atlas (GSE246791)

The published atlas contains 2.3 million nuclei and 1,482 cell populations. The exact target interval is absent from the thresholded union candidate-enhancer catalog, but raw 500-bp tile/insertion matrices contain signal in hypothalamic dissections.

Seven currently provenance-complete hypothalamic samples contributed:

- `CEMBA200312_6H` / GSM7877104
- `CEMBA200319_6H` / GSM7877105
- `CEMBA200305_7J` / GSM7877102
- `CEMBA200520_7J` / GSM7877103
- `CEMBA200312_8K` / GSM7877106
- `CEMBA200319_8K` / GSM7877107
- `CEMBA200305_9L` / GSM7876880

The intended eighth sample is `GSM7876882` (`9L,rep2`). A prior workflow guessed its processed filename and failed to download it. A repair workflow now resolves the H5AD URL from the official GSM SOFT record instead of guessing. Until that repair result is incorporated, all subclass numbers below are explicitly seven-sample results.

Across the seven samples:

- total atlas nuclei in the official metadata table: **73,007**
- target-tile positive nuclei: **133**
- overall target-tile rate: **0.182%**
- official `Sample + Barcode -> Subclass` join: **133/133 matched**

### Subclass structure

Pooled subclass enrichment is descriptive: subclass target-tile rate divided by the seven-sample overall target-tile rate. Biological samples are preserved rather than treated as independent cells.

Selected results:

| subclass | target / total nuclei | pooled enrichment | samples with target signal |
|---|---:|---:|---:|
| `127 DMH-LHA Vgll2 Glut` | 4 / 322 | 6.82x | 3 / 7 |
| `090 BST-MPN Six3 Nrgn Gaba` | 8 / 1,192 | 3.68x | 3 / 7 |
| `131 LHA-AHN-PVH Otp Trh Glut` | 9 / 1,800 | **2.74x** | **5 / 7** |
| `104 TU-ARH Otp Six6 Gaba` | 7 / 1,603 | 2.40x | 4 / 7 |
| `132 AHN-RCH-LHA Otp Fezf1 Glut` | 6 / 1,515 | 2.17x | 4 / 7 |
| `106 PVpo-VMPO-MPN Hmx2 Gaba` | 8 / 2,831 | 1.55x | **6 / 7** |
| `318 Astro-NT` | 9 / 11,640 | 0.42x | — |
| `327 Oligo` | 7 / 11,060 | 0.35x | — |

The largest point estimate (`DMH-LHA Vgll2`) is based on only four target nuclei and is therefore fragile. The more defensible reproducibility signal is `LHA-AHN-PVH Otp Trh Glut`: nine target nuclei, pooled enrichment ~2.74x, and target signal in five of seven available samples. `PVpo-VMPO-MPN Hmx2 Gaba` has broader sample coverage but weaker enrichment.

This changes the interpretation of the raw signal: it is not merely a glial/ubiquitous accessibility artifact. The exact-locus-containing tile has a measurable neuronal hypothalamic structure.

## Evidence line 2 — independent paired Multiome (GSE226277)

GSE226277 provides matched RNA/ATAC Multiome hypothalamus pairs from male WT mice. It is a recurrent-hypoglycaemia study, not an rs1421085 genotype experiment, so it is used only as an orthogonal same-nucleus feasibility dataset.

Four independently processed hypothalamus pairs were verified from official GEO metadata:

- AH replicate 1
- AH replicate 2
- RH replicate 1
- RH replicate 2

A streaming fragment reader queried `mm10 chr8:91,374,000-91,374,500` without requiring a whole-file download and joined fragment barcodes to the paired filtered RNA matrix.

Result:

- **2 / 4 pairs** contain at least one filtered nucleus with both an ATAC fragment overlapping the target interval and `Irx3` RNA.
- total such nuclei: **2**.
- the other two pairs contain zero `locus+ / Irx3+` nuclei.

The two nuclei are ultra-sparse observations and are not evidence of statistical enrichment by themselves. They are important because they provide an independent modality/dataset showing that the target interval and `Irx3` transcription can coexist in the same adult hypothalamic nucleus.

Candidate transcript co-detection in the two `locus+ / Irx3+` nuclei:

- ARID5B: 2/2
- CUX1: 2/2
- TET1: 2/2
- KDM2B: 2/2
- MECP2: 1/2
- DNMT1: 0/2
- CXXC1: 0/2

Because `n=2`, these counts must **not** be used to rank mediators.

## Relation to prior evidence

The open-data result sits between established literature facts:

- the obesity-associated FTO interval can contact/regulate IRX3/IRX5 in brain/developmental contexts;
- the exact rs1421085 T>C edit increases Irx3 in adult male posterior hypothalamus and an IRX3-positive PH neuronal circuit can alter feeding/body weight;
- the missing edge is the adult-PH regulatory state connecting the exact substitution to increased Irx3.

The present pass does not fill that edge. It narrows the adult cell-state surface on which the missing mechanism should be tested.

## Novelty search result

A targeted primary-literature search was performed for combinations of:

- rs1421085 + adult hypothalamus + ATAC/chromatin accessibility;
- rs1421085 + single nucleus/single cell + brain;
- FTO obesity interval + adult hypothalamic accessibility;
- rs1421085 + Irx3 + same-nucleus Multiome.

The key prior papers establish adipocyte ARID5B circuitry, long-range FTO-IRX3 regulation, temporally restricted hypothalamic/developmental enhancer effects, the adult rs1421085 mouse phenotype/IRX3 circuit, and the underlying single-cell brain atlas resources. The search did not identify a prior primary paper reporting the specific combination reconstructed here: exact-ortholog adult hypothalamic raw accessibility resolved to subclasses plus independent same-nucleus target-interval ATAC / `Irx3` RNA co-detection.

This is a bounded literature-search statement, not proof that no unindexed or unpublished analysis exists.

## Publication-gate decision

### Gate A — novelty

**PASS for a narrow Level-1 claim**, subject to normal peer-review novelty checking.

### Gate B — provenance / exact context

**PASS.** Exact coordinate, assemblies, GEO accessions, barcodes and author-provided subclass metadata are preserved.

### Gate C — replication / orthogonal support

**PASS for Level 1.** GSE246791 provides cell-resolved snATAC structure; GSE226277 independently provides same-nucleus ATAC+RNA feasibility.

### Gate D — discrimination

**PASS for Level 1.** The result rejects the oversimplified alternatives that the adult locus is simply closed everywhere or that raw target signal is explained primarily by abundant glial populations. It prioritizes specific hypothalamic neuronal contexts for direct perturbation/occupancy experiments.

### Gate E — claim strength

- Level 1 computational discovery candidate: **PASS**
- Level 2 mechanistic computational evidence: **NOT PASSED**
- Level 3 causal molecular mechanism: **NOT PASSED**

## What would promote the result

Level 2 needs, in the relevant neuronal subclasses, at least one of:

1. allele-specific accessibility in exact T/C material;
2. direct TF occupancy at the exact site;
3. cell-state-specific enhancer-to-Irx3 contact;
4. methylation/epigenomic state that discriminates the live TF hypotheses.

Level 3 still requires exact-allele and mediator perturbation/rescue in the relevant adult PH context, or equivalently strong causal identification.

## Reproducibility anchors

- CAUSAL-DNA draft PR #1 remains open/unmerged.
- GSE226277 same-nucleus workflow run: `32237141430` — four pair jobs + aggregate SUCCESS.
- GSE246791 official whole-brain Supplementary Table 2 was used for the barcode/subclass join.
- the current 9L replicate-2 repair is provenance-first and resolves the file from the official GSM record.

## Bottom line

The discovery is **not** “we found the rs1421085 mechanism.”

The defensible discovery is narrower:

> The adult mouse hypothalamic interval containing the exact rs1421085 ortholog exhibits rare, subthreshold chromatin accessibility with reproducible neuronal subclass structure, and an independent paired Multiome dataset contains rare nuclei in which accessibility at the same interval co-occurs with Irx3 transcription.

That result is now strong enough to preserve as a Level-1 computational discovery article while GAP-001 remains open.
