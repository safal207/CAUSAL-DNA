# DISCOVERY-SPRINT-001 — resolve the adult posterior-hypothalamus mediator

## Discovery target

Resolve at least one previously unresolved, independently reproducible step in:

```text
rs1421085 T>C
  -> M_adult_male_PH ?
  -> Irx3 up
```

The sprint is **evidence-first**. Its purpose is to discover a new fact, not to improve the planning machinery.

## What counts as progress

| Level | Meaning | Can close GAP-001? |
|---|---|---|
| D0 | rediscovery of a published result | no |
| D1 | new descriptive pattern in one public dataset | no |
| D2 | same new pattern reproduced in an independent dataset or modality | no |
| D3 | variant-region regulatory state is linked to Irx3 in the relevant adult PH cell context and competing mechanisms are discriminated | candidate causal edge only |
| D4 | exact variant or mediator is perturbed/rescued in the relevant adult PH context with independent verification | potentially, subject to the causal-proof protocol |

No computational result may set `cause_found=true`.

## Load-bearing question

> Which adult posterior-hypothalamic cell state and regulatory mechanism connects the conserved rs1421085 element to increased Irx3 expression, particularly in males?

## Work packages

### WP1 — cell-of-action localization

Use adult/juvenile posterior-hypothalamus transcriptomic resources to identify reproducible Irx3+ cell populations and quantify sex/context structure.

Primary resources:

- GSE146692 — ventral posterior hypothalamus scRNA-seq, male and female samples.
- GSE208355 — HypoMap, integrated murine hypothalamus atlas.

Exit criterion: a cell population definition reproducible across at least two sources or an explicit `INSUFFICIENT_RESOLUTION` result.

### WP2 — exact conserved locus mapping

Map the published mouse ortholog of human rs1421085 onto a current mouse assembly using the published conserved sequence/CRISPR guide rather than assuming coordinates.

Published anchors:

- Laber et al. 2021: orthologous CRM, guide `TAATCAATACGATGCCTT`, donor sequence, and the rs1421085 enhancer-region fosmid at mm9 chr8:93,911,541-93,950,923.
- Wang et al. 2023: homologous mouse point mutation reported as `Fto c.45+60752 T>C`.

Exit criterion: exact assembly + strand + coordinate + local sequence, with provenance and a reproducible mapping script.

### WP3 — adult hypothalamic chromatin accessibility

Interrogate the mapped locus in adult brain/hypothalamus chromatin resources.

Primary resources:

- GSE246791 — adult mouse whole-brain snATAC atlas, ~2.3M cells, 117 anatomical dissections.
- GSE226277 — hypothalamus single-nucleus Multiome RNA+ATAC.
- GSE288514 — male mouse hypothalamus Multiome RNA+ATAC in diet-induced obesity / incretin perturbation context; context-mismatched but useful as a metabolic-state replication lane.

Exit criterion: determine whether the orthologous rs1421085 element is accessible in Irx3-relevant hypothalamic cell classes and whether that pattern reproduces across resources.

### WP4 — regulatory interpreter discrimination

Within cells where the locus is accessible, compare candidate mechanisms:

```text
H1 ARID5B reuse
H2 CUX1 exact-site interpreter
H3 accessibility-state mechanism
H4 enhancer-to-Irx3 contact mechanism
```

Presence is not occupancy. Motif matches are not binding. Binding is not causal mediation.

Exit criterion: at least one hypothesis is materially demoted or one candidate is promoted to a sharper, falsifiable D3 candidate.

### WP5 — 3D contact

Test whether the adult hypothalamic cell class carrying the accessible locus also supports Fto-intron -> Irx3 regulatory contact using cell-type-resolved contact/enhancer-link resources where available.

Exit criterion: explicit positive, negative, or unresolved contact result with cell-context provenance.

## First execution order

```text
WP2 exact locus coordinate
   -> WP3 locus accessibility
   -> WP1 exact Irx3+ cell class alignment
   -> WP4 TF/interpreter discrimination
   -> WP5 contact integration
```

The coordinate gate goes first because a locus-specific ATAC claim without a verified orthologous coordinate is not admissible evidence.

## Discovery claim template

A candidate discovery must be written as:

```text
OBSERVED:
  what was measured

CONTEXT:
  species / assembly / region / cell type / sex / state

NOVELTY CHECK:
  what prior primary literature already showed and what it did not show

REPLICATION:
  independent dataset / modality

ALTERNATIVES:
  competing explanations still compatible

FALSIFICATION:
  experiment or analysis that would refute the candidate

STATUS:
  D0 / D1 / D2 / D3 / D4
```

## Stop rule

If public data cannot resolve the mechanism, the sprint must terminate with the **minimum decisive wet-lab experiment** rather than filling the gap with inference.
