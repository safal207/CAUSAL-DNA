# GAP-001 Decisive Regulatory Matrix

Status: proposed experiment specification
Target: `rs1421085 T>C -> ??? -> Irx3 up in adult male posterior hypothalamus`

## Objective

Discriminate among local TF-occupancy, chromatin-accessibility, 3D-contact, cell-context, and developmental-state explanations using the smallest matched experiment that measures the complete DNA-to-transcription bridge.

## Experimental factors

| Factor | Minimum levels | Why load-bearing |
|---|---|---|
| genotype | T/T, C/C; add T/C if feasible | establishes exact-allele contrast and enables allele-specific measurements |
| sex | male, female | adult PH Irx3 genotype response is sex-dependent in current mouse evidence |
| region | posterior hypothalamus | avoids dilution by whole/bulk hypothalamus |
| cell identity | resolved empirically | native SNP cell of action is not established |
| metabolic state | matched | diet/temperature/state can gate locus effects |
| age/development | matched adult; developmental comparator optional | interval regulation is temporally gated |

## Matched readout panel

| Layer | Readout | Candidate hypothesis tested | Minimal interpretation |
|---|---|---|---|
| target RNA | `Irx3` | boundary condition | confirm cell-resolved genotype effect |
| nearby RNA | `Fto`, `Rpgrip1l`, `Irx5` | alternative target routing | detect multi-target enhancer grammar |
| accessibility | rs1421085-element ATAC/snATAC | H3 | allele/state-dependent element accessibility |
| TF occupancy | ARID5B CUT&RUN/CUT&Tag/ChIP | H1 | test adipocyte-repressor reuse in PH |
| TF occupancy | CUX1 or suitable state proxy | H2a | test exact-site alternative/combinatorial factor model |
| 3D contact | element -> `Irx3` promoter Capture-C/4C | H4 | test allele-dependent target contact |
| lineage/state | cell-type + developmental/state markers | H7/H9/H10 | locate effect context and gating |

## Required contrasts

### Contrast A — allele effect in male PH

```text
C/C male PH vs T/T male PH
```

Must reproduce `Irx3` increase in the same material used for mediator measurements.

### Contrast B — sex interaction

```text
(C/C - T/T) male
vs
(C/C - T/T) female
```

This tests the **allele-by-sex interaction**, not merely baseline male/female differences.

### Contrast C — within-cell-context mediator

For each responding cell population:

```text
allele -> mediator M -> Irx3
```

Do not pool nonresponding populations into the same estimate.

## Decision table

| Observation | Update |
|---|---|
| C allele lowers ARID5B occupancy; accessibility/contact/Irx3 rise | H1 strongly rises |
| C allele changes CUX1 occupancy/state; ARID5B stable | H2a rises |
| C allele changes accessibility but neither candidate TF occupancy changes | H3 rises; broaden TF search |
| C allele changes enhancer-Irx3 contact with little local occupancy/accessibility change | H4 rises |
| effect exists only in one cell subtype | cell identity becomes part of causal edge |
| effect is male-specific at mediator layer | sex-dependent cofactor/state becomes load-bearing |
| adult mediator absent but developmental state differs persistently | H9 developmental memory rises |
| no local mediator differs despite Irx3 change | H5 indirect mediator / composition / distal regulation rises |

## Intervention tier

Observation is not enough for the final edge. After prioritizing mediator `M`:

1. perturb `M` without changing rs1421085;
2. measure `Irx3` and the relevant chromatin/contact state;
3. rescue `M` in the risk-allele background when technically possible;
4. preserve negative outcomes.

### Strong causal pattern

```text
C allele -> M changes -> Irx3 up
M rescue/normalization -> Irx3 returns toward T-allele state
```

### Strong falsification pattern

```text
C allele -> M changes
but manipulating M does not alter Irx3
```

Then `M` is correlated with the allele but is not load-bearing for the `Irx3` edge under the tested context.

## Anti-confounding requirements

- Randomize processing order across genotype and sex.
- Use biological replicates; cells are nested within animals and are not independent animals.
- Blind genotype during primary image/peak/contact quantification where practical.
- Predefine cell-population inclusion rules.
- Report failed/negative mediator tests.
- Record genome assembly, exact orthologous variant coordinates, reagents, batches and environmental conditions.
- Do not compare a deletion perturbation to a point mutation as if they were equivalent interventions.

## Success criterion

The experiment succeeds scientifically even if every current candidate mediator fails, provided it produces a verified statement such as:

> In adult male PH cell population X, rs1421085 C increases Irx3, but no allele-dependent change was detected in ARID5B occupancy, CUX1 occupancy, local accessibility, or enhancer-Irx3 contact under condition Y.

That negative result removes major causal paths and defines the next search space.
