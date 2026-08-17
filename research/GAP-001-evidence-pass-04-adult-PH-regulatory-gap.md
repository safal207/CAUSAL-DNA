# GAP-001 Evidence Pass 04 — the adult posterior-hypothalamus regulatory gap

Date: 2026-08-17
Status: unresolved molecular edge isolated; decisive measurement specified

## The gap is now narrower than `rs1421085 -> Irx3`

The available evidence supports several boundary facts:

1. Exact `rs1421085 T>C` is causal for an ARID5B-dependent enhancer mechanism in human adipocyte progenitors (Claussnitzer et al., 2015; PMID 26287746).
2. Deleting the local rs1421085 CRM can alter enhancer-to-`Irx3/Irx5` contact in an adipocyte differentiation state, while having no detectable effect on bulk adult hypothalamic `Irx3/Irx5` expression (Laber et al., 2021; PMID 34290091).
3. The broader FTO obesity interval can contact `IRX3/IRX5` in brain/hypothalamic-neuron contexts and shows strong developmental/cell-state gating (Smemo et al., 2014; Sobreira et al., 2021).
4. The exact homologous `T>C` substitution increases `Irx3` in adult male posterior hypothalamus and increased PH neuronal IRX3 is sufficient to drive feeding/body-weight effects (Sullivan et al., 2025; PMID 40835181).

What remains missing is the molecular bridge inside the adult posterior-hypothalamic context:

```text
rs1421085 T>C
  -> [adult PH allele-dependent regulatory state ???]
  -> Irx3 transcription up
```

## Search result for the missing bridge

In the primary studies reviewed for GAP-001, we have **not identified direct adult-PH, exact-allele measurements** of all or any of the following at the rs1421085-containing element:

- chromatin accessibility (`T` vs `C`);
- ARID5B occupancy (`T` vs `C`);
- CUX1 occupancy / functional factor state (`T` vs `C`);
- enhancer-to-`Irx3` contact (`T` vs `C`);
- a cell-resolved mediation chain connecting one of those molecular changes to `Irx3` expression.

This statement is deliberately scoped to the sources/searches reviewed by CAUSAL-DNA; it is **not** a claim that no such dataset can exist anywhere.

## Why adipocyte contact evidence cannot fill this edge

Laber et al. showed that the rs1421085-DEL82 perturbation reduced enhancer interactions with `Irx3/Irx5` in an induced adipocyte state. The interaction effect itself was state-dependent.

Therefore:

```text
adipocyte enhancer-contact effect
  != proof of
adult PH exact-SNP enhancer-contact effect
```

The context is part of the causal claim.

## Why Sullivan 2025 does not fill this edge

Sullivan et al. provides a particularly valuable exact-SNP boundary condition:

```text
exact T>C genotype
  -> Irx3 increases in male PH
```

and a downstream sufficiency chain:

```text
PH neuronal Irx3 increase
  -> neuronal excitability/activity down
  -> feeding/body weight up
```

But these two experimentally supported boundaries still leave a molecular mediator gap between the DNA base and `Irx3` transcription.

Downstream neuronal sufficiency also does not by itself identify the native cell type in which the SNP first changes regulatory state.

## Updated causal graph

```text
                        ADIPOCYTE BRANCH
rs1421085 T>C -> ARID5B repression loss -> enhancer activity up -> IRX3/IRX5 up

                        CNS CONTEXT
FTO obesity interval -> long-range IRX3/IRX5 regulation [development/state gated]

                        ADULT MALE PH — GAP-001
rs1421085 T>C
   |
   +--> accessibility change? --------------------+
   |                                              |
   +--> ARID5B occupancy change? -----------------+
   |                                              |
   +--> CUX1/combinatorial occupancy change? -----+--> enhancer state/contact ? --> Irx3 up
   |                                              |
   +--> developmental-memory/state gate? ---------+
   |                                              |
   +--> cell-composition / indirect mediator? ----+

Irx3 up -> PH IRX3+ neuronal activity down -> food intake/body weight up
```

## Current hypothesis ranking

### H3 — allele-dependent accessibility

Status: `HIGH_VALUE_UNTESTED_MEDIATOR`.

Reason: exact point mutation can preserve the enhancer while changing its state; deletion-vs-substitution discordance is compatible with this model.

### H4 — allele-dependent 3D contact

Status: `LOCUS_CONTACT_SUPPORTED; EXACT_ALLELE_ADULT_PH_EFFECT_UNTESTED`.

Long-range FTO-locus/IRX3 contact is supported in brain/hypothalamic-neuron contexts, but the exact `T>C -> contact change` edge in adult PH is unresolved.

### H1 — ARID5B reuse

Status: `CELL_FEASIBLE_BUT_UNPROVEN_WITH_NEGATIVE_TENSION`.

`Arid5b` co-detection with `Irx3` is common enough in VPH cells to keep feasibility open, but DEL82 adult/bulk hypothalamus is negative and no adult-PH allele-specific ARID5B occupancy evidence has been identified.

### H2a — CUX1 exact-site role

Status: `EXACT_SITE_AND_CELL_FEASIBLE_WITH_TARGET_GENE_CONFLICT`.

CUX1 has exact-site allele-sensitive evidence and substantial VPH co-detection with `Irx3`, but existing neuronal CUX1 work points to an `FTO/RPGRIP1L` program rather than a general `IRX3` program.

### H9 — temporal/developmental gate

Status: `SUPPORTED_ARCHITECTURAL_MODIFIER`.

Developmental hypothalamus is much more sensitive to interval deletion than adult bulk hypothalamus. This can explain why target routing changes with state, but does not yet identify the adult PH molecular mediator.

## The smallest decisive experiment

The most informative next study is not another broad expression association. It is a **same-cell-context, exact-allele mediation experiment**.

### Factorial design

```text
allele: T/T vs C/C (heterozygous T/C useful for allele-specific assays)
sex: male vs female
region: posterior hypothalamus
cell population: cell-resolved; do not assume neuron before measurement
state: matched age, diet, temperature, hormonal/metabolic context
```

### Measurements from matched biological material

1. `Irx3` expression — single-cell/single-nucleus or spatial RNA.
2. Accessibility at the rs1421085 element — targeted ATAC or snATAC.
3. ARID5B occupancy — CUT&RUN / CUT&Tag / ChIP as technically appropriate.
4. CUX1 occupancy or factor-state proxy — same element.
5. Enhancer-to-`Irx3` contact — targeted Capture-C/4C or compatible single-cell/multiome approach.
6. Nearby alternative targets — `Fto`, `Rpgrip1l`, `Irx5` where biologically measurable.

### Causal discriminator

A useful mediation logic is:

```text
T>C
 -> molecular mediator M changes
 -> Irx3 changes
```

Then perturb `M` while preserving the SNP-containing element.

Examples:

- If `C` changes ARID5B occupancy and restoring ARID5B occupancy normalizes `Irx3`, H1 rises to causal.
- If `C` changes contact with `Irx3` without an ARID5B/CUX1 occupancy change, H4 rises and local single-TF stories fall.
- If the molecular change exists only in one PH cell subtype or only in males, that context becomes part of the causal edge.
- If no local regulatory measurement changes while `Irx3` changes, indirect mediator / developmental-memory hypotheses rise.

## Definition of a CAUSAL-DNA win

GAP-001 is not solved by finding another correlation around the FTO locus.

A meaningful resolution is:

```text
exact allele
  -> measured molecular state change
  -> target-gene change
```

in the relevant adult PH context, plus an intervention that breaks or rescues the proposed mediator edge.

Until then the correct status remains:

`GAP-001 = OPEN; MOLECULAR MEDIATOR ISOLATED AS THE NEXT LOAD-BEARING EDGE`.
