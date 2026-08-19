# DNA-RNA Exchange Process Graph, Causal Orientation Kernel, and Boundary

## Purpose

CAUSAL-DNA treats the relationship between DNA state and RNA state as a network of exchange processes rather than a single arrow or an unspecified field.

The model separates three questions:

1. **EPG — Exchange Process Graph:** which molecular, chemical, electrical, mechanical, spatial, and feedback processes can connect a DNA state to an RNA state?
2. **COK — Causal Orientation Kernel:** which smallest experimentally supported intermediate set intercepts every currently validated causal path from the source DNA state to the target RNA state?
3. **COB — Causal Orientation Boundary:** which multivariate state transition precedes the RNA trajectory, can be experimentally crossed, can be reversed/rescued, and survives independent verification?

These are scientific contracts, not claims that a universal orientation centre has already been discovered.

## Exchange Process Graph

An EPG edge has one of four meanings:

- `CAUSES` — intervention-supported causal change;
- `ENABLES` — a state is required or permissive for a downstream transition;
- `MODULATES` — a state changes the magnitude/probability of a downstream transition;
- `FEEDBACK` — a downstream state modifies a prior regulatory state.

Each edge is separately labelled:

- `HYPOTHESIS`
- `OBSERVATIONAL`
- `INTERVENTION_SUPPORTED`
- `RESCUE_SUPPORTED`
- `VERIFIED`

Only intervention-supported or stronger forward edges enter COK path computation. Observational centrality is never enough.

## Causal Orientation Kernel

For source `D` and target RNA state `R`, let `P(D,R)` be the set of currently validated causal paths. A candidate kernel `K` must intersect every validated path:

`for every p in P(D,R): K intersects interior(p)`

The detector returns the smallest such sets.

If a direct validated `D -> R` edge exists and bypasses every candidate mediator, no intermediate COK can be claimed. This is deliberate: a known genotype-to-expression boundary effect does not identify the missing mechanism.

Independent verification is stricter still: every edge on every currently relevant causal path must be `VERIFIED` before the kernel certificate becomes `VERIFIED_KERNEL`.

## Causal Orientation Boundary

COK identifies a causal cut; COB asks whether a state transition behaves like an experimentally real boundary.

A boundary progresses only through this evidence chain:

```text
PROJECTIVE
  -> TEMPORALLY_SUPPORTED
  -> INTERVENTION_SUPPORTED
  -> REVERSIBLE
  -> VERIFIED
```

Requirements:

1. **Temporal precedence** — the candidate state changes before the target RNA trajectory.
2. **Boundary-crossing intervention** — forcing the state across the proposed boundary changes the target RNA trajectory.
3. **Reversal/rescue** — reversing the state in the causal background moves the RNA trajectory back as predicted.
4. **Independent verification** — the transition replicates independently in the same declared context.

Prediction, entropy reduction, graph centrality, or a visually sharp threshold can prioritize a boundary experiment, but cannot satisfy these requirements.

## CDNA-001 / GAP-001

Current source and target:

```text
exact rs1421085-ortholog T>C
  -> ???
  -> Irx3 transcription increase
```

The current EPG records:

- an intervention-supported exact-edit boundary effect from the variant to increased `Irx3`;
- hypothesis-level local sequence grammar and ARID5B/CUX1/reader occupancy routes;
- observational adult-hypothalamic exact-locus accessibility;
- observational wild-type structural compatibility with the `Irx3` promoter;
- explicit placeholders for chemical, ionic, electrical, or mechanical signalling without pretending that rs1421085-specific evidence exists for them.

Because the known exact-edit source-to-target effect still bypasses every unverified mechanistic candidate, the present COK result is:

```text
kernel_found = false
```

Three COB candidates are retained:

- `COB_A_OCCUPANCY` — exact-site occupancy switch;
- `COB_B_CHROMATIN_CONTACT` — joint accessibility/contact switch;
- `COB_C_DISTRIBUTED` — distributed occupancy + chromatin + contact switch.

All remain `PROJECTIVE` because public datasets do not provide the required temporal, crossing-intervention, reversal, and independent-verification chain.

## Relation to CMG and NSG

```text
CMG
  identifies whether a mediator can cross from association toward causal materialization
    ↓
NSG
  classifies necessity, sufficiency, partial mediation, redundancy, and residual paths
    ↓
EPG
  assembles validated and unresolved exchange processes into one typed graph
    ↓
COK
  searches for the smallest causal cut across validated DNA -> RNA paths
    ↓
COB
  tests whether the proposed orienting state is an actual reversible transition boundary
```

## Scientific claim ceiling

The framework can currently say that a candidate orientation kernel or boundary is **defined and testable**.

It cannot yet say that a universal DNA-RNA orientation centre exists.

A broad principle would require successful COK/COB identification across multiple independent variants, genes, cell contexts, and perturbation systems, followed by out-of-sample prediction and independent replication.
