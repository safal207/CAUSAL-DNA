# 4D Causal Lattice

CAUSAL-DNA extends the three-space graph with a fourth, orthogonal structure: time, while also making evidence level and observer explicit coordinates of every state.

The lattice is epistemic, not metaphysical.

## Axes

A state is represented as:

```text
State = Space × Time × Evidence × Observer
```

### Space

- `projective` — hypotheses, planned experiments, possible causal structures;
- `bardo` — unresolved transition states where competing explanations are preserved;
- `material` — observations or interventions backed by explicit evidence provenance.

### Time

- `past` — already observed or historically established evidence;
- `present` — the currently unresolved research boundary;
- `future` — planned tests and possible resolutions.

Future states cannot be represented as material facts.

### Evidence

- `none`
- `supportive`
- `perturbed`
- `rescued`
- `established`

Strong evidence levels require provenance.

### Observer

- `model` — proposes or updates hypotheses;
- `experiment` — produces measurements or interventions;
- `independent_verifier` — evaluates whether the materialization criteria were actually met.

Observer separation is intentional: the model that proposes a causal route is not allowed to silently become the verifier that declares it established.

## CDNA-001 contour

```text
PAST / MATERIAL
rs1421085 exact intervention
        |
        v
Irx3 PH increase
        |
        +------> neuronal activity decrease
        |
        v
PRESENT / PROJECTIVE
What molecular mediator connects T>C to Irx3 up?
        |
        v
PRESENT / BARDO
ARID5B | CUX1 | accessibility | 3D contact | ...
        |
        v
FUTURE / PROJECTIVE
exact-allele, sex-stratified, cell-resolved decisive experiment
        |
        v
FUTURE / BARDO
independent resolution gate
        |
        X
No future material state exists until new evidence is produced.
```

The final materialization step is intentionally absent from the current lattice because `GAP-001` remains open and `cause_found=false`.

## Integrity rules

1. `projective -> material` bypass is forbidden.
2. Future states cannot be material facts.
3. Material states require evidence provenance.
4. `materializes` transitions must be `bardo -> material` and carry evidence references.
5. Strong evidence levels require provenance.
6. An independent-verifier state must carry the protocol/evidence basis it is expected to evaluate.
7. Historical material evidence may update the present model, but forward testing transitions cannot reverse temporal order.

## Why this matters

The ordinary causal graph asks:

> What is connected to what?

The three-space graph adds:

> Is this a possibility, an unresolved transition, or an observed fact?

The 4D lattice adds:

> When did this knowledge exist, how strong was the evidence, and from which observer role was the state asserted?

This prevents three common errors in causal research:

- treating a future hypothesis as if it were a current fact;
- treating historical evidence as automatically valid in a new biological context;
- allowing the proposing model to self-certify its own causal conclusion.

## Implementation

- engine: `causal_dna/temporal_lattice.py`
- schema: `schemas/temporal-lattice.schema.json`
- CDNA-001 instance: `cases/CDNA-001-rs1421085.temporal-lattice.json`
- tests: `tests/test_temporal_lattice.py`
- repository gate: `analysis/validate_repository.py`
