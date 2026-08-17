# Three-Space Causal Graph

CAUSAL-DNA models a second dimension in addition to the ordinary causal chain:
**the state of causal knowledge itself**.

The three spaces are formal graph semantics. They are useful for reasoning about
how a possible mechanism becomes an evidence-backed mechanism without silently
turning a hypothesis into a fact.

## 1. Projective space — possibility and design

Projective nodes encode things we can coherently propose before they are proven:

- hypotheses;
- candidate mediators;
- predicted transitions;
- planned experiments;
- alternative causal programs.

A projective node may be elegant, biologically plausible and well motivated. It
is still not a material fact.

## 2. Bardo — transition and competition

Bardo is the explicit liminal layer between a proposed mechanism and an
observed mechanism. It stores:

- competing causal explanations;
- partial support;
- conflicting evidence;
- unresolved cell-of-action questions;
- experiment gates;
- rejected and superseded paths.

Every Bardo node has a resolution state:

- `open`
- `materialized`
- `rejected`
- `superseded`

Rejected paths are preserved rather than deleted. This makes the graph an audit
trail of causal search, not only a picture of the winning story.

## 3. Material space — evidence-backed state

Material nodes represent observations or interventions with provenance:

- measured molecular changes;
- perturbations;
- rescues;
- circuit states;
- organismal phenotypes.

A material node **must** contain `evidence_refs`.

An edge entering material space **must** also contain `evidence_refs`.

The graph engine forbids a direct `projective -> material` edge. A hypothesis
cannot become a fact without an explicit Bardo transition.

## CDNA-001 contour

The current case has a material backbone:

```text
engineered rs1421085 allele
  -> observed PH Irx3 increase
  -> PH neuronal activity decrease
  -> food intake increase
  -> body-weight increase
```

But the molecular mediator between the first two material states is still open.
That gap is represented as a second contour:

```text
MATERIAL BOUNDARY
rs1421085 exact edit
        |
        v
PROJECTIVE
GAP-001 question
  |-- ARID5B reuse
  |-- CUX1 exact-site state
  |-- chromatin accessibility
  |-- enhancer-Irx3 contact
  |-- composition/state shift
  |-- sex/context modifier
  |-- intact enhancer grammar
  `-- temporal gating
        |
        v
BARDO
competing evidence-bearing states
        |
        v
B_INTEGRATION_GATE
        |
        | exact allele x sex x cell-state experiment
        | + mediator perturbation/rescue
        v
MATERIAL BOUNDARY
observed Irx3 increase
```

The lower material boundary does **not** mean the mediator has already
materialized. It constrains the search: any winning Bardo path must explain an
already observed total genotype-to-expression effect.

## Executable API

Load and inspect the graph:

```bash
python -m causal_dna.space_graph cases/CDNA-001-rs1421085.space-graph.json
```

Ask for a full three-space contour:

```bash
python -m causal_dna.space_graph \
  cases/CDNA-001-rs1421085.space-graph.json \
  --start P_H3_ACCESS \
  --end M_IRX3_UP \
  --contours
```

Generate Mermaid from the same source-of-truth JSON:

```bash
python -m causal_dna.space_graph \
  cases/CDNA-001-rs1421085.space-graph.json \
  --mermaid-out /tmp/cdna-001.mmd
```

No hand-maintained diagram is required. The renderer and verifier consume the
same graph document.

## Integrity rules

The implementation currently enforces:

1. unique node and edge IDs;
2. every edge references existing nodes;
3. all material nodes have provenance;
4. every edge entering material space has provenance;
5. `projective -> material` bypass is forbidden;
6. `materializes_as` is valid only for `bardo -> material`;
7. rejected/superseded Bardo states require a reason;
8. materialized Bardo states require `resolved_by` evidence;
9. GAP-001 cannot be open while its integration gate is resolved;
10. each encoded GAP-001 projective hypothesis must have a full
    `projective -> bardo -> material` contour to the downstream material boundary.

These rules are checked in CI by `analysis/validate_repository.py` and
`tests/test_space_graph.py`.

## Why this matters

A conventional causal graph represents a proposed world:

```text
A -> B -> C
```

The three-space graph also represents **epistemic motion**:

```text
possible A->B
    -> contested/tested A->B
        -> observed/rescued A->B
```

This lets CAUSAL-DNA reason not only about biological causality, but about the
controlled transition from hypothesis to evidence.
