# CAUSAL-DNA

**CAUSAL-DNA** is an open research project for turning DNA variant associations into explicit, testable causal chains.

The core question is not only:

> Which variants correlate with a phenotype?

but:

> Through which molecular and cellular transitions does a specific DNA variant cause a measurable phenotype, and where does the evidence stop?

## Core model

```text
Variant
  -> regulatory / molecular effect
  -> target gene or transcript effect
  -> cellular state transition
  -> tissue / circuit effect
  -> organismal phenotype
  -> independent verification
```

Each arrow is a **causal edge**. A causal edge is never promoted to `established` only because two nodes are associated.

Every edge should carry:

- `claim`
- `from` / `to`
- `context` (species, tissue, cell type, sex, developmental state)
- `evidence_type`
- `intervention`
- `source`
- `confidence`
- `alternative_explanations`
- `falsification_test`
- `status`

See [`protocol/causal-proof-v0.1.md`](protocol/causal-proof-v0.1.md).

## Three-space causal graph

CAUSAL-DNA models the **transition of causal knowledge** through three graph spaces:

```text
PROJECTIVE
possible mechanism / hypothesis / planned experiment
        |
        v
BARDO
competing, partial, conflicting or unresolved transition state
        |
        v
MATERIAL
observation / perturbation / rescue with provenance
```

This is implemented as an executable graph, not a hand-maintained diagram:

- [`causal_dna/space_graph.py`](causal_dna/space_graph.py) — graph engine, traversal, materialization frontier and Mermaid renderer;
- [`schemas/space-graph.schema.json`](schemas/space-graph.schema.json) — machine-readable contract;
- [`cases/CDNA-001-rs1421085.space-graph.json`](cases/CDNA-001-rs1421085.space-graph.json) — first real three-space case graph;
- [`docs/three-space-causal-graph.md`](docs/three-space-causal-graph.md) — semantics and integrity rules.

A projective hypothesis is forbidden from jumping directly into material space. It must cross an explicit Bardo state, and every material node/transition requires evidence provenance. Rejected and superseded Bardo paths are preserved rather than deleted.

Example:

```bash
python -m causal_dna.space_graph \
  cases/CDNA-001-rs1421085.space-graph.json \
  --start P_H3_ACCESS \
  --end M_IRX3_UP \
  --contours
```

Mermaid can be generated from the same JSON source of truth:

```bash
python -m causal_dna.space_graph \
  cases/CDNA-001-rs1421085.space-graph.json \
  --mermaid-out /tmp/cdna-001.mmd
```

## 4D causal lattice

The three-space graph is extended with explicit **time, evidence strength, and observer role**:

```text
State = Space × Time × Evidence × Observer
```

The temporal axis distinguishes:

```text
past    = already observed evidence
present = current unresolved causal boundary
future  = planned tests / possible resolutions
```

A future hypothesis cannot be represented as a material fact. The current CDNA-001 lattice therefore stops at a future independent Bardo resolution gate while `GAP-001` remains open.

Implementation:

- [`causal_dna/temporal_lattice.py`](causal_dna/temporal_lattice.py) — 4D lattice engine;
- [`schemas/temporal-lattice.schema.json`](schemas/temporal-lattice.schema.json) — lattice contract;
- [`cases/CDNA-001-rs1421085.temporal-lattice.json`](cases/CDNA-001-rs1421085.temporal-lattice.json) — first lattice instance;
- [`docs/4d-causal-lattice.md`](docs/4d-causal-lattice.md) — semantics and invariants.

The observer axis separates `model`, `experiment`, and `independent_verifier`, so the same component that proposes a mechanism cannot silently certify it as established.

## Self-updating causal processor

The authoritative research history is append-only:

```text
Evidence -> Event -> Generation -> Rebuildable Projection
```

`cases/CDNA-001.events.json` is replayed by [`causal_dna/processor.py`](causal_dna/processor.py). New generations may add observations, perturbations, experiment outcomes, model updates, hypothesis changes, or independent verification, but historical events cannot be rewritten.

The current model projection is therefore disposable and reproducible: if a projection is lost, it can be rebuilt from the event history.

See [`docs/self-updating-causal-processor.md`](docs/self-updating-causal-processor.md).

## Active experiment selection

CAUSAL-DNA ranks candidate experiments by expected reduction in planning uncertainty:

```text
E* = argmax_E E[ H(before) - H(after | E) ]
```

The planning weights are **not biological evidence**. They exist only to prioritize experiments.

The first CDNA-001 plan compares ARID5B occupancy, CUX1 occupancy, allele-specific chromatin accessibility, enhancer-to-Irx3 3D contact, and a matched multi-omic design.

Implementation:

- [`causal_dna/experiment_selector.py`](causal_dna/experiment_selector.py)
- [`schemas/experiment-plan.schema.json`](schemas/experiment-plan.schema.json)
- [`cases/CDNA-001.experiment-plan.json`](cases/CDNA-001.experiment-plan.json)
- [`docs/active-causal-discovery.md`](docs/active-causal-discovery.md)

## Adaptive replanning

The selector now closes the sequential loop:

```text
rank experiments
      -> experiment outcome
      -> generation N+1
      -> update planning weights
      -> recompute entropy / information gain
      -> choose next experiment
```

There are two separate modes:

- **simulation mode** — `simulation_only=true` scenarios test how ranking would change under hypothetical outcomes without touching evidence history;
- **real event mode** — an evidence-backed `experiment_outcome` from the `experiment` observer is appended to the next processor generation and triggers replanning.

A real experiment outcome can change planning priority, but it still cannot set `cause_found`, materialize a Bardo path, or create a verified causal edge. Those transitions remain behind independent verification.

Implementation:

- [`causal_dna/adaptive_replanner.py`](causal_dna/adaptive_replanner.py)
- [`schemas/adaptive-replan.schema.json`](schemas/adaptive-replan.schema.json)
- [`schemas/causal-event.schema.json`](schemas/causal-event.schema.json)
- [`cases/CDNA-001.replanning-scenarios.json`](cases/CDNA-001.replanning-scenarios.json)
- [`docs/adaptive-causal-loop.md`](docs/adaptive-causal-loop.md)

## Multi-step research strategy

CAUSAL-DNA can also search a short **adaptive policy tree** under a declared budget:

```text
experiment E1
   |
   +-- outcome A -> experiment E2
   |
   +-- outcome B -> experiment E4
```

Instead of asking only “what should we test next?”, this asks:

> Which contingent 2–3 step strategy is expected to reduce the most uncertainty without exceeding the research budget?

The planner tracks expected terminal entropy, expected information gain, expected cost, maximum path cost, and gain per cost. It performs exact short-horizon enumeration with Pareto pruning. The horizon is intentionally capped at three steps so planning assumptions do not compound indefinitely before new evidence arrives.

This layer is still planning-only. It cannot write `cause_found`, `causal_status`, `materialized`, or `verification_status`.

Implementation:

- [`causal_dna/strategy_planner.py`](causal_dna/strategy_planner.py)
- [`schemas/strategy-config.schema.json`](schemas/strategy-config.schema.json)
- [`cases/CDNA-001.strategy-configs.json`](cases/CDNA-001.strategy-configs.json)
- [`docs/multi-step-research-strategy.md`](docs/multi-step-research-strategy.md)

## First case: CDNA-001 — rs1421085

The obesity-associated non-coding variant **rs1421085 T>C** is a useful benchmark because one branch is unusually well supported experimentally, while another branch remains incomplete.

### Adipocyte branch — strong mechanistic evidence

```text
rs1421085 T>C
  -> disruption of an ARID5B repressor motif
  -> IRX3 / IRX5 derepression in adipocyte progenitors
  -> shift away from thermogenic beige adipocytes
  -> mitochondrial thermogenesis decreases
  -> lipid storage increases
```

The 2015 study used directed perturbations and endogenous CRISPR-Cas9 editing to test this pathway.

### CNS branch — open causal gap

A 2025 mouse study reported that the obesity-associated rs1421085 allele increases brain `Irx3`, including in the posterior hypothalamus (PH), and that increasing IRX3 in the PH is sufficient to increase food intake and body weight. Increased IRX3 also reduced excitability/activity of IRX3-positive PH neurons.

The unresolved edge we start with is:

```text
rs1421085
  -> ???
  -> IRX3 up in posterior hypothalamus
  -> PH neuronal activity down
  -> food intake up
  -> body weight up
```

The project will first test whether the known adipocyte regulatory mechanism generalizes to the relevant CNS cell types, or whether a different regulatory mechanism explains this branch.

See [`cases/CDNA-001-rs1421085.md`](cases/CDNA-001-rs1421085.md).

## Evidence policy

We distinguish:

- **association** — statistical co-occurrence only;
- **supportive** — mechanistically compatible observational evidence;
- **perturbed** — intervention changes the proposed mediator or outcome;
- **rescued** — reversing/restoring the proposed causal factor reverses the effect;
- **established** — convergent evidence with context and alternatives explicitly tested.

A causal claim must preserve negative, conflicting, and superseded evidence instead of deleting it.

## Initial sources

1. Claussnitzer M, et al. *FTO Obesity Variant Circuitry and Adipocyte Browning in Humans.* N Engl J Med. 2015. PMID: 26287746. https://pubmed.ncbi.nlm.nih.gov/26287746/
2. Laber S, et al. *Linking the FTO obesity rs1421085 variant circuitry to cellular, metabolic, and organismal phenotypes in vivo.* Sci Adv. 2021. PMID: 34290091. https://pubmed.ncbi.nlm.nih.gov/34290091/
3. Sullivan AI, et al. *Mice harboring the obesity-associated SNP rs1421085 exhibit increased body weight and reveal an IRX3 neuronal circuit regulating body weight.* Mol Metab. 2025. PMID: 40835181. https://pubmed.ncbi.nlm.nih.gov/40835181/

## Research discipline

CAUSAL-DNA is a computational and evidence-mapping research project. It does **not** treat exploratory hypotheses as medical advice, diagnostic conclusions, or validated therapeutic targets.

---

**Status:** bootstrap / causal proof protocol + three-space graph + 4D lattice + append-only processor + active experiment selection + adaptive replanning + budgeted multi-step strategy planning
