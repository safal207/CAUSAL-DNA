# Active Causal Discovery

CAUSAL-DNA can rank planned experiments by how much uncertainty they are expected to remove from the current hypothesis space.

This layer is a **planning engine**, not an evidence engine.

## Objective

For open hypotheses `H`, define planning entropy:

```text
Entropy(H) = -sum p(h) log2 p(h)
```

For candidate experiment `E` with possible outcomes `o`:

```text
ExpectedPosteriorEntropy(E)
  = sum P_plan(o | E) * Entropy(H | o)
```

and expected information gain:

```text
IG(E) = Entropy(H_before) - ExpectedPosteriorEntropy(E)
```

The selector can optimize either:

```text
max IG(E)
```

or a simple resource-aware objective:

```text
max IG(E) / cost(E)
```

## Critical epistemic boundary

The probabilities and hypothesis weights in `cases/CDNA-001.experiment-plan.json` are explicit **planning priors only**.

They are not:

- measured biological probabilities;
- posterior causal confidence;
- evidence for any causal edge;
- permission to materialize a Bardo state;
- a substitute for perturbation, rescue, or independent verification.

A selector score can choose the next test. It cannot change `cause_found`, causal-edge status, the append-only material evidence log, or independent-verifier state.

## CDNA-001 candidate tests

The first plan compares:

- exact-allele ARID5B occupancy;
- exact-allele CUX1 occupancy;
- sex-stratified cell-resolved chromatin accessibility;
- enhancer-to-Irx3 contact;
- a matched multi-omic experiment measuring RNA, accessibility, TF occupancy, and 3D contact together.

The declared priors are deliberately transparent and editable. Changing them changes the **planning recommendation**, not the scientific record.

## Integration with the causal processor

```text
append-only evidence history
          |
          v
current generation / projection
          |
          v
open hypothesis set
          |
          v
active experiment selector
          |
          v
ranked planned tests
          |
          X
(no causal-status mutation)

real experiment
          |
          v
new evidence event
          |
          v
independent verification
          |
          v
new generation
```

This keeps active learning downstream of the evidence projection and upstream of future experiments, never inside the materialization boundary.

## Files

- engine: `causal_dna/experiment_selector.py`
- schema: `schemas/experiment-plan.schema.json`
- CDNA-001 plan: `cases/CDNA-001.experiment-plan.json`
- tests: `tests/test_experiment_selector.py`
- repository gate: `analysis/validate_repository.py`
