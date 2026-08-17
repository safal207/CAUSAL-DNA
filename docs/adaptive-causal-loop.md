# Adaptive Causal Replanning Loop

CAUSAL-DNA now separates three operations that are often accidentally mixed:

1. **evidence update** — append-only material/verification events;
2. **planning update** — change hypothesis weights used to choose the next test;
3. **causal promotion** — change the status of a causal edge only after evidence and independent verification.

Adaptive replanning belongs only to operation 2.

## Loop

```text
current planning weights
        |
        v
rank candidate experiments
        |
        v
perform / simulate an outcome
        |
        v
apply outcome likelihood multipliers
        |
        v
new planning generation
        |
        v
recompute entropy + information gain
        |
        v
select next experiment
```

The planning state can change while the authoritative causal state remains unchanged.

## Weight semantics

Each experiment outcome contains `posterior_weights`. In adaptive mode these are interpreted as likelihood-like planning multipliers:

```text
new_weight[h] = current_weight[h] * outcome_multiplier[h]
```

The weights are normalized only when entropy is computed. They are not probabilities of biological truth and must not be cited as evidence.

With the initial equal planning weights, this interpretation reproduces the original information-gain ranking. After an outcome, the same experiment model can be replayed against the updated planning state.

## Generations

Replanning outcomes advance monotonically:

```text
generation 1 = initial planning state
generation 2 = after first outcome
generation 3 = after second outcome
...
```

A completed experiment is excluded from immediate reselection in the same replan history.

## Safety / integrity boundary

Adaptive replanning is forbidden from setting any of the following:

- `cause_found`
- `causal_status`
- `edge_status`
- `materialized`
- `verification_status`

Those belong to the evidence processor and independent verification path.

The repository includes explicit simulation scenarios in `cases/CDNA-001.replanning-scenarios.json`. They are marked `simulation_only=true`; their purpose is to test the planner's behavior under opposite hypothetical outcomes without pretending either result occurred biologically.

## CDNA-001 example

The first simulations branch on the ARID5B occupancy assay:

```text
E1_ARID5B_OCCUPANCY
        |
        +--> allele_specific_occupancy
        |        -> increase ARID5B planning weight
        |        -> recompute next experiment
        |
        +--> no_allele_specific_occupancy
                 -> decrease ARID5B planning weight
                 -> recompute next experiment
```

The resulting next-test choice is printed by the repository validator and is therefore CI-visible and reproducible.

## Implementation

- `causal_dna/adaptive_replanner.py`
- `causal_dna/experiment_selector.py`
- `schemas/adaptive-replan.schema.json`
- `cases/CDNA-001.replanning-scenarios.json`
- `tests/test_adaptive_replanner.py`
- `analysis/validate_repository.py`
