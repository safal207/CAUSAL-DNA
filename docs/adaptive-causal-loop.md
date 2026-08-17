# Adaptive Causal Replanning Loop

CAUSAL-DNA separates three operations that are often accidentally mixed:

1. **evidence update** — append-only material/experimental/verification events;
2. **planning update** — change hypothesis weights used to choose the next test;
3. **causal promotion** — change the status of a causal edge only after evidence and independent verification.

Adaptive replanning belongs to operation 2, but it can now be triggered by a real operation-1 `experiment_outcome` event without collapsing the two layers.

## Loop

```text
current planning weights
        |
        v
rank candidate experiments
        |
        v
perform experiment
        |
        v
EXPERIMENT_OUTCOME event + provenance
        |
        +----------------------> append-only processor generation N+1
        |
        v
apply outcome planning multipliers
        |
        v
recompute entropy + information gain
        |
        v
select next experiment
        |
        X
no causal promotion without independent verification
```

The planning state can change while the authoritative causal projection remains unchanged.

## Two execution modes

### Simulation mode

`cases/CDNA-001.replanning-scenarios.json` contains explicit `simulation_only=true` branches. They test how the planner would behave under opposite hypothetical outcomes. They are never appended to `cases/CDNA-001.events.json` and therefore never become biological evidence.

### Real event mode

A real result must arrive as a causal-processor event:

```json
{
  "event_type": "experiment_outcome",
  "observer": "experiment",
  "generation": 2,
  "evidence_refs": ["..."],
  "payload": {
    "experiment_id": "E1_ARID5B_OCCUPANCY",
    "outcome_name": "no_allele_specific_occupancy"
  }
}
```

`AdaptiveReplanner.replan_from_event(...)` requires:

- the event targets exactly the next processor generation;
- `observer=experiment`;
- non-empty evidence provenance;
- a known experiment and declared outcome;
- the experiment has not already been consumed in the same planning history.

The event is appended to the processor, but `experiment_outcome` is deliberately ignored by the causal-status projection. The implementation checks that the sets of open, rejected, superseded, and verified hypotheses are identical before and after the replan.

## Weight semantics

Each experiment outcome contains `posterior_weights`. In adaptive mode these are interpreted as likelihood-like planning multipliers:

```text
new_weight[h] = current_weight[h] * outcome_multiplier[h]
```

The weights are normalized only when entropy is computed. They are not probabilities of biological truth and must not be cited as evidence.

With the initial equal planning weights, this interpretation reproduces the original information-gain ranking. After an outcome, the same experiment model is rescored against the updated planning state.

## Generations

Replanning advances monotonically:

```text
generation 1 = initial planning / processor projection
generation 2 = after first real outcome
generation 3 = after second outcome
...
```

For real event mode, planning generation and authoritative processor generation must match before the next event is accepted.

## Safety / integrity boundary

Adaptive replanning is forbidden from setting any of the following:

- `cause_found`
- `causal_status`
- `edge_status`
- `materialized`
- `verification_status`

A high planning weight means only “prioritize tests around this route.” It never means “this route is causal.”

The evidence path remains:

```text
experiment_outcome
      -> evidence review
      -> Bardo resolution candidate
      -> independent_verifier
      -> verification event
      -> only then causal promotion may be considered
```

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

The same declared outcome vocabulary is used by real `experiment_outcome` events, so simulation and execution share one planning model while remaining different epistemic states.

## Implementation

- `causal_dna/adaptive_replanner.py`
- `causal_dna/experiment_selector.py`
- `causal_dna/processor.py`
- `schemas/adaptive-replan.schema.json`
- `schemas/causal-event.schema.json`
- `cases/CDNA-001.replanning-scenarios.json`
- `tests/test_adaptive_replanner.py`
- `tests/test_processor.py`
- `analysis/validate_repository.py`
