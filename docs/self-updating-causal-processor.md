# Self-Updating Causal Processor

CAUSAL-DNA treats causal knowledge as an append-only event stream plus rebuildable projections.

```text
material outcome
    -> independent verification
    -> append causal event
    -> generation N+1
    -> rebuild projection
    -> open/reject/supersede hypotheses
    -> new projective space
```

The event stream is history. Projections are disposable views.

## Invariants

1. Existing event IDs cannot be rewritten.
2. New events cannot be appended into an older generation.
3. Observation, perturbation and verification events require provenance.
4. Verification must be asserted by `independent_verifier`.
5. Model updates create a new generation rather than mutating historical material evidence.
6. Rejected and superseded hypotheses remain in the event history.
7. `GAP-001` remains open until a real verification event closes the missing edge.

## Generation semantics

- generation 0: initial material observations / interventions;
- generation 1: first model projection and competing hypotheses;
- generation N+1: created only after new evidence or a model update caused by evidence.

A projection can always be reconstructed by replaying `cases/CDNA-001.events.json` through `causal_dna.processor.CausalProcessor`.

This follows the same architectural principle used in reliable state systems: authoritative append-only records, rebuildable derived projections.
