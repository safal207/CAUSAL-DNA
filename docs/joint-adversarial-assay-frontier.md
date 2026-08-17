# Joint adversarial assay frontier

CAUSAL-DNA now supports a planning layer where the cost and reliability of several experiments can move at the same time.

This is deliberately different from a single-target phase map. A single-target map answers:

> What happens if one assay becomes cheaper, more expensive, stronger, or weaker while competitors stay fixed?

The joint frontier asks:

> Which first experiment minimizes worst-case planning regret when several assays may simultaneously have different cost and reliability states?

## Current CDNA-001 joint surface

The first joint pass varies:

- `E1_ARID5B_OCCUPANCY`
- `E3_ALLELE_SNATAC`
- `E5_MATCHED_MULTIOMIC`

Each receives an independent state from:

```text
cost multiplier: 0.75, 1.00, 1.25
reliability:     0.70, 0.85, 1.00
```

That is nine states per assay and therefore:

```text
9^3 = 729 joint assay configurations
```

Each configuration is evaluated across a deterministic 165-point simplex of planning priors. Outcome probabilities are also allowed a bounded mass transfer of `0.10`.

`E2_CUX1_OCCUPANCY`, `E4_CAPTURE_C`, and `STOP` remain explicit competitors at baseline settings rather than being silently removed.

## Decision rule

For every prior and assay configuration, each experiment receives a conservative utility interval:

```text
utility = information gain / effective cost
```

The interval comes from the admissible outcome-probability shifts. Reliability shrinks the outcome posterior multipliers toward neutral `1.0`, so a less reliable assay is modeled as less discriminating.

For action `a` at a prior state:

```text
interval_regret(a) = best_possible_upper_utility - lower_utility(a)
```

The joint-cell recommendation is the action minimizing maximum interval regret across the prior grid. Ties are resolved deterministically.

This is a robust planning recommendation only. It is not a biological probability or evidence grade.

## Joint boundaries

The analyzer records adjacent cells where the winning first action changes. Because the surface is six-dimensional — cost and reliability for each of three assays — boundaries are reported both as action-pair counts and by the coordinate axis along which the switch occurs.

This lets the repository answer questions such as:

- Is the E1/E5 switch mostly controlled by E5 cost or by E1 reliability?
- Does E3 appear only when accessibility measurement becomes relatively cheap and reliable?
- Which assay uncertainty contributes most to recommendation instability?

## Policy Pareto frontier

A separate view uses a declared reference assay profile and the existing short-horizon strategy planner.

For horizons `1`, `2`, and `3`, candidate adaptive policies are annotated with:

- expected information gain;
- expected cost;
- maximum path cost;
- first-step maximum robust regret;
- planning horizon.

A policy is removed when another policy is no worse in all three optimization dimensions:

```text
higher expected information gain
lower expected cost
lower first-step robust regret
```

Planning horizon is retained as an explanatory dimension rather than treated as automatically better or worse. This is important because the earlier CDNA-001 analyses already showed that horizon can change the preferred first experiment.

## Scientific boundary

The joint frontier cannot write or modify:

- `cause_found`
- `causal_status`
- `edge_status`
- `materialized`
- `verification_status`

No number of planning-grid evaluations can turn a hypothesis into evidence. New causal status still requires material evidence and the independent-verification path enforced by the causal processor.
