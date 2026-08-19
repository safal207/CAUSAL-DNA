# Joint adversarial assay frontier

CAUSAL-DNA now supports a planning layer where the cost and reliability of several experiment families can move at the same time.

This is deliberately different from a single-target phase map. A single-target map asks:

> What happens if one assay becomes cheaper, more expensive, stronger, or weaker while competitors stay fixed?

The joint frontier asks:

> Which first experiment minimizes worst-case planning regret when several assay families may simultaneously have different cost and reliability states?

## Current CDNA-001 joint surface

The first joint pass exposes three independent stress axes:

- `E1_ARID5B_OCCUPANCY`
- `E3_ALLELE_SNATAC`
- `E5_MATCHED_MULTIOMIC`

Each axis receives an independent state from:

```text
cost multiplier: 0.75, 1.00, 1.25
reliability:     0.70, 0.85, 1.00
```

That is nine states per axis and therefore:

```text
9^3 = 729 joint assay configurations
```

Each configuration is evaluated across a deterministic 165-point simplex of planning priors. Outcome probabilities are also allowed a bounded mass transfer of `0.10`.

## Linked assay-family settings

An initial unconstrained version stressed E1 while leaving the structurally symmetric E2 planning assay at baseline, and likewise stressed E3 while leaving E4 fixed. That created an unfair decision artifact: the unstressed twin could dominate simply because its cost/reliability assumptions were protected.

The production joint pass therefore declares:

```text
E2_CUX1_OCCUPANCY <- shares E1_ARID5B_OCCUPANCY cost/reliability axis
E4_CAPTURE_C      <- shares E3_ALLELE_SNATAC cost/reliability axis
```

This linking is about the **current synthetic planning model only**. It does not state that ARID5B and CUX1 are biologically equivalent, or that snATAC and Capture-C measure the same mechanism. The hypotheses and actions remain separate. The shared axis only removes an artificial advantage caused by stress-testing one member of a currently symmetric cost/likelihood pair while freezing the other.

`STOP` also remains an explicit action.

## Decision rule

For every prior and joint assay configuration, each experiment receives a conservative utility interval:

```text
utility = information gain / effective cost
```

The interval comes from admissible outcome-probability shifts. Reliability shrinks the outcome posterior multipliers toward neutral `1.0`, so a less reliable assay is modeled as less discriminating.

For action `a` at a prior state:

```text
interval_regret(a) = best_possible_upper_utility - lower_utility(a)
```

The joint-cell recommendation is the action minimizing maximum interval regret across the prior grid. Ties are resolved deterministically.

This is a robust planning recommendation only. It is not a biological probability or evidence grade.

## Current result

After the linked-family correction, the declared 729-cell box is one connected robust region:

```text
E5_MATCHED_MULTIOMIC = 729 / 729 cells
boundary pairs       = none
boundary axes        = none
```

The calculation covers:

```text
729 cells x 165 priors = 120,285 prior-cell evaluations
```

The absence of an internal boundary is itself a result. It means the current coarse uncertainty box does not reach a planning state where another first action has lower maximum interval regret than E5.

It does **not** imply that E5 is universally optimal. The nearest observed cells have only about `0.007654 bits/cost-unit` of regret margin over E1 or E4. A recommendation switch may therefore lie outside the scanned cost/reliability box, between coarser grid points, or appear when another uncertainty dimension is widened.

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

At the current reference profile, **62 policies** remain on the Pareto frontier:

```text
horizon 1:  3
horizon 2: 13
horizon 3: 46
```

First-action representation on that frontier is:

```text
E1_ARID5B_OCCUPANCY = 15
E2_CUX1_OCCUPANCY   = 15
E3_ALLELE_SNATAC    = 12
E4_CAPTURE_C        = 12
E5_MATCHED_MULTIOMIC=  8
```

So E5 can own the full **one-step minimax-regret region** while E1/E2/E3/E4 still remain non-dominated once expected information, expected cost, robust first-step regret, and future branching are considered together.

This is exactly why planning horizon remains load-bearing: the earlier depth-3 maximum-information policy starts with E3 even though E3 does not win the current one-step minimax box.

## Scientific boundary

The joint frontier cannot write or modify:

- `cause_found`
- `causal_status`
- `edge_status`
- `materialized`
- `verification_status`

No number of planning-grid evaluations can turn a hypothesis into evidence. New causal status still requires material evidence and the independent-verification path enforced by the causal processor.
