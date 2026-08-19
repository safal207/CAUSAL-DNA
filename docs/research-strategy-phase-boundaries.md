# Robust research strategy phase boundaries

CAUSAL-DNA now maps where a planning decision changes as assumptions about a target experiment change.

This is **not** a biological phase transition. The word `phase` is used only for a decision region in a declared planning-parameter space.

## Question

For the current GAP-001 experiment set, when does the robust first-step choice switch as the matched multi-omic experiment changes in:

- assay reliability;
- cost;
- outcome-probability uncertainty?

The target is `E5_MATCHED_MULTIOMIC`. Competitor experiments remain at their baseline reliability and baseline cost.

## Decision rule

For every phase cell:

1. vary E5 reliability and E5 cost according to the grid;
2. allow bounded pairwise transfers of outcome probability mass for every experiment;
3. sweep a deterministic simplex grid of planning priors;
4. compute lower/upper information-gain-per-cost utility intervals;
5. calculate interval regret against the best possible competitor;
6. choose the action with the smallest worst-case interval regret.

Conceptually:

```text
cell = (E5 reliability, E5 cost multiplier, probability-transfer radius)

for prior in prior_grid:
    utility_interval(action | cell, prior)
    regret(action) = best_possible_utility - action_lower_utility

winner(cell) = argmin_action max_prior regret(action)
```

`STOP` remains an explicit action.

## Grid

The first map uses:

```text
E5 reliability:       0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0
E5 cost multiplier:   0.5, 0.75, 1.0, 1.25, 1.5, 2.0
probability transfer: 0.0, 0.1, 0.2
```

This gives 126 phase cells.

Each cell uses a 165-point deterministic prior simplex (`12` total units, minimum `1` per hypothesis), for 20,790 prior-cell evaluations before per-experiment outcome scenarios are expanded.

The smaller prior grid is intentional. The full 1,140-prior sensitivity layer remains separate. This phase scan spends that computational budget on the additional reliability/cost/probability dimensions.

## Boundary semantics

A boundary is recorded only when two **adjacent declared grid cells** choose different winners.

Therefore a boundary such as:

```text
cost 1.25x -> winner E5
cost 1.50x -> winner E1
```

means only:

> under the declared model, a decision switch was detected somewhere inside the bracket `[1.25x, 1.50x]`.

It does **not** justify reporting a fabricated continuous threshold such as `1.37x`.

## Integrity boundary

Phase mapping is planning-only. It cannot write or imply:

- `cause_found`;
- `causal_status`;
- `edge_status`;
- `materialized`;
- `verification_status`.

No number of planning-grid evaluations can promote a biological mechanism. Only new material evidence and the independent-verification path can change causal status.

## Implementation

- `causal_dna/phase_boundary.py`
- `schemas/phase-boundary.schema.json`
- `cases/CDNA-001.phase-boundary.json`
- `tests/test_phase_boundary.py`
- `analysis/validate_phase_boundary.py`

The CI report prints winner counts, transition-pair counts, E5 cost ceilings on the declared grid, and the cells with the smallest winner-versus-runner-up regret margin.
