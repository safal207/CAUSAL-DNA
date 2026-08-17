# Multi-step Research Strategy Planner

CAUSAL-DNA can now rank not only the next experiment, but a short adaptive
research policy under a declared budget.

This is a planning layer only. It does not create biological evidence and it
cannot promote a causal edge.

## Objective

For current planning weights `H` and a candidate adaptive policy `pi`:

```text
expected information gain(pi)
  = entropy(current weights)
  - E[entropy(weights after executing pi)]
```

The planner also reports:

```text
information efficiency
  = expected information gain / expected total cost
```

Two objectives are supported:

- `information_gain` — maximize expected reduction in hypothesis uncertainty;
- `gain_per_cost` — maximize expected information gain per expected cost.

## Adaptive policy, not fixed sequence

A strategy is a tree:

```text
experiment E1
   |
   +-- outcome A --> choose E2
   |
   +-- outcome B --> choose E4
```

Different outcomes may therefore produce different next experiments. This is
stronger than pre-committing to a fixed sequence because the strategy can react
to what was learned.

## Budget semantics

`budget` is a maximum cumulative cost along any realized path. The planner
never returns a policy whose `max_path_cost` exceeds that value.

Expected cost is also tracked separately because different outcome branches may
stop at different times.

## Exact short-horizon search

The implementation enumerates adaptive policy branches up to depth 3 and
Pareto-prunes policies that are simultaneously no more informative and no less
costly than another policy.

The depth bound is deliberate. Planning probabilities and likelihood-like
multipliers are explicit assumptions, so very deep trees would amplify model
assumptions faster than evidence. A short horizon keeps the result inspectable
and suitable for re-planning after each real experimental outcome.

## Integrity boundary

The output contains only planning quantities:

- first experiment;
- contingent next experiments by outcome;
- expected terminal entropy;
- expected information gain;
- expected and maximum-path cost;
- information gain per cost.

It has no field for `cause_found`, `causal_status`, `materialized`, or
`verification_status`.

A real experiment still follows the authoritative path:

```text
strategy recommendation
        -> experiment
        -> experiment_outcome event + provenance
        -> adaptive replan
        -> independent verification if causal promotion is proposed
```

## CDNA-001 configurations

`cases/CDNA-001.strategy-configs.json` contains three reproducible planning
queries:

- two cost units, depth 2, maximize information gain;
- three cost units, depth 3, maximize gain per cost;
- 4.5 cost units, depth 3, maximize information gain.

The repository validator executes all three so the selected strategies are
visible in CI.

## Implementation

- `causal_dna/strategy_planner.py`
- `schemas/strategy-config.schema.json`
- `cases/CDNA-001.strategy-configs.json`
- `tests/test_strategy_planner.py`
- `analysis/validate_repository.py`
