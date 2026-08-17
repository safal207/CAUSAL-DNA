# Distributionally robust experiment planning

CAUSAL-DNA now separates two robustness questions:

1. **prior robustness** — does the recommendation survive changes in hypothesis weights?
2. **model robustness** — does the recommendation survive changes in the experiment model itself?

The distributional layer addresses the second question.

## Uncertainty set

For every candidate first experiment, the planner varies four planning dimensions:

```text
Hypothesis prior
× outcome probabilities
× assay reliability
× experiment cost
```

The current CDNA-001 uncertainty set is deterministic and rectangular:

- all 1,140 prior states from `CDNA-001-SENS-001`;
- up to 0.15 probability mass transferred between outcome categories;
- assay reliability in `[0.6, 1.0]`;
- cost multiplier in `[0.75, 1.25]`.

These are planning stress parameters, not biological confidence intervals.

## Reliability semantics

Outcome multipliers are likelihood-like planning weights. Reliability shrinks them toward the neutral multiplier `1.0`:

```text
m_stressed = 1 + reliability * (m - 1)
```

Thus a less reliable assay becomes less discriminating without pretending that its biological result has changed.

## Utility interval

For each experiment and each prior state, the uncertainty set induces an interval:

```text
U(E, p) = [lower information-gain-per-cost, upper information-gain-per-cost]
```

Because uncertainty is rectangular across actions, the conservative interval regret is:

```text
Regret(E, p) = max_j U_upper(j, p) - U_lower(E, p)
```

The distributionally robust first action minimizes the largest such regret across all prior states.

## Why cost enters here

The earlier minimax-regret pass used information gain in bits. Once cost itself becomes uncertain, raw information gain no longer captures the decision problem. The distributional layer therefore optimizes information gain per cost unit.

## Scientific boundary

This entire layer is advisory planning infrastructure. It cannot:

- set `cause_found=true`;
- promote a causal edge;
- materialize a Bardo state;
- create verification evidence;
- overwrite the append-only event stream.

Only new material evidence plus the independent-verification path can change scientific causal status.

## Implementation

- `causal_dna/distributional_robust.py`
- `schemas/distributional-robustness.schema.json`
- `cases/CDNA-001.distributional-robustness.json`
- `tests/test_distributional_robust.py`
- `analysis/validate_distributional_robustness.py`
