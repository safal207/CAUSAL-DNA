# Minimax-regret robust planning

Sensitivity analysis asks how often a planning policy wins across alternative priors. That is useful, but frequency alone does not answer a second question:

> If the planning priors are wrong, which first experiment is the least damaging choice?

CAUSAL-DNA therefore adds a separate robust first-action layer.

## Decision rule

For each prior state `p` in the deterministic sensitivity grid and each affordable first action `a`, compute one-step expected information gain:

```text
IG(a, p)
```

Let the best available action at that prior be:

```text
IG*(p) = max_a IG(a, p)
```

The regret of choosing action `a` is:

```text
R(a, p) = IG*(p) - IG(a, p)
```

The minimax-regret action is:

```text
argmin_a max_p R(a, p)
```

A second conservative criterion is maximin information gain:

```text
argmax_a min_p IG(a, p)
```

The two criteria are reported separately because they answer different questions.

## Why this is first-step only

The depth-3 strategy planner evaluates contingent policy trees. The minimax-regret layer intentionally evaluates the **first action only**. This prevents long-horizon planning assumptions from being presented as a robustness guarantee.

The layers therefore answer different questions:

- sensitivity: how often does a strategy win across priors?
- depth-3 strategy: what contingent policy is best under the declared model?
- minimax regret: what first experiment has the smallest worst-case loss relative to the locally best first experiment?
- maximin gain: what first experiment preserves the largest worst-case information gain?

## STOP is explicit

`STOP` is included as a legal action with zero information gain. This matters because some highly concentrated planning priors can make every modeled experiment yield zero expected entropy reduction under the declared outcome model.

## Scientific boundary

Minimax regret and maximin gain are **planning robustness metrics only**. They do not:

- promote a causal edge;
- materialize a Bardo state;
- set `cause_found=true`;
- create biological evidence;
- replace independent verification.

The analysis varies planning priors, not biology.

## Implementation

- `causal_dna/robust_decision.py` — robust first-step analyzer;
- `cases/CDNA-001.robust-decision.json` — analysis configuration;
- `schemas/robust-decision.schema.json` — machine-readable configuration contract;
- `analysis/validate_robust_decision.py` — CI gate and report;
- `tests/test_robust_decision.py` — invariants and minimax/maximin consistency tests.
