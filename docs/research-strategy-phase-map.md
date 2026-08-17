# Research strategy phase map

CAUSAL-DNA distinguishes a **recommendation phase boundary** from a biological phase transition.

A recommendation phase boundary is simply the point where the preferred next experiment changes as declared planning assumptions vary.

For this first map, the target experiment is `E5_MATCHED_MULTIOMIC`. The axes are:

- target experiment cost multiplier;
- target assay reliability;
- bounded outcome-probability uncertainty;
- planning prior over the open mechanistic hypotheses.

The objective is conservative one-step efficiency:

```text
utility(E) = min admissible outcome-probability shift EIG(E) / effective_cost(E)
```

The map reports two views at every point:

1. **baseline-prior winner** — the best action under the experiment plan's declared equal planning weights;
2. **dominant-prior winner** — the action that wins on the largest fraction of the same 1,140 deterministic prior states used by the sensitivity analysis.

Neither view is biological probability or causal confidence.

## Break-even price

For the baseline prior, the target experiment's break-even cost is solved directly:

```text
conservative_gain(target) / break_even_cost
    = utility(best competitor)
```

Therefore:

```text
break_even_cost = conservative_gain(target) / utility(best competitor)
```

and the break-even multiplier is that cost divided by the target experiment's declared baseline cost.

This provides an interpretable statement such as:

> At reliability R and probability-transfer bound D, E5 is more efficient than the best competitor only below cost multiplier X.

## Why this is separate from minimax regret

The distributionally robust layer asks which action minimizes worst-case regret. This phase map asks a different question: where does **conservative information-per-cost** switch winner?

Those answers are allowed to differ. In fact, disagreement is informative:

- minimax regret may prefer a broad expensive assay because it is safer against model misspecification;
- information-per-cost may prefer a narrow assay because it buys more discrimination per unit cost.

CAUSAL-DNA preserves both rather than collapsing them into one score.

## Integrity rules

- phase-map parameters are planning assumptions only;
- recommendation switches cannot set `cause_found`;
- a phase boundary cannot materialize a Bardo state;
- no point in the map is evidence;
- STOP remains an explicit action;
- the same deterministic 1,140-prior grid is reused for comparability.

Implementation:

- `causal_dna/phase_map.py`
- `schemas/phase-map.schema.json`
- `cases/CDNA-001.phase-map.json`
- `analysis/validate_phase_map.py`
- `tests/test_phase_map.py`
