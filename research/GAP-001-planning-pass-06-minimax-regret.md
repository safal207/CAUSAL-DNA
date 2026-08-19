# GAP-001 planning pass 06 — minimax regret

## Question

The depth-3 sensitivity pass showed that the baseline `E3_ALLELE_SNATAC -> (E5/E1)` policy is modal but prior-sensitive. This pass asks a different question:

> If the planning priors are wrong, which **first experiment** minimizes the worst loss in information gain relative to the locally best first experiment?

This is a planning-robustness question, not a biological causality claim.

## Method

Reuse the exact deterministic 1,140-point simplex prior grid from `CDNA-001-SENS-001`.

For every prior state `p` and first action `a`:

```text
R(a,p) = max_b IG(b,p) - IG(a,p)
```

Then compute:

```text
minimax regret = argmin_a max_p R(a,p)
maximin gain   = argmax_a min_p IG(a,p)
```

`STOP` is kept as an explicit zero-information action.

The robust layer is intentionally **first-step only**. It does not claim robustness for a depth-3 policy tree.

## CI result

Workflow run `32035902947` / job `95406142201`: PASS.

- 67/67 unit tests passed.
- repository evidence validator passed.
- robust-decision validator passed.
- `GAP-001` remains `OPEN`.
- `cause_found=false`.

### Action table

| action | worst-case IG | mean IG | max regret | mean regret | optimal prior points |
|---|---:|---:|---:|---:|---:|
| E1_ARID5B_OCCUPANCY | 0.0000 | 0.1603 | 0.4544 | 0.1554 | 125 |
| E2_CUX1_OCCUPANCY | 0.0000 | 0.1603 | 0.4544 | 0.1554 | 125 |
| E3_ALLELE_SNATAC | 0.0000 | 0.1762 | 0.3528 | 0.1395 | 131 |
| E4_CAPTURE_C | 0.0000 | 0.1762 | 0.3528 | 0.1395 | 131 |
| **E5_MATCHED_MULTIOMIC** | **0.0000** | **0.2839** | **0.3124** | **0.0318** | **692** |
| STOP | 0.0000 | 0.0000 | 0.4904 | 0.3157 | 14 |

## Result

Both robust criteria select:

```text
E5_MATCHED_MULTIOMIC
```

- minimax regret: `0.3124 bits` maximum regret;
- maximin information gain: tied at `0` worst-case gain across actions, then the deterministic tie-break selects E5 because it also has the smallest max regret and strongest mean gain.

The zero worst-case gain is itself informative: under some highly concentrated planning priors, every modeled experiment can have zero expected entropy reduction. Therefore **maximin gain alone is not discriminative on this grid**. Minimax regret is the more useful conservative criterion here.

## Interpretation

The planning layers now disagree in an interpretable way:

```text
one-step nominal IG          -> E5
one-step gain/cost           -> E1
three-step nominal strategy  -> E3 then conditional E5/E1
three-step modal sensitivity -> E3, but not robust majority
minimax-regret first action  -> E5
```

This is not a contradiction. Each objective answers a different decision problem.

The robust conclusion is limited to the declared planning model:

> If we want a first experiment that loses the least information relative to the locally best first experiment when our hypothesis priors are wrong, E5 is the safest modeled choice.

It does **not** mean E5 is biologically more causal, clinically preferable, or experimentally feasible in the real world without further cost/assay constraints.
