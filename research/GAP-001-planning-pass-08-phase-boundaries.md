# GAP-001 planning pass 08 — research strategy phase boundaries

## Question

Where does the preferred first experiment switch as the cost and reliability of the matched multi-omic assay change?

This pass maps **planning recommendation boundaries**, not biological phase transitions.

Objective:

```text
conservative utility = worst admissible one-step information gain / effective cost
```

The surface varies:

- E5 matched-multiomic cost multiplier: 0.4 to 1.6;
- E5 reliability: 0.5 to 1.0;
- outcome-probability mass-transfer bound: 0.0, 0.05, 0.10, 0.15;
- the same 1,140 deterministic planning priors used by the sensitivity analysis.

Competitor reliability is held at 1.0 for this phase-map pass so the E5 cost/reliability boundary remains interpretable.

## Result

The 312-point phase surface contains only two recommendation phases under this objective:

```text
baseline equal-prior phase occupancy:
E1_ARID5B_OCCUPANCY = 291 / 312
E5_MATCHED_MULTIOMIC = 21 / 312

dominant winner across the 1,140-prior grid:
E1_ARID5B_OCCUPANCY = 308 / 312
E5_MATCHED_MULTIOMIC =   4 / 312
```

`E2`, `E3`, `E4`, and `STOP` do not become phase winners in this scanned surface.

This is not in conflict with the earlier depth-3 result where E3 starts the best adaptive policy: that result optimizes multi-step expected information gain under a budget. The present map optimizes conservative **one-step information per cost**.

## Exact baseline-prior E5 ↔ E1 break-even cost

With no outcome-probability shift:

| E5 reliability | break-even multiplier | break-even absolute cost |
|---:|---:|---:|
| 0.50 | 0.2824 | 0.8472 |
| 0.60 | 0.3537 | 1.0611 |
| 0.70 | 0.4235 | 1.2705 |
| 0.80 | 0.4912 | 1.4737 |
| 0.90 | 0.5567 | 1.6702 |
| 1.00 | 0.6199 | 1.8597 |

With outcome-probability mass-transfer uncertainty of 0.15:

| E5 reliability | break-even multiplier | break-even absolute cost |
|---:|---:|---:|
| 0.50 | 0.2326 | 0.6979 |
| 0.60 | 0.2911 | 0.8733 |
| 0.70 | 0.3482 | 1.0446 |
| 0.80 | 0.4035 | 1.2105 |
| 0.90 | 0.4569 | 1.3707 |
| 1.00 | 0.5083 | 1.5249 |

Interpretation under the declared planning model:

- at reliability 1.0 and no outcome uncertainty, E5 beats E1 on conservative information-per-cost only below cost 1.8597 (61.99% of its declared baseline cost 3.0);
- at reliability 1.0 and probability-transfer bound 0.15, that ceiling falls to 1.5249 (50.83% of baseline);
- at reliability 0.8, the corresponding ceilings are 1.4737 and 1.2105;
- lower assay reliability and wider outcome uncertainty both shrink the economically rational E5 region.

## Grid slices

At transfer 0.0, reliability 1.0:

```text
E5 phase: cost multiplier 0.40–0.60
E1 phase: cost multiplier 0.70–1.60
```

At transfer 0.15, reliability 1.0:

```text
E5 phase: cost multiplier 0.40–0.50
E1 phase: cost multiplier 0.60–1.60
```

Across the full prior grid the E5 region is substantially smaller. For example, at transfer 0.0 and reliability 1.0, E5 is the dominant prior-grid winner only at multiplier 0.40; at 0.50 and above E1 dominates.

## Minimax-regret boundary layer

A second layer asks a different question: which first action minimizes the largest interval regret across alternative planning priors?

`causal_dna/phase_boundary.py` scans:

```text
E5 reliability:       0.40 .. 1.00
E5 cost multiplier:   0.50, 0.75, 1.00, 1.25, 1.50, 2.00
probability transfer: 0.00, 0.10, 0.20
priors per cell:       165
phase cells:           126
prior-cell states:     20,790
```

The robust scan again contains only two winner phases:

```text
E1_ARID5B_OCCUPANCY  = 61 / 126 cells
E5_MATCHED_MULTIOMIC = 65 / 126 cells
```

All **73 adjacent decision-boundary edges** are `E1 <-> E5` switches. E3 and E4 never become minimax-regret winners on this one-step surface.

### E5 winning cost ceiling on the sampled minimax grid

`none` means E5 never wins at a sampled cost. `2.00x+` means E5 still wins at the largest sampled cost, so the true switch — if one exists — lies beyond the scan.

At transfer `0.00`:

```text
r=0.40 none
r=0.50 none
r=0.60 none
r=0.70 none
r=0.80 none
r=0.90 0.50x
r=1.00 1.50x
```

At transfer `0.10`:

```text
r=0.40 1.50x
r=0.50 2.00x+
r=0.60 2.00x+
r=0.70 2.00x+
r=0.80 2.00x+
r=0.90 2.00x+
r=1.00 2.00x+
```

At transfer `0.20`:

```text
r=0.40 none
r=0.50 none
r=0.60 none
r=0.70 none
r=0.80 2.00x+
r=0.90 2.00x+
r=1.00 2.00x+
```

The non-monotonic pattern across transfer radii is a property of this discrete interval-regret uncertainty model. It must not be extrapolated into a smooth continuous law without a denser uncertainty-set analysis.

## Cross-objective interpretation

CAUSAL-DNA now preserves several intentionally different decision views:

```text
maximum raw one-step information gain      -> E5
maximum one-step information gain / cost   -> E1
best depth-3 adaptive policy               -> E3 -> E5 / E1
minimax regret over planning priors        -> E5
distributional minimax interval regret     -> E5 (narrow margin)
conservative efficiency phase map          -> E1 except low-cost/high-reliability E5 region
minimax-regret phase map                    -> E1 / E5 boundary depends strongly on reliability and uncertainty radius
```

There is therefore no single context-free "best experiment". The recommended first action is a function of:

```text
BestAction = f(
  objective,
  planning horizon,
  prior,
  outcome uncertainty,
  assay reliability,
  cost
)
```

A particularly useful negative result is that **E3 is the best depth-3 first move but never wins either one-step phase surface**. Planning horizon is therefore itself load-bearing context and must be preserved in any recommendation record.

## Scientific status

None of these phase boundaries are biological evidence. They do not resolve the adult posterior-hypothalamus mediator and do not alter causal status.

```text
GAP-001 = OPEN
cause_found = false
```
