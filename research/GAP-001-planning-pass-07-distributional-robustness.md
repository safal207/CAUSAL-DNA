# GAP-001 planning pass 07 — distributional robustness

## Question

Does the earlier minimax recommendation `E5_MATCHED_MULTIOMIC` remain the safest first action when uncertainty is expanded beyond hypothesis priors to include the experiment model itself?

## Declared planning uncertainty set

This pass uses the same 1,140 deterministic prior states as `CDNA-001-SENS-001` and additionally stresses:

```text
outcome probability mass transfer <= 0.15
assay reliability                  = 0.6 .. 1.0
cost multiplier                    = 0.75 .. 1.25
```

The objective is information gain per cost unit.

These ranges are explicit planning stress assumptions. They are **not** biologically calibrated confidence intervals and are not evidence.

## Result

The distributional minimax-interval-regret action remains:

```text
E5_MATCHED_MULTIOMIC
```

with:

```text
max interval regret = 0.5585 bits/cost-unit
mean interval regret = 0.4190
```

Other maximum interval regrets:

```text
E1_ARID5B_OCCUPANCY = 0.5736
E2_CUX1_OCCUPANCY   = 0.5736
E3_ALLELE_SNATAC    = 0.5736
E4_CAPTURE_C        = 0.5736
STOP                = 0.5736
```

Thus the minimax margin over the alternatives is only:

```text
0.5736 - 0.5585 = 0.0151 bits/cost-unit
```

## Important negative result

`robust_optimal_points = 0` for **every** action.

This means no experiment's worst plausible utility dominates every competitor's best plausible utility at any one of the 1,140 prior states under the declared rectangular uncertainty set.

So the correct interpretation is **not**:

> E5 is universally robust.

It is:

> E5 has the smallest worst-case interval regret among the available first actions under this declared uncertainty set, but the decision is fragile and uncertainty intervals overlap strongly.

## Mean-vs-worst-case tension

The broad E5 assay wins the minimax criterion, but cheaper focused assays have lower mean interval regret:

```text
E1/E2 mean interval regret = 0.4018
E3/E4                     = 0.4080
E5                        = 0.4190
```

This exposes a genuine decision tradeoff:

- E5 protects against the worst planning-model misspecification;
- E1/E2 can be preferable under average-regret criteria;
- the recommendation therefore depends on risk attitude, not biological truth.

## Scientific boundary

This planning pass does not change the causal state:

```text
GAP-001 = OPEN
cause_found = false
```

No new material biological evidence was introduced.

## Next useful question

The next robustness target is a **break-even / phase-boundary map**: identify the combinations of cost uncertainty, assay reliability, and outcome-probability uncertainty at which the minimax first action switches between E1/E2, E3/E4, and E5.
