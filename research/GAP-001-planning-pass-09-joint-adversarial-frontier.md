# GAP-001 planning pass 09 — joint adversarial assay frontier

## Question

Does the robust first-experiment recommendation survive when cost and reliability move for multiple assay families at the same time, rather than perturbing only E5?

This pass concerns **planning robustness only**. It does not add biological evidence or change causal status.

## Joint uncertainty box

Three independent planning axes are varied:

```text
E1_ARID5B_OCCUPANCY
E3_ALLELE_SNATAC
E5_MATCHED_MULTIOMIC
```

For each axis:

```text
cost multiplier = {0.75, 1.00, 1.25}
reliability     = {0.70, 0.85, 1.00}
```

Therefore:

```text
9 states per axis
9^3 = 729 joint configurations
```

Every cell is evaluated on a deterministic 165-state prior simplex with outcome-probability mass-transfer uncertainty `0.10`.

Total evaluation surface:

```text
729 x 165 = 120,285 prior-cell states
```

The cell decision rule is minimax interval regret over one-step information gain per cost.

## Methodological correction discovered during the pass

The first implementation independently stressed E1, E3, and E5 while leaving E2 and E4 at pristine baseline settings.

That run produced:

```text
E1 = 162 / 729
E2 = 414 / 729
E5 = 153 / 729
```

This was a useful diagnostic, not a biological result. In the current planning model E1/E2 have symmetric cost and discrimination structure, as do E3/E4. Stressing only one member of each pair gave the fixed twin an artificial advantage.

The production analysis therefore links planning settings:

```text
E2_CUX1_OCCUPANCY <- E1_ARID5B_OCCUPANCY axis
E4_CAPTURE_C      <- E3_ALLELE_SNATAC axis
```

The hypotheses remain distinct. This linkage does not assert biological equivalence. It only states that, for this synthetic stress test, currently symmetric assay assumptions must be stressed symmetrically.

The repository now tests this explicitly.

## Corrected result

After linking the symmetric planning families:

```text
E5_MATCHED_MULTIOMIC = 729 / 729 joint cells
boundary pairs       = none
boundary axes        = none
```

Thus the entire declared coarse uncertainty box is one one-step minimax-regret region for E5.

This is stronger than the earlier single-target result inside this particular box, but it is not a universal-optimality claim.

The nearest observed decision cells have only approximately:

```text
regret margin = 0.007654 bits/cost-unit
```

between E5 and E1 or E4. Therefore a switch may be nearby outside the box or between coarser parameter levels.

The correct conclusion is:

> E5 is robust throughout the currently declared linked-family box, but the measured safety margin near the edge is small enough that the robustness radius remains unresolved.

## Pareto policy frontier

A separate policy view keeps planning horizon visible and evaluates short adaptive strategies under the reference assay profile and budget `4.5`.

After pruning policies dominated in expected information gain, expected cost, and first-step robust regret, **62 policies** remain:

```text
horizon 1 =  3
horizon 2 = 13
horizon 3 = 46
```

First-action representation:

```text
E1_ARID5B_OCCUPANCY = 15
E2_CUX1_OCCUPANCY   = 15
E3_ALLELE_SNATAC    = 12
E4_CAPTURE_C        = 12
E5_MATCHED_MULTIOMIC=  8
```

This matters because the one-step robust box and the multi-step frontier answer different questions.

Current decision views now include:

```text
maximum raw one-step information gain        -> E5
maximum one-step information gain / cost     -> E1
best depth-3 expected-information policy     -> E3 -> E5 / E1
prior-only minimax regret                     -> E5
distributional minimax interval regret        -> E5, narrow margin
single-target conservative efficiency map     -> mostly E1
linked-family joint one-step minimax box       -> E5 in 729 / 729 cells
multi-objective short-horizon Pareto frontier  -> E1/E2/E3/E4/E5 all represented
```

There is no contradiction. The recommendation is a function of objective, uncertainty model, budget, and horizon.

## Engineering result

The final CI run passed:

```text
102 unit tests PASS
CAUSAL-DNA VALIDATION PASS
minimax-regret validation PASS
distributional robustness PASS
phase-map validation PASS
robust phase-boundary validation PASS
joint assay frontier validation PASS
three-space graph validation PASS
```

The joint analyzer caches information-gain intervals by prior, experiment, and reliability, allowing the 729-cell surface to remain practical in CI.

## Next causal-planning question

Because no corrected joint boundary lies inside the current box, brute-force expansion is less informative than directly searching for the nearest counterexample.

Define a baseline assumption state `theta_0` and search for:

```text
d* = min ||theta - theta_0||_W
     subject to winner(theta) != E5
```

This **adversarial robustness radius** would report:

- the smallest assumption change that dislodges E5;
- which action replaces it;
- which cost/reliability/outcome dimension carries the switch;
- a bracket around the nearest decision boundary;
- whether the counterexample persists under multiple priors.

That is the natural next pass.

## Scientific status

Nothing in this analysis resolves the adult posterior-hypothalamic mediator.

```text
GAP-001 = OPEN
cause_found = false
```
