# GAP-001 Planning Pass 05 — Strategy robustness across planning priors

## Question

Is the budget-4.5, depth-3 strategy selected under equal planning weights —

```text
E3_ALLELE_SNATAC
  accessibility_shift -> E5_MATCHED_MULTIOMIC
  no_accessibility_shift -> E1_ARID5B_OCCUPANCY
```

— robust to reasonable changes in the hypothesis planning priors?

## Method

This is a **planning sensitivity analysis**, not biological evidence.

We enumerated the deterministic positive simplex grid

```text
H1 + H2 + H3 + H4 = 21
Hi >= 1
```

for the four currently open planner hypotheses:

- H1_ARID5B_REUSE
- H2_CUX1
- H3_ACCESSIBILITY
- H4_3D_CONTACT

The grid contains exactly **1,140 distinct prior states**. For every state we reran the same budgeted depth-3 strategy search with:

- budget = 4.5 cost units;
- objective = expected information gain;
- experiment costs held fixed;
- outcome planning probabilities held fixed;
- outcome likelihood-like multipliers held fixed.

The sensitivity grid therefore isolates **prior-weight sensitivity** only.

## Result

The equal-weight baseline still selects:

```text
E3_ALLELE_SNATAC
  accessibility_shift -> E5_MATCHED_MULTIOMIC
  no_accessibility_shift -> E1_ARID5B_OCCUPANCY
```

Across all 1,140 priors, first-step winners were:

| First action | Count | Share |
|---|---:|---:|
| E3_ALLELE_SNATAC | 346 | 30.35% |
| E4_CAPTURE_C | 266 | 23.33% |
| E1_ARID5B_OCCUPANCY | 220 | 19.30% |
| E2_CUX1_OCCUPANCY | 164 | 14.39% |
| E5_MATCHED_MULTIOMIC | 142 | 12.46% |
| STOP | 2 | 0.18% |

The baseline first action `E3_ALLELE_SNATAC` is also the **modal** first action, but its support is only **30.35%**. Under the declared operational robustness threshold of 60%, this is **not first-step robust**.

The exact baseline contingent policy is also the **modal policy signature**, but it appears in only **10.00%** of the prior grid. Under the declared 25% policy threshold, it is **not policy-signature robust**.

Machine validator output:

```text
priors=1140
baseline_first=E3_ALLELE_SNATAC
baseline_first_support=0.3035
dominant_first=E3_ALLELE_SNATAC
dominant_first_support=0.3035
baseline_policy_support=0.1000
first_step_robust=false
policy_robust=false
```

## Interpretation

The important finding is **not** that E3 is a bad experiment. E3 remains the single most frequently optimal first action over the tested prior region and is the equal-prior optimum.

The finding is that there is no dominant planning policy across this region. The recommendation landscape is fragmented among chromatin accessibility, 3D contact, and TF-occupancy assays.

Therefore:

> `E3 -> E5/E1` should be treated as the current equal-prior optimum, not as a prior-insensitive research strategy.

This result also exposes a second planning dimension: horizon matters. At depth 2, a preliminary sweep favored `E5_MATCHED_MULTIOMIC` much more strongly, while the original E3 policy emerges at depth 3. Strategy choice is therefore sensitive both to priors and to lookahead horizon.

## What this does not change

Nothing in this pass changes the biological causal record:

```text
GAP-001 = OPEN
cause_found = false
```

No hypothesis is promoted, rejected, materialized, or independently verified by this analysis.

## Next decision problem

The next planning layer should not merely choose the modal strategy. It should explicitly optimize **robustness under model uncertainty**, for example by comparing:

1. minimax regret across the prior grid;
2. worst-case information gain;
3. distributionally robust expected gain;
4. value of obtaining better prior/outcome-model calibration before committing to an expensive experiment.

That would distinguish "best under one planning model" from "safest research decision across plausible planning models."
