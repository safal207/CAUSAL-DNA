# GAP-001 planning pass 10 — adversarial robustness radius

## Question

How far must the declared planning assumptions move from the reference profile before the current one-step minimax-regret recommendation (`E5_MATCHED_MULTIOMIC`) first changes?

This pass measures **decision-model stability**, not biological confidence.

## Metric and search

The declared distance is weighted L-infinity:

```text
d(theta, theta0) = max_i |theta_i - theta0_i| / scale_i
```

Current scales:

- assay-family cost axis: `0.25` cost-multiplier units per radius unit;
- E5 cost axis: `0.25` cost-multiplier units per radius unit;
- E5 reliability axis: `0.10` reliability units per radius unit.

The four adverse coordinates are:

1. TF-family cost down (`E1`, with `E2` linked);
2. chromatin-family cost down (`E3`, with `E4` linked);
3. E5 cost up;
4. E5 reliability down.

The search does **not** assume that moving every adverse coordinate together is globally worst. It evaluates a normalized ray family with per-axis fractions:

```text
0.0, 0.5, 1.0
```

and retains profiles with `max(fraction)=1.0`, producing **65 rays**.

Each ray uses:

```text
coarse scan step = 0.05
-> first non-E5 bracket
-> binary refinement tolerance = 0.002
```

This first-switch scan is necessary because some rays show re-entry (`E5 -> competitor -> E5`) at larger perturbations.

## Baseline

At the reference profile:

```text
winner = E5_MATCHED_MULTIOMIC
regret margin to runner-up = 0.014509
```

## Nearest detected ray-family boundary

The minimum first-switch radius on the declared 65-ray family is approximately:

```text
d* = 1.8016
```

Ray fractions:

```text
TF_FAMILY_COST_DOWN       = 0.0
CHROMATIN_FAMILY_COST_DOWN = 1.0
E5_COST_UP                = 1.0
E5_RELIABILITY_DOWN       = 1.0
```

Boundary settings are approximately:

```text
TF-family cost multiplier        = 1.0000
chromatin-family cost multiplier = 0.5496
E5 cost multiplier               = 1.4504
E5 reliability                   = 0.8198
```

The deterministic winner at the first refined non-E5 point is `E4_CAPTURE_C`, with `E3_ALLELE_SNATAC` effectively tied at model precision because the present planning surface is symmetric under the H3/H4 exchange. Therefore the meaningful interpretation is:

> the recommendation switches from E5 to the **chromatin-family assays**, not specifically to a biologically preferred Capture-C mechanism.

The target-versus-winner max-regret gap immediately after the refined boundary is only about `3.7e-05`, as expected for a first-switch point.

## Single-axis thresholds

Only one single-axis switch is found inside the declared limits:

```text
CHROMATIN_FAMILY_COST_DOWN
normalized distance ~= 2.8156
cost multiplier ~= 0.2961
winner -> chromatin family (deterministic E4 tie-break)
```

No switch is found within the declared single-axis limits for:

```text
TF family cost down to 0.25x
E5 cost up to 2.00x
E5 reliability down to 0.40
```

This does **not** mean those directions are infinitely robust; only that the declared limits do not cross a first-action boundary when moved alone.

## All-axis coordinated ray

If all four adverse coordinates move together at fraction `1.0`, the first switch occurs later:

```text
d ~= 2.9203
winner -> E1_ARID5B_OCCUPANCY
```

Approximate settings:

```text
TF-family cost multiplier        = 0.2699
chromatin-family cost multiplier = 0.2699
E5 cost multiplier               = 1.7301
E5 reliability                   = 0.7080
```

The fact that this is farther than the minimum ray is exactly why the all-axis path must not define the robustness radius by itself.

## Re-entry regression

One declared ray with:

```text
TF fraction        = 0.5
chromatin fraction = 1.0
E5 cost fraction   = 0.0
E5 reliability     = 0.0
```

shows an early chromatin-family win but returns to E5 at the far endpoint. The test suite now requires the analyzer to preserve the **first switch** rather than incorrectly reporting “no boundary” from the endpoint alone.

## Interpretation

The current planning model says:

```text
single-factor E5 degradation        -> relatively stable
chromatin cost reduction alone      -> switch only after a large reduction
chromatin cheaper + E5 worse        -> substantially closer boundary
```

So the nearest observed vulnerability of the E5 recommendation is a **relative assay-economics / reliability interaction**, not one isolated parameter.

This is a planning result only. It does not promote H3/H4, E3/E4, or any causal edge.

## Scope limitation

`d* ~= 1.8016` is the minimum over the declared finite 65-ray family with a 0.05 first-switch scan and 0.002 local refinement. It is **not** claimed as the exact global minimum over every continuous direction in parameter space.

A future pass can perform adaptive continuous boundary optimization around this detected chromatin/E5 transition.

## Scientific status

```text
GAP-001 = OPEN
cause_found = false
```
