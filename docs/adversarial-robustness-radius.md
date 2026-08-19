# Adversarial robustness radius

The adversarial robustness radius asks a stricter question than a bounded phase map:

> How far must declared planning assumptions move from the reference profile before the current minimax-regret first action changes?

For CDNA-001 the current target action is `E5_MATCHED_MULTIOMIC`.

## Metric

The implementation uses a declared weighted L-infinity distance:

```text
d(theta, theta0) = max_i |theta_i - theta0_i| / scale_i
```

The scales are explicit configuration, not hidden calibration. The validator therefore reports both the normalized distance and raw boundary values.

Current adverse directions are:

- TF-family assay cost decreases (`E1`, with `E2` linked to the same stress axis);
- chromatin-family assay cost decreases (`E3`, with `E4` linked to the same stress axis);
- `E5` cost increases;
- `E5` reliability decreases.

## Why the all-axis corner is not enough

Minimax interval regret is a relative decision criterion. Moving another apparently adverse coordinate can change the regret of several actions at once. Therefore the path where all coordinates move together is **not assumed to be globally worst**.

CAUSAL-DNA keeps the all-axis coordinated ray as a diagnostic, but the primary radius is the minimum first-switch distance across a declared family of normalized rays.

For the current configuration, every axis receives a fraction from:

```text
0.0, 0.5, 1.0
```

and only profiles with `max(fraction)=1.0` are retained. With four axes this yields 65 normalized rays.

This is still a finite ray family, not a proof over every possible continuous direction in parameter space. That limitation is explicit.

## First-switch detection and re-entry

A second subtlety is that a ray can show:

```text
E5 -> competitor -> E5
```

as the perturbation grows. Looking only at the ray endpoint would miss the earlier decision boundary.

Each ray therefore uses:

```text
coarse scan -> first non-E5 bracket -> binary refinement
```

The coarse scan resolution and binary tolerance are both declared in configuration. The reported radius is a **scan-refined first-switch radius** over the declared ray family.

## Boundary views

### Minimum ray-family radius

The primary result is the smallest first-switch radius found across all normalized fractional rays. The result records:

- the normalized distance;
- the ray fractions;
- the replacement winner;
- the target-versus-winner regret gap;
- the raw cost/reliability settings at the boundary.

### Single-axis thresholds

Each adverse axis is also moved alone until the first-action winner changes or the declared search limit is reached.

A missing threshold means only:

> no switch was found inside that declared limit.

It does not imply infinite robustness.

### All-axis coordinated radius

All adverse axes move together at the same normalized radius. This remains useful for comparison but is no longer treated as the definition of the adversarial radius.

## Family-linked stress settings

`E2_CUX1_OCCUPANCY` inherits the TF-family cost/reliability settings from `E1_ARID5B_OCCUPANCY`, and `E4_CAPTURE_C` inherits the chromatin-family settings from `E3_ALLELE_SNATAC`.

This prevents a symmetric baseline competitor from receiving an artificial advantage merely because it was excluded from the stress coordinates. It does **not** claim that ARID5B and CUX1 assays, or snATAC and Capture-C assays, are biologically equivalent.

## Scientific boundary

This layer cannot write or modify:

- `cause_found`
- `causal_status`
- `edge_status`
- `materialized`
- `verification_status`

Changing the planning winner cannot resolve GAP-001. Causal promotion still requires material evidence and independent verification.
