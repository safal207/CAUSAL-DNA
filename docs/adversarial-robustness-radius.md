# Adversarial robustness radius

The adversarial robustness radius asks a stricter question than a bounded phase map:

> How far must declared planning assumptions move from the reference profile before the current minimax-regret first action changes?

For CDNA-001 the current target action is `E5_MATCHED_MULTIOMIC`.

## Metric

The first implementation uses a declared weighted L-infinity distance:

```text
d(theta, theta0) = max_i |theta_i - theta0_i| / scale_i
```

The scales are explicit configuration, not hidden calibration. The validator therefore reports both the normalized distance and raw boundary values.

Current adverse directions are:

- TF-family assay cost decreases (`E1`, with `E2` linked to the same stress axis);
- chromatin-family assay cost decreases (`E3`, with `E4` linked to the same stress axis);
- `E5` cost increases;
- `E5` reliability decreases.

These directions worsen E5 relative to its competitors under the declared one-step information-gain-per-cost minimax-regret model.

## Two boundary views

### Single-axis threshold

Each adverse axis is moved alone until the first-action winner changes or the declared search limit is reached.

A missing threshold means only:

> no switch was found inside that declared limit.

It does not imply infinite robustness.

### Coordinated worst-corner radius

All adverse axes move together at the same normalized radius. An axis saturates when it reaches its declared limit.

The first radius at which `E5` is no longer the minimax-regret action is refined by binary search to the configured tolerance.

This coordinated radius is a robustness property of the declared rectangular planning model. It is not a biological distance, confidence interval, or probability.

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
